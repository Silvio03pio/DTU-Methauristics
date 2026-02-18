# ===========================
# CTSP.jl  (SOP solver)
# ===========================

# richiamo il file instance reader decidendo quale istanza leggere
include(joinpath(@__DIR__, "InstanceReader.jl"))

# ===========================
# DATI PER SOP 
# ===========================

struct CTSPInstance
    name::String
    ub::Int
    n::Int
    cost::Matrix{Int}                 # cost[i,j] = costo i -> j (NO -1)
    preds::Vector{Vector{Int}}        # preds[i] = lista predecessori che devono venire PRIMA di i
    start::Int                        # tipicamente 1
    goal::Int                         # tipicamente n
end

mutable struct CTSPSolution
    route::Vector{Int}                # path es. [1, ..., n]
    objective::Int                    # costo totale del path
    visited::AbstractVector{Bool}     # può essere BitVector o Vector{Bool}
end

# Costruttore comodo per una soluzione vuota
function CTSPSolution(n::Int)
    return CTSPSolution(Int[], 0, falses(n))
end

# perturbazione della initial solution
function perturb_route!(route::Vector{Int}, preds::Vector{Vector{Int}}; K::Int=5)
    n = length(route)
    
    for _ in 1:K
        # Random insertion: pick node, move to random feasible position
        i = rand(2:n-1)  # Don't move start/goal
        node = route[i]
        
        # Find new feasible position
        new_pos = rand(2:n)
        if try_insert_feasible(route, preds, i, new_pos)
            # Apply insertion
            deleteat!(route, i)
            insert!(route, new_pos, node)
        end
    end
end

# ===========================
# 1) COSTRUISTO matrice PRECEDENZE + COSTI 
# ===========================

# raw_cost viene dal file .sop:
#  - raw_cost[i,j] >= 0  -> costo vero
#  - raw_cost[i,j] == -1 -> vincolo di precedenza (j deve venire prima di i)
function build_preds_and_cost(raw_cost::AbstractMatrix{<:Integer})
    n = size(raw_cost, 1)

    preds = [Int[] for _ in 1:n]
    cost  = Matrix{Int}(raw_cost)             # normalizzo a Int

    BIG = typemax(Int) ÷ 4                    # infinito pratico

    for i in 1:n, j in 1:n
        if raw_cost[i, j] == -1
            # SOP convention: -1 in (i,j) => j precede i
            push!(preds[i], j)
            cost[i, j] = BIG                  # non è un costo, lo rendo "non scelto"
        end
    end

    return preds, cost
end



# ===========================
# 2) FEASIBILITY: posso visitare v adesso?
# ===========================

@inline function feasible(v::Int, visited::AbstractVector{Bool}, preds::Vector{Vector{Int}})
    visited[v] && return false
    for p in preds[v]
        visited[p] || return false
    end
    return true
end


# ===========================
# 3) NEAREST NEIGHBOR FEASIBLE
# ===========================

# trova il nearest neighbor tra i candidati AMMISSIBILI
function get_nearest_feasible_neighbor(cost::AbstractMatrix{<:Integer},
                                       preds::Vector{Vector{Int}},
                                       visited::AbstractVector{Bool},
                                       city::Int;
                                       forbid_goal_until_end::Bool=true,
                                       goal::Int=size(cost,1),
                                       remaining::Int=0)

    best = 0
    best_cost = typemax(Int)

    n = size(cost, 1)

    for v in 1:n

        # opzionale: non scegliere il goal fino all'ultimo step
        if forbid_goal_until_end && v == goal && remaining > 1
            continue
        end

        if feasible(v, visited, preds)
            c = Int(cost[city, v])
            if c < best_cost
                best_cost = c
                best = v
            end
        end
    end

    return best   # 0 se non esiste candidato feasible
end


# costruzione NN per SOP
function nearest_neighbor_ctsp(inst::CTSPInstance)
    sol = CTSPSolution(inst.n)

    # start fisso
    push!(sol.route, inst.start)
    sol.visited[inst.start] = true

    # costruisco fino a visitare tutti i nodi
    for step in 2:inst.n
        current = sol.route[end]
        remaining = inst.n - length(sol.route)    # <-- CORRETTO: quanti nodi mancano

        nxt = get_nearest_feasible_neighbor(inst.cost, inst.preds, sol.visited, current;
                                            forbid_goal_until_end=true,
                                            goal=inst.goal,
                                            remaining=remaining)

        if nxt == 0
            # debug utile
            println("ROUTE finora: ", sol.route)
            println("Nodi non visitati: ", findall(!, sol.visited))
            error("Bloccato: nessun nodo feasible disponibile.")
        end

        sol.objective += inst.cost[current, nxt]
        push!(sol.route, nxt)
        sol.visited[nxt] = true
    end

    return sol
end


# ===========================
# 4) LOCAL SEARCH (INSERTION) CON PRECEDENZE
# ===========================

# calcolo objective del PATH (non ciclo)
function compute_objective(cost::AbstractMatrix{<:Integer}, route::Vector{Int})
    obj = 0
    for i in 1:length(route)-1
        obj += Int(cost[route[i], route[i+1]])
    end
    return obj
end

# controllo se una route rispetta tutte le precedenze
function is_feasible_route(route::Vector{Int}, preds::Vector{Vector{Int}})
    pos = Dict{Int,Int}()
    for (i, v) in enumerate(route)
        pos[v] = i
    end

    for v in route
        for p in preds[v]
            if pos[p] > pos[v]
                return false
            end
        end
    end

    return true
end

# Local Search: reinserisco un nodo in un'altra posizione (solo se migliora + resta feasible)
function insertion_local_search(inst::CTSPInstance, route::Vector{Int})
    best_route = copy(route)
    best_obj = compute_objective(inst.cost, best_route)

    improved = true
    while improved
        improved = false
        n = length(best_route)

        # non muovo start (pos 1) e goal (pos n)
        for i in 2:n-1
            for j in 2:n-1
                i == j && continue

                candidate = copy(best_route)
                node = candidate[i]
                deleteat!(candidate, i)
                insert!(candidate, j, node)

                # controlla feasibility e objective
                if is_feasible_route(candidate, inst.preds)
                    cand_obj = compute_objective(inst.cost, candidate)

                    if cand_obj < best_obj
                        best_obj = cand_obj
                        best_route = candidate
                        improved = true
                        break
                    end
                end
            end
            improved && break
        end
    end

    return best_route, best_obj
end

# ===========================
# 6) 2-opt CON PRECEDENZE
# ===========================
function two_opt_local_search(inst::CTSPInstance, route::Vector{Int})
    best_route = copy(route)
    best_obj = compute_objective(inst.cost, best_route)

    improved = true
    while improved
        improved = false
        n = length(best_route)

        # non tocchiamo start (pos 1) e goal (pos n)
        for i in 2:n-2
            for k in i+1:n-1
                candidate = copy(best_route)

                # reverse del segmento [i:k]
                candidate[i:k] = reverse(candidate[i:k])

                if is_feasible_route(candidate, inst.preds)
                    cand_obj = compute_objective(inst.cost, candidate)
                    if cand_obj < best_obj
                        best_obj = cand_obj
                        best_route = candidate
                        improved = true
                        break
                    end
                end
            end
            improved && break
        end
    end

    return best_route, best_obj
end


# ===========================
# MAIN
# ===========================

function main()

    # 0) scelgo istanza
    inst_path = joinpath(@__DIR__, "Instances", "Instances", "br17.10.sop")

    # 1) leggo istanza
    name, ub, n, raw_cost = read_instance(inst_path)

    # 2) costruisco preds e cost pulita
    preds, cost = build_preds_and_cost(raw_cost)

    # 3) creo oggetto istanza CTSP
    inst = CTSPInstance(name, ub, n, cost, preds, 1, n)

    println("Instance: ", inst.name, " | n=", inst.n, " | UB=", inst.ub)

    # 4) costruzione NN feasible
    sol = nearest_neighbor_ctsp(inst)

    println("\n--- Nearest Neighbor ---")
    println("Objective: ", sol.objective)

    route_best = copy(sol.route)
    obj_best = sol.objective

    

    # 5) Local search: prima 2-opt (macro), poi insertion (fine tuning)
    improved = true
    while improved
        improved = false

        # 5.1) 2-opt
        route2, obj2 = two_opt_local_search(inst, route_best)
        if obj2 < obj_best
            route_best = route2
            obj_best = obj2
            improved = true
        end
        #5.2) insertion
        route3, obj3 = insertion_local_search(inst, route_best)
        if obj3 < obj_best
            route_best = route3
            obj_best = obj3
            improved = true
        end 
    end

    println("\n--- Final solution after LS ---")
    println("Route: ", route_best)
    println("Objective: ", obj_best)
    println("Gap vs UB: ", obj_best - inst.ub)
    #println("cost:", cost)
    #println("preds:", preds)

end

main()


