# ===========================
# s243130.jl  (CTSP / SOP solver)
# Iterated Local Search: kick=insertion, LS=2-opt + insertion
# CLI format:
#   julia s243130.jl <instance_path> <solution_path> <time_limit_seconds>
# ===========================

using Random

include(joinpath(@__DIR__, "InstanceReader.jl"))

# ===========================
# DATA STRUCTURES
# ===========================

struct CTSPInstance
    name::String
    ub::Int
    n::Int
    cost::Matrix{Int}                 # cost[i,j] = cost i -> j (NO -1)
    preds::Vector{Vector{Int}}        # preds[i] = nodes that must come before i
    start::Int
    goal::Int
end

mutable struct CTSPSolution
    route::Vector{Int}
    objective::Int
    visited::AbstractVector{Bool}
end

function CTSPSolution(n::Int)
    return CTSPSolution(Int[], 0, falses(n))
end

# ===========================
# 1) BUILD PRECEDENCES + CLEAN COST
# ===========================

# raw_cost:
#  - raw_cost[i,j] >= 0  -> actual travel cost
#  - raw_cost[i,j] == -1 -> precedence (j must come before i)
function build_preds_and_cost(raw_cost::AbstractMatrix{<:Integer})
    n = size(raw_cost, 1)

    preds = [Int[] for _ in 1:n]
    cost  = Matrix{Int}(raw_cost)

    BIG = typemax(Int) ÷ 4

    @inbounds for i in 1:n, j in 1:n
        if raw_cost[i, j] == -1
            push!(preds[i], j)
            cost[i, j] = BIG
        end
    end

    return preds, cost
end

# ===========================
# 2) FEASIBILITY (online construction)
# ===========================

@inline function feasible(v::Int, visited::AbstractVector{Bool}, preds::Vector{Vector{Int}})
    visited[v] && return false
    @inbounds for p in preds[v]
        visited[p] || return false
    end
    return true
end

# ===========================
# 3) NEAREST NEIGHBOR FEASIBLE
# ===========================

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

    @inbounds for v in 1:n
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

    return best
end

function nearest_neighbor_ctsp(inst::CTSPInstance)
    sol = CTSPSolution(inst.n)

    push!(sol.route, inst.start)
    sol.visited[inst.start] = true

    for step in 2:inst.n
        current = sol.route[end]
        remaining = inst.n - length(sol.route)

        nxt = get_nearest_feasible_neighbor(inst.cost, inst.preds, sol.visited, current;
                                            forbid_goal_until_end=true,
                                            goal=inst.goal,
                                            remaining=remaining)

        if nxt == 0
            # fallback: prova a prendere QUALUNQUE nodo feasible (anche goal) se sei bloccato
            @inbounds for v in 1:inst.n
                if feasible(v, sol.visited, inst.preds)
                    nxt = v
                    break
                end
            end
            nxt == 0 && error("Bloccato: nessun nodo feasible disponibile.")
        end

        sol.objective += inst.cost[current, nxt]
        push!(sol.route, nxt)
        sol.visited[nxt] = true
    end

    return sol
end

# ===========================
# 4) OBJECTIVE + ROUTE FEASIBILITY
# ===========================

function compute_objective(cost::AbstractMatrix{<:Integer}, route::Vector{Int})
    obj = 0
    @inbounds for i in 1:length(route)-1
        obj += Int(cost[route[i], route[i+1]])
    end
    return obj
end

function is_feasible_route(route::Vector{Int}, preds::Vector{Vector{Int}})
    pos = Dict{Int,Int}()
    @inbounds for (i, v) in enumerate(route)
        pos[v] = i
    end
    @inbounds for v in route
        for p in preds[v]
            if pos[p] > pos[v]
                return false
            end
        end
    end
    return true
end

# ===========================
# 5) LOCAL SEARCH: INSERTION (bounded)
# ===========================

function insertion_local_search(inst::CTSPInstance, route::Vector{Int}; deadline::Float64=Inf)
    best_route = copy(route)
    best_obj = compute_objective(inst.cost, best_route)

    improved = true
    while improved && time() < deadline
        improved = false
        n = length(best_route)

        for i in 2:n-1
            time() >= deadline && break
            for j in 2:n-1
                time() >= deadline && break
                i == j && continue

                candidate = copy(best_route)
                node = candidate[i]
                deleteat!(candidate, i)
                insert!(candidate, j, node)

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
# 6) LOCAL SEARCH: 2-OPT (bounded)
# ===========================

function two_opt_local_search(inst::CTSPInstance, route::Vector{Int}; deadline::Float64=Inf)
    best_route = copy(route)
    best_obj = compute_objective(inst.cost, best_route)

    improved = true
    while improved && time() < deadline
        improved = false
        n = length(best_route)

        for i in 2:n-2
            time() >= deadline && break
            for k in i+1:n-1
                time() >= deadline && break

                candidate = copy(best_route)
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
# 7) PERTURBATION: RANDOM FEASIBLE INSERTIONS (K accepted moves)
# ===========================

function perturbation_insertion(route::Vector{Int},
                                preds::Vector{Vector{Int}};
                                Kaccepted::Int=5,
                                max_trials::Int=500)

    pert = copy(route)
    n = length(pert)
    n <= 3 && return pert

    accepted = 0
    trials = 0

    while accepted < Kaccepted && trials < max_trials
        trials += 1

        i = rand(2:n-1)
        j = rand(2:n-1)
        i == j && continue

        cand = copy(pert)
        node = cand[i]
        deleteat!(cand, i)
        insert!(cand, j, node)

        if is_feasible_route(cand, preds)
            pert = cand
            accepted += 1
        end
    end

    return pert
end

# ===========================
# 8) ITERATED LOCAL SEARCH (bounded by time limit)
# ===========================

function iterated_local_search(inst::CTSPInstance,
                               init_route::Vector{Int};
                               time_limit::Float64=60.0,
                               K0::Int=3,
                               Kmin::Int=1,
                               Kmax::Int=50)

    deadline = time() + time_limit

    current_route, current_obj = two_opt_local_search(inst, init_route; deadline=deadline)
    current_route, current_obj = insertion_local_search(inst, current_route; deadline=deadline)

    best_route = copy(current_route)
    best_obj   = current_obj

    K = K0

    while time() < deadline
        kicked = perturbation_insertion(current_route, inst.preds; Kaccepted=K, max_trials=500)

        new_route, new_obj = two_opt_local_search(inst, kicked; deadline=deadline)
        new_route, new_obj = insertion_local_search(inst, new_route; deadline=deadline)

        if new_obj < current_obj
            current_route = new_route
            current_obj   = new_obj
            K = max(Kmin, K - 1)
        else
            K = min(Kmax, K + 1)
        end

        if current_obj < best_obj
            best_route = copy(current_route)
            best_obj   = current_obj
        end
    end

    return best_route, best_obj
end

# ===========================
# 9) OUTPUT WRITER (assignment format)
# ===========================

function write_solution(sol_path::String, route::Vector{Int})
    mkpath(dirname(sol_path))
    open(sol_path, "w") do io
        # one line, nodes separated by a space
        println(io, join(route, " "))
    end
end

# ===========================
# MAIN (CLI)
# ===========================

function main()
    if length(ARGS) != 3
        println("Usage:")
        println("  julia s243130.jl <instance_path> <solution_path> <time_limit_seconds>")
        exit(1)
    end

    inst_path = ARGS[1]
    sol_path  = ARGS[2]
    time_limit = parse(Float64, ARGS[3])

    name, ub, n, raw_cost = read_instance(inst_path)
    preds, cost = build_preds_and_cost(raw_cost)

    inst = CTSPInstance(name, ub, n, cost, preds, 1, n)

    # Build initial solution
    sol0 = nearest_neighbor_ctsp(inst)

    # Improve with ILS within time limit
    route_best, obj_best = iterated_local_search(inst, sol0.route; time_limit=time_limit)

    # Safety: ensure feasibility (should already be true)
    if !is_feasible_route(route_best, inst.preds)
        # fallback to initial feasible route if something goes wrong
        route_best = sol0.route
        obj_best = sol0.objective
    end

    # Write solution file in required format
    write_solution(sol_path, route_best)

    # (optional) minimal stderr/info for you; comment out if you want silent
    # println("Solved: ", inst.name, " | obj=", obj_best, " | UB=", inst.ub, " | gap=", obj_best - inst.ub)
end

main()