#da cambiare con quello del professore una volta fuori


############################################################
# TSP – Reactive GRASP + 2-opt Local Search + Path Relinking
# (PR with same neighborhood as LS = 2-opt) vs (PR with SWAP)
############################################################

using Random
using Printf

############################
# Instance reader
############################
function readInstance(filename)
    file = open(filename)
    name = split(readline(file))[2]
    readline(file); readline(file)
    dim = parse(Int32, split(readline(file))[2])
    readline(file); readline(file)

    coord = zeros(Float32, dim, 2)
    for i in 1:dim
        data = parse.(Float32, split(readline(file)))
        coord[i, :] = data[2:3]
    end
    close(file)
    return name, coord, dim
end

############################
# Distance matrix
############################
function getDistanceMatrix(coord::Array{Float32,2}, dim::Int32)
    dist = zeros(Float32, dim, dim)
    for i in 1:dim
        for j in 1:dim
            if i != j
                dist[i,j] = round(sqrt((coord[i,1]-coord[j,1])^2 + (coord[i,2]-coord[j,2])^2), digits=2)
            end
        end
    end
    return dist
end

############################
# Solver / Solution structs
############################
struct TSPSolver
    coord::Array{Float32,2}
    dim::Int32
    dist::Array{Float32,2}
    startTime::UInt64
    timelimit::UInt32
    TSPSolver(coord, dim, dist, timelimit) = new(coord, dim, dist, time_ns(), timelimit)
end

mutable struct TSPSolution
    route::Vector{Int32}
    objective::Float32
    TSPSolution(dim::Int32) = new(zeros(Int32, dim), 0.0f0)
end

############################
# Time
############################
function elapsed_time(m::TSPSolver)
    return (time_ns() - m.startTime) / 1e9  # seconds (Float64)
end

############################
# Basic utilities
############################
@inline edgekey(a::Int32, b::Int32) = (min(a,b), max(a,b))

function tour_cost(m::TSPSolver, route::Vector{Int32})
    n = length(route)
    c = 0.0f0
    for i in 1:n-1
        c += m.dist[route[i], route[i+1]]
    end
    c += m.dist[route[n], route[1]]
    return c
end

############################
# Nearest Neighbor
############################
function get_nearest_neighbor(dist, visited, city)
    bestd = Inf
    next = 0
    for i in 1:length(dist[city,:])
        if !visited[i] && bestd > dist[city,i]
            bestd = dist[city,i]
            next = i
        end
    end
    return next
end

function nearest_neighbor_heuristic(m::TSPSolver)
    sol = TSPSolution(m.dim)
    visited = falses(m.dim)

    sol.route[1] = 1
    visited[1] = true
    sol.objective = 0.0f0

    for i in 1:m.dim-1
        j = get_nearest_neighbor(m.dist, visited, sol.route[i])
        sol.route[i+1] = j
        visited[j] = true
        sol.objective += m.dist[sol.route[i], sol.route[i+1]]
    end

    # Close tour (FIXED: last -> first)
    sol.objective += m.dist[sol.route[end], sol.route[1]]
    return sol
end

############################
# 2-opt delta evaluation
############################
function eval_two_opt_step(m::TSPSolver, sol::TSPSolution, i::Int, j::Int)
    # assumption i < j
    n = length(sol.route)

    if i == 1 && j == n
        return sol.objective
    elseif i == 1
        return sol.objective -
               m.dist[sol.route[i], sol.route[n]] -
               m.dist[sol.route[j], sol.route[j+1]] +
               m.dist[sol.route[j], sol.route[n]] +
               m.dist[sol.route[i], sol.route[j+1]]
    elseif j == n
        return sol.objective -
               m.dist[sol.route[i-1], sol.route[i]] -
               m.dist[sol.route[j], sol.route[1]] +
               m.dist[sol.route[i-1], sol.route[j]] +
               m.dist[sol.route[i], sol.route[1]]
    else
        return sol.objective -
               m.dist[sol.route[i-1], sol.route[i]] -
               m.dist[sol.route[j], sol.route[j+1]] +
               m.dist[sol.route[i-1], sol.route[j]] +
               m.dist[sol.route[i], sol.route[j+1]]
    end
end

function best_two_opt(m::TSPSolver, sol::TSPSolution)
    iBest = 0
    jBest = 0
    bestObj = Inf32

    for i in 2:m.dim-1
        for j in i+1:m.dim
            obj = eval_two_opt_step(m, sol, i, j)
            if obj < bestObj
                bestObj = obj
                iBest = i
                jBest = j
            end
        end
    end

    return bestObj, iBest, jBest
end

function local_search!(m::TSPSolver, sol::TSPSolution)
    improved = true
    while improved && elapsed_time(m) < m.timelimit
        val, i, j = best_two_opt(m, sol)
        if i == 0
            break
        end
        if val < sol.objective
            reverse!(sol.route, i, j)
            sol.objective = val
            improved = true
        else
            improved = false
        end
    end
    return sol
end

function local_search(m::TSPSolver)
    sol = nearest_neighbor_heuristic(m)
    return local_search!(m, sol)
end

############################
# GRASP construction (alpha-RCL)
############################
function grasp_construction(m::TSPSolver; alpha::Float32=0.5f0)
    sol = TSPSolution(m.dim)
    visited = falses(m.dim)

    current = 1
    sol.route[1] = current
    visited[current] = true
    sol.objective = 0.0f0

    for pos in 2:m.dim
        c = Float32[]
        cand_ids = Int32[]
        for j in 1:m.dim
            if !visited[j]
                push!(c, m.dist[current, j])
                push!(cand_ids, Int32(j))
            end
        end

        c_min = minimum(c)
        c_max = maximum(c)
        threshold = c_min + alpha*(c_max - c_min)

        rcl = Int32[]
        for k in 1:length(c)
            if c[k] <= threshold
                push!(rcl, cand_ids[k])
            end
        end
        if isempty(rcl)
            # fallback: best
            idx = argmin(c)
            rcl = [cand_ids[idx]]
        end

        next_city = rand(rcl)
        sol.route[pos] = next_city
        visited[next_city] = true
        sol.objective += m.dist[current, next_city]
        current = next_city
    end

    sol.objective += m.dist[sol.route[end], sol.route[1]]
    return sol
end

############################
# Reactive alpha stats
############################
mutable struct AlphaStats
    alphas::Vector{Float32}
    all_costs::Vector{Vector{Float32}}
    counts::Vector{Int}
    window::Int
end

function AlphaStats(alphas::Vector{Float32}, window::Int)
    AlphaStats(alphas, [Float32[] for _ in 1:length(alphas)], zeros(Int, length(alphas)), window)
end

function update_stats!(stats::AlphaStats, alpha_idx::Int, cost::Float32)
    push!(stats.all_costs[alpha_idx], cost)
    if length(stats.all_costs[alpha_idx]) > stats.window
        popfirst!(stats.all_costs[alpha_idx])
    end
end

function select_alpha!(stats::AlphaStats)
    m = length(stats.alphas)

    s_star = Inf32
    has_data = false
    for costs in stats.all_costs
        if !isempty(costs)
            s_star = min(s_star, minimum(costs))
            has_data = true
        end
    end

    if !has_data
        p = fill(1.0f0/m, m)
        idx = argmax(cumsum(p) .>= rand(Float32))
        stats.counts[idx] += 1
        return stats.alphas[idx], idx
    end

    q = Float32[]
    for i in 1:m
        costs_i = stats.all_costs[i]
        if isempty(costs_i)
            push!(q, 1.0f0)
        else
            A_i = sum(costs_i) / length(costs_i)
            push!(q, s_star / A_i)
        end
    end
    p = q / sum(q)

    idx = argmax(cumsum(p) .>= rand(Float32))
    stats.counts[idx] += 1
    return stats.alphas[idx], idx
end

############################
# Edge sets for similarity
############################
function edgeset(route::Vector{Int32})
    n = length(route)
    E = Set{Tuple{Int32,Int32}}()
    for i in 1:n-1
        push!(E, edgekey(route[i], route[i+1]))
    end
    push!(E, edgekey(route[n], route[1]))
    return E
end

common_edges(routeA::Vector{Int32}, routeB::Vector{Int32}) = length(intersect(edgeset(routeA), edgeset(routeB)))

############################
# PR with SAME neighborhood as LS: 2-opt
# Choose 2-opt moves that increase overlap with guiding edges
############################
function delta_common_two_opt(route::Vector{Int32}, i::Int, j::Int, Eg::Set{Tuple{Int32,Int32}})
    n = length(route)

    if i == 1 && j == n
        return 0
    elseif i == 1
        removed1 = edgekey(route[1], route[n])
        removed2 = edgekey(route[j], route[j+1])
        added1   = edgekey(route[j], route[n])
        added2   = edgekey(route[1], route[j+1])
    elseif j == n
        removed1 = edgekey(route[i-1], route[i])
        removed2 = edgekey(route[n], route[1])
        added1   = edgekey(route[i-1], route[n])
        added2   = edgekey(route[i], route[1])
    else
        removed1 = edgekey(route[i-1], route[i])
        removed2 = edgekey(route[j], route[j+1])
        added1   = edgekey(route[i-1], route[j])
        added2   = edgekey(route[i], route[j+1])
    end

    inc = (added1 in Eg ? 1 : 0) + (added2 in Eg ? 1 : 0)
    dec = (removed1 in Eg ? 1 : 0) + (removed2 in Eg ? 1 : 0)
    return inc - dec
end

function path_relinking_two_opt(m::TSPSolver, sol_start::TSPSolution, sol_guide::TSPSolution;
                                max_steps::Int=200, do_ls::Bool=true)
    cur = TSPSolution(m.dim)
    cur.route .= sol_start.route
    cur.objective = sol_start.objective

    best = TSPSolution(m.dim)
    best.route .= cur.route
    best.objective = cur.objective

    Eg = edgeset(sol_guide.route)
    cur_common = common_edges(cur.route, sol_guide.route)

    steps = 0
    while steps < max_steps
        steps += 1
        if cur_common == m.dim
            break
        end

        best_i = 0
        best_j = 0
        best_dcom = -10^9
        best_obj = Inf32

        for i in 2:m.dim-1
            for j in i+1:m.dim
                dcom = delta_common_two_opt(cur.route, i, j, Eg)
                if dcom < 0
                    continue
                end
                obj = eval_two_opt_step(m, cur, i, j)

                if (dcom > best_dcom) || (dcom == best_dcom && obj < best_obj)
                    best_dcom = dcom
                    best_obj = obj
                    best_i = i
                    best_j = j
                end
            end
        end

        if best_i == 0
            break
        end

        reverse!(cur.route, best_i, best_j)
        cur.objective = best_obj
        cur_common += best_dcom

        if cur.objective < best.objective
            best.route .= cur.route
            best.objective = cur.objective
        end
    end

    if do_ls
        local_search!(m, best)
    end
    return best
end

############################
# Alternative PR neighborhood: SWAP two cities (positions)
# (Your “one of your choice” neighborhood)
############################
function eval_swap_step(m::TSPSolver, route::Vector{Int32}, obj::Float32, i::Int, j::Int)
    # swap positions i and j in the tour. i<j.
    # We keep city 1 fixed by only calling with i,j >= 2.
    n = length(route)

    # helper to get route positions with wrap
    prevpos(p) = (p == 1 ? n : p-1)
    nextpos(p) = (p == n ? 1 : p+1)

    a = route[i]
    b = route[j]

    # Collect affected directed edges (u->v) in current tour
    affected = Set{Tuple{Int32,Int32}}()

    function add_edge(u::Int32, v::Int32)
        push!(affected, (u,v))
    end

    # edges around i and j (before swap)
    ip = prevpos(i); inx = nextpos(i)
    jp = prevpos(j); jnx = nextpos(j)

    add_edge(route[ip], route[i])
    add_edge(route[i], route[inx])
    add_edge(route[jp], route[j])
    add_edge(route[j], route[jnx])

    # If adjacent, some edges overlap; Set handles duplicates

    # compute old sum on affected edges
    oldsum = 0.0f0
    for (u,v) in affected
        oldsum += m.dist[u,v]
    end

    # simulate swap in-place on a copy of just the needed neighborhoods
    # easiest: read neighbors after swap by logic
    # We'll create a tiny function that returns city at position p after swap:
    city_after(p) = (p == i ? b : (p == j ? a : route[p]))

    # recompute new affected edges with swapped cities
    affected2 = Set{Tuple{Int32,Int32}}()
    add2(u,v) = push!(affected2, (u,v))

    add2(city_after(ip), city_after(i))
    add2(city_after(i), city_after(inx))
    add2(city_after(jp), city_after(j))
    add2(city_after(j), city_after(jnx))

    newsum = 0.0f0
    for (u,v) in affected2
        newsum += m.dist[u,v]
    end

    return obj - oldsum + newsum
end

function delta_common_swap(route::Vector{Int32}, i::Int, j::Int, Eg::Set{Tuple{Int32,Int32}})
    # approximate exact delta common-edges by checking only changed undirected edges around swapped positions
    n = length(route)

    prevpos(p) = (p == 1 ? n : p-1)
    nextpos(p) = (p == n ? 1 : p+1)

    a = route[i]
    b = route[j]

    ip = prevpos(i); inx = nextpos(i)
    jp = prevpos(j); jnx = nextpos(j)

    # edges BEFORE (undirected)
    before = Set{Tuple{Int32,Int32}}()
    push!(before, edgekey(route[ip], route[i]))
    push!(before, edgekey(route[i], route[inx]))
    push!(before, edgekey(route[jp], route[j]))
    push!(before, edgekey(route[j], route[jnx]))

    city_after(p) = (p == i ? b : (p == j ? a : route[p]))

    after = Set{Tuple{Int32,Int32}}()
    push!(after, edgekey(city_after(ip), city_after(i)))
    push!(after, edgekey(city_after(i), city_after(inx)))
    push!(after, edgekey(city_after(jp), city_after(j)))
    push!(after, edgekey(city_after(j), city_after(jnx)))

    # delta common = (#after∩Eg) - (#before∩Eg) restricted to changed edges
    inc = 0
    dec = 0
    for e in after
        inc += (e in Eg) ? 1 : 0
    end
    for e in before
        dec += (e in Eg) ? 1 : 0
    end
    return inc - dec
end

function path_relinking_swap(m::TSPSolver, sol_start::TSPSolution, sol_guide::TSPSolution;
                             max_steps::Int=200, do_ls::Bool=true)
    cur = TSPSolution(m.dim)
    cur.route .= sol_start.route
    cur.objective = sol_start.objective

    best = TSPSolution(m.dim)
    best.route .= cur.route
    best.objective = cur.objective

    Eg = edgeset(sol_guide.route)
    cur_common = common_edges(cur.route, sol_guide.route)

    steps = 0
    while steps < max_steps
        steps += 1
        if cur_common == m.dim
            break
        end

        best_i = 0
        best_j = 0
        best_dcom = -10^9
        best_obj = Inf32

        # swap neighborhood (restrict: keep position 1 fixed => i,j >=2)
        for i in 2:m.dim-1
            for j in i+1:m.dim
                dcom = delta_common_swap(cur.route, i, j, Eg)
                if dcom < 0
                    continue
                end
                obj = eval_swap_step(m, cur.route, cur.objective, i, j)

                if (dcom > best_dcom) || (dcom == best_dcom && obj < best_obj)
                    best_dcom = dcom
                    best_obj = obj
                    best_i = i
                    best_j = j
                end
            end
        end

        if best_i == 0
            break
        end

        # apply swap
        tmp = cur.route[best_i]
        cur.route[best_i] = cur.route[best_j]
        cur.route[best_j] = tmp
        cur.objective = best_obj
        cur_common += best_dcom

        if cur.objective < best.objective
            best.route .= cur.route
            best.objective = cur.objective
        end
    end

    if do_ls
        local_search!(m, best)
    end
    return best
end

############################
# Elite set (simple)
############################
function clone_solution(sol::TSPSolution)
    c = TSPSolution(Int32(length(sol.route)))
    c.route .= sol.route
    c.objective = sol.objective
    return c
end

function update_elite!(elite::Vector{TSPSolution}, sol::TSPSolution; elite_size::Int=10)
    push!(elite, clone_solution(sol))
    sort!(elite, by = s -> s.objective)
    if length(elite) > elite_size
        resize!(elite, elite_size)
    end
end

function pick_guiding(elite::Vector{TSPSolution}, best_sol::TSPSolution)
    # guiding: sometimes best, sometimes random elite for diversity
    if isempty(elite)
        return best_sol
    end
    if rand() < 0.7
        return best_sol
    else
        return elite[rand(1:length(elite))]
    end
end

############################
# Reactive GRASP + PR
############################
"""
pr_mode:
  :none     -> no path relinking
  :twoopt   -> PR using same neighborhood as LS (2-opt)
  :swap     -> PR using swap neighborhood
"""
function reactive_grasp(m::TSPSolver;
                        max_iters::Int=50,
                        alphas::Vector{Float32}=collect(Float32, 0.1:0.1:0.9),
                        window::Int=20,
                        elite_size::Int=10,
                        pr_mode::Symbol=:twoopt)

    stats = AlphaStats(alphas, window)
    best_sol = TSPSolution(m.dim)
    best_sol.objective = Inf32
    elite = TSPSolution[]  # store best locals

    for iter in 1:max_iters
        if elapsed_time(m) >= m.timelimit
            break
        end

        alpha, alpha_idx = select_alpha!(stats)

        # 1) construct
        sol = grasp_construction(m; alpha=alpha)

        # 2) local search
        local_search!(m, sol)
        final_cost = sol.objective

        # 3) Path Relinking (intensification)
        if pr_mode != :none && isfinite(best_sol.objective)
            guiding = pick_guiding(elite, best_sol)

            pr_best = nothing
            if pr_mode == :twoopt
                pr_best = path_relinking_two_opt(m, sol, guiding; max_steps=200, do_ls=true)
            elseif pr_mode == :swap
                pr_best = path_relinking_swap(m, sol, guiding; max_steps=200, do_ls=true)
            end

            if pr_best !== nothing && pr_best.objective < sol.objective
                sol.route .= pr_best.route
                sol.objective = pr_best.objective
                final_cost = sol.objective
            end
        end

        # 4) update stats
        update_stats!(stats, alpha_idx, final_cost)

        # 5) keep global best + elite
        if final_cost < best_sol.objective
            best_sol.route .= sol.route
            best_sol.objective = final_cost
        end
        update_elite!(elite, sol; elite_size=elite_size)

        @printf("Iter %3d | α=%.2f | cost=%.3f | best=%.3f\n", iter, alpha, final_cost, best_sol.objective)
    end

    println("Alpha selection counts: ", stats.counts)
    return best_sol
end

############################
# Main: run comparison
############################
function main()
    # Read instance
    name, coord, dim = readInstance("tsp_fun.tsp")
    dist = getDistanceMatrix(coord, dim)

    println("Instance: $name | n=$(dim)")

    # Baseline: NN + LS
    m0 = TSPSolver(coord, dim, dist, 10)
    sol_nn = nearest_neighbor_heuristic(m0)
    println("NN objective: ", sol_nn.objective)
    sol_ls = local_search(m0)
    println("LS (2-opt) objective: ", sol_ls.objective)

    # Reactive GRASP + PR(2-opt)
    println("\n=== Reactive GRASP + Path Relinking (2-opt) ===")
    m1 = TSPSolver(coord, dim, dist, 10)
    alphas = collect(Float32, 0.1:0.1:0.9)
    best_pr_2opt = reactive_grasp(m1; max_iters=40, alphas=alphas, pr_mode=:twoopt, elite_size=10)
    println("Best (PR 2-opt): ", best_pr_2opt.objective)

    # Reactive GRASP + PR(SWAP)
    println("\n=== Reactive GRASP + Path Relinking (swap) ===")
    m2 = TSPSolver(coord, dim, dist, 10)
    best_pr_swap = reactive_grasp(m2; max_iters=40, alphas=alphas, pr_mode=:swap, elite_size=10)
    println("Best (PR swap): ", best_pr_swap.objective)

    println("\nDONE.")
end

main()