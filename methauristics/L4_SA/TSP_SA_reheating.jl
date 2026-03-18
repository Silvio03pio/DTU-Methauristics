#***** Instance reader *********************************************************
using Random
using CSV
using DataFrames
using Plots


#=
#------------------------------------------------------------------------------
# RE-HEATING STRATEGY IN SIMULATED ANNEALING
#
# In Simulated Annealing, the temperature gradually decreases in order to move
# the algorithm from an exploratory phase (where even worse solutions can be
# accepted) to a more exploitative phase similar to a hill-climbing procedure.
# When the temperature becomes very small, the probability of accepting worse
# solutions approaches zero and the algorithm may become stuck in a local
# optimum. In that situation the search may continue for many iterations
# without producing any improvement.
#
# To mitigate this behaviour we introduce a re-heating mechanism. The algorithm
# keeps track of the number of iterations since the last improvement of the
# best solution found so far. If this number exceeds a predefined threshold
# (stagnation_limit), we assume that the search is currently trapped in a
# local optimum.
#
# Instead of restarting the algorithm from scratch, the temperature is reset
# to the initial temperature T0. This temporarily increases the probability
# of accepting worse solutions again, allowing the search process to escape
# from the local optimum and explore new regions of the solution space.
#
# Importantly, the best solution found so far is always preserved. Therefore,
# re-heating only affects the exploration dynamics and does not discard any
# valuable information obtained during the search.
#
# This mechanism improves the balance between exploration and exploitation
# and helps prevent the algorithm from wasting many iterations at extremely
# low temperatures where almost no moves are accepted.
#------------------------------------------------------------------------------
=#

# Arguments:
#     filename::String  The full path of the instance file
# Returns:
#     name::String  The name of the TSP instance
#     coord::Array{Float32,2} An array of coordinate pairs (x,y)
#     dim::Int32  The dimension of the coord array
#*******************************************************************************
function readInstance(filename)
    file = open(filename)

    name = split(readline(file))[2]

    readline(file)
    readline(file)

    dim = parse(Int32, split(readline(file))[2])

    readline(file)
    readline(file)

    coord = zeros(Float32, dim, 2)

    for i in 1:dim
        data = parse.(Float32, split(readline(file)))
        coord[i, :] = data[2:3]
    end

    close(file)

    return name, coord, dim
end

#***** Creates a distance matrix ***********************************************
# Arguments:
#     coord::Array{Float32,2}  An array of coordinate pairs (x,y)
#     dim::Int32  The dimension of the coord array
# Returns:
#     dist::Array{Float32,2} Distance matrix based on the straight line distance
#*******************************************************************************
function getDistanceMatrix(coord::Array{Float32,2}, dim::Int32)
    dist = zeros(Float32, dim, dim)
    for i in 1:dim
        for j in 1:dim
            if i != j
                dist[i, j] = round(
                    sqrt((coord[i, 1] - coord[j, 1])^2 + (coord[i, 2] - coord[j, 2])^2),
                    digits = 2
                )
            end
        end
    end
    return dist
end

#***** TSP solution structure **************************************************
mutable struct TSPSolution
    route::Vector{Int32}
    objective::Float32

    TSPSolution(dim) = new(zeros(Int32, dim), 0.0f0)
end

#***** Objective computation ***************************************************
# Computes the total tour length, including the closing edge from last to first
function compute_objective(route::Vector{Int32}, dist)
    n = length(route)
    obj = 0.0f0
    for k in 1:(n - 1)
        obj += dist[route[k], route[k + 1]]
    end
    obj += dist[route[end], route[1]]
    return obj
end

#***** Random initial solution *************************************************
# Builds a completely random tour (permutation of 1:dim) and computes its cost
function random_initial_solution(dist, dim::Int32, rng::AbstractRNG)
    sol = TSPSolution(dim)
    sol.route .= Int32.(randperm(rng, Int(dim)))
    sol.objective = compute_objective(sol.route, dist)
    return sol
end

#***** Nearest neighbor helper *************************************************
# Kept for completeness, not used by SA
function get_nearest_neighbor(dist, visited, city)
    distance = Inf
    next = 0
    for i in 1:length(dist[city, :])
        if !visited[i] && distance > dist[city, i]
            distance = dist[city, i]
            next = i
        end
    end
    return next
end

#***** Nearest neighbor heuristic **********************************************
# Kept for completeness, not used by SA
function nearest_neighbor_heuristic(dist, dim)
    sol = TSPSolution(dim)
    visited = zeros(Bool, dim)

    sol.route[1] = 1
    visited[1] = true

    for i in 1:dim-1
        j = get_nearest_neighbor(dist, visited, sol.route[i])
        sol.route[i + 1] = j
        visited[j] = true
        sol.objective += dist[sol.route[i], sol.route[i + 1]]
    end

    sol.objective += dist[sol.route[1], sol.route[end]]
    return sol
end

#***** Best 2-opt move *********************************************************
# Kept from local search version, not used as main SA neighborhood
function best_2opt_move(sol::TSPSolution, dist)
    route = sol.route
    n = length(route)
    best_delta = 0.0f0
    best_i, best_j = 0, 0

    for i in 1:(n - 3)
        for j in (i + 2):(n - 1)
            a = route[i]
            b = route[i + 1]
            c = route[j]
            d = route[j + 1]

            delta = (dist[a, c] + dist[b, d]) - (dist[a, b] + dist[c, d])

            if delta < best_delta
                best_delta = delta
                best_i, best_j = i, j
            end
        end
    end

    return (best_i, best_j), best_delta
end

#***** First improving 2-opt move **********************************************
# Kept from local search version, not used as main SA neighborhood
function first_improving_2opt_move(sol::TSPSolution, dist)
    route = sol.route
    n = length(route)

    for i in 1:(n - 3)
        for j in (i + 2):(n - 1)
            a = route[i]
            b = route[i + 1]
            c = route[j]
            d = route[j + 1]

            delta = (dist[a, c] + dist[b, d]) - (dist[a, b] + dist[c, d])

            if delta < 0
                return (i, j), delta
            end
        end
    end

    return (0, 0), 0.0f0
end

#***** Full 2-opt move *********************************************************
# Kept as exact version for debugging/checking
function switch_2opt(sol::TSPSolution, i::Int, j::Int, dist)
    new_sol = deepcopy(sol)
    reverse!(new_sol.route, i + 1, j)
    new_sol.objective = compute_objective(new_sol.route, dist)
    return new_sol
end

#***** 2-opt delta evaluation **************************************************
# Computes the objective variation of a 2-opt move in O(1)
function eval_2opt_delta(sol::TSPSolution, i::Int, j::Int, dist)
    route = sol.route

    a = route[i]
    b = route[i + 1]
    c = route[j]
    d = route[j + 1]

    delta = (dist[a, c] + dist[b, d]) - (dist[a, b] + dist[c, d])
    return Float32(delta)
end

#***** Apply 2-opt move in place ***********************************************
# Applies the move to the route and updates the objective with the delta
function apply_2opt_move!(sol::TSPSolution, i::Int, j::Int, delta::Float32)
    reverse!(sol.route, i + 1, j)
    sol.objective += delta
    return sol
end

#***** Simulated Annealing *****************************************************
# Random initial solution + random 2-opt neighborhood + probabilistic accept
# Uses delta evaluation instead of full objective recomputation
function simulated_annealing(dist, dim::Int32;
                             T0::Float32 = 100.0f0,
                             alpha::Float32 = 0.995f0,
                             max_iter::Int = 10_000,
                             time_limit::Float64 = 10.0,
                             stagnation_limit::Int = 5000,
                             seed::Int = 42)

    n = Int(dim)
    rng = MersenneTwister(seed)

    # --- RANDOM INITIAL SOLUTION ---
    current_sol = random_initial_solution(dist, dim, rng)
    initial_sol = deepcopy(current_sol)
    best_sol = deepcopy(current_sol)

    # --- TEMPERATURE INITIALIZATION ---
    T = T0

    # --- ITERATION COUNTERS ---
    iter = 0
    accepted_moves = 0
    reheat_count = 0
    iters_since_improvement = 0

    # --- DATA LOGGING ---
    history = DataFrame(
        iteration = Int[],
        temperature = Float64[],
        current_obj = Float64[],
        candidate_obj = Float64[],
        delta = Float64[],
        accepted = Int[],
        improving = Int[],
        best_obj = Float64[],
        reheated = Int[]
    )

    # --- EXECUTION TIME MEASUREMENT ---
    start_time = time()

    while (time() - start_time) < time_limit 
        iter += 1

        if n <= 3
            break
        end

        # --- RANDOM 2-OPT MOVE ---
        i = rand(rng, 1:(n - 3))
        j = rand(rng, (i + 2):(n - 1))

        cur_obj = current_sol.objective

        # --- DELTA EVALUATION ---
        delta = eval_2opt_delta(current_sol, i, j, dist)
        new_obj = cur_obj + delta

        accepted = 0
        improving = 0
        reheated = 0

        # --- PROBABILISTIC ACCEPTANCE ---
        if new_obj < cur_obj
            apply_2opt_move!(current_sol, i, j, delta)
            accepted = 1
            improving = 1
            accepted_moves += 1
        else
            if rand(rng) < exp(-delta / T)
                apply_2opt_move!(current_sol, i, j, delta)
                accepted = 1
                accepted_moves += 1
            end
        end

        # --- UPDATE GLOBAL BEST ---
        if current_sol.objective < best_sol.objective
            best_sol = deepcopy(current_sol)
            iters_since_improvement = 0
        else
            iters_since_improvement += 1
        end

        # --- RE-HEATING AFTER STAGNATION ---
        if iters_since_improvement >= stagnation_limit
            T = T0
            iters_since_improvement = 0
            reheat_count += 1
            reheated = 1
        end

        # --- LOGGING ---
        push!(history, (
            iter,
            Float64(T),
            Float64(cur_obj),
            Float64(new_obj),
            Float64(delta),
            accepted,
            improving,
            Float64(best_sol.objective),
            reheated
        ))

        # --- TEMPERATURE COOLING ---
        T *= alpha
    end

    elapsed_time = time() - start_time

    return initial_sol, best_sol, iter, accepted_moves, reheat_count, elapsed_time, history
end

#=
#***** Plot utilities **********************************************************
function plot_sa_history(history::DataFrame)
    xvals = Int.(history.iteration)
    maxit = maximum(xvals)

    if maxit <= 100
        step = 10
    elseif maxit <= 1000
        step = 100
    elseif maxit <= 10000
        step = 500
    elseif maxit <= 100000
        step = 5000
    else
        step = 50000
    end

    ticks = collect(1:step:maxit)
    if ticks[end] != maxit
        push!(ticks, maxit)
    end

    # 1. Best objective over iterations
    p1 = plot(
        xvals,
        history.best_obj,
        xlabel = "Iteration",
        ylabel = "Best objective",
        title = "Best Objective over Iterations",
        label = "Best objective",
        lw = 2,
        xlims = (1, maxit),
        xticks = ticks,
        xformatter = x -> string(Int(round(x)))
    )
    display(p1)
    savefig(p1, "sa_best_objective.png")

    # 2. Current and candidate objective
    p2 = plot(
        xvals,
        history.current_obj,
        xlabel = "Iteration",
        ylabel = "Objective value",
        title = "Current vs Candidate Objective",
        label = "Current objective",
        lw = 2,
        xlims = (1, maxit),
        xticks = ticks,
        xformatter = x -> string(Int(round(x)))
    )
    plot!(
        p2,
        xvals,
        history.candidate_obj,
        label = "Candidate objective",
        lw = 2
    )
    display(p2)
    savefig(p2, "sa_current_candidate.png")

    # 3. Accepted / rejected
    p3 = scatter(
        xvals,
        history.accepted,
        xlabel = "Iteration",
        ylabel = "Accepted",
        title = "Accepted Moves over Iterations",
        label = "Accepted = 1 / Rejected = 0",
        markersize = 2,
        xlims = (1, maxit),
        xticks = ticks,
        xformatter = x -> string(Int(round(x)))
    )
    display(p3)
    savefig(p3, "sa_acceptance.png")
end

=#

function print_move(route, i, j)
    a = route[i]
    b = route[i + 1]
    c = route[j]
    d = route[j + 1]
    println("Cut indices: i=$i, j=$j")
    println("Remove edges: ($a,$b) and ($c,$d)")
    println("Add edges:    ($a,$c) and ($b,$d)")
    println("Segment reversed route[", i + 1, ":", j, "] = ", route[(i + 1):j])
end

#***** Main ********************************************************************
function main()
    name, coord, dim = readInstance("berlin52_7542.tsp")
    dist = getDistanceMatrix(coord, dim)

    initial_sol, best_sol, iters, accepted_moves, reheat_count, elapsed_time, history =
        simulated_annealing(dist, dim; time_limit = 10.0, stagnation_limit = 5000)

    println("=== SIMULATED ANNEALING TSP ===")
    println("Instance name: ", name)
    println()

    println("Initial random route: ", initial_sol.route)
    println("Initial random objective: ", initial_sol.objective)
    println()

    println("Best route found by SA: ", best_sol.route)
    println("Best objective found by SA: ", best_sol.objective)
    println()

    println("Iterations performed: ", iters)
    println("Accepted moves: ", accepted_moves)
    println("Re-heating operations: ", reheat_count)
    println("Execution time (seconds): ", elapsed_time)

    #CSV.write("sa_history.csv", history)
    #println("History saved to sa_history.csv")

    println("Min iteration in history: ", minimum(history.iteration))
    println("Max iteration in history: ", maximum(history.iteration))

    #plot_sa_history(history)
end

main()