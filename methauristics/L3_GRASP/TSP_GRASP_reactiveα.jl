
using Plots
function readInstance(filename)
    #open file for reading
    file = open(filename)
    #read the name of the instance
    name = split(readline(file))[2]
    #The next two lines are not interesting for us. Skip them
    readline(file);readline(file)
    #Read the size of the instance (the number of cities)
    dim = parse(Int32,split(readline(file))[2])
    #The next two lines are not interesting for us. Skip them
    readline(file);readline(file)
    #Create a Matrix (dim ⋅ 2) to hold the coordinates
    coord = zeros(Float32,dim,2)
    #Read coordinates
    for i in 1:dim
        data = parse.(Float32,split(readline(file)))
        coord[i,:]=data[2:3]
    end
    #Close the file
    close(file)
    #return the data we need
    return name,coord,dim
end

#** Creates a distance matrix ****************
# Arguments:
#     coord::Array{Float32,2}  An array of coordinate pairs (x,y)
#     dim::Int32  The dimention of the coord array
# Returns:
#     dist::Array{Float32,2} Distance matrix based on the straight line distance
#***************************
function getDistanceMatrix(coord::Array{Float32,2},dim::Int32)
    dist = zeros(Float32,dim,dim)
    for i in 1:dim
       for j in 1:dim
            if i!=j
                dist[i,j]=round(sqrt((coord[i,1]-coord[j,1])^2+(coord[i,2]-coord[j,2])^2),digits=2)
            end
        end
    end
    return dist
end


# struct representing the solver
struct TSPSolver
    coord::Array{Float32,2}
    dim::Int32
    dist::Array{Float32,2}
    startTime::UInt64
    timelimit::UInt32
    TSPSolver(coord, dim, dist, timelimit) = new(coord,dim,dist,time_ns(),timelimit)
end
# struct representing a solution
mutable struct TSPSolution
    # we represent a solution as a list of cities.
    route::Array{Int32,1}
    # placeholder for the objective value
    objective::Float32

    # solution default constructor
    TSPSolution(dim) = new(zeros(Int32,dim),0)
end



# given a distance matrix, and a list of visited cities
# returns the city closed to "city"
function get_nearest_neighbor(dist,visited,city)
    # variable to keep track of the smallest distance 
    distance = Inf
    # variable to keep track of the next city
    next =0
    # run through all the cities
    for i in 1:length(dist[city,:])
        # if city is not visited and the distance is smaller
        # than what we have seen
        if !visited[i] && distance>dist[city,i]
            # update the best distance and the next city
            distance = dist[city,i]
            next = i
        end
    end

    #return the next city
    return next
end



# Nearest neighbor construction heuristic
function nearest_neighbor_heuristic(m::TSPSolver)
    # We initialize the empty solution
    sol = TSPSolution(m.dim)

    # keep track of whcih cities have been visited
    visited = zeros(Bool,m.dim)

    # we start by assigning the first city and flag it as visited
    sol.route[1] = 1
    visited[1] = true

    # we make a loop where i is the index of the last city inserted
    # since we will be assigning the next (i+1) city, we need to stop at dim-1
    for i in 1:m.dim-1
        # get the ID of the neighor closest to city in position i (the last one we added)
        j = get_nearest_neighbor(m.dist, visited, sol.route[i])
        # assign the neighbor to the solution and flag it as visited
        sol.route[i+1] = j
        visited[j]=true

        # update the objective value
        sol.objective += m.dist[sol.route[i],sol.route[i+1]]
    end
    # we need to remember to add cost between the last ans the first city
    # since we are building a cycle
    sol.objective+=m.dist[sol.route[1],sol.route[end]]
    
    return sol
end



# Given the solver, which stores the start time, this function
# returns the elapsed time in seconds.
function elapsed_time(m::TSPSolver)
   return round((time_ns()-m.startTime)/1e9,digits=3)
end

# To save memory and time I have made a delta-evaluation of a 2-opt
# operation. Given two cities, i and j, the function calculates the
# objective value of performin the 2-opt without actually changing 
# or copying the solution.
function eval_two_opt_step(m::TSPSolver,sol::TSPSolution,i,j)
    #assumption: i<j
    #case 1: i=1,j=end, simple reverse no cost chage
    if i==1 && j==length(sol.route)
        return sol.objective
    elseif i==1 #case 1: i=1, j in the middle
        return sol.objective -
               m.dist[sol.route[i],sol.route[end]] -
               m.dist[sol.route[j],sol.route[j+1]] +
               m.dist[sol.route[j],sol.route[end]] +
               m.dist[sol.route[i],sol.route[j+1]]
    elseif j==m.dim #case 3, i in the middle, j=end
        return sol.objective -
            m.dist[sol.route[i-1],sol.route[i]] -
            m.dist[sol.route[j],sol.route[1]] +
            m.dist[sol.route[i-1],sol.route[j]] +
            m.dist[sol.route[i],sol.route[1]]
    else #case i and j in the middle
        return sol.objective -
            m.dist[sol.route[i-1],sol.route[i]] -
            m.dist[sol.route[j],sol.route[j+1]] +
            m.dist[sol.route[i-1],sol.route[j]] +
            m.dist[sol.route[i],sol.route[j+1]]
       end
end

# this is the neighborhood operator, which includes the
# step function best improved.
function best_two_opt(m::TSPSolver, sol::TSPSolution)
    # I initialize the variables that will be holding
    # the i and j city over which a 2-opt is to be perfomed.
    iBest = 0
    jBest = 0
    bestObj = Inf
    # I go though all the possible 2-opt (i,j) intervals
    for i in 2:m.dim-1
        for j in i+1:m.dim
            # I calculated the objtive function of the 2-opt via delta-evaluation
            obj = eval_two_opt_step(m,sol,i,j) #delta evaluation
            # If the objective is better than the current best, make this the best
            if obj < bestObj
                iBest = i
                jBest = j
                bestObj = obj
            end
        end
    end
    # return the potential objective value and the 2-opt interval cities
    return bestObj,iBest, jBest
end

# This is the local search function. Localt search usually includes, as one
# of the first steps, the generation of an initial solution. Since we will be 
# using this fuction for GRASP, I have decided to pass the initial solution as
# a parameter.
function local_search!(m::TSPSolver, sol::TSPSolution)
    improved=true
    # Terminate if no improvement found or timeout
    while(improved && elapsed_time(m)<m.timelimit)
        # Neighborhood + Step function
        # this function is returning the 2-opt to perform rather than
        # an entire solution. I do that to save memory and the extra
        # runtime cost of having to copy the soltution
        val, i, j = best_two_opt(m,sol)
        # Acceptance criteria
        if val < sol.objective
            # implement solution change
            # since the neighborhood function is returning the 2-opt
            # and not a solution, accepting the solution means
            # applying the 2-op to the current solution.
            # Which we do by reversing the cities between i and j.
            reverse!(sol.route,i,j)
            sol.objective = val
            improved = true
        else
            improved = false
        end
    end
    return sol
end


# This function is just a simplified version of the calculated
# for local serach where the initial solution is already included
function local_search(m::TSPSolver)
    sol = nearest_neighbor_heuristic(m)
    return local_search!(m,sol)
end


"""
Greedy Randomized Construction (GRASP) for TSP with alpha-based RCL
- alpha ∈ [0,1]: greediness parameter (0=pure greedy, 1=fully random)
"""
function grasp_construction(m::TSPSolver; alpha::Float32=0.5f0)
    sol = TSPSolution(m.dim)
    visited = falses(m.dim)

    # start always at 1
    current = 1
    sol.route[1] = current
    visited[current] = true
    sol.objective = 0.0f0

    for pos in 2:m.dim
        # build candidate list: unvisited cities with incremental costs
        c = Float32[]  # c(e)
        cand_ids = Int[]
        for j in 1:m.dim
            if !visited[j]
                push!(c, m.dist[current, j])
                push!(cand_ids, j)
            end
        end

        if isempty(c)
            error("No candidates left")
        end

        # find min and max incremental costs
        c_min = minimum(c)
        c_max = maximum(c)

        # RCL: candidates where c(e) ∈ [c_min, c_min + α*(c_max - c_min)]
        threshold = c_min + alpha * (c_max - c_min)
        rcl_mask = c .<= threshold  #rcl include c if the value of c is minor of treshold 
        #creates a boolean array where each element is true if the corresponding cost in c is ≤ the threshold, false otherwise.
        rcl = cand_ids[rcl_mask]

        if isempty(rcl)
            # fallback: use at least the min
            min_idx = argmin(c)
            rcl = [cand_ids[min_idx]]
        end

        # pick one uniformly at random from RCL
        next_city = rand(rcl)

        # add to route
        sol.route[pos] = next_city
        visited[next_city] = true
        sol.objective += m.dist[current, next_city]
        current = next_city
    end

    # close the tour
    sol.objective += m.dist[sol.route[end], sol.route[1]]

    return sol
end

"""
GRASP metaheuristic
- max_iters: number of GRASP iterations
- alpha: RCL parameter in construction
"""
function grasp(m::TSPSolver; max_iters::Int=50, alpha::Float32=0.5f0)
    best_sol = TSPSolution(m.dim)
    best_sol.objective = Inf32

    iter = 0
    while iter < max_iters && elapsed_time(m) < m.timelimit
        iter += 1

        # 1) Greedy randomized construction
        sol = grasp_construction(m; alpha=alpha)

        # 2) Local search (2-opt)
        local_search!(m, sol)

        # 3) Keep best
        if sol.objective < best_sol.objective
            best_sol.route .= sol.route
            best_sol.objective = sol.objective
        end
    end

    return best_sol
end

function experiment_alpha_curve()
    # read instance
    name, coord, dim = readInstance("tsp_toy50.tsp")
    dist = getDistanceMatrix(coord, dim)

    alphas = 0.0:0.1:1.0
    best_vals = zeros(Float64, length(alphas))

    for (idx, alpha) in enumerate(alphas)
        m = TSPSolver(coord, dim, dist, 10)  # 10s time limit
        sol_alpha = grasp(m; max_iters=30, alpha=Float32(alpha))
        best_vals[idx] = sol_alpha.objective
        println("alpha=$alpha -> objective = $(best_vals[idx])")
    end

    # plot
    p = plot(alphas, best_vals,
             xlabel = "α (RCL greediness)",
             ylabel = "Best objective",
             title = "GRASP performance vs α",
             marker = :o,
             legend = false)

    display(p)
    savefig(p, "grasp_alpha_curve.png")
end


# Track stats for each alpha
mutable struct AlphaStats
    alphas::Vector{Float32}
    all_costs::Vector{Vector{Float32}}
    counts::Vector{Int}
    window::Int
end

function AlphaStats(alphas::Vector{Float32}, window::Int)
    AlphaStats(alphas,
               [Float32[] for _ in 1:length(alphas)],
               zeros(Int, length(alphas)),
               window)
end

function select_alpha!(stats::AlphaStats)
    m = length(stats.alphas)
    
    # s* = best solution EVER
    s_star = Inf32
    has_data = false
    for costs in stats.all_costs
        if !isempty(costs)
            s_star = min(s_star, minimum(costs))
            has_data = true
        end
    end
    
    if !has_data
        println("No data: uniform probs [0.111...]")
        p = fill(1.0f0/m, m)
    else
        q = Float32[]
        for i in 1:m
            costs_i = stats.all_costs[i]
            if isempty(costs_i)
                # NO SOLUTIONS with this αᵢ → LOW quality (small qᵢ)
                push!(q, 1.0f0)  # small constant
            else
                # qᵢ = s*/Aᵢ
                A_i = sum(costs_i) / length(costs_i)
                push!(q, s_star / A_i)
            end
        end
        p = q / sum(q)
    end
    
    println("s*=$(round(s_star,digits=1)) | Probs=$(round.(p,digits=3))")
    
    idx = argmax(cumsum(p) .>= rand(Float32))
    stats.counts[idx] += 1
    return stats.alphas[idx], idx
end



function update_stats!(stats::AlphaStats, alpha_idx::Int, cost::Float32)
    push!(stats.all_costs[alpha_idx], cost)

    # Mantieni solo gli ultimi "window" valori
    if length(stats.all_costs[alpha_idx]) > stats.window
        popfirst!(stats.all_costs[alpha_idx])
    end
end


function reactive_grasp(m::TSPSolver; max_iters::Int=50, alphas::Vector{Float32}=0.1f0:0.1f0:0.9f0, window::Int=20)
    stats = AlphaStats(alphas, window)
    best_sol = TSPSolution(m.dim)
    best_sol.objective = Inf32

    iter = 0
    #while iter < max_iters && elapsed_time(m) < m.timelimit
    while iter < max_iters 
        iter += 1

        # 1) Select alpha probabilistically
        alpha, alpha_idx = select_alpha!(stats)

        # 2) Greedy randomized construction
        sol = grasp_construction(m; alpha=alpha)

        # 3) Local search
        local_search!(m, sol)
        final_cost = sol.objective

        # 4) Update stats with this cost
        update_stats!(stats, alpha_idx, final_cost)

        # 5) Keep global best
        if final_cost < best_sol.objective
            best_sol.route .= sol.route
            best_sol.objective = final_cost
        end

        println("Iter $iter: α=$(alpha), cost=$(final_cost), best=$(best_sol.objective)")
    end

    println("Alpha selection counts: ", stats.counts)
    return best_sol
end

# main function
function main()
    # read the instance
    name, coord, dim = readInstance("tsp_toy.tsp")
    # get the distance mastrix
    dist = getDistanceMatrix(coord,dim)
    # create the solver
    m = TSPSolver(coord,dim,dist,10)
    # call local search
    sol = nearest_neighbor_heuristic(m)
    println("RouteNeighrest Neighbor: ",sol.route)
    println("Objective: ",sol.objective)

    sol = local_search(m)
    # print the solution
    println("Route Local Search: ",sol.route)
    println("Objective: ",sol.objective)

    # 4) run GRASP
   # alpha = 0.5f0  # alpha parameter (instead of K)
    #iters = 50
    #best_sol = grasp(m; max_iters=iters, alpha=alpha)
   # println("Best route (GRASP, α=$alpha): ", best_sol.route)
    #println("Best objective: ", best_sol.objective)
    #experiment_alpha_curve()
    
    # Reactive GRASP
    alphas = collect(Float32, 0.1:0.1:0.9)
    best_reactive = reactive_grasp(m; max_iters=40, alphas=alphas)
    println("Reactive GRASP best: ", best_reactive.objective)
end

main()