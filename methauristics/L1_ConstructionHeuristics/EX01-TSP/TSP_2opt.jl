
#***** Instance reader *********************************************************
# Arguments:
#     filename::String  The full path of the instance file
# Returns:
#     name::String  The name of the TSP instance
#     coord::Array{Float32,2} An array of coordinate pairs (x,y)
#     dim::Int32  The dimention of the coord array
#*******************************************************************************
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

#***** Creates a distance matrix ***********************************************
# Arguments:
#     coord::Array{Float32,2}  An array of coordinate pairs (x,y)
#     dim::Int32  The dimention of the coord array
# Returns:
#     dist::Array{Float32,2} Distance matrix based on the straight line distance
#*******************************************************************************
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
function nearest_neighbor_heuristic(dist,dim)
    # We initialize the empty solution
    sol = TSPSolution(dim)

    # keep track of whcih cities have been visited
    visited = zeros(Bool,dim)

    # we start by assigning the first city and flag it as visited
    sol.route[1] = 1
    visited[1] = true

    # we make a loop where i is the index of the last city inserted
    # since we will be assigning the next (i+1) city, we need to stop at dim-1
    for i in 1:dim-1
        # get the ID of the neighor closest to city in position i (the last one we added)
        j = get_nearest_neighbor(dist, visited, sol.route[i])
        # assign the neighbor to the solution and flag it as visited
        sol.route[i+1] = j
        visited[j]=true

        # update the objective value
        sol.objective += dist[sol.route[i],sol.route[i+1]]
    end
    # we need to remember to add cost between the last ans the first city
    # since we are building a cycle
    sol.objective+=dist[sol.route[1],sol.route[end]]
    
    return sol
end

"""
Trova il miglior 2-opt move (i,j) e il relativo Δ (negativo = miglioramento).

Assume che `sol` sia un tour "chiuso" con l'ultimo nodo uguale al primo,
cioè length(sol) == n+1 e sol[end] == sol[1].
"""

function best_2opt_move(sol::TSPSolution, dist)
    route = sol.route
    n = length(route)              # ✅ route non chiuso
    best_delta = 0.0f0
    best_i, best_j = 0, 0

    for i in 1:(n-3)
        for j in (i+2):(n-1)
            a = route[i]
            b = route[i+1]
            c = route[j]
            d = route[j+1]

            delta = (dist[a,c] + dist[b,d]) - (dist[a,b] + dist[c,d])

            if delta < best_delta
                best_delta = delta
                best_i, best_j = i, j
            end
        end
    end

    return (best_i, best_j), best_delta
end





function first_improving_2opt_move(sol::TSPSolution, dist)
    route = sol.route
    n = length(route)

    for i in 1:(n-3)
        for j in (i+2):(n-1)
            a = route[i]; b = route[i+1]
            c = route[j]; d = route[j+1]

            delta = (dist[a,c] + dist[b,d]) - (dist[a,b] + dist[c,d])

            if delta < 0
                return (i, j), delta
            end
        end
    end

    return (0, 0), 0.0f0
end



function compute_objective(route::Vector{Int32}, dist)
    n = length(route)
    obj = 0.0f0
    for k in 1:(n-1)
        obj += dist[route[k], route[k+1]]
    end
    obj += dist[route[end], route[1]]   # chiusura ciclo
    return obj
end


function switch_2opt(sol::TSPSolution, i::Int, j::Int, dist)
    new_sol = deepcopy(sol)
    reverse!(new_sol.route, i+1, j)                 # 2-opt vero
    new_sol.objective = compute_objective(new_sol.route, dist)  # ricalcolo completo
    return new_sol
end


#=
function debug_moves(sol::TSPSolution, dist)
    route = sol.route
    n = length(route)
    best = (0,0); bestΔ = 0.0f0
    first = (0,0); firstΔ = 0.0f0; found=false

    for i in 1:(n-3), j in (i+2):(n-1)
        a,b,c,d = route[i],route[i+1],route[j],route[j+1]
        Δ = (dist[a,c]+dist[b,d]) - (dist[a,b]+dist[c,d])
        if Δ < 0 && !found
            first = (i,j); firstΔ = Δ; found=true
        end
        if Δ < bestΔ
            bestΔ = Δ; best = (i,j)
        end
    end

    println("first = ", first, " Δ=", firstΔ)
    println("best  = ", best,  " Δ=", bestΔ)
end

=#



function print_move(route, i, j)
    a = route[i]; b = route[i+1]
    c = route[j]; d = route[j+1]
    println("Cut indices: i=$i, j=$j")
    println("Remove edges: ($a,$b) and ($c,$d)")
    println("Add edges:    ($a,$c) and ($b,$d)")
    println("Segment reversed route[", i+1, ":", j, "] = ", route[(i+1):j])
end

function main()
    name, coord, dim = readInstance("tsp_toy.tsp")
    dist = getDistanceMatrix(coord, dim)

    sol = nearest_neighbor_heuristic(dist, dim)
    println("=== BASE SOLUTION (Nearest Neighbor) ===")
    println("Route: ", sol.route)
    println("Objective: ", sol.objective)

    # ----- FIRST -----
    (m_first, d_first) = first_improving_2opt_move(sol, dist)
    println("\n=== FIRST-IMPROVEMENT 2-OPT ===")
    println("Move (i,j): ", m_first, "  Δ=", d_first)
    if d_first < 0
        i, j = m_first
        print_move(sol.route, i, j)
        sol_first = switch_2opt(sol, i, j, dist)
        println("New route: ", sol_first.route)
        println("New objective: ", sol_first.objective)
    else
        println("No improving move found.")
    end

    # ----- BEST -----
    (m_best, d_best) = best_2opt_move(sol, dist)
    println("\n=== BEST-IMPROVEMENT 2-OPT ===")
    println("Move (i,j): ", m_best, "  Δ=", d_best)
    if d_best < 0
        i, j = m_best
        print_move(sol.route, i, j)
        sol_best = switch_2opt(sol, i, j, dist)
        println("New route: ", sol_best.route)
        println("New objective: ", sol_best.objective)
    else
        println("No improving move found.")
    end
end

main()

