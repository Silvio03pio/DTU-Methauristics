
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



function main()
    name, coord, dim = readInstance("tsp_toy.tsp")
    dist = getDistanceMatrix(coord,dim)

    sol = nearest_neighbor_heuristic(dist,dim)
    println("Route: ",sol.route)
    println("Objective: ",sol.objective)

end

main()