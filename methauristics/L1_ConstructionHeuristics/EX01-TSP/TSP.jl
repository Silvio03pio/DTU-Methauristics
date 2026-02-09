include("TSP_Solver.jl") 

function main() #  per ottenere stessi nomi come nel solver
    name, coord, dim = readInstance(joinpath(@__DIR__, "toy_slides.tsp"))
    dist = getDistanceMatrix(coord, dim)

    n = size(dist, 1) #sapere le dimensioni di un array/matrice

    tot = 0.0                 # costo totale
    visited = falses(n)        # vettore Bool: false = non visitata
    tour = Int[]               # qui salviamo l’ordine delle città visitate

    start = 1  # nodo di partenza
    current = start

    visited[current] = true
    push!(tour, current) #aggiungi un elemento alla fine di un array 

    for step in 2:n

        best_city = 0
        best_d = Inf # inizializzo a infinito 

        for cand in 1:n
            if !visited[cand] && dist[current, cand] < best_d
                best_d = dist[current, cand]
                best_city = cand
            end
        end

        tot += best_d
        current = best_city
        visited[current] = true
        push!(tour, current)

    end

    tot += dist[current, start]
    push!(tour, start)

    println("Tour: ", tour)
    println("Total cost: ", tot)


end

main()


