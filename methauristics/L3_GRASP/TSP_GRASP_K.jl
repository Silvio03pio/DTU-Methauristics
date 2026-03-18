using Random

# struct representing a solution (dal file originale)
mutable struct TSPSolution
    # we represent a solution as a list of cities.
    route::Array{Int32,1}
    # placeholder for the objective value
    objective::Float32

    # solution default constructor
    TSPSolution(dim) = new(zeros(Int32,dim),0)
end

#***** Instance reader *********************************************************
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

# Struct per tenere traccia delle statistiche
mutable struct SearchStats
    iterations::Int
    improvements::Int
    time::Float64
    objective_values::Vector{Float32}
    
    SearchStats() = new(0, 0, 0.0, Float32[])
end

# Calcolo obiettivo con delta evaluation 
function compute_delta(route::Vector{Int32}, dist, i::Int, j::Int)
    n = length(route)
    a = route[i]
    b = route[i+1]
    c = route[j]
    d = route[j+1]
    
    # Delta = (nuovi archi) - (vecchi archi)
    return (dist[a,c] + dist[b,d]) - (dist[a,b] + dist[c,d])
end

# 2-opt move con aggiornamento incrementale dell'obiettivo 
function apply_2opt_move!(sol::TSPSolution, i::Int, j::Int, dist)
    # PRIMA calcola il delta usando il route originale
    delta = compute_delta(sol.route, dist, i, j)
    
    # POI inverti il segmento
    reverse!(sol.route, i+1, j)
    
    # AGGIORNA l'obiettivo
    sol.objective += delta
    
    return sol
end

# Best improvement 2-opt 
function best_improvement_2opt!(sol::TSPSolution, dist, stats::SearchStats)
    improved = true
    n = length(sol.route)
    
    while improved
        improved = false
        best_delta = 0.0f0
        best_i, best_j = 0, 0
        
        # Cerca il miglior miglioramento
        for i in 1:(n-3)
            for j in (i+2):(n-1)
                stats.iterations += 1
                delta = compute_delta(sol.route, dist, i, j)
                
                if delta < best_delta
                    best_delta = delta
                    best_i, best_j = i, j
                    improved = true
                end
            end
        end
        
        # Applica il miglior miglioramento trovato (se migliora)
        if improved && best_delta < 0
            apply_2opt_move!(sol, best_i, best_j, dist)
            stats.improvements += 1
            push!(stats.objective_values, sol.objective)
        end
    end
    
    return sol
end

# First improvement 2-opt 
function first_improvement_2opt!(sol::TSPSolution, dist, stats::SearchStats)
    improved = true
    n = length(sol.route)
    
    while improved
        improved = false
        
        for i in 1:(n-3)
            for j in (i+2):(n-1)
                stats.iterations += 1
                delta = compute_delta(sol.route, dist, i, j)
                
                if delta < 0  # Primo miglioramento trovato
                    apply_2opt_move!(sol, i, j, dist)
                    improved = true
                    stats.improvements += 1
                    push!(stats.objective_values, sol.objective)
                    break
                end
            end
            if improved break end
        end
    end
    
    return sol
end

# Genera soluzione casuale
function random_solution(dim::Int32)
    sol = TSPSolution(dim)
    sol.route = shuffle!([i for i in 1:dim])
    return sol
end

# Calcola obiettivo per soluzione random
function compute_route_objective!(sol::TSPSolution, dist)
    n = length(sol.route)
    obj = 0.0f0
    for i in 1:(n-1)
        obj += dist[sol.route[i], sol.route[i+1]]
    end
    obj += dist[sol.route[end], sol.route[1]]
    sol.objective = obj
    return sol
end

# Funzione di perturbazione per ILS (double bridge move)
function double_bridge_perturbation!(sol::TSPSolution, dist)
    n = length(sol.route)
    
    # Scegli 4 punti di taglio casuali
    a = rand(2:div(n,4))
    b = rand(a+1:div(n,3))
    c = rand(b+1:div(n,2))
    d = rand(c+1:n-1)
    
    # Double bridge move: A B C D → A C B D
    route = sol.route
    segment1 = route[1:a]
    segment2 = route[a+1:b]
    segment3 = route[b+1:c]
    segment4 = route[c+1:d]
    segment5 = route[d+1:end]
    
    new_route = vcat(segment1, segment4, segment3, segment2, segment5)
    sol.route = new_route
    
    # Ricalcola obiettivo
    compute_route_objective!(sol, dist)
    
    return sol
end

# Iterated Local Search
function iterated_local_search(dist, dim, max_iterations::Int=100)
    stats = SearchStats()
    
    # Genera soluzione iniziale (NN)
    sol = nearest_neighbor_heuristic(dist, dim)
    best_sol = deepcopy(sol)
    
    println("Soluzione iniziale NN: ", round(sol.objective, digits=2))
    
    for iter in 1:max_iterations
        # Perturbazione
        perturbed = deepcopy(sol)
        double_bridge_perturbation!(perturbed, dist)
        
        # Ricerca locale (first improvement)
        first_improvement_2opt!(perturbed, dist, stats)
        
        # Criterio di accettazione (miglioramento o probabilità)
        if perturbed.objective < sol.objective
            sol = deepcopy(perturbed)
            if sol.objective < best_sol.objective
                best_sol = deepcopy(sol)
                println("Iterazione $iter: nuovo best = ", round(best_sol.objective, digits=2))
            end
        else
            # Accetta con probabilità basata sul deterioramento
            acceptance_prob = exp(-(perturbed.objective - sol.objective) / sol.objective)
            if rand() < acceptance_prob
                sol = deepcopy(perturbed)
            end
        end
    end
    
    return best_sol, stats
end

# ===============================
# GRASP CON K FISSO
# ===============================

# GRASP con K fisso - Fase di costruzione
function grasp_constructive_fixedK(dist, dim, K::Int)
    # K è il numero fisso di candidati nella RCL
    
    sol = TSPSolution(dim)
    visited = zeros(Bool, dim)
    
    # Scegli città iniziale casuale
    start_city = rand(1:dim)
    sol.route[1] = start_city
    visited[start_city] = true
    
    for position in 1:dim-1
        current_city = sol.route[position]
        
        # Trova tutte le città non visitate
        unvisited = [c for c in 1:dim if !visited[c]]
        
        # Calcola distanze per ogni città non visitata
        distances = [(city, dist[current_city, city]) for city in unvisited]
        
        # Ordina per distanza crescente
        sort!(distances, by = x -> x[2])
        
        # Prendi le prime K città (RCL a dimensione fissa)
        # Se ci sono meno di K città non visitate, prendi tutte
        rcl_size = min(K, length(unvisited))
        rcl = [distances[i][1] for i in 1:rcl_size]
        
        # Scegli casualmente una città dalla RCL
        next_city = rcl[rand(1:length(rcl))]
        
        # Assegna alla soluzione
        sol.route[position+1] = next_city
        visited[next_city] = true
        
        # Aggiorna obiettivo
        sol.objective += dist[current_city, next_city]
    end
    
    # Aggiungi costo del ritorno alla città iniziale
    sol.objective += dist[sol.route[end], sol.route[1]]
    
    return sol
end

# GRASP principale con K fisso
function grasp_fixedK(dist, dim, max_iterations::Int=100, K::Int=3; verbose=true)
    best_sol = nothing
    best_obj = Inf
    
    # Statistiche
    stats = SearchStats()
    
    if verbose
        println("\n" * "="^50)
        println("ESECUZIONE GRASP CON K FISSO = $K")
        println("="^50)
    end
    
    for iter in 1:max_iterations
        # FASE 1: Costruzione randomizzata con K fisso
        sol = grasp_constructive_fixedK(dist, dim, K)
        
        # FASE 2: Miglioramento con 2-opt (uso first improvement)
        first_improvement_2opt!(sol, dist, stats)
        
        # Tieni traccia della migliore soluzione
        if sol.objective < best_obj
            best_obj = sol.objective
            best_sol = deepcopy(sol)
            if verbose
                println("Iterazione $iter: nuovo best = ", round(best_obj, digits=2))
            end
        end
        
        # Stampa progresso ogni 10 iterazioni
        if verbose && iter % 10 == 0
            println("  Progresso: $iter/$max_iterations - Best corrente: ", round(best_obj, digits=2))
        end
    end
    
    if verbose
        println("\n" * "-"^40)
        println("RISULTATI GRASP (K=$K):")
        println("  Best obiettivo: ", round(best_obj, digits=2))
        println("  Iterazioni totali 2-opt: ", stats.iterations)
        println("  Miglioramenti applicati: ", stats.improvements)
    end
    
    return best_sol, stats
end

# Test GRASP con diversi valori di K
function test_grasp_different_K(dist, dim, max_iterations::Int=50)
    println("\n" * "="^60)
    println("TEST GRASP CON DIVERSI VALORI DI K")
    println("="^60)
    
    # Test diversi valori di K
    K_values = [1, 2, 3, 4, 5, 7, 10]
    results = []
    
    # Nearest Neighbor normale per confronto
    sol_nn = nearest_neighbor_heuristic(dist, dim)
    println("\nNearest Neighbor (greedy puro): ", round(sol_nn.objective, digits=2))
    
    for K in K_values
        print("\nTest con K = $K... ")
        sol, stats = grasp_fixedK(dist, dim, max_iterations, K, verbose=false)
        println("done! Best = ", round(sol.objective, digits=2))
        
        push!(results, (K, sol.objective, stats.iterations))
    end
    
    # Tabella comparativa
    println("\n" * "-"^60)
    println("TABELLA COMPARATIVA GRASP")
    println("-"^60)
    println("K     | Obiettivo | Miglioramento % | Iterazioni 2-opt")
    println("-"^60)
    
    for (K, obj, iter) in results
        improvement = 100 * (1 - obj/sol_nn.objective)
        println(rpad(K, 6), "| ", rpad(round(obj, digits=2), 10), 
                "| ", rpad(round(improvement, digits=2), 16), 
                "| ", iter)
    end
    
    # Trova il miglior K
    best_idx = argmin([r[2] for r in results])
    best_K, best_obj, best_iter = results[best_idx]
    println("\n🏆 Miglior K = $best_K con obiettivo ", round(best_obj, digits=2))
    
    return results
end

# Funzione mean per comodità
mean(v) = sum(v)/length(v)

# Funzione per eseguire e confrontare gli esperimenti
function run_experiments(filename::String)
    println("="^70)
    println("ESPERIMENTI COMPLETI: NN vs ILS vs GRASP (K fisso)")
    println("="^70)
    
    # Carica istanza
    name, coord, dim = readInstance(filename)
    dist = getDistanceMatrix(coord, dim)
    
    println("\n📊 Istanza: $name, dimensioni: $dim città")
    
    # 1. Nearest Neighbor + 2-opt
    println("\n" * "-"^50)
    println("1. NEAREST NEIGHBOR + 2-OPT")
    println("-"^50)
    
    sol_nn = nearest_neighbor_heuristic(dist, dim)
    stats_nn_best = SearchStats()
    stats_nn_first = SearchStats()
    
    sol_nn_best = deepcopy(sol_nn)
    sol_nn_first = deepcopy(sol_nn)
    
    time_best = @elapsed best_improvement_2opt!(sol_nn_best, dist, stats_nn_best)
    time_first = @elapsed first_improvement_2opt!(sol_nn_first, dist, stats_nn_first)
    
    println("  Soluzione iniziale NN: ", round(sol_nn.objective, digits=2))
    println("  Best improvement: ", round(sol_nn_best.objective, digits=2), 
            " (iter: ", stats_nn_best.iterations, ", miglioramenti: ", stats_nn_best.improvements, ")")
    println("  First improvement: ", round(sol_nn_first.objective, digits=2), 
            " (iter: ", stats_nn_first.iterations, ", miglioramenti: ", stats_nn_first.improvements, ")")
    
    # 2. ILS
    println("\n" * "-"^50)
    println("2. ITERATED LOCAL SEARCH (50 iterazioni)")
    println("-"^50)
    
    sol_ils, stats_ils = iterated_local_search(dist, dim, 50)
    println("  Obiettivo finale ILS: ", round(sol_ils.objective, digits=2))
    println("  Iterazioni 2-opt totali: ", stats_ils.iterations)
    
    # 3. GRASP con diversi K
    println("\n" * "-"^50)
    println("3. GRASP CON K FISSO")
    println("-"^50)
    
    grasp_results = test_grasp_different_K(dist, dim, 30)  # 30 iterazioni per K
    
    # 4. GRASP con K ottimale (più iterazioni)
    println("\n" * "-"^50)
    println("4. GRASP CON MIGLIOR K (100 iterazioni)")
    println("-"^50)
    
    # Trova il miglior K dai test precedenti
    best_idx = argmin([r[2] for r in grasp_results])
    best_K, _, _ = grasp_results[best_idx]
    
    println("  Eseguo GRASP con K=$best_K per 100 iterazioni...")
    sol_grasp_best, stats_grasp_best = grasp_fixedK(dist, dim, 100, best_K, verbose=true)
    
    # Tabella comparativa finale
    println("\n" * "="^70)
    println("🏁 TABELLA COMPARATIVA FINALE")
    println("="^70)
    println("Metodo                          Obiettivo    Iterazioni 2-opt    Miglioramento %")
    println("-"^70)
    
    impr_nn_best = 100 * (1 - sol_nn_best.objective/sol_nn.objective)
    impr_nn_first = 100 * (1 - sol_nn_first.objective/sol_nn.objective)
    impr_ils = 100 * (1 - sol_ils.objective/sol_nn.objective)
    impr_grasp = 100 * (1 - sol_grasp_best.objective/sol_nn.objective)
    
    println("Nearest Neighbor (iniziale)     ", rpad(round(sol_nn.objective, digits=2), 13), 
            rpad("-", 19), "-")
    println("NN + Best 2-opt                  ", rpad(round(sol_nn_best.objective, digits=2), 13), 
            rpad(stats_nn_best.iterations, 19), 
            rpad(round(impr_nn_best, digits=2), 16))
    println("NN + First 2-opt                  ", rpad(round(sol_nn_first.objective, digits=2), 13), 
            rpad(stats_nn_first.iterations, 19), 
            rpad(round(impr_nn_first, digits=2), 16))
    println("ILS (50 iter)                    ", rpad(round(sol_ils.objective, digits=2), 13), 
            rpad(stats_ils.iterations, 19), 
            rpad(round(impr_ils, digits=2), 16))
    println("GRASP (K=$best_K, 100 iter)        ", rpad(round(sol_grasp_best.objective, digits=2), 13), 
            rpad(stats_grasp_best.iterations, 19), 
            rpad(round(impr_grasp, digits=2), 16))
    
    return (sol_nn, sol_nn_best, sol_nn_first, sol_ils, sol_grasp_best, grasp_results)
end

# ===============================
# ESECUZIONE
# ===============================

println("\n" * "🚀"^30)
println("🚀 AVVIO ESPERIMENTI GRASP SU TSP")
println("🚀"^30 * "\n")

# Esegui gli esperimenti
results = run_experiments("tsp_toy50.tsp")

println("\n" * "✅"^30)
println("✅ ESPERIMENTI COMPLETATI")
println("✅"^30)