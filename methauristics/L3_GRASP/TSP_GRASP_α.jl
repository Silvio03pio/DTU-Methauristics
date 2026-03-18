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

# ============================================================
# GRASP CON PARAMETRO ALPHA (CORRETTO CON FLOAT32)
# ============================================================

"""
    grasp_constructive_alpha(dist, dim, alpha)

Costruisce una soluzione TSP usando GRASP con parametro alpha.
Alpha = 0 -> puramente greedy (solo la città migliore)
Alpha = 1 -> completamente casuale (tutte le città)
"""
function grasp_constructive_alpha(dist, dim, alpha::Float32)
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
        
        # Calcola distanze
        distances = [dist[current_city, c] for c in unvisited]
        
        # Trova min e max
        c_min = minimum(distances)
        c_max = maximum(distances)
        
        # Calcola soglia: c_min + alpha * (c_max - c_min)
        threshold = c_min + alpha * (c_max - c_min)
        
        # Crea RCL: città con distanza <= threshold
        rcl = []
        for (idx, c) in enumerate(unvisited)
            if distances[idx] <= threshold
                push!(rcl, c)
            end
        end
        
        # Se RCL è vuota (caso limite), prendi la città migliore
        if isempty(rcl)
            rcl = [unvisited[argmin(distances)]]
        end
        
        # Scegli casualmente dalla RCL
        next_city = rcl[rand(1:length(rcl))]
        
        # Assegna alla soluzione
        sol.route[position+1] = next_city
        visited[next_city] = true
        
        # Aggiorna obiettivo
        sol.objective += dist[current_city, next_city]
    end
    
    # Aggiungi costo del ritorno
    sol.objective += dist[sol.route[end], sol.route[1]]
    
    return sol
end

# GRASP principale con alpha (ORA ACCETTA FLOAT64 E CONVERTE)
function grasp_alpha(dist, dim, max_iterations::Int=100, alpha::Float64=0.3; verbose=true)
    # Converti alpha in Float32 per la funzione di costruzione
    alpha32 = Float32(alpha)
    best_sol = nothing
    best_obj = Inf
    
    # Statistiche
    stats = SearchStats()
    
    if verbose
        println("\n" * "="^50)
        println("ESECUZIONE GRASP CON ALPHA = $alpha")
        println("="^50)
    end
    
    for iter in 1:max_iterations
        # FASE 1: Costruzione randomizzata con alpha
        sol = grasp_constructive_alpha(dist, dim, alpha32)
        
        # FASE 2: Miglioramento con 2-opt
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
        println("RISULTATI GRASP (alpha=$alpha):")
        println("  Best obiettivo: ", round(best_obj, digits=2))
        println("  Iterazioni totali 2-opt: ", stats.iterations)
        println("  Miglioramenti applicati: ", stats.improvements)
    end
    
    return best_sol, stats
end

# ============================================================
# GRASP CON K FISSO
# ============================================================

function grasp_constructive_fixedK(dist, dim, K::Int)
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
        
        # Prendi le prime K città
        rcl_size = min(K, length(unvisited))
        rcl = [distances[i][1] for i in 1:rcl_size]
        
        # Scegli casualmente dalla RCL
        next_city = rcl[rand(1:length(rcl))]
        
        # Assegna alla soluzione
        sol.route[position+1] = next_city
        visited[next_city] = true
        
        # Aggiorna obiettivo
        sol.objective += dist[current_city, next_city]
    end
    
    # Aggiungi costo del ritorno
    sol.objective += dist[sol.route[end], sol.route[1]]
    
    return sol
end

function grasp_fixedK(dist, dim, max_iterations::Int=100, K::Int=3; verbose=true)
    best_sol = nothing
    best_obj = Inf
    stats = SearchStats()
    
    if verbose
        println("\n" * "="^50)
        println("ESECUZIONE GRASP CON K FISSO = $K")
        println("="^50)
    end
    
    for iter in 1:max_iterations
        sol = grasp_constructive_fixedK(dist, dim, K)
        first_improvement_2opt!(sol, dist, stats)
        
        if sol.objective < best_obj
            best_obj = sol.objective
            best_sol = deepcopy(sol)
            if verbose
                println("Iterazione $iter: nuovo best = ", round(best_obj, digits=2))
            end
        end
        
        if verbose && iter % 10 == 0
            println("  Progresso: $iter/$max_iterations - Best corrente: ", round(best_obj, digits=2))
        end
    end
    
    if verbose
        println("\n" * "-"^40)
        println("RISULTATI GRASP (K=$K):")
        println("  Best obiettivo: ", round(best_obj, digits=2))
        println("  Iterazioni totali 2-opt: ", stats.iterations)
    end
    
    return best_sol, stats
end

# ============================================================
# FUNZIONI DI TEST E CONFRONTO
# ============================================================

# Test GRASP con diversi valori di alpha
function test_grasp_different_alpha(dist, dim, max_iterations::Int=30)
    println("\n" * "="^70)
    println("📊 TEST GRASP CON DIVERSI VALORI DI ALPHA")
    println("="^70)
    
    # Test diversi valori di alpha (come Float64)
    alpha_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    results = []
    
    # Nearest Neighbor per confronto
    sol_nn = nearest_neighbor_heuristic(dist, dim)
    println("\n📌 Nearest Neighbor (greedy puro): ", round(sol_nn.objective, digits=2))
    println("\n" * "-"^70)
    println("Alpha   | Obiettivo | Miglioramento % | Iterazioni 2-opt")
    println("-"^70)
    
    for alpha in alpha_values
        print("Test alpha = $alpha... ")
        sol, stats = grasp_alpha(dist, dim, max_iterations, alpha, verbose=false)
        println("done! → ", round(sol.objective, digits=2))
        
        improvement = 100 * (1 - sol.objective/sol_nn.objective)
        push!(results, (alpha, sol.objective, improvement, stats.iterations))
        
        println(rpad(alpha, 8), "| ", rpad(round(sol.objective, digits=2), 9), 
                "| ", rpad(round(improvement, digits=2), 16), 
                "| ", stats.iterations)
    end
    
    # Trova il miglior alpha
    best_idx = argmin([r[2] for r in results])
    best_alpha, best_obj, best_impr, best_iter = results[best_idx]
    println("\n🏆 MIGLIOR ALPHA = $best_alpha con obiettivo ", round(best_obj, digits=2))
    
    return results
end

# Test GRASP con diversi K (per confronto)
function test_grasp_different_K(dist, dim, max_iterations::Int=30)
    println("\n" * "="^70)
    println("📊 TEST GRASP CON DIVERSI VALORI DI K")
    println("="^70)
    
    K_values = [1, 2, 3, 4, 5, 7, 10, 15]
    results = []
    
    sol_nn = nearest_neighbor_heuristic(dist, dim)
    println("\n📌 Nearest Neighbor: ", round(sol_nn.objective, digits=2))
    println("\n" * "-"^70)
    println("K       | Obiettivo | Miglioramento % | Iterazioni 2-opt")
    println("-"^70)
    
    for K in K_values
        print("Test K = $K... ")
        sol, stats = grasp_fixedK(dist, dim, max_iterations, K, verbose=false)
        println("done! → ", round(sol.objective, digits=2))
        
        improvement = 100 * (1 - sol.objective/sol_nn.objective)
        push!(results, (K, sol.objective, improvement, stats.iterations))
        
        println(rpad(K, 8), "| ", rpad(round(sol.objective, digits=2), 9), 
                "| ", rpad(round(improvement, digits=2), 16), 
                "| ", stats.iterations)
    end
    
    best_idx = argmin([r[2] for r in results])
    best_K, best_obj, best_impr, best_iter = results[best_idx]
    println("\n🏆 MIGLIOR K = $best_K con obiettivo ", round(best_obj, digits=2))
    
    return results
end

# Funzione per confrontare direttamente alpha e K
function compare_alpha_vs_K(filename::String, max_iter::Int=30)
    println("\n" * "🚀"^60)
    println("🚀 CONFRONTO GRASP: ALPHA vs K FISSO")
    println("🚀"^60)
    
    # Carica istanza
    name, coord, dim = readInstance(filename)
    dist = getDistanceMatrix(coord, dim)
    
    println("\n📊 Istanza: $name, dimensioni: $dim città")
    
    # Test alpha
    println("\n" * "🔵"^40)
    println("🔵 TEST CON ALPHA")
    println("🔵"^40)
    alpha_results = test_grasp_different_alpha(dist, dim, max_iter)
    
    # Test K
    println("\n" * "🟢"^40)
    println("🟢 TEST CON K FISSO")
    println("🟢"^40)
    k_results = test_grasp_different_K(dist, dim, max_iter)
    
    # Estrai migliori
    best_alpha_idx = argmin([r[2] for r in alpha_results])
    best_alpha_val, best_alpha_obj, best_alpha_impr, _ = alpha_results[best_alpha_idx]
    
    best_k_idx = argmin([r[2] for r in k_results])
    best_k_val, best_k_obj, best_k_impr, _ = k_results[best_k_idx]
    
    # Tabella comparativa
    println("\n" * "="^70)
    println("🏁 CONFRONTO DIRETTO: MIGLIOR ALPHA vs MIGLIOR K")
    println("="^70)
    println("Metodo          | Valore | Obiettivo | Miglioramento %")
    println("-"^70)
    println("Miglior Alpha   | ", rpad(best_alpha_val, 7), "| ", rpad(round(best_alpha_obj, digits=2), 9), 
            "| ", round(best_alpha_impr, digits=2))
    println("Miglior K       | ", rpad(best_k_val, 7), "| ", rpad(round(best_k_obj, digits=2), 9), 
            "| ", round(best_k_impr, digits=2))
    
    # Calcola differenza percentuale
    if best_alpha_obj < best_k_obj
        diff = 100 * (1 - best_alpha_obj/best_k_obj)
        println("\n✅ Alpha batte K di ", round(diff, digits=2), "%")
    else
        diff = 100 * (1 - best_k_obj/best_alpha_obj)
        println("\n✅ K batte Alpha di ", round(diff, digits=2), "%")
    end
    
    return (alpha_results, k_results)
end

# Funzione per eseguire GRASP con alpha ottimale (più iterazioni)
function run_optimal_grasp(filename::String, max_iter::Int=100)
    println("\n" * "🌟"^50)
    println("🌟 GRASP CON PARAMETRI OTTIMALI")
    println("🌟"^50)
    
    name, coord, dim = readInstance(filename)
    dist = getDistanceMatrix(coord, dim)
    
    println("\nIstanza: $name, dimensioni: $dim città")
    
    # Trova miglior alpha con test rapido
    println("\n🔍 Ricerca miglior alpha (test rapido)...")
    alpha_test = test_grasp_different_alpha(dist, dim, 20)
    best_alpha_idx = argmin([r[2] for r in alpha_test])
    best_alpha, best_alpha_test_obj, _, _ = alpha_test[best_alpha_idx]
    println("  Miglior alpha trovato: $best_alpha (obiettivo test: ", round(best_alpha_test_obj, digits=2), ")")
    
    # Trova miglior K con test rapido
    println("\n🔍 Ricerca miglior K (test rapido)...")
    k_test = test_grasp_different_K(dist, dim, 20)
    best_k_idx = argmin([r[2] for r in k_test])
    best_k, best_k_test_obj, _, _ = k_test[best_k_idx]
    println("  Miglior K trovato: $best_k (obiettivo test: ", round(best_k_test_obj, digits=2), ")")
    
    # Esegui GRASP con alpha ottimale
    println("\n" * "-"^50)
    println("🚀 ESECUZIONE GRASP CON ALPHA=$best_alpha ($max_iter iterazioni)")
    println("-"^50)
    sol_alpha, stats_alpha = grasp_alpha(dist, dim, max_iter, best_alpha, verbose=true)
    
    # Esegui GRASP con K ottimale
    println("\n" * "-"^50)
    println("🚀 ESECUZIONE GRASP CON K=$best_k ($max_iter iterazioni)")
    println("-"^50)
    sol_k, stats_k = grasp_fixedK(dist, dim, max_iter, best_k, verbose=true)
    
    # Confronto finale
    println("\n" * "="^70)
    println("🏆 CONFRONTO FINALE (con $max_iter iterazioni)")
    println("="^70)
    println("Metodo          | Obiettivo | Iterazioni 2-opt | Miglioramento %")
    println("-"^70)
    
    sol_nn = nearest_neighbor_heuristic(dist, dim)
    impr_alpha = 100 * (1 - sol_alpha.objective/sol_nn.objective)
    impr_k = 100 * (1 - sol_k.objective/sol_nn.objective)
    
    println("Nearest Neighbor | ", rpad(round(sol_nn.objective, digits=2), 9), 
            "| -               | -")
    println("GRASP Alpha=$best_alpha | ", rpad(round(sol_alpha.objective, digits=2), 9), 
            "| ", rpad(stats_alpha.iterations, 15), "| ", round(impr_alpha, digits=2))
    println("GRASP K=$best_k     | ", rpad(round(sol_k.objective, digits=2), 9), 
            "| ", rpad(stats_k.iterations, 15), "| ", round(impr_k, digits=2))
    
    return (sol_alpha, stats_alpha, sol_k, stats_k, best_alpha, best_k)
end

# Funzione mean per comodità
mean(v) = sum(v)/length(v)

# ============================================================
# ESECUZIONE PRINCIPALE
# ============================================================

println("\n" * "🎯"^40)
println("🎯 GRASP PER TSP: CONFRONTO ALPHA vs K FISSO")
println("🎯"^40)

# Scegli quale test eseguire:

# 1. Confronto rapido alpha vs K
println("\n" * "⚡"^40)
println("⚡ TEST RAPIDO (30 iterazioni per parametro)")
println("⚡"^40)
compare_alpha_vs_K("tsp_toy.tsp", 30)

# 2. Esecuzione con parametri ottimali (più iterazioni)
println("\n" * "🔥"^40)
println("🔥 TEST APPROFONDITO (100 iterazioni con migliori parametri)")
println("🔥"^40)
run_optimal_grasp("tsp_toy.tsp", 100)

println("\n" * "✅"^40)
println("✅ TEST COMPLETATI")
println("✅"^40)