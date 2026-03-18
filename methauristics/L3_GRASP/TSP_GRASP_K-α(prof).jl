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

# ======================================================================
# GRASP IMPLEMENTATION (CON CITTÀ INIZIALE FISSA = 1)
# ======================================================================

# Struct per statistiche
mutable struct GRASPStats
    iterations::Int
    improvements::Int
    construction_time::Float64
    search_time::Float64
    objective_values::Vector{Float32}
    
    GRASPStats() = new(0, 0, 0.0, 0.0, Float32[])
end

# GRASP costruttivo con K fisso (città iniziale = 1)
function grasp_constructive_fixedK(m::TSPSolver, K::Int)
    sol = TSPSolution(m.dim)
    visited = zeros(Bool, m.dim)
    
    # Città iniziale FISSA = 1 (come Nearest Neighbor)
    sol.route[1] = 1
    visited[1] = true
    
    for position in 1:m.dim-1
        current_city = sol.route[position]
        
        # Trova tutte le città non visitate
        unvisited = [c for c in 1:m.dim if !visited[c]]
        
        # Calcola distanze per ogni città non visitata
        distances = [(city, m.dist[current_city, city]) for city in unvisited]
        
        # Ordina per distanza crescente
        sort!(distances, by = x -> x[2])
        
        # Prendi le prime K città
        rcl_size = min(K, length(unvisited))
        rcl = [distances[i][1] for i in 1:rcl_size]
        
        # Scegli casualmente una città dalla RCL
        next_city = rcl[rand(1:length(rcl))]
        
        # Assegna alla soluzione
        sol.route[position+1] = next_city
        visited[next_city] = true
        
        # Aggiorna obiettivo
        sol.objective += m.dist[current_city, next_city]
    end
    
    # Aggiungi costo del ritorno
    sol.objective += m.dist[sol.route[end], sol.route[1]]
    
    return sol
end

# GRASP costruttivo con alpha (città iniziale = 1)
function grasp_constructive_alpha(m::TSPSolver, alpha::Float64)
    alpha32 = Float32(alpha)
    sol = TSPSolution(m.dim)
    visited = zeros(Bool, m.dim)
    
    # Città iniziale FISSA = 1 (come Nearest Neighbor)
    sol.route[1] = 1
    visited[1] = true
    
    for position in 1:m.dim-1
        current_city = sol.route[position]
        
        # Trova tutte le città non visitate
        unvisited = [c for c in 1:m.dim if !visited[c]]
        
        # Calcola distanze
        distances = [m.dist[current_city, c] for c in unvisited]
        
        # Trova min e max
        c_min = minimum(distances)
        c_max = maximum(distances)
        
        # Calcola soglia: c_min + alpha * (c_max - c_min)
        threshold = c_min + alpha32 * (c_max - c_min)
        
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
        sol.objective += m.dist[current_city, next_city]
    end
    
    # Aggiungi costo del ritorno
    sol.objective += m.dist[sol.route[end], sol.route[1]]
    
    return sol
end

# GRASP principale con K fisso
function grasp_fixedK(m::TSPSolver, max_iterations::Int, K::Int; verbose=true)
    best_sol = nothing
    best_obj = Inf
    stats = GRASPStats()
    
    if verbose
        println("\n" * "="^60)
        println("GRASP CON K FISSO = $K ($max_iterations iterazioni)")
        println("="^60)
    end
    
    for iter in 1:max_iterations
        # Controllo time limit
        if elapsed_time(m) >= m.timelimit
            if verbose
                println("\n⚠️ Time limit raggiunto dopo $iter iterazioni")
            end
            break
        end
        
        # Fase costruttiva
        constr_time = @elapsed sol = grasp_constructive_fixedK(m, K)
        stats.construction_time += constr_time
        
        # Fase di ricerca locale
        search_time = @elapsed local_search!(m, sol)
        stats.search_time += search_time
        stats.iterations += 1
        
        # Aggiorna migliore soluzione
        if sol.objective < best_obj
            best_obj = sol.objective
            best_sol = deepcopy(sol)
            stats.improvements += 1
            push!(stats.objective_values, best_obj)
            if verbose
                println("Iterazione $iter: NUOVO BEST = ", round(best_obj, digits=2))
            end
        end
        
        # Stampa progresso
        if verbose && iter % 10 == 0
            println("  Progresso: $iter/$max_iterations - Best corrente: ", round(best_obj, digits=2))
        end
    end
    
    if verbose
        println("\n" * "-"^60)
        println("RISULTATI GRASP K=$K:")
        println("  Best obiettivo: ", round(best_obj, digits=2))
        println("  Iterazioni completate: ", stats.iterations)
        println("  Miglioramenti trovati: ", stats.improvements)
        println("  Tempo costruzione: ", round(stats.construction_time, digits=3), "s")
        println("  Tempo ricerca locale: ", round(stats.search_time, digits=3), "s")
        println("  Tempo totale: ", round(stats.construction_time + stats.search_time, digits=3), "s")
    end
    
    return best_sol, stats
end

# GRASP principale con alpha
function grasp_alpha(m::TSPSolver, max_iterations::Int, alpha::Float64; verbose=true)
    best_sol = nothing
    best_obj = Inf
    stats = GRASPStats()
    
    if verbose
        println("\n" * "="^60)
        println("GRASP CON ALPHA = $alpha ($max_iterations iterazioni)")
        println("="^60)
    end
    
    for iter in 1:max_iterations
        # Controllo time limit
        if elapsed_time(m) >= m.timelimit
            if verbose
                println("\n⚠️ Time limit raggiunto dopo $iter iterazioni")
            end
            break
        end
        
        # Fase costruttiva
        constr_time = @elapsed sol = grasp_constructive_alpha(m, alpha)
        stats.construction_time += constr_time
        
        # Fase di ricerca locale
        search_time = @elapsed local_search!(m, sol)
        stats.search_time += search_time
        stats.iterations += 1
        
        # Aggiorna migliore soluzione
        if sol.objective < best_obj
            best_obj = sol.objective
            best_sol = deepcopy(sol)
            stats.improvements += 1
            push!(stats.objective_values, best_obj)
            if verbose
                println("Iterazione $iter: NUOVO BEST = ", round(best_obj, digits=2))
            end
        end
        
        # Stampa progresso
        if verbose && iter % 10 == 0
            println("  Progresso: $iter/$max_iterations - Best corrente: ", round(best_obj, digits=2))
        end
    end
    
    if verbose
        println("\n" * "-"^60)
        println("RISULTATI GRASP alpha=$alpha:")
        println("  Best obiettivo: ", round(best_obj, digits=2))
        println("  Iterazioni completate: ", stats.iterations)
        println("  Miglioramenti trovati: ", stats.improvements)
        println("  Tempo costruzione: ", round(stats.construction_time, digits=3), "s")
        println("  Tempo ricerca locale: ", round(stats.search_time, digits=3), "s")
        println("  Tempo totale: ", round(stats.construction_time + stats.search_time, digits=3), "s")
    end
    
    return best_sol, stats
end

# Test per verificare che K=1 dia lo stesso risultato di NN
function test_k1_equals_nn(m::TSPSolver)
    println("\n" * "🔍"^50)
    println("🔍 TEST: K=1 vs NEAREST NEIGHBOR")
    println("🔍"^50)
    
    # Nearest Neighbor + Local Search
    m_nn = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_nn = local_search(m_nn)
    println("\n📌 Nearest Neighbor + Local Search:")
    println("  Obiettivo: ", round(sol_nn.objective, digits=2))
    println("  Route: ", sol_nn.route)
    
    # GRASP K=1 + Local Search (una sola iterazione)
    m_k1 = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_k1, stats = grasp_fixedK(m_k1, 1, 1, verbose=false)
    println("\n📌 GRASP K=1 (1 iterazione):")
    println("  Obiettivo: ", round(sol_k1.objective, digits=2))
    println("  Route: ", sol_k1.route)
    
    # Verifica uguaglianza
    println("\n" * "-"^50)
    if sol_nn.objective == sol_k1.objective && sol_nn.route == sol_k1.route
        println("✅ SUCCESSO: K=1 produce lo STESSO risultato di NN!")
    else
        println("❌ ERRORE: K=1 produce un risultato DIVERSO da NN")
        println("\nDifferenza obiettivo: ", abs(sol_nn.objective - sol_k1.objective))
    end
end

# Test GRASP con diversi K (con città iniziale fissa)
function test_grasp_different_K(m::TSPSolver, max_iterations::Int=30)
    println("\n" * "📊"^40)
    println("📊 TEST GRASP CON DIVERSI VALORI DI K")
    println("📊"^40)
    
    K_values = [1, 2, 3, 4, 5, 7, 10, 15]
    results = []
    
    # Nearest Neighbor per confronto
    m_nn = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_nn = local_search(m_nn)
    println("\n📌 Nearest Neighbor: ", round(sol_nn.objective, digits=2))
    
    println("\n" * "-"^70)
    println("K       | Obiettivo | Miglioramento % | Iterazioni | Tempo (s)")
    println("-"^70)
    
    for K in K_values
        print("Test K = $K... ")
        # Resetta il timer per ogni test
        m_copy = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
        sol, stats = grasp_fixedK(m_copy, max_iterations, K, verbose=false)
        improvement = 100 * (1 - sol.objective/sol_nn.objective)
        total_time = stats.construction_time + stats.search_time
        
        push!(results, (K, sol.objective, improvement, stats.iterations, total_time))
        
        println("done! → ", round(sol.objective, digits=2))
        println(rpad(K, 8), "| ", rpad(round(sol.objective, digits=2), 9), 
                "| ", rpad(round(improvement, digits=2), 16), 
                "| ", rpad(stats.iterations, 11),
                "| ", round(total_time, digits=3))
    end
    
    # Trova il miglior K
    best_idx = argmin([r[2] for r in results])
    best_K, best_obj, best_impr, best_iter, best_time = results[best_idx]
    println("\n🏆 MIGLIOR K = $best_K con obiettivo ", round(best_obj, digits=2))
    
    return results
end

# Test GRASP con diversi alpha (con città iniziale fissa)
function test_grasp_different_alpha(m::TSPSolver, max_iterations::Int=30)
    println("\n" * "📊"^40)
    println("📊 TEST GRASP CON DIVERSI VALORI DI ALPHA")
    println("📊"^40)
    
    alpha_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    results = []
    
    # Nearest Neighbor per confronto
    m_nn = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_nn = local_search(m_nn)
    println("\n📌 Nearest Neighbor: ", round(sol_nn.objective, digits=2))
    
    println("\n" * "-"^70)
    println("Alpha   | Obiettivo | Miglioramento % | Iterazioni | Tempo (s)")
    println("-"^70)
    
    for alpha in alpha_values
        print("Test alpha = $alpha... ")
        # Resetta il timer per ogni test
        m_copy = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
        sol, stats = grasp_alpha(m_copy, max_iterations, alpha, verbose=false)
        improvement = 100 * (1 - sol.objective/sol_nn.objective)
        total_time = stats.construction_time + stats.search_time
        
        push!(results, (alpha, sol.objective, improvement, stats.iterations, total_time))
        
        println("done! → ", round(sol.objective, digits=2))
        println(rpad(alpha, 8), "| ", rpad(round(sol.objective, digits=2), 9), 
                "| ", rpad(round(improvement, digits=2), 16), 
                "| ", rpad(stats.iterations, 11),
                "| ", round(total_time, digits=3))
    end
    
    # Trova il miglior alpha
    best_idx = argmin([r[2] for r in results])
    best_alpha, best_obj, best_impr, best_iter, best_time = results[best_idx]
    println("\n🏆 MIGLIOR ALPHA = $best_alpha con obiettivo ", round(best_obj, digits=2))
    
    return results
end

# Confronto diretto tra alpha e K
function compare_alpha_vs_K(m::TSPSolver, max_iter::Int=30)
    println("\n" * "🚀"^70)
    println("🚀 CONFRONTO GRASP: ALPHA vs K FISSO")
    println("🚀"^70)
    
    # Test alpha
    println("\n" * "🔵"^50)
    println("🔵 TEST CON ALPHA")
    println("🔵"^50)
    alpha_results = test_grasp_different_alpha(m, max_iter)
    
    # Test K
    println("\n" * "🟢"^50)
    println("🟢 TEST CON K FISSO")
    println("🟢"^50)
    k_results = test_grasp_different_K(m, max_iter)
    
    # Estrai migliori
    best_alpha_idx = argmin([r[2] for r in alpha_results])
    best_alpha_val, best_alpha_obj, best_alpha_impr, _, _ = alpha_results[best_alpha_idx]
    
    best_k_idx = argmin([r[2] for r in k_results])
    best_k_val, best_k_obj, best_k_impr, _, _ = k_results[best_k_idx]
    
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
    
    # Calcola differenza
    if best_alpha_obj < best_k_obj
        diff = 100 * (1 - best_alpha_obj/best_k_obj)
        println("\n✅ Alpha batte K di ", round(diff, digits=2), "%")
    else
        diff = 100 * (1 - best_k_obj/best_alpha_obj)
        println("\n✅ K batte Alpha di ", round(diff, digits=2), "%")
    end
    
    return (alpha_results, k_results)
end

# GRASP con parametri ottimali
function run_optimal_grasp(m::TSPSolver, max_iter::Int=100)
    println("\n" * "🌟"^60)
    println("🌟 GRASP CON PARAMETRI OTTIMALI")
    println("🌟"^60)
    
    # Trova miglior alpha con test rapido
    println("\n🔍 Ricerca miglior alpha (test rapido)...")
    m_alpha = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    alpha_test = test_grasp_different_alpha(m_alpha, 20)
    best_alpha_idx = argmin([r[2] for r in alpha_test])
    best_alpha, best_alpha_test_obj, _, _, _ = alpha_test[best_alpha_idx]
    println("  Miglior alpha trovato: $best_alpha (obiettivo test: ", round(best_alpha_test_obj, digits=2), ")")
    
    # Trova miglior K con test rapido
    println("\n🔍 Ricerca miglior K (test rapido)...")
    m_k = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    k_test = test_grasp_different_K(m_k, 20)
    best_k_idx = argmin([r[2] for r in k_test])
    best_k, best_k_test_obj, _, _, _ = k_test[best_k_idx]
    println("  Miglior K trovato: $best_k (obiettivo test: ", round(best_k_test_obj, digits=2), ")")
    
    # Esegui GRASP con alpha ottimale
    println("\n" * "-"^60)
    println("🚀 ESECUZIONE GRASP CON ALPHA=$best_alpha ($max_iter iterazioni)")
    println("-"^60)
    m_alpha_opt = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_alpha, stats_alpha = grasp_alpha(m_alpha_opt, max_iter, best_alpha, verbose=true)
    
    # Esegui GRASP con K ottimale
    println("\n" * "-"^60)
    println("🚀 ESECUZIONE GRASP CON K=$best_k ($max_iter iterazioni)")
    println("-"^60)
    m_k_opt = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_k, stats_k = grasp_fixedK(m_k_opt, max_iter, best_k, verbose=true)
    
    # Confronto finale
    println("\n" * "="^70)
    println("🏆 CONFRONTO FINALE (con $max_iter iterazioni)")
    println("="^70)
    println("Metodo          | Obiettivo | Iterazioni | Miglioramento % | Tempo (s)")
    println("-"^70)
    
    m_nn = TSPSolver(m.coord, m.dim, m.dist, m.timelimit)
    sol_nn = local_search(m_nn)
    impr_alpha = 100 * (1 - sol_alpha.objective/sol_nn.objective)
    impr_k = 100 * (1 - sol_k.objective/sol_nn.objective)
    
    println("Nearest Neighbor | ", rpad(round(sol_nn.objective, digits=2), 9), 
            "| -         | -               | -")
    println("GRASP Alpha=$best_alpha | ", rpad(round(sol_alpha.objective, digits=2), 9), 
            "| ", rpad(stats_alpha.iterations, 10), 
            "| ", rpad(round(impr_alpha, digits=2), 15),
            "| ", round(stats_alpha.construction_time + stats_alpha.search_time, digits=3))
    println("GRASP K=$best_k     | ", rpad(round(sol_k.objective, digits=2), 9), 
            "| ", rpad(stats_k.iterations, 10), 
            "| ", rpad(round(impr_k, digits=2), 15),
            "| ", round(stats_k.construction_time + stats_k.search_time, digits=3))
    
    return (sol_alpha, stats_alpha, sol_k, stats_k, best_alpha, best_k)

end 
# ======================================================================
# MAIN FUNCTION WITH GRASP
# ======================================================================

function main_with_grasp()
    println("\n" * "🎯"^50)
    println("🎯 TSP SOLVER WITH GRASP (K FIXED & ALPHA)")
    println("🎯"^50)
    
    # read the instance
    name, coord, dim = readInstance("tsp_toy50.tsp")
    # get the distance matrix
    dist = getDistanceMatrix(coord, dim)
    # create the solver (time limit 30 secondi per tutti i test)
    m = TSPSolver(coord, dim, dist, 30)
    
    println("\n📊 Istanza: $name, dimensioni: $dim città, time limit: 30s")
    
    # 1. TEST K=1 vs NN (per verificare che siano uguali)
    test_k1_equals_nn(m)
    
    # 2. Nearest Neighbor + Local Search
    println("\n" * "-"^50)
    println("1. NEAREST NEIGHBOR + LOCAL SEARCH")
    println("-"^50)
    m_nn = TSPSolver(coord, dim, dist, 30)
    start_time = time_ns()
    sol_nn = local_search(m_nn)
    nn_time = (time_ns() - start_time) / 1e9
    println("  Obiettivo: ", round(sol_nn.objective, digits=2))
    println("  Tempo: ", round(nn_time, digits=3), "s")
    
    # 3. Confronto GRASP Alpha vs K
    compare_alpha_vs_K(m, 30)
    
    # 4. GRASP con parametri ottimali
    m_opt = TSPSolver(coord, dim, dist, 30)
    run_optimal_grasp(m_opt, 50)
    
    println("\n" * "✅"^50)
    println("✅ ESPERIMENTI COMPLETATI")
    println("✅"^50)
end

# Esegui la funzione principale
main_with_grasp()