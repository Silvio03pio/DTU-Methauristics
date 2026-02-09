struct Instance
    G::Int32
    T::Int32
    S::Int32
    shared_interests::Array{Int32,2}
    partner::Array{Int32,1}
    gender::Array{Int32,1}
    relation::Array{Bool,2}
end

function Instance(filename)
    f = open(filename)
    readline(f)#header line comment

    #read instance size data
    data = split(readline(f))
    G = parse(Int32,data[1]) # no. guests
    T = parse(Int32,data[2]) # no. tables
    S = parse(Int32,data[3]) # no. sits per table

    #read common interests data
    readline(f);readline(f) #empty line + header
    shared_interests = zeros(Int32,G,G)
    for g1 in 1:G
        data = parse.(Int32,split(readline(f)))
        for g2 in 1:G
            shared_interests[g1,g2] = data[g2]
        end
    end

    #read guests' data
    readline(f);readline(f) #empty line + header
    partner = zeros(Int32,G)
    gender = zeros(Int32,G)
    relation = zeros(Bool,G,G)

    for g in 1:G
        data = split(readline(f))#read guest data
        partner[g] = parse(Int32,data[1])#read partner
        gender[g]=data[2]=="M" ? -1 : 1 #read gender Male -1, Female 1 #NOTE: modified from original
        #read the list of guest that guest g knows.
        for gg in parse.(Int32,data[3:end])
            relation[g,gg] = true
        end
    end
    close(f)

    return Instance(G, T, S, shared_interests, partner, gender, relation)

end


#= 
    I am creating a solver structure to keep all the needed parameters of the solver.
    This is not necessary for this simple exercise, but as we increase the complexity
    of the methods, more information might need to be stored.
=#
mutable struct WPSolver
    #Instance
    inst::Instance

end

#= 
    This stuct represents the solution and all the
    helper variables I used in the heuristic
=#
mutable struct WPSolution
    # short list of guests. This is composed of the first guest in a couple
    # and all the single guests. I simply assume that a gues who is in a
    # couple occupies 2 seats
    guests::Array{Int32,1}

    # number of guests in the short list
    G::Int32

    # number of seats occupied by a guest = 1 for single 2 for couple
    seats::Array{Int32,1}

    # Keeps an updated table capacity (the number of seats available in the table)
    tableCap::Array{Int32,1}

    # solution representation
    # an array of size G where the guest -> table assignment is represented
    table::Array{Int32,1}

    # dual solution representation
    # tables are represented with sets of guests, this makes the
    # calculation of the objective function easier to calculate
    sittingPlan::Array{Set{Int32},1}

    # storing the objective value
    objective::Int32
end

# Constructor for an empty solution.
function WPSolution(m::WPSolver)
    G, guests = get_short_guest_list(m)
    WPSolution(guests, 
               G, 
               [m.inst.partner[guests[g]]>0 ? 2 : 1 for g in 1:G],
               [m.inst.S for t in 1:m.inst.T],
               zeros(Int32,G),
               [Set{Int32}() for t in 1:m.inst.T],0)
end

# Generates a list of guests where all the couples are first, but only the first person 
# of a couple is present. All the single guest come after. It returns the sizi and the list.
function get_short_guest_list(m::WPSolver)
    couples = [(i,m.inst.partner[i]) for i in 1:m.inst.G if m.inst.partner[i]>0 && m.inst.partner[i]>i]
    C = length(couples)
    GG = m.inst.G-C
    guests = zeros(Int32,GG)
    # first all the couples
    for i in 1:C
        guests[i]=couples[i][1]
    end
    i = 1
    for g in 1:m.inst.G
        # then all the singles
        if m.inst.partner[g]==0
            guests[C+i]=g
            i+=1
        end
    end
    return GG, guests
end

# Simple construction heuristic that does not look at cost, but 
# finds a feasible solution.
#
# NOTE: This could have been implementde in different ways, some probably better than this
#
# ASSUMPTION: All guest couples are first in the guest array
function initial_solution(m::WPSolver)
    # create solution struct
    sol = WPSolution(m)
    t = 1 # table ID
    g = 1 # guest ID
    
    # the heuristic goes two rounds. First it assigns all the couples
    # and all the singles it can. One round is enough fo even number fo seats.
    # is the seats are odd, there will be empty seats and the heuristic
    # will fill those in the second round. If nmore than one round is needed,
    # the instance is infeasible as it would require the couples to sit in different tables
    round = 1
    # so long as we have guests to assign
    while g<=sol.G
        # check if the table has enough capacity
        if sol.tableCap[t]>=sol.seats[g]
            # assign the guest
            assign_seat(m,sol,t,g)
            # go to the next guest
            g+=1
        else
            # if there is not enough capacity, go to the next table
            t+=1
        end
        # if there are no more tables go to round 2
        if t>m.inst.T
            t=1
            if round == 1
                round = 2
            else
                # if we have been twice in round 2 we have an infeasible instance
                throw("The tables were checked twice and we could not fir the guests. Infeasible instance!")
            end
        end
    end
    # calculate the objetive value
    sol.objective = objective(m,sol)
    return sol
end

# assigns a guest to a table
function assign_seat(m::WPSolver, sol::WPSolution,t,g)
    # assign table t to the guest
    sol.table[g]=t
    # update the table seat capacity
    sol.tableCap[t]-=sol.seats[g]
    #add the guest to the dual solution representation
    push!(sol.sittingPlan[t],g)
end

# solution checker
function check_solution(m::WPSolver, sol::WPSolution)
    T = m.inst.T
    G = sol.G
    
    # keep track of how many seats we are using per table
    seats = zeros(Int32,T)

    # go though every guest
    for g in 1:G
        # check that a table has been assinged
        if sol.table[g]==0
            throw("Guest $(sol.guests[g]) has not been assigned a table!")
        # check that the table ID is valid
        elseif sol.table[g]>T|
            throw("Guest $(sol.guests[g]) is seated at table $(sol.table[g]) which does not exist!")
        else
            # update the number of available seats
            seats[sol.table[g]]+=1
        end
        if m.inst.partner[sol.guests[g]]>0
            # update the number of available seats in case of a couple
            seats[sol.table[g]]+=1
        end
    end

    # check that we do not exceed the number of seats per table
    if maximum(seats)>m.inst.S
        throw("One of the table has more guests than there are chair!")
    end   
end


# calculate the value of a single table, where the table is represented by a set of guests
function tableValue(m::WPSolver, sol::WPSolution, table::Set{Int32})
    genderBalance = 0
    interests = 0
    relations = 0

    # create a table set with the original guest IDs
    TT = Set{Int32}()
    # for each guest
    for g in table
        #add the orginal guest ID to the set
        push!(TT,sol.guests[g])
        if sol.seats[g]>1
            # add the partner guest ID to the set
            push!(TT,m.inst.partner[sol.guests[g]])
        end
    end
    # for each guest
    for g in TT
        # updated the gender balance. Note that I have changed gender values to be 1 for females and -1 for man
        genderBalance+=m.inst.gender[g]
        for gg in TT 
            # for each other guest in the table (only in half the matrix)
            if g<gg
                # update relation and interest value
                relations += m.inst.relation[g,gg] ? 1 : 0
                interests += m.inst.shared_interests[g,gg]
            end
        end
    end
    # calculatde gender balance objective over 2
    genderBalance = abs(genderBalance)
    genderCost = genderBalance>2 ? genderBalance-2 : 0

    #return full objective for the table
    return interests + relations - genderCost

end

# function to calculate the complete objective value
function objective(m::WPSolver, sol::WPSolution)
    obj = 0
    for t in 1:m.inst.T
        # sum the value of each table
        obj+= tableValue(m,sol,sol.sittingPlan[t])
    end
    return obj
end

# print solution on the terminal
function print_solution(m::WPSolver, sol::WPSolution)
    for t in 1:m.inst.T
        print("Table $t [")
        for g in sol.sittingPlan[t]
            print(" $(sol.guests[g])")
            if sol.seats[g]>1
                print(" $(m.inst.partner[sol.guests[g]])")
            end
        end
        println(" ]")
    end
    println("Obj: ",objective(m,sol))
end

function main()
    m = WPSolver(Instance("WeddingData_100_10_50_40.dat"))
    sol = initial_solution(m)
    print_solution(m,sol)
end

main()