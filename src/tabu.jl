#=
tabu:
- Julia version: 
- Author: sujin
- Date: 2022-11-03
=#

include("tabu_movs.jl")

#=

procède à une recherche tabou

=#

function tabu(obj::Int64, sol::Vector{Vector{Int64}},Q::Int64,tenure::Int64,k::Float64,
                c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
                d::Array{Float64,2})
    cpt_pas_amelioration::Int64 = 0
    nb_iterations = 1
    best_solution = deepcopy(sol)
    best_candidate = deepcopy(sol)
    tabu_list = [sol for i in 1:tenure]
    start = time()
    timelim = 1
    while (time()<start + timelim)

        neighborhood = getNeighbors_swap(best_candidate,b,s)
        #neighborhood = getNeighbors_shift(best_candidate,Q,b,s)
        best_candidate = neighborhood[1]
        found_amelio = false
        i = 1

        #performs the first non-tabu move that leads to an improvement of the current solution
        # or a tabu move that reaches the aspiration move (< best solution)
        while !found_amelio && i < length(neighborhood)
            
            if !(neighborhood[i] in tabu_list)
                #display(neighborhood[i][2])
                if evaluate_solution(obj,neighborhood[i],d,c,b,s) < evaluate_solution(obj,best_candidate,d,c,b,s)
                #if gain_swap(obj,best_candidate,neighborhood[i][2],neighborhood[i][1][1][1],neighborhood[i][1][1][2],c,b,s,d) < 0
                found_amelio = true
                    best_candidate = neighborhood[i]
                end
            else

                if evaluate_solution(obj,neighborhood[i],d,c,b,s) < evaluate_solution(obj,best_solution,d,c,b,s)
                    best_candidate = neighborhood[i]
                    found_amelio = true
                end    
            end
            
            i += 1
        end

        if !found_amelio
            best_candidate = tabu_list[findmin(x -> evaluate_solution(obj,x,d,c,b,s),tabu_list)[2]]
            cpt_pas_amelioration += 1  

        else 
            if evaluate_solution(obj,best_candidate,d,c,b,s) < evaluate_solution(obj,best_solution,d,c,b,s)
                best_solution = best_candidate
                cpt_pas_amelioration = 0
            end
        end

        tabu_list[(nb_iterations+1)%tenure + 1] = best_candidate
        nb_iterations += 1
        #=
        println(cpt_pas_amelioration)
        if (nb_iterations % 10 == 0)
            println(evaluate_solution(best_solution,d,c,b,s))
        end    
        =#
    end

    return best_solution
end

function getNeighbors_shift(sol::Vector{Vector{Int64}},Q::Int64,
                            b::Array{Float64,2},s::Vector{Float64})::Vector{Vector{Vector{Vector{Int64}}}}
    
    neighbors = Array{Vector{Vector{Vector{Int64}}},2}(undef,(length(sol[1]),length(sol[2])))

    for i in 1:length(sol[1])
        for j in 1:length(sol[2])
            neighbors[i,j] = [[[i,j]],shift(sol,Q,i,sol[1][i],j,b,s)]
        end
    end
    
    return vec(neighbors)

end    


function getNeighbors_swap(sol::Vector{Vector{Int64}},
    b::Array{Float64,2},s::Vector{Float64})::Vector{Vector{Vector{Int64}}}

    neighbors = Array{Vector{Vector{Int64}},2}(undef,(length(sol[3]),length(setdiff(1:length(sol[2]),sol[3]))))

    for i in 1:length(sol[3])
            for j in 1:length(setdiff(1:length(sol[2]),sol[3]))
                neighbors[i,j] = swap(sol,sol[3][i],j,b,s)
            end
            
    end

    return vec(neighbors)

end    


        



