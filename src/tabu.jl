#=
tabu:
- Julia version: 
- Author: sujin
- Date: 2022-11-03
=#

include("tabu_movs.jl")

#=

tabu selon l'obj en paramètre

mouvement : swap (fermer un cclvl1 et ouvrir un autre, et transférer les terminaux)

critère d'arrêt : k*nb de cclvl1 itérations sans améliorer best_solutions (k=0.5 pour l'instant, à tuner)

critère d'aspiration : < best solution


=#

function tabu(obj::Int64, sol::Vector{Vector{Int64}},Q::Int64,tenure::Int64,k::Float64,
                c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
                d::Array{Float64,2})
    cpt_pas_amelioration::Int64 = 0
    nb_iterations = 1
    best_solution = deepcopy(sol)
    best_candidate = deepcopy(sol)
    tabu_list = []
    while cpt_pas_amelioration < k*length(sol[2])

        found_amelio = false

        #performs the first non-tabu move that leads to an improvement of the current solution
        # or a tabu move that reaches the aspiration move (< best solution)

        # VOISINAGE = SWAP
        
        i = 1
        j = 1

        
        while !found_amelio && i < length(best_candidate[3])
            
            while !found_amelio && j < length(setdiff(1:length(sol[2]),sol[3]))

                neighbor = swap(best_candidate,best_candidate[3][i],setdiff(1:length(best_candidate[2]),best_candidate[3])[j],b,s)

            
                if !([i,j] in tabu_list)
                    #display(neighborhood[i][2])
                    if evaluate_solution(obj,neighbor,d,c,b,s) < evaluate_solution(obj,best_candidate,d,c,b,s)
                    #if gain_swap(obj,best_candidate,neighborhood[i][2],neighborhood[i][1][1][1],neighborhood[i][1][1][2],c,b,s,d) < 0
                        found_amelio = true
                        best_candidate = neighbor
                    end
                else
                    #critère d'aspiration
                    if evaluate_solution(obj,neighbor,d,c,b,s) < evaluate_solution(obj,best_solution,d,c,b,s)
                        best_candidate = neighbor
                        found_amelio = true
                    end    
                end

                j += 1
            end    
            
            i += 1
        end
        
        

        # VOISINAGE = Shift

        #=

        conc_possibles = [j for j in eachindex(sol[2]) if length(findall(x -> x == j,sol[1])) != Q ]

        i = 1
        j = conc_possibles[1]

        while !found_amelio && i < length(best_candidate[1])
            
            while !found_amelio && j < length(conc_possibles)

                neighbor = shift(best_candidate,Q, i,best_candidate[1][i],conc_possibles[j], b,s)

            
                if !([i,j] in tabu_list)
                    #display(neighborhood[i][2])
                    if evaluate_solution(obj,neighbor,d,c,b,s) < evaluate_solution(obj,best_candidate,d,c,b,s)
                    #if gain_swap(obj,best_candidate,neighborhood[i][2],neighborhood[i][1][1][1],neighborhood[i][1][1][2],c,b,s,d) < 0
                        found_amelio = true
                        best_candidate = neighbor
                    end
                else
                    #critère d'aspiration
                    if evaluate_solution(obj,neighbor,d,c,b,s) < evaluate_solution(obj,best_solution,d,c,b,s)
                        best_candidate = neighbor
                        found_amelio = true
                    end    
                end

                j += 1
            end    
            
            i += 1
        end
        =#
        

        if !found_amelio

            if !isempty(tabu_list)
            
            
                least_worst_move = tabu_list[findmin(move -> 
                                    evaluate_solution(obj,swap(best_candidate,move[1],move[2],b,s),d,c,b,s),tabu_list)[2]]
                #best_candidate = tabu_list[findmin(x -> evaluate_solution(obj,x,d,c,b,s),tabu_list)[2]]
                best_candidate =swap(best_candidate,least_worst_move[1],least_worst_move[2],b,s)
            
                if length(tabu_list) < tenure
                    push!(tabu_list,least_worst_move)
                else    
                    tabu_list[(nb_iterations+1)%tenure + 1] = least_worst_move
            
                end

            end 
            cpt_pas_amelioration += 1



        else 
            if evaluate_solution(obj,best_candidate,d,c,b,s) < evaluate_solution(obj,best_solution,d,c,b,s)
                best_solution = best_candidate
                cpt_pas_amelioration = 0
            else
                cpt_pas_amelioration += 1

            end
            if length(tabu_list) < tenure
                push!(tabu_list,[best_candidate[3][i],setdiff(1:length(best_candidate[2]),best_candidate[3])[j]])
            else    
                tabu_list[(nb_iterations+1)%tenure + 1] = [best_candidate[3][i],setdiff(1:length(best_candidate[2]),best_candidate[3])[j]]
        
            end
        end    
        #tabu_list[(nb_iterations+1)%tenure + 1] = best_candidate
        nb_iterations += 1

        
        
        
    end

    return best_solution
end

function getNeighbors_shift(sol::Vector{Vector{Int64}},Q::Int64,
                            b::Array{Float64,2},s::Vector{Float64})::Vector{Vector{Vector{Int64}}}
    
    conc_possibles = [j for j in eachindex(sol[2]) if length(findall(x -> x == j,sol[1])) != Q ]

    neighbors = Array{Vector{Vector{Int64}},2}(undef,(length(sol[1]),length(conc_possibles)))

    for i in 1:length(sol[1])
        for j in eachindex(conc_possibles)
            neighbors[i,j] = shift(sol,Q,i,sol[1][i],conc_possibles[j],b,s)
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


        



