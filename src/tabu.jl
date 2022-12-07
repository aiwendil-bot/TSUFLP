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

function tabu(obj::Int64, sol::Solution,Q::Int64,tenure::Int64,k::Float64,
                c::Array{Int64,2},b::Array{Int64,2},s::Vector{Int64},
                d::Array{Int64,2})::Solution

    cpt_pas_amelioration::Int64 = 0
    nb_iterations = 1
    best_solution = deepcopy(sol)
    best_candidate = deepcopy(sol)
    tabu_list = Int64[]
    while cpt_pas_amelioration < k*length(sol.assign_conclvl1)

        found_amelio = false

        #performs the first non-tabu move that leads to an improvement of the current solution
        # or a tabu move that reaches the aspiration move (< best solution)

        # VOISINAGE = SWAP
        
        i = 1
        j = 1

        conc_depart = 1
        conc_arrivée = 1

        
        while !found_amelio && i <= length(best_candidate.conclvl1_ouverts)
            
            while !found_amelio && j <= length(setdiff(1:length(best_candidate.assign_conclvl1),best_candidate.conclvl1_ouverts))
                conc_depart = best_candidate.conclvl1_ouverts[i]
                conc_arrivée = setdiff(1:length(best_candidate.assign_conclvl1),best_candidate.conclvl1_ouverts)[j]

                neighbor = swap(best_candidate,conc_depart,conc_arrivée,b,s)

            
                if !(conc_arrivée in tabu_list)
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
              

        if !found_amelio

            if !isempty(tabu_list)
            
                moves_potentiels = [[i,j] for i in best_candidate.conclvl1_ouverts for j in tabu_list] 
                least_worst_move = moves_potentiels[findmin(move -> 
                                    evaluate_solution(obj,swap(best_candidate,move[1],move[2],b,s),d,c,b,s),moves_potentiels)[2]]
                #best_candidate = tabu_list[findmin(x -> evaluate_solution(obj,x,d,c,b,s),tabu_list)[2]]
                best_candidate =swap(best_candidate,least_worst_move[1],least_worst_move[2],b,s)
            
                if length(tabu_list) < tenure
                    push!(tabu_list,least_worst_move[1])
                else    
                    tabu_list[(nb_iterations+1)%tenure + 1] = least_worst_move[1]
            
                end
                
            else

                moves_potentiels = [[i,j] for i in best_candidate.conclvl1_ouverts for j in setdiff(1:length(sol.assign_conclvl1),best_candidate.conclvl1_ouverts)]
                least_worst_move = moves_potentiels[findmin(move -> 
                                    evaluate_solution(obj,swap(best_candidate,move[1],move[2],b,s),d,c,b,s),moves_potentiels)[2]]
                best_candidate =swap(best_candidate,least_worst_move[1],least_worst_move[2],b,s)
                push!(tabu_list,least_worst_move[1])
                

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
                push!(tabu_list,conc_depart)
            else    
                tabu_list[(nb_iterations+1)%tenure + 1] = conc_depart
        
            end
        end    
        #tabu_list[(nb_iterations+1)%tenure + 1] = best_candidate
        nb_iterations += 1
 
    end

    return best_solution
end
#=
function getNeighbors_shift(sol::Vector{Vector{Int64}},Q::Int64,
                            b::Array{Float64,2},s::Vector{Float64})::Vector{Vector{Vector{Int64}}}
    
    conc_possibles = [j for j in eachindex(sol.assign_conclvl1) if length(findall(x -> x == j,sol.assign_term)) != Q ]

    neighbors = Array{Vector{Vector{Int64}},2}(undef,(length(sol.assign_term),length(conc_possibles)))

    for i in 1:length(sol.assign_term)
        for j in eachindex(conc_possibles)
            neighbors[i,j] = shift(sol,Q,i,sol.assign_term[i],conc_possibles[j],b,s)
        end
    end
    
    return vec(neighbors)

end    


function getNeighbors_swap(sol::Vector{Vector{Int64}},
    b::Array{Float64,2},s::Vector{Float64})::Vector{Vector{Vector{Int64}}}

    neighbors = Array{Vector{Vector{Int64}},2}(undef,(length(sol.conclvl1_ouverts),length(setdiff(1:length(sol.assign_conclvl1),sol.conclvl1_ouverts))))

    for i in 1:length(sol.conclvl1_ouverts)
            for j in 1:length(setdiff(1:length(sol.assign_conclvl1),sol.conclvl1_ouverts))
                neighbors[i,j] = swap(sol,sol.conclvl1_ouverts[i],j,b,s)
            end
            
    end

    return vec(neighbors)

end    
=#
        



