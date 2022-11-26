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
    tabu_list = Vector{Vector{Vector{Int64}}}(undef,tenure)

    while cpt_pas_amelioration <= k*length(sol[2])
        
        neighborhood = getNeighbors_shift(best_candidate,Q,b,s)
        best_candidate = neighborhood[1]

        for neighbor in neighborhood
            if !(neighbor in tabu_list) && gain_shift(obj,best_candidate,)
    end
end

function getNeighbors_shift(sol::Vector{Vector{Int64}},Q::Int64,
                            b::Array{Float64,2},s::Vector{Float64})::Array{Vector{Vector{Int64}},2}
    
    neighbors = Array{Vector{Vector{Int64}},2}(undef,(length(sol[1]),length(sol[2])))

    for i in 1:length(sol[1])
        for j in 1:length(sol[2])
            neighbors[i,j] = shift(sol,Q,i,sol[1][i],j,b,s)
        end
    end
    
    return neighbors

end    


        



