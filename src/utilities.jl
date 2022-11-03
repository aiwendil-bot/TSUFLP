#=
utilities:
- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-31
=#

include("nearest_neighbors.jl")

# fonction d'utilité liée à la fonction de coût entre terminaux libres et le CLVL1 candidat
# candidat clvl1
# c : matrice de coûts entre terminaux et clvl1
# b : matrice de coûts entre clvl1 et clvl2
# d matrice de distances
# free_terminals

# calcul le coût d'insertion du clvl1 candidat à la solution, prenant compte de :
# coût d'assignation des Q plus proches terminaux
# coût d'ouverture du clvl1 et de son assignation (au mieux) s'il n'est pas déjà ouvert dans la solution

function g1(candidat::Int64,c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
    d::Vector{Vector{Int64}},free_terminals::Vector{Int64},Q::Int64,clvl1_ouverts::Vector{Int64},
    clvl2_ouverts::Vector{Int64})::Float64

    if length(free_terminals) <= Q
        return 1/(sum([c[candidat][k] for k in free_terminals])
               + minimum(b[candidat,:]) + (argmin(s[argmin(b[candidat,:])]) in clvl2_ouverts ? 0 : minimum(s[argmin(b[candidat,:])])))
    else
        nearest = nearest_neighbors(candidat,d,Q)
        return 1/(sum([b[candidat][k] for k in nearest])
        + minimum(b[candidat,:]) + (argmin(s[argmin(b[candidat,:])]) in clvl2_ouverts ? 0 : minimum(s[argmin(b[candidat,:])])))
    end
end

# calcul le coût d'insertion du clvl1 candidat à la solution, prenant compte de :
# distance max entre les dépôts ouverts et les terminaux
# coût d'ouverture du clvl1 et de son assignation (au mieux) s'il n'est pas déjà ouvert dans la solution

function g2(candidat::Int64,c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
    d::Vector{Vector{Int64}},free_terminals::Vector{Int64},Q::Int64,clvl1_ouverts::Vector{Int64},
    clvl2_ouverts::Vector{Int64})::Float64

    min_ctri = zeros(Int64,length(clvl1_ouverts)+1)
    for i in 1:length(vcat(candidat,clvl1_ouverts))
        min_ctri[i] = minimum(d[vcat(candidat,clvl1_ouverts)[i]][j] for j in free_terminals)
    end
    return 1/maximum(min_ctri)
    #return maximum([d[i][j] for i in vcat(candidat,clvl1_ouverts) for j in free_terminals ])

end