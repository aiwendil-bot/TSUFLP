#=
utilities:
- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-31
=#
using DelimitedFiles
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

function g1(candidat::Int64, c::Array{Float64,2}, b::Array{Float64,2}, s::Vector{Float64},
    terminaux_tries_couts::Array{Int64,2}, free_terminals::Vector{Int64}, Q::Int64,
    clvl2_ouverts::Vector{Int64})::Float64

    cout_affectation::Float64 = 0.0

    cpt = 0
    k = 1

    while cpt < Q && k <= size(c,1)
        if terminaux_tries_couts[k,candidat] in free_terminals
            cout_affectation += c[k,candidat]
            cpt +=1
        end
        k += 1
    end        

    return 1 / (cout_affectation
               + minimum(b[candidat, :]) + (argmin(s[argmin(b[candidat, :])]) in clvl2_ouverts ? 0 : minimum(s[argmin(b[candidat, :])])))

end

# calcul le coût d'insertion du clvl1 candidat à la solution, prenant compte de :
# distance max entre les dépôts ouverts et les terminaux
# coût d'ouverture du clvl1 et de son assignation (au mieux) s'il n'est pas déjà ouvert dans la solution

function g2(candidat::Int64, d::Matrix{Float64}, free_terminals::Vector{Int64}, preferences::Array{Int64,2})::Float64

    for i in 1:size(d,1)
        if preferences[i,candidat] in free_terminals
            return 1/d[preferences[i,candidat],candidat]
        end    

    end    
end

#calcul le coût d'affectation du terminal candid au clvl1 conc au regard des objectifs
# avec un jeu de poids λ

function cout_affectation_term(candid::Int64,conc::Int64,λ::Float64,c::Array{Float64,2},
    d::Matrix{Float64},clvl1_pris::Vector{Int64}, terminaux_pris::Vector{Int64})::Float64

    current_distance = d[candid][conc]
    max_distance = length(terminaux_pris) == 0 ? current_distance : maximum([d[i][j] for i in clvl1_pris for j in terminaux_pris])

    cout_obj_2 = current_distance <= max_distance ? 0 : current_distance - max_distance

    return λ * c[conc][candid] + (1-λ) * cout_obj_2
end

function distance(point1::Vector{Float64},point2::Vector{Float64})::Float64
    lat1,long1 = point1
    lat2,long2 = point2

    x = (long2-long1) * cos((lat1+lat2)/2)
    y = lat2 - lat1
    return sqrt(x^2 + y^2)*1.852*60
end

#affecte les Q "meilleurs" terminaux au candidat selon la matrice de cout "matrice"

function affect_rapide!(candidat::Int64,matrice::Array{Int64,2},free_terminals::Vector{Int64},Q::Int64,sol::Vector{Vector{Int64}})

        #compteurs pour contrôler l'affectation des terminaux
        cpt_affectation::Int64 = 0
        k = 1
        #on lui affecte les Q terminaux les plus éloignés
        while cpt_affectation < Q && k <= length(sol[1])
            if matrice[k,candidat] in free_terminals
                sol[1][matrice[k,candidat]] = candidat
                cpt_affectation += 1
            end
            k += 1
        end    

end

#=
solution : vecteur composé de :

[1] : vecteur de nI éléments, valeurs = concentrateurs de lvl 1 assignés
[2] : vecteur de nJ éléments, valeurs = concentrateurs de lvl 2 assignés

[3] concentrateurs de lvl 1 ouverts
[4] concentrateurs de lvl 2 ouverts
=#

function evaluate_solution(sol::Vector{Vector{Int64}}, d::Array{Float64,2}, 
                            c::Array{Float64}, b::Array{Float64},s::Vector{Float64})

    z1 = sum([c[i,sol[1][i]] for i in 1:length(sol[1])]) + 
         sum([b[j,sol[2][j]] for j in 1:length(sol[2]) if sol[2][j] != 0])+
         sum([s[k] for k in sol[4]])
    
    z2 = maximum([d[i,sol[1][i]] for i in 1:length(sol[1])])

    return [z1, z2]
    
end

function gain_shift(obj::Int64, solInit::Vector{Vector{Int64}},solModif::Vector{Vector{Int64}},
                    terminal::Int64,c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
                    d::Array{Float64,2})

    gain = 0                
    if obj == 1

        clvl1_depart = solInit[1][terminal]
        clvl1_arrivee = solModif[1][terminal]

        clvl2_depart = solInit[2][clvl1_depart]
        clvl2_arrivee = solModif[2][clvl1_arrivee]

        gain += c[terminal,clvl1_arrivee] - c[terminal,clvl1_depart]
        gain += b[clvl1_arrivee,clvl2_arrivee] - b[clvl1_depart,clvl2_depart]
        gain += s[clvl2_arrivee] - s[clvl2_depart]
        
    else

        gain += maximum([d[i,solModif[1][i]] for i in 1:length(solModif[1])]) - 
                maximum([d[i,solInit[1][i]] for i in 1:length(solInit[1])])
        
    end    

    return gain

end                    