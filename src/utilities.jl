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

function g1(candidat::Int64, c::Array{Int64,2}, b::Array{Int64,2}, s::Vector{Int64},
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
    
    if (cout_affectation
        + minimum(b[candidat, :]) + (argmin(s[argmin(b[candidat, :])]) in clvl2_ouverts ? 0 : minimum(s[argmin(b[candidat, :])]))) == 0
        return 1/0.000001

    end
    return 1 / (cout_affectation
               + minimum(b[candidat, :]) + (argmin(s[argmin(b[candidat, :])]) in clvl2_ouverts ? 0 : minimum(s[argmin(b[candidat, :])])))

end

# calcul le coût d'insertion du clvl1 candidat à la solution, prenant compte de :
# distance max entre les dépôts ouverts et les terminaux
# coût d'ouverture du clvl1 et de son assignation (au mieux) s'il n'est pas déjà ouvert dans la solution

function g2(candidat::Int64, d::Matrix{Int64}, free_terminals::Vector{Int64}, preferences::Array{Int64,2})::Float64

    for i in 1:size(d,1)
           
            return d[preferences[i,candidat],candidat] == 0 ? 1/0.000001 : 1/d[preferences[i,candidat],candidat]    

    end    
end

#calcul le coût d'affectation du terminal candid au clvl1 conc au regard des objectifs
# avec un jeu de poids λ

function cout_affectation_term(candid::Int64,conc::Int64,λ::Float64,c::Array{Int64,2},
    d::Matrix{Int64},clvl1_pris::Vector{Int64}, terminaux_pris::Vector{Int64})::Float64

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

function affect_rapide!(candidat::Int64,matrice::Array{Int64,2},free_terminals::Vector{Int64},Q::Int64,sol::Solution)

        #compteurs pour contrôler l'affectation des terminaux
        cpt_affectation::Int64 = 0
        k = 1
        #on lui affecte les Q terminaux les plus éloignés
        while cpt_affectation < Q && k <= length(sol.assign_term)
            if matrice[k,candidat] in free_terminals
                sol.assign_term[matrice[k,candidat]] = candidat
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


function evaluate_solution(obj::Int64,sol::Solution, d::Array{Int64,2}, 
                            c::Array{Int64,2}, b::Array{Int64,2},s::Vector{Int64})

    if obj == 1

    return sum([c[i,sol.assign_term[i]] for i in eachindex(sol.assign_term)]) + 
         sum([b[j,sol.assign_conclvl1[j]] for j in eachindex(sol.assign_conclvl1) if sol.assign_conclvl1[j] != 0])+
         sum([s[k] for k in sol.conclvl2_ouverts])
    end

    if obj == 2

    return maximum([d[i,sol.assign_term[i]] for i in 1:length(sol.assign_term)])
    end
    if obj == 3
        
        return [sum([c[i,sol.assign_term[i]] for i in 1:length(sol.assign_term)]) + 
         sum([b[j,sol.assign_conclvl1[j]] for j in 1:length(sol.assign_conclvl1) if sol.assign_conclvl1[j] != 0])+
         sum([s[k] for k in sol.conclvl2_ouverts]), maximum([d[i,sol.assign_term[i]] for i in 1:length(sol.assign_term)])]
    end
end
#=
function gain_shift(obj::Int64, solInit::Vector{Vector{Int64}},solModif::Vector{Vector{Int64}},
                    terminal::Int64,c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
                    d::Array{Float64,2})

    gain = 0                
    if obj == 1

        clvl1_depart = solInit[1][terminal]
        clvl1_arrivee = solModif[1][terminal]

        clvl2_depart = solInit[2][clvl1_depart]
        clvl2_arrivee = solModif[2][clvl1_arrivee]
        #=
        println(clvl1_depart)
        println(clvl1_arrivee)
        println(clvl2_depart)
        println(clvl2_arrivee)
        =#
        gain += c[terminal,clvl1_arrivee] - c[terminal,clvl1_depart]

        
        if clvl2_depart == 0
            gain += b[clvl1_arrivee,clvl2_arrivee] + s[clvl2_arrivee]

        else    

            gain += b[clvl1_arrivee,clvl2_arrivee] - b[clvl1_depart,clvl2_depart]
            gain += s[clvl2_arrivee] - s[clvl2_depart]
        end
        
    else

        gain += maximum([d[i,solModif[1][i]] for i in 1:length(solModif[1])]) - 
                maximum([d[i,solInit[1][i]] for i in 1:length(solInit[1])])
        
    end    

    return gain

end              

function gain_swap(obj::Int64, solInit::Vector{Vector{Int64}},solModif::Vector{Vector{Int64}},
    clvl1_depart::Int64,clvl1_arrivee::Int64,c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
    d::Array{Float64,2})

    gain = 0               
    clvl2_depart = solInit[2][clvl1_depart]
        clvl2_arrivee = solModif[2][clvl1_arrivee] 
    if obj == 1

        
        #=
        println(clvl1_depart)
        println(clvl1_arrivee)
        println(clvl2_depart)
        println(clvl2_arrivee)
        =#
        gain += sum([c[i,clvl1_arrivee] for i in eachindex(solInit.assign_term) if solModif.assign_term[i] == clvl1_arrivee])
                - sum([c[i,clvl1_depart] for i in eachindex(solInit.assign_term) if solInit.assign_term[i] == clvl1_depart])

                if clvl2_depart == 0
                    gain += b[clvl1_arrivee,clvl2_arrivee] + s[clvl2_arrivee]
        
                else    
        
                    gain += b[clvl1_arrivee,clvl2_arrivee] - b[clvl1_depart,clvl2_depart]
                    gain += s[clvl2_arrivee] - s[clvl2_depart]
                end     
    else


        gain += maximum([d[i,clvl1_arrivee] for i in eachindex(solModif.assign_term) if solModif.assign_term[i] == clvl1_arrivee ]) - 
                maximum([d[i,clvl1_depart] for i in eachindex(solInit.assign_term) if solInit.assign_term[i] == clvl1_depart])
        
    end    

    return gain

end    
=#

function distance_solutions(sol1::Solution,sol2::Solution)::Int64

    return length(setdiff(sol1.conclvl1_ouverts,sol2.conclvl1_ouverts)) + length(setdiff(sol2.conclvl1_ouverts,sol1.conclvl1_ouverts)) +
           length(setdiff(sol1.conclvl2_ouverts,sol2.conclvl2_ouverts)) + length(setdiff(sol2.conclvl2_ouverts,sol1.conclvl2_ouverts))

end