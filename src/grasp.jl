"""
    grasp()

- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-30

# Arguments

I : terminaux
J : CLVL1 = concentrateurs de niveau 1
K : CLVL2 = concentrateurs de niveau 2
Q : capacité max des concentrateurs
b,c,s : arrays de coûts (vérifier si les couts sont float ou int)
d : arrays de distances entre terminaux et CLVL1
a : alpha de grasp
P : nb de solutions à générer, supposé pair


"""
#=
solution : vecteur composé de :

[1] : vecteur de nI éléments, valeurs = concentrateurs de lvl 1 assignés
[2] : vecteur de nJ éléments, valeurs = concentrateurs de lvl 2 assignés

[3] concentrateurs de lvl 1 ouverts
[4] concentrateurs de lvl 2 ouverts
=#

#=

    Entrée : données du pb (taille, couts, nb de solutions)

    -> générer P/2 premières solutions :
        choisir un CLVL1 au hasard
        l'affecter aux Q terminaux les plus proches
        tant que les terminaux ne sont pas tous affectés:
            -> calculer la restricted list des CLVL1 restants (via g1)
                -> on calcule tous les 1/f1 pour les CLVL1 restants
                -> on prend ceux > threshold
            -> en choisir un au hasard parmi la RL
            -> l'affecter aux Q plus proches terminaux non affectés, et le supprimer de la CL

        tant que les CLVL1 ouverts ne sont pas tous affectés:
            on ouvre Q' fois moins de CLVL 2 que de CLVL 1
             -> calculer la restricted list des CLVL2 restants (via g1)
                -> on calcule tous les 1/f1 pour les CLVL2 restants
                -> on prend ceux > threshold
            -> en choisir un au hasard parmi la RL
            -> l'affecter aux Q' plus proches CLVL1 non affectés, et le supprimer de la CL

    choix des CLVL2 dépendant que de l'objectif 1 => ensemble des solutions biaisé
    pondérer les objectifs ?

    Sortie : P solutions (P/2 guidées par la fonction de coûts, P/2 par pull

=#

using Random
include("nearest_neighbors.jl")
include("utilities.jl")

function grasp(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Int64,2},
    c::Array{Int64,2},s::Array{Int64,2},d::Array{Float64,2},a::Float64,P::Int64)::Vector{Vector{Int64},3}

    #construit une solution
    function genSol(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Int64,2},
    c::Array{Int64,2},s::Array{Int64,2},d::Array{Float64,2},a::Float64,g::Function)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)),zeros(Int64,length(J)),Int64[],Int64[]]
        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)
        #on le supprime de la CL
        deleteat!(CL_CLVL1, find(x->x==random_first,CL_CLVL1))
        #terminaux eligibles
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]
        # assignation des terminaux libres les plus proches
        if length(free_terminals) <= Q
            for i in free_terminals
                sol[1][i] = random_first
            end
        else
            for i in nearest_neighbors(d,random_first,Q,free_terminals)
                sol[1][i] = random_first
            end
        end
        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g1, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) > 1
            gmin = minimum([g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])
            gmax = maximum([g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])
            RL = filter(x->g1(x,c,b,s,d,free_terminals,Q,sol[3],sol[4]) > (gmax - a*(gmax-gmin)),CL_CLVL1)
            chosen_one = rand(RL)
            deleteat!(CL_CLVL1, find(x->x==chosen_one,CL_CLVL1))
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]
            if length(free_terminals) <= Q
                for i in free_terminals
                    sol[1][i] = chosen_one
                end
            else
                for i in nearest_neighbors(d,chosen_one,Q,free_terminals)
                    sol[1][i] = chosen_one
                end
            end
            #on ouvre le CLVL1 choisi
            push!(sol[3],chosen_one)

            #on l'affecte au clvl2 le moins couteux
            sol[2][chosen_one] = argmin(s[chosen_one])
            #et on ouvre si ce n'est pas déjà le cas ce clvl2
            if !(argmin(s[argmin(b[chosen_one,:])]) in sol[4])
                push!(sol[4],argmin(s[argmin(b[chosen_one,:])]))
            end
        end

        return sol
    end



end