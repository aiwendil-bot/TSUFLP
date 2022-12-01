"""
    grasp()

- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-30

# Arguments

I : terminaux
J : CLVL1 = concentrateurs de niveau 1 (couts = 0.8 à 1.2 fois Q)
K : CLVL2 = concentrateurs de niveau 2
Q : capacité max des concentrateurs (entre 5 et 9, 7 si tous la même)
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

    Sortie : P solutions (P/3 guidées par la fonction de coûts, P/3 par pull, P/3 moitié moitié

=#

using Random
using Dates

include("nearest_neighbors.jl")
include("utilities.jl")

function grasp(I::Vector{Int64}, J::Vector{Int64}, K::Vector{Int64}, Q::Int64, b::Array{Float64,2},
    c::Array{Float64,2}, s::Vector{Float64}, d::Array{Float64,2}, a::Float64, λ::Float64, P::Int64)::Vector{Vector{Vector{Int64}}}

    #construit une array des terminaux triés par distance décroissante de chaque CLVL1

    terminaux_tries_distance = Array{Int64,2}(undef,(length(I),length(J)))
    for j in 1:length(J)
        terminaux_tries_distance[:,j] = sort(I,by= term -> d[term,j],rev=true)
    end

    
    #construit une array des terminaux triés par cout croissant de chaque CLVL1

    terminaux_tries_couts = Array{Int64,2}(undef,(length(I),length(J)))
    for j in 1:length(J)
        terminaux_tries_couts[:,j] = sort(I,by= term -> c[term,j])
    end
    
    #construit une solution de type 1 (objectif minimiser couts)
    function genSol1(I::Vector{Int64}, J::Vector{Int64}, K::Vector{Int64}, Q::Int64, b::Array{Float64,2},
    c::Array{Float64,2}, s::Vector{Float64}, a::Float64)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)), zeros(Int64,length(J)), Int64[], Int64[]]
        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)
        #on le supprime de la CL
        deleteat!(CL_CLVL1, findfirst(x->x == random_first,CL_CLVL1))
        # on ajoute le candidat a la liste des concetrateurs ouverts
        push!(sol[3],random_first)
        #terminaux eligibles
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]
        # assignation des terminaux libres les plus proches
        affect_rapide!(random_first,terminaux_tries_couts,free_terminals,Q,sol)

        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g1, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) >= 1
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

            evals = [g1(k,c,b,s,terminaux_tries_couts,free_terminals,Q,sol[4]) for k in CL_CLVL1]
            gmin = minimum(evals)
            gmax = maximum(evals)
            #RL = filter(x->g1(x,c,b,s,d,free_terminals,Q,sol[3],sol[4]) > (gmax - a*(gmax-gmin)),CL_CLVL1)
            RL_index = filter(i->evals[i] > (gmax - a*(gmax-gmin)),[i for i in 1:length(evals)])
            RL = [CL_CLVL1[i] for i in RL_index]
            chosen_one = rand(RL)
            deleteat!(CL_CLVL1, findfirst(x->x==chosen_one,CL_CLVL1))
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

            affect_rapide!(chosen_one,terminaux_tries_couts,free_terminals,Q,sol)

            #on ouvre le CLVL1 choisi
            push!(sol[3],chosen_one)

            #on l'affecte au clvl2 le moins couteux
            clvl2_moins_cher = argmin([b[chosen_one,k] + k in sol[4] ? 0 : s[k] for k in K])
            sol[2][chosen_one] = clvl2_moins_cher
                
            #et on ouvre si ce n'est pas déjà le cas ce clvl2
            if !(clvl2_moins_cher in sol[4])
                push!(sol[4],clvl2_moins_cher)
            end
        end
        return sol
    end


    #construit une solution de type 2 (objectif minimiser distance max des concentrateurs niv1 et minimiser couts concentrateurs niv 2)
    function genSol2(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Float64,2},
    s::Vector{Float64},d::Array{Float64,2},a::Float64)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)),zeros(Int64,length(J)),Int64[],Int64[]]
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]
        compteur_affectation_conc = zeros(Int64,length(J))


        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)
        push!(sol[3], random_first)
        
        # on ajoute le candidat a la liste des concetrateurs ouverts
        #=
        push!(sol[3],random_first)
        compteur_affectation_conc[random_first] += 1
        sol[1][terminaux_tries_distance[end,random_first]] = random_first
        =#
        deleteat!(CL_CLVL1, findfirst(x->x==random_first,CL_CLVL1))

        affect_rapide!(random_first,terminaux_tries_distance,free_terminals,Q,sol)





        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) >= 1
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

            evals = [g2(k,d,free_terminals,terminaux_tries_distance) for k in CL_CLVL1]
            gmin = minimum(evals)
            gmax = maximum(evals)
            # filtre en fonction du threshold
            RL_index = filter(i->evals[i] > (gmax - a*(gmax-gmin)),[i for i in 1:length(evals)])
            RL = [CL_CLVL1[i] for i in RL_index]
            chosen_one = rand(RL)

            #on ouvre le CLVL1 choisi
            #version "rapide"
            
            push!(sol[3],chosen_one)
            deleteat!(CL_CLVL1, findfirst(x->x==chosen_one,CL_CLVL1))
            
            affect_rapide!(chosen_one,terminaux_tries_distance,free_terminals,Q,sol)
            
            #version "lente"
            #=
            if !(chosen_one in sol[3])
                push!(sol[3], chosen_one)
            end
            
            k = 1

            while k <= length(sol[1])
                if terminaux_tries_distance[k,chosen_one] in free_terminals
                    sol[1][terminaux_tries_distance[k,chosen_one]] = chosen_one
                    compteur_affectation_conc[chosen_one] += 1
                    break
                end
                k += 1
            end
            
            if compteur_affectation_conc[chosen_one] == Q
                deleteat!(CL_CLVL1, findfirst(x->x==chosen_one,CL_CLVL1))
            end
            =#
        end

        # on connecte chaque concentrateur de niveau 1 ouvert a leur concentrateur de niveau 2 le moins couteux
        for ctr_ouvert in sol[3]
            #on l'affecte au clvl2 le moins couteux
            clvl2_moins_cher = argmin([b[ctr_ouvert,k] + k in sol[4] ? 0 : s[k] for k in K])
            sol[2][ctr_ouvert] = clvl2_moins_cher
                
            #et on ouvre si ce n'est pas déjà le cas ce clvl2
            if !(clvl2_moins_cher in sol[4])
                push!(sol[4],clvl2_moins_cher)
            end
        end
        return sol
    end



    #construit une solution de type 3 (0.5-0.5 entre
    #minimiser distance max des concentrateurs niv1 et les coûts puis minimiser couts concentrateurs niv 2)
    function genSol_compromis(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Float64,2},
    c::Array{Float64,2},s::Vector{Float64},d::Array{Float64,2},a::Float64, λ::Float64)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)),zeros(Int64,length(J)),Int64[],Int64[]]
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]


        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)
        
        # on ajoute le candidat a la liste des concentrateurs ouverts

        push!(sol[3],random_first)

        deleteat!(CL_CLVL1, findfirst(x->x==random_first,CL_CLVL1))

        #dans 50% des cas : affecte les Q plus éloignés
        #dans 50% des cas : affecte les Q moins couteux
        if rand() < 0.5
            affect_rapide!(random_first,terminaux_tries_distance,free_terminals,Q,sol)
        else
            affect_rapide!(random_first,terminaux_tries_couts,free_terminals,Q,sol)
        end    

        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) >= 1

            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]


            evals = [λ* g1(k,c,b,s,terminaux_tries_couts,free_terminals,Q,sol[4])+
            (1-λ)*g2(k,d,free_terminals,terminaux_tries_distance) for k in CL_CLVL1]

            gmin = minimum(evals)
            gmax = maximum(evals)

            # filtre en fonction du threshold
            RL_index = filter(i->evals[i] > (gmax - a*(gmax-gmin)),[i for i in 1:length(evals)])
            RL = [CL_CLVL1[i] for i in RL_index]
            chosen_one = rand(RL)

            #on ouvre le CLVL1 choisi
            push!(sol[3],chosen_one)
            deleteat!(CL_CLVL1, findfirst(x->x==chosen_one,CL_CLVL1))


            if rand() < 0.5
                affect_rapide!(chosen_one,terminaux_tries_distance,free_terminals,Q,sol)
            else
                affect_rapide!(chosen_one,terminaux_tries_couts,free_terminals,Q,sol)
            end 

            
        end

        # on connecte chaque concentrateur de niveau 1 ouvert a leur concentrateur de niveau 2 le moins couteux
        for ctr_ouvert in sol[3]
            cout_ctr_ouvert = [b[ctr_ouvert,k] + (s[k] in sol[4] ? 0 : s[k]) for k in K]
            CL_CLVL2 = argmin(cout_ctr_ouvert)
            sol[2][ctr_ouvert] = CL_CLVL2
            if !(CL_CLVL2 in sol[4])
                push!(sol[4],CL_CLVL2)
            end
        end
        return sol
    end

    solutions = Vector{Vector{Vector{Int64}}}(undef, P)
    #=
    for k in 1:Int(P/3)
        solutions[k] = genSol1(I, J, K, Q, b, c, s, a)
    end
    
    for k in Int(P/3)+1:Int(2*P/3)
        solutions[k] = genSol2(I, J, K, Q, b, s, d, a)
    end

    for k in Int(2*P/3)+1:P
        solutions[k] = genSol_compromis(I,J,K,Q,b,c,s,d,a, λ)
    end

    =#

    for k in 1:Int(P/5)
        solutions[k] = genSol1(I, J, K, Q, b, c, s, a)
    end
    
    for k in Int(P/5)+1:Int(2*P/5)
        solutions[k] = genSol_compromis(I,J,K,Q,b,c,s,d,a, 0.75)
    end

    for k in Int(2*P/5)+1:Int(3*P/5)
        solutions[k] = genSol_compromis(I,J,K,Q,b,c,s,d,a, 0.5)
    end

    for k in Int(3*P/5)+1:Int(4*P/5)
        solutions[k] = genSol_compromis(I,J,K,Q,b,c,s,d,a, 0.25)
    end

    for k in Int(4*P/5)+1:P
        solutions[k] = genSol2(I, J, K, Q, b, s, d, a)
    end
    
    return solutions

end
