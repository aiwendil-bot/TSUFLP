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

    Sortie : P solutions (P/2 guidées par la fonction de coûts, P/2 par pull

=#

using Random
include("nearest_neighbors.jl")
include("utilities.jl")

function grasp(I::Vector{Int64}, J::Vector{Int64}, K::Vector{Int64}, Q::Int64, b::Array{Float64,2},
    c::Array{Float64,2}, s::Vector{Float64}, d::Array{Float64,2}, a::Float64, P::Int64)::Vector{Vector{Int64}}

    #construit une solution de type 1 (objectif minimiser couts)
    function genSol1(I::Vector{Int64}, J::Vector{Int64}, K::Vector{Int64}, Q::Int64, b::Array{Float64,2},
    c::Array{Float64,2}, s::Vector{Float64}, d::Array{Float64,2}, a::Float64)::Vector{Vector{Int64}}

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
        while length(filter(x->x==0,sol[1])) >= 1
            gmin = minimum([g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])
            gmax = maximum([g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])
            RL = filter(x->g1(x,c,b,s,d,free_terminals,Q,sol[3],sol[4]) > (gmax - a*(gmax-gmin)),CL_CLVL1)
            println(sol[1])
            chosen_one = rand(RL)
            deleteat!(CL_CLVL1, findfirst(x->x==chosen_one,CL_CLVL1))
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]
            println(free_terminals)
            if length(free_terminals) <= Q
                for i in free_terminals
                    sol[1][i] = chosen_one
                end
            else
                println(nearest_neighbors(d,chosen_one,Q,free_terminals))
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

    return genSol1(I, J, K, Q, b, c, s, d, a)

    #construit une solution de type 2 (objectif minimiser distance max des concentrateurs niv1 et minimiser couts concentrateurs niv 2)
    function genSol2(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Float64,2},
    c::Array{Float64,2},s::Vector{Float64},d::Array{Float64,2},a::Float64)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)),zeros(Int64,length(J)),Int64[],Int64[]]
        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)

        #on ne doit pas supprimer de la CL ?
        #deleteat!(CL_CLVL1, find(x->x==random_first,CL_CLVL1))


        #on lui affecte son terminal le plus éloigné
        farthest = argmax(d[chosen_one][j] for j in sol[1])
        sol[1][farthest] = chosen_one
        cpt_affectation_ctr = zeros(Int64,length(J))
        cpt_affectation_ctr[chosen_one] = 1

        #terminaux eligibles
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

        # on ajoute le candidat a la liste des concetrateurs ouverts
        push!(sol[3],random_first)

        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) >= 1
            gmin = minimum([g2(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])
            gmax = maximum([g2(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])

            # filtre en fonction du threshold
            RL = filter(x->g2(x,c,b,s,d,free_terminals,Q,sol[3],sol[4]) > (gmax - a*(gmax-gmin)),CL_CLVL1)
            chosen_one = rand(RL)
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

            #on lui affecte son terminal le plus éloigné
            farthest = argmax(d[chosen_one][j] for j in free_terminals) #correction (avant in sol[1])
            sol[1][farthest] = chosen_one
            cpt_affectation_ctr[chosen_one] += 1

            #si le concentrateur n'est pas encore ouvert,
            if !(chosen_one in sol[3])
                #on ouvre le CLVL1 choisi
                push!(sol[3],chosen_one)
            end

            #si le chosen one est saturé il ne peut plus etre choisi
            if cpt_affectation_ctr[chosen_one] == Q
                deleteat!(CL_CLVL1, find(x->x==chosen_one,CL_CLVL1))
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

    #construit une solution de type 3 (0.5-0.5 entre
    #minimiser distance max des concentrateurs niv1 et les coûts puis minimiser couts concentrateurs niv 2)
    function genSol_compromis(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Float64,2},
    c::Array{Float64,2},s::Vector{Float64},d::Array{Float64,2},a::Float64, λ::Float64)::Vector{Vector{Int64}}

        #Initialisation
        #terminaux = 0 => non affectés
        #CLVL1 = 0 => non ouvert

        sol = [zeros(Int64,length(I)),zeros(Int64,length(J)),Int64[],Int64[]]
        CL_CLVL1 = deepcopy(J) #on initialise la candidate list à J
        #choisir un CLVL1 au hasard
        random_first = rand(J)

        #on ne doit pas supprimer de la CL ?
        #deleteat!(CL_CLVL1, find(x->x==random_first,CL_CLVL1))

        #on lui affecte son terminal avec le meilleur compromis

        terminal_compromis = argmin([cout_affectation_term(terminal,random_first,0.5,c,
        d,sol[3], Int64[]) for terminal in 1:length(sol[1])])

        sol[1][terminal_compromis] = chosen_one
        cpt_affectation_ctr = zeros(Int64,length(J))
        cpt_affectation_ctr[chosen_one] = 1

        #terminaux eligibles
        free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

        # on ajoute le candidat a la liste des concetrateurs ouverts
        push!(sol[3],random_first)

        #=tant qu'il reste des terminaux à assigner :
        on calcule la RCL selon g, on choisit un CLVL1 au hasard dedans
        on affecte les Q terminaux libres les plus proches à ce chosen one
        et on affecte le clvl2 le moins couteux au chosen one

        =#
        while length(filter(x->x==0,sol[1])) >= 1
            gmin = minimum([λ* g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4])+
            (1-λ)*g2(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])

            gmax = maximum([λ* g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4])+
            (1-λ)*g2(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) for k in CL_CLVL1])

            # filtre en fonction du threshold
            RL = filter(x->λ* g1(k,c,b,s,d,free_terminals,Q,sol[3],sol[4])+
            (1-λ)*g2(k,c,b,s,d,free_terminals,Q,sol[3],sol[4]) > (gmax - a*(gmax-gmin)),CL_CLVL1)
            chosen_one = rand(RL)
            free_terminals = [i for i in 1:length(sol[1]) if sol[1][i] == 0]

            #on lui affecte le terminal compromis
            compromis = argmin([cout_affectation_term(terminal,chosen_one,λ,c,
            d,sol[3], Int64[]) for terminal in free_terminals])

            sol[1][compromis] = chosen_one
            cpt_affectation_ctr[chosen_one] += 1

            #si le concentrateur n'est pas encore ouvert,
            if !(chosen_one in sol[3])
                #on ouvre le CLVL1 choisi
                push!(sol[3],chosen_one)
            end

            #si le chosen one est saturé il ne peut plus etre choisi
            if cpt_affectation_ctr[chosen_one] == Q
                deleteat!(CL_CLVL1, find(x->x==chosen_one,CL_CLVL1))
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

end
