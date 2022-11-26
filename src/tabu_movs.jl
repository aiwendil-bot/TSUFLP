

#Shift : Réassigne un terminal à un autre concentrateur (ouvert ou non)

function shift(solInit::Vector{Vector{Int64}},Q::Int64, terminal::Int64,conc_depart::Int64,conc_arrivee::Int64,
                b::Array{Float64,2},s::Vector{Float64})

    #si le conc_arrivee est déjà saturé on ne fait rien
    sol = deepcopy(solInit)
    if length(findall(x -> x == conc_arrivee,sol[1])) == Q
        return sol
    end    

    #si c était le dernier terminal affecté à ce conc_depart, on le ferme

    if length(findall(x -> x == conc_depart,sol[1])) == 1
        deleteat!(sol[3], findfirst(x->x==conc_depart,sol[3]))

        #si c était le dernier cclvl1 affecté à un cclvl2 on ferme celui-ci

        cclvl2 = sol[2][conc_depart]
        if length(findall(x -> x == cclvl2,sol[2])) == 1
            deleteat!(sol[4], findfirst(x->x==cclvl2,sol[4]))
        end    

        sol[2][conc_depart] = 0
    end
    
    #on affecte le nouveau conc au terminal

    sol[1][terminal] = conc_arrivee

    #si le conc_arrivee n'était pas ouvert on l'ouvre

    if !(conc_arrivee in sol[3])
        push!(sol[3],conc_arrivee)
        #et on l'affecte au clvl2 le moins cher
        clvl2_moins_cher = argmin([b[conc_arrivee,k] + k in sol[4] ? 0 : s[k] for k in sol[4]])
        sol[2][conc_arrivee] = clvl2_moins_cher
            
        #et on ouvre si ce n'est pas déjà le cas ce clvl2
        if !(clvl2_moins_cher in sol[4])
            push!(sol[4],clvl2_moins_cher)
        end
    end

    return sol

end 

#=
Swap: Interchange les terminaux assignés à 2 concentrateurs différents
 -> on peut swaper avec un CLVL1 fermé:
   - ouvrir le CLVL1 d'arrivé
   - transférer les terminaux affectés
   - fermer le CLVL1 de départ

=#

function swap!(sol::Vector{Vector{Int64}},conc_depart::Int64,
    conc_arrivee::Int64,b::Array{Float64,2},s::Vector{Float64})
    
    #si conc_arrivee est déjà ouvert on ne fait rien

    if conc_arrivee in sol[3]
        return nothing
    end

    #on ferme conc_depart

    deleteat!(sol[3], findfirst(x->x==conc_depart,sol[3]))

    cclvl2 = sol[2][conc_depart]

    if length(findall(x -> x == cclvl2,sol[2])) == 1
        deleteat!(sol[4], findfirst(x->x==cclvl2,sol[4]))
    end    

    sol[2][conc_depart] = 0

    #on reaffecte les terminaux

    for i in 1:length(sol[1])
        if sol[1][i] == conc_depart
            sol[1][i] = conc_arrivee
        end
    end
    
    # on ouvre conc_arrivee

    push!(sol[3],conc_arrivee)

    #et on l'affecte au clvl2 le moins cher
    clvl2_moins_cher = argmin([b[conc_arrivee,k] + k in sol[4] ? 0 : s[k] for k in sol[4]])
    sol[2][conc_arrivee] = clvl2_moins_cher
                
    #et on ouvre si ce n'est pas déjà le cas ce clvl2
    if !(clvl2_moins_cher in sol[4])
        push!(sol[4],clvl2_moins_cher)
    end

    return nothing

end