

#Shift : Réassigne un terminal à un autre concentrateur (ouvert ou non)

function shift(solInit::Solution,Q::Int64, terminal::Int64,conc_depart::Int64,conc_arrivee::Int64,
                b::Array{Float64,2},s::Vector{Float64})::Solution

    

    #si le conc_arrivee est déjà saturé on ne fait rien
    sol = deepcopy(solInit)
    if length(findall(x -> x == conc_arrivee,sol.assign_term)) == Q
        return sol
    end    

    #si c était le dernier terminal affecté à ce conc_depart, on le ferme
    
    if length(findall(x -> x == conc_depart,sol.assign_term)) == 1
        deleteat!(sol.conclvl1_ouverts, findfirst(x->x==conc_depart,sol.conclvl1_ouverts))

        #si c était le dernier cclvl1 affecté à un cclvl2 on ferme celui-ci

        cclvl2 = sol.assign_conclvl1[conc_depart]
        if length(findall(x -> x == cclvl2,sol.assign_conclvl1)) == 1
            deleteat!(sol.conclvl2_ouverts, findfirst(x->x==cclvl2,sol.conclvl2_ouverts))
        end    

        sol.assign_conclvl1[conc_depart] = 0
    end
    
    #on affecte le nouveau conc au terminal

    sol.assign_term[terminal] = conc_arrivee

    #si le conc_arrivee n'était pas ouvert on l'ouvre

    if !(conc_arrivee in sol.conclvl1_ouverts)
        push!(sol.conclvl1_ouverts,conc_arrivee)
        #et on l'affecte au clvl2 le moins cher
        clvl2_moins_cher = argmin([b[conc_arrivee,k] + k in sol.conclvl2_ouverts ? 0 : s[k] for k in sol.conclvl2_ouverts])
        sol.assign_conclvl1[conc_arrivee] = clvl2_moins_cher
            
        #et on ouvre si ce n'est pas déjà le cas ce clvl2
        if !(clvl2_moins_cher in sol.conclvl2_ouverts)
            push!(sol.conclvl2_ouverts,clvl2_moins_cher)
        end
    end

    return sol

end 

#=
Swap: Interchange les terminaux assignés à 2 concentrateurs différents
 -> on peut swaper avec un CLVL1 fermé:
   - ouvrir le CLVL1 d'arrivée
   - transférer les terminaux affectés
   - fermer le CLVL1 de départ

=#

function swap(solInit::Solution,conc_depart::Int64,
    conc_arrivee::Int64,b::Array{Int64,2},s::Vector{Int64})::Solution
    
    sol = deepcopy(solInit)

    #si le conc départ n'est pas ouvert on inverse les concentrateurs

   
    if !(conc_depart in sol.conclvl1_ouverts)
        conc_depart, conc_arrivee = conc_arrivee, conc_depart
    end    

    #si conc_arrivee est déjà ouvert on ne fait rien

    if conc_arrivee in sol.conclvl1_ouverts || ( !(conc_depart in sol.conclvl1_ouverts) && !(conc_arrivee in sol.conclvl1_ouverts))
        return sol
    end

    #on ferme conc_depart


    deleteat!(sol.conclvl1_ouverts, findfirst(x->x==conc_depart,sol.conclvl1_ouverts))

    cclvl2 = sol.assign_conclvl1[conc_depart]

    if length(findall(x -> x == cclvl2,sol.assign_conclvl1)) == 1
        deleteat!(sol.conclvl2_ouverts, findfirst(x->x==cclvl2,sol.conclvl2_ouverts))
    end    

    sol.assign_conclvl1[conc_depart] = 0

    #on reaffecte les terminaux

    for i in 1:length(sol.assign_term)
        if sol.assign_term[i] == conc_depart
            sol.assign_term[i] = conc_arrivee
        end
    end
    
    # on ouvre conc_arrivee

    push!(sol.conclvl1_ouverts,conc_arrivee)

    #et on l'affecte au clvl2 le moins cher
    clvl2_moins_cher = argmin([b[conc_arrivee,k] + k in sol.conclvl2_ouverts ? 0 : s[k] for k in eachindex(s)])
    sol.assign_conclvl1[conc_arrivee] = clvl2_moins_cher
                
    #et on ouvre si ce n'est pas déjà le cas ce clvl2
    if !(clvl2_moins_cher in sol.conclvl2_ouverts)
        push!(sol.conclvl2_ouverts,clvl2_moins_cher)
    end

    return sol

end

function swap2(solInit::Solution,conc_in::Int64,conc_out::Int64)::Solution

    sol = deepcopy(solInit)
    if !(conc_in in sol.conclvl2_ouverts)
    push!(sol.conclvl2_ouverts,conc_in)
    end
    deleteat!(sol.conclvl2_ouverts, findfirst(x -> x==conc_out,sol.conclvl2_ouverts))

    for k in eachindex(sol.assign_conclvl1)
        if sol.assign_conclvl1[k] == conc_out
            sol.assign_conclvl1[k] = conc_in
        end
    end
    
    return sol

    
end