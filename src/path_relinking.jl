

#=

#path relinking autorisant des solutions intermédiaires non réalisables

-> evalue les solutions entre solInit et solFin

    -> tant qu'il y a des conc lvl 1 différents, en ajouter un et en retirer un (swap en gros)
        -> si pas même nombre d'ouverts :
            -> si un ouvert de trop dans solInit, le fermer et affecter les terminaux comme dans la solFin
            -> si un ouvert de trop dans solFin : l'ouvrir et affecter les terminaux comme dans la solFin
            
            
    -> tant qu'il y a des conc lvl 2 différents, en ajouter un et en retirer un        
        -> de même qu'avant

    -> tant qu'il y a des terminaux mal attribués, les swap (encore pareil)

-> pour chaque move, vérifier si la solution obtenue est non dominée

    -> celles qui sont non dominées sont stockées dans la skip list

-> s'il y a au moins un move qui conduit à une solution non dominée, en choisir un pour avancer

-> sinon, en choisir un au hasard

=#


function path_relinking!(solInit::Solution,solFin::Solution,skiplist::SkipList,Q::Int64,
                        c::Array{Int64,2},b::Array{Int64,2},d::Array{Int64,2},s::Vector{Int64})::Vector{Solution}

    sols_intermediaires = Solution[]                        
    sol_transit = deepcopy(solInit)
    
    # réconcilier les conc lvl1

    conc1_in = setdiff(solFin.conclvl1_ouverts,solInit.conclvl1_ouverts)
    conc1_out = setdiff(solInit.conclvl1_ouverts,solFin.conclvl1_ouverts)

    while length(conc1_in) > 0 && length(conc1_out) > 0
        
        moves_dominated = []
        moves_nondominated = []

        for i in eachindex(conc1_in)
            for j in eachindex(conc1_out)
                candidate = swap(sol_transit,conc1_in[i],conc1_out[j],b,s)
                evals = evaluate_solution(3,candidate,d,c,b,s)
                elem = Elem(Point(evals[1],evals[2]),candidate)
                if SL_insert!(skiplist,elem,0.5) #si efficace
                    push!(moves_nondominated,[i,j])
                else
                    push!(moves_dominated,[i,j])
                end
            
            end
        end

        move = !(isempty(moves_nondominated)) ? rand(moves_nondominated) : rand(moves_dominated)

        sol_transit = swap(sol_transit,conc1_in[move[1]],conc1_out[move[2]],b,s)
        push!(sols_intermediaires,sol_transit)

        deleteat!(conc1_in,move[1])
        deleteat!(conc1_out,move[2])

    end



    if length(conc1_in) > 0

        for i in eachindex(conc1_in)
            push!(sol_transit[3], conc1_in[i])
            
            sol_transit.assign_conclvl1[conc_in[i]] = solFin.assign_conclvl1[conc_in[i]] 

            for k in eachindex(sol_transit.assign_term)
                if solFin.assign_term[k] == conc1_in[i]
                    sol_transit.assign_term[k] = conc1_in[i]
                end
            end
            
            if isFeasible(sol_transit,Q)
                evals = evaluate_solution(3,sol_transit,d,c,b,s)
                elem = Elem(Point(evals[1],evals[2]),sol_transit)
                SL_insert!(skiplist,elem,0.5)
            end    
            push!(sols_intermediaires,sol_transit)

        end

    end              

    if length(conc1_out) > 0

        for j in eachindex(conc1_out)
            deleteat!(sol_transit[3],findfirst(x -> x==conc1_out[j]))
            sol_transit.assign_conclvl1[conc1_out[j]] = 0

            for k in eachindex(sol_transit[1])
                if sol_transit.assign_term[k] == conc1_out[j]
                    sol_transit.assign_term[k] == solFin.assign_term[k]
                end
            end
            
            if isFeasible(sol_transit,Q)
                #inserer(skiplist,sol_transit)
                push!(sols_intermediaires,sol_transit)

            end
            
        end
    end


    # réconcilier les cclvl2

    conc2_in = setdiff(solFin.conclvl2_ouverts,sol_transit.conclvl2_ouverts)
    conc2_out = setdiff(sol_transit.conclvl2_ouverts,solFin.conclvl2_ouverts)

    while length(conc2_in) > 0 && length(conc2_out) > 0

        moves_dominated = []
        moves_nondominated = []

        for i in eachindex(conc2_in)
            for j in eachindex(conc2_out)

                candidate = swap2(sol_transit,conc2_in[i],conc2_out[j])
                
                
                if isFeasible(candidate, Q)
                    evals = evaluate_solution(3,candidate,d,c,b,s)
                    elem = Elem(Point(evals[1],evals[2]),candidate)
                    if SL_insert!(skiplist,elem,0.5) #si efficace
                        push!(moves_nondominated,[i,j])
                    else
                        push!(moves_dominated,[i,j])
                    end
                end
            end
        end
        
        move = !(isempty(moves_nondominated)) ? rand(moves_nondominated) : rand(moves_dominated)

        sol_transit = swap2(sol_transit,conc2_in[move[1]],conc2_out[move[2]])
        push!(sols_intermediaires,sol_transit)

        deleteat!(conc2_in,move[1])
        deleteat!(conc2_out,move[2])
    end    

    if length(conc2_in) > 0
        for i in eachindex(conc2_in)

            if !(conc2_in[i] in sol_transit.conclvl2_ouverts)
                push!(sol_transit.conclvl2_ouverts, conc2_in[i])
            end
            
            for k in eachindex(sol_transit.assign_conclvl1)
                
                if solFin.assign_conclvl1[k] == conc2_in[i]
                    sol_transit.assign_conclvl1[k] = conc2_in[i]
                end

            end    
            if isFeasible(sol_transit, Q)
                evals = evaluate_solution(3,sol_transit,d,c,b,s)
                elem = Elem(Point(evals[1],evals[2]),sol_transit)
                SL_insert!(skiplist,elem,0.5)
            end
                push!(sols_intermediaires,sol_transit)

        end
            
    end

    if length(conc2_out) > 0
        for j in eachindex(conc2_out)
            deleteat!(sol_transit.conclvl2_ouverts,findfirst(x -> x==conc2_out[j],sol_transit.conclvl2_ouverts))

            for k in eachindex(sol_transit.assign_conclvl1)
                if sol_transit.assign_conclvl1[k] == conc2_out[j]
                    sol_transit.assign_conclvl1[k] = solFin.assign_conclvl1[k]
                end
            end
            
            if isFeasible(sol_transit, Q)
                
                evals = evaluate_solution(3,sol_transit,d,c,b,s)
                elem = Elem(Point(evals[1],evals[2]),sol_transit)
                SL_insert!(skiplist,elem,0.5)

            end
            push!(sols_intermediaires,sol_transit)


        end
    end

    #=
    on change un à un chaque terminal
        pour chaque changement on vérifie que le conc d'arrivée n'est pas saturé
        si c est le cas on regarde si la solution est dominée
        
    =#
    
    wrong_terminals = [k for k in eachindex(sol_transit.assign_term) if sol_transit.assign_term[k] != solFin.assign_term[k]]
    
    while length(wrong_terminals) > 0

        i = rand(wrong_terminals)

        sol_transit.assign_term[i] = solFin.assign_term[i]


        if length(wrong_terminals) % Int(floor(length(sol_transit.assign_term)/2)) == 0    
            push!(sols_intermediaires,sol_transit)
        end

        if isFeasible(sol_transit, Q)
                evals = evaluate_solution(3,sol_transit,d,c,b,s)
                elem = Elem(Point(evals[1],evals[2]),sol_transit)
                SL_insert!(skiplist,elem,0.5)
        end

        #end
        
        
        deleteat!(wrong_terminals, findfirst(x -> x == i, wrong_terminals))
    end

    wrong_conc = [k for k in eachindex(sol_transit.assign_conclvl1) if sol_transit.assign_conclvl1[k] != solFin.assign_conclvl1[k]]


    while length(wrong_conc) > 0

        i = rand(wrong_conc)

        sol_transit.assign_conclvl1[i] = solFin.assign_conclvl1[i]

        if isFeasible(sol_transit, Q)
            evals = evaluate_solution(3,sol_transit,d,c,b,s)
            elem = Elem(Point(evals[1],evals[2]),sol_transit)
            SL_insert!(skiplist,elem,0.5)
        end
        push!(sols_intermediaires,sol_transit)

        
        deleteat!(wrong_conc, findfirst(x -> x == i, wrong_conc))
    end

    return sols_intermediaires
end    

