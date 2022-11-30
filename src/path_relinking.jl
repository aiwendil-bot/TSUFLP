

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

include("tabu_movs.jl")

function path_relinking!(solInit::Vector{Vector{Int64}},solFin::Vector{Vector{Int64}},skiplist,Q::Int64,
                        c::Array{Float64,2},b::Array{Float64,2},d::Array{Float64,2},s::Vector{Float64})

    sols_intermediaires = Vector{Vector{Int64}}[]                        
    sol_transit = deepcopy(solInit)
    
    # réconcilier les conc lvl1

    conc1_in = setdiff(solFin[3],solInit[3])
    conc1_out = setdiff(solInit[3],solFin[3])

    while length(conc1_in) > 0 && length(conc1_out) > 0
        
        moves_dominated = []
        moves_nondominated = []

        for i in eachindex(conc1_in)
            for j in eachindex(conc1_out)
                candidate = swap(sol_transit,conc1_in[i],conc1_out[j],b,s)
                #=
                if inserer(skiplist,candidate) #si efficace
                    push!(moves_nondominated,[i,j])
                else
                    push!(moves_dominated,[i,j])
                end
                =#
                push!(moves_nondominated,[i,j])
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
            
            sol_transit[2][conc_in[i]] = solFin[2][conc_in[i]] 

            for k in eachindex(sol_transit[1])
                if solFin[1][k] == conc1_in[i]
                    sol_transit[1][k] = conc1_in[i]
                end
            end
            
            #inserer(skiplist,sol_transit)
            push!(sols_intermediaires,sol_transit)

        end

    end              

    if length(conc1_out) > 0

        for j in eachindex(conc1_out)
            deleteat!(sol_transit[3],findfirst(x -> x==conc1_out[j]))
            sol_transit[2][conc1_out[j]] = 0

            for k in eachindex(sol_transit[1])
                if sol_transit[1][k] == conc1_out[j]
                    sol_transit[1][k] == solFin[1][k]
                end
            end
            
            if isFeasible(sol_transit,Q)
                #inserer(skiplist,sol_transit)
                push!(sols_intermediaires,sol_transit)

            end
            
        end
    end

    # réconcilier les cclvl2

    conc2_in = setdiff(solFin[4],solInit[4])
    conc2_out = setdiff(solInit[4],solFin[4])

    while length(conc2_in) > 0 && length(conc2_out) > 0

        moves_dominated = []
        moves_nondominated = []

        for i in eachindex(conc2_in)
            for j in eachindex(conc2_out)
                candidate = swap2(sol_transit,conc2_in[i],conc2_out[j])
                push!(moves_nondominated,[i,j])
                #=
                if feasible
                if inserer(skiplist,candidate) #si efficace
                    push!(moves_nondominated,[i,j])
                else
                    push!(moves_dominated,[i,j])
                end
                =#
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

            if !(conc2_in[i] in sol_transit[4])
                push!(sol_transit[4], conc2_in[i])
            end
            
            for k in eachindex(sol_transit[2])
                
                if solFin[2][k] == conc2_in[i]
                    sol_transit[2][k] = conc2_in[i]
                end

            end    
                #inserer(skiplist,sol_transit)
                push!(sols_intermediaires,sol_transit)

        end
            
    end

    if length(conc2_out) > 0
        for j in eachindex(conc2_out)
            deleteat!(sol_transit[4],findfirst(x -> x==conc2_out[j],sol_transit[4]))

            for k in eachindex(sol_transit[2])
                if sol_transit[2][k] == conc2_out[j]
                    sol_transit[2][k] = solFin[2][k]
                end
            end
            
            #inserer(skiplist,sol_transit)
            push!(sols_intermediaires,sol_transit)


        end
    end
    #=
    on change un à un chaque terminal
        pour chaque changement on vérifie que le conc d'arrivée n'est pas saturé
        si c est le cas on regarde si la solution est dominée
        
    =#
    
    wrong_terminals = [k for k in eachindex(sol_transit[1]) if sol_transit[1][k] != solFin[1][k]]

    while length(wrong_terminals) > 0

        i = rand(wrong_terminals)

        sol_transit[1][i] = solFin[1][i]

        #relâcher contrainte saturation ?

        #if length(findall(x -> x == solFin[1][i],solFin[1])) <= Q
            #inserer(skiplist, sol_transit)
        push!(sols_intermediaires,sol_transit)

        #end
        
        deleteat!(wrong_terminals, findfirst(x -> x == i, wrong_terminals))
    end

    wrong_conc = [k for k in eachindex(sol_transit[2]) if sol_transit[2][k] != solFin[2][k]]


    while length(wrong_conc) > 0

        i = rand(wrong_conc)

        sol_transit[2][i] = solFin[2][i]

        #inserer(skiplist, sol_transit)
        push!(sols_intermediaires,sol_transit)

        
        deleteat!(wrong_conc, findfirst(x -> x == i, wrong_conc))
    end

    return sols_intermediaires
end    

