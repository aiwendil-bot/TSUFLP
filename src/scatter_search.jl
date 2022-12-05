function scatter_search(instance_name::String, c::Array{Int64,2},b::Array{Int64,2},s::Vector{Int64},d::Array{Int64,2},
    Q::Int64, a::Float64, P::Int64, tenure::Int64, k::Float64, β::Int64)

    nondominated_solutions = creer_SL(5) 

    I = [i for i in 1:size(c,1)]
    J = [i for i in 1:size(c,2)]
    K = [i for i in 1:length(s)]
    
    # Création de la population initiale
    println("generation pop de taille $P w/ grasp")

    pop_init = @time grasp(I, J, K, Q, b, c, s, d, a, P)

    for i in eachindex(pop_init)
        if !isFeasible(pop_init[i],Q)
            println(i)
            println("faux")
        end
    end

    plot_grasp(pop_init,instance_name,c,b,s,d)

    # Première amélioration (Tabu Search)

    println("tabu 1st population, k = $k, tenure = $tenure")

    first_improvment = @time vcat([tabu(1,pop_init[i],Q,tenure,k,c,b,s,d) for i in 1:Int(P/2)],
                            [tabu(2,pop_init[i],Q,tenure,k,c,b,s,d) for i in 1+Int(P/2):P])

    plot_tabu(first_improvment,instance_name, c,b,s,d)

    for sol in first_improvment
        if !isFeasible(sol,Q)
            println("faux")
        end
    end

    
    # Création des RefSets de taille β (=6)

    println("calcul refsets, β = $β")
    refSets = @time create_refset(first_improvment,β, b,c, s, d,Q)

    # solutions qui ont déjà été dans la même paire

    paires_interdites = Vector{Solution}[]
    
    # stop criterion : nouvelle solution dans le refset

    stop_criterion::Bool = false

    cpt_iteration = 1

    while !stop_criterion

        #subset generation method

        trial_solutions = Solution[] 

        paires_iteration = [[i,j] for i in refSets[1] for j in refSets[2] 
                            if !([i,j] in paires_interdites) && !([j,i] in paires_interdites)]
        
        #path relinking
        n_paire = 1
        for paire in paires_iteration

            push!(paires_interdites,[paire[1],paire[2]])
            println("path relinking et intensification sur la paire $n_paire / $(length(paires_iteration))...")
            sols_intermediaires = path_relinking!(paire[1],paire[2],nondominated_solutions,Q,c,b,d,s)
        
        # intensification sur chaque solution du path relinking
            #println("intensification sur $(length(sols_intermediaires)) solutions")
            for sol in sols_intermediaires
                tabu1 = tabu(1,sol,Q,tenure,k,c,b,s,d)
                tabu2 = tabu(2,sol,Q,tenure,k,c,b,s,d)
                push!(trial_solutions,tabu1)
                push!(trial_solutions,tabu2)
                if isFeasible(tabu1, Q)
                    evals = evaluate_solution(3,tabu1,d,c,b,s)
                    elem = Elem(Point(evals[1],evals[2]),tabu1)
                    SL_insert!(nondominated_solutions,elem,0.5)
                end

                if isFeasible(tabu2, Q)
                    evals = evaluate_solution(3,tabu2,d,c,b,s)
                    elem = Elem(Point(evals[1],evals[2]),tabu2)
                    SL_insert!(nondominated_solutions,elem,0.5)
                end

            end  

            n_paire += 1
        end      
#=
        for k in sort(refSets[1], by = x -> evaluate_solution(1,x,d,c,b,s))
            println(evaluate_solution(1,k,d,c,b,s))
        end
        
        for k in sort(refSets[2], by = x -> evaluate_solution(2,x,d,c,b,s))
            println(evaluate_solution(2,k,d,c,b,s))
        end
        
        println("-----------")
        =#
        stop_criterion = update_refset!(trial_solutions, refSets, b,c, s, d, Q)
        #=
        for k in sort(refSets[1], by = x -> evaluate_solution(1,x,d,c,b,s))
            println(evaluate_solution(1,k,d,c,b,s))
        end
        
        for k in sort(refSets[2], by = x -> evaluate_solution(2,x,d,c,b,s))
            println(evaluate_solution(2,k,d,c,b,s))
        end
        =#

        new_skiplist = creer_SL(5)

        for elem in get_elems(nondominated_solutions)
            SL_insert!(new_skiplist,elem,0.5)
        end    

        plot_pareto(new_skiplist,instance_name,cpt_iteration)

        cpt_iteration += 1



    end
    new_skiplist = creer_SL(5)

    for elem in get_elems(nondominated_solutions)
        SL_insert!(new_skiplist,elem,0.5)
    end 

    return [elem.point for elem in get_elems(new_skiplist)]

end

