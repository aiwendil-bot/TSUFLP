include("SkipList.jl")


function scatter_search(c::Array{Float64,2},b::Array{Float64,2},s::Vector{Float64},
    Q::Int64, a::Float64, λ::Float64, P::Int64, tenure::Int64, k::Float64, β::Int64)

    nondominated_solutions = SkipList(5) 
    
    # Création de la population initiale
    println("generation pop de taille $P w/ grasp")

    pop_init = @time grasp(I, J, K, Q, b, c, s, distances, a, λ, P)

    # Première amélioration (Tabu Search)

    println("tabu 1st population, k = $k, tenure = $tenure")

    first_improvment = @time vcat([tabu(1,pop_init[i],Q,tenure,k,c,b,s,distances) for i in 1:Int(P/2)],
                            [tabu(2,pop_init[i],Q,tenure,k,c,b,s,distances) for i in 1+Int(P/2):P])

    # Création des RefSets de taille β (=6)

    println("calcul refsets, β = $β")
    refSets = @time create_refset(first_improvment,β, b,c, s, distances)

    # solutions qui ont déjà été dans la même paire

    paires_interdites = Vector{Vector{Vector{Int64}}}[]
    
    # stop criterion : nouvelle solution dans le refset

    stop_criterion::Bool = false

    cpt_iteration = 1

    while !stop_criterion

        #subset generation method

        trial_solutions = Vector{Vector{Int64}}[] 

        paires_iteration = [[i,j] for i in refSets[1] for j in refSets[2] 
                            if !([i,j] in paires_interdites) && !([j,i] in paires_interdites)]
        
        #path relinking

        for paire in paires_iteration

            push!(paires_interdites,[paire[1],paire[2]])

            sols_intermediaires = path_relinking!(paire[1],paire[2],skiplist,Q,c,b,d,s)
        
        # intensification sur chaque solution du path relinking

            for sol in sols_intermediaires
                tabu1 = tabu(1,sol,Q,tenure,k,c,b,s,d)
                tabu2 = tabu(2,sol,Q,tenure,k,c,b,s,d)
                push!(trial_solutions,tabu1)
                push!(trial_solutions,tabu2)
                if isFeasible(tabu1, Q)
                    evals = evaluate_solution(3,tabu1,d,c,b,s)
                    elem = Elem(Point(evals[1],evals[2]),Solution(tabu1[1],tabu1[2],tabu1[3]))
                    SL_insert!(skiplist,elem,0.5)
                end

                if isFeasible(tabu2, Q)
                    evals = evaluate_solution(3,tabu2,d,c,b,s)
                    elem = Elem(Point(evals[1],evals[2]),Solution(tabu2[1],tabu2[2],tabu2[3]))
                    SL_insert!(skiplist,elem,0.5)
                end

            end  
        end      

        stop_criterion = update_refset!(trial_solutions, refSets, β, b,c, s, d,Q)



        cpt_iteration += 1
        

    end

    return nondominated_solutions

end
