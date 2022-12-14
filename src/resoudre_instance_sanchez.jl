include("Solution.jl")
include("SkipList.jl")
include("generation_couts.jl")
include("recup_data.jl")
include("grasp.jl")
include("tabu.jl")
include("refset.jl")
include("path_relinking.jl")
include("isFeasible.jl")
include("scatter_search.jl")
include("parsing_instance.jl")
include("generation_couts_sanchez.jl")
include("plot_sols.jl")
include("crossover.jl")
include("vopt_model.jl")


function resoudre_instance_sanchez(filename::String, Q::Int64=7 , 
                    a::Float64=0.4, P::Int64=300,  k::Float64=0.5 , β::Int64=10, crossover_on::Bool=false, tabu_cross_on::Bool=false)

    nameinstance = String(split(Vector(split(filename, '/'))[end],'.')[1])               

    terms, conclvl1, conclvl2 = parsing_instance(filename)

    # génération des matrices de distance, c , b, s
    
    distances = generation_distance_sanchez(terms, conclvl1)

    c = generation_c_sanchez(length(terms),length(conclvl1),minimum(distances),maximum(distances))

    b = generation_distance_sanchez(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in eachindex(conclvl2)]

    starting_time = time()

    YN, XE = scatter_search(nameinstance,c,b,s,distances,Q,a,P,k,β,crossover_on,tabu_cross_on)

    ending_time = round(time()- starting_time, digits=4)

    if crossover_on

        if tabu_cross_on

            writedlm("out/$nameinstance/YN_scatter_tabu_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))
            writedlm("out/$nameinstance/XE_scatter_tabu_cross_$k.txt",XE)
            


        else

            writedlm("out/$nameinstance/YN_scatter_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))
            writedlm("out/$nameinstance/XE_scatter_cross_$k.txt",XE)

        end

    else    

        writedlm("out/$nameinstance/YN_scatter_$k.txt",vcat([[ending_time,length(YN)]],YN))
        writedlm("out/$nameinstance/XE_scatter_$k.txt",XE)

    end

    println("fichiers de résultats (plots, Y_N, X_E) écrits dans out/$nameinstance")

end
