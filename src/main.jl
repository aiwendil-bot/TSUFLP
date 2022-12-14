#=
main:
- Julia version: 1.8.0
- Author: sujin
- Date: 2022-11-02
=#
using Random, Plots, LaTeXStrings, StatsBase
Random.seed!(1235)
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
include("resoudre_instance_sanchez.jl")
include("resoudre_instance_via_dataset.jl")
include("plot_sols.jl")
include("crossover.jl")
include("vopt_model.jl")


function main()



    Q::Int64 = 7

    a::Float64 = 0.4

    P = 300

    k = 0.1

    β = 10

    resoudre_instance_sanchez("data/files/small2.txt", Q,a,P,k,β,false,false)

    for k in [0.1,0.5,0.9]
        resoudre_instance_via_dataset("angers", "data/terms_49.txt",[47.40904,47.50549,-0.65266,-0.42212],"data/clvl1.csv","data/clvl2.csv",[100,25,4], Q, 
    a, P,  k , β, false, false)
        resoudre_instance_via_dataset("angers", "data/terms_49.txt",[47.40904,47.50549,-0.65266,-0.42212],"data/clvl1.csv","data/clvl2.csv",[100,25,4], Q, 
    a, P,  k , β, true, false)
        resoudre_instance_via_dataset("angers", "data/terms_49.txt",[47.40904,47.50549,-0.65266,-0.42212],"data/clvl1.csv","data/clvl2.csv",[100,25,4], Q, 
    a, P,  k , β, true, false)
    end

    for file in readdir("data/instances_test",join=true)
        nameinstance = String(split(Vector(split(file, '/'))[end],'.')[1])   
        for k in [0.1,0.5,0.9]

        resoudre_instance_sanchez(file, Q,a,P,k,β,false,false)
        resoudre_instance_sanchez(file, Q,a,P,k,β,true,false)
        resoudre_instance_sanchez(file, Q,a,P,k,β,true,true)
        comparaison_scatter_vOpt(nameinstance,k)

        end
    end

    
    
    

end

main()