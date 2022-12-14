
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
include("plot_sols.jl")
include("crossover.jl")
include("vopt_model.jl")

#=

nameinstance : nom de sortie (fichiers dans out/nameinstance)
data_terms : csv où trouver les terminaux
zone : [latitude_inf,latitude_sup,longitude_inf,longitude_sup]
csvconc1, csvconc2 csv des coordonnées des concentrateurs niveau 1 et 2 (via overpass-turbo)
taille [nb_terms, nb_clvl1, nb_clvl2]

=#


function resoudre_instance_via_dataset(taille::Vector{Int64}, Q::Int64=7 , 
    a::Float64=0.4, P::Int64=300,  k::Float64=0.5 , β::Int64=10, crossover_on::Bool=false, tabu_cross_on::Bool=false)

    nameinstance = "angers"

    terms, conclvl1, conclvl2 = readdlm("data/angers/terminaux.txt")[1:taille[1],:], readdlm("data/angers/clvl1.txt")[1:taille[2],:], 
    readdlm("data/angers/clvl2.txt")[1:taille[3],:]
    # génération des matrices de distance, c , b, s
    
    distances = generation_matrice_distance(terms, conclvl1)

    c = generation_c_sanchez(size(terms,1),size(conclvl1,1),minimum(distances),maximum(distances))

    b = generation_matrice_distance(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in 1:size(conclvl2,1)]

    starting_time = time()

    YN, XE = scatter_search(nameinstance,c,b,s,distances,Q,a,P,k,β,crossover_on,tabu_cross_on)

    ending_time = round(time()- starting_time, digits=4)

    if crossover_on

        if tabu_cross_on

            writedlm("out/angers/YN_scatter_tabu_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))
            writedlm("out/angers/XE_scatter_tabu_cross_$k.txt",XE)
            


        else

            writedlm("out/angers/YN_scatter_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))
            writedlm("out/angers/XE_scatter_cross_$k.txt",XE)

        end

    else    

        writedlm("out/angers/YN_scatter_$k.txt",vcat([[ending_time,length(YN)]],YN))
        writedlm("out/angers/XE_scatter_$k.txt",XE)

    end

    println("fichiers de résultats (plots, Y_N, X_E) écrits dans out/angers")
    
end

