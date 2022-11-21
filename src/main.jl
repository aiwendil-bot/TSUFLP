#=
main:
- Julia version: 1.8.0
- Author: sujin
- Date: 2022-11-02
=#
using Random
Random.seed!(1234)

include("generation_couts.jl")
include("recup_data.jl")
include("grasp.jl")

function main()

    terminaux = readdlm("data/terminaux.txt")
    coord_terminaux = terminaux[rand(1:size(terminaux,1),200),:]
    coord_clvl1 = readdlm("data/clvl1.txt")
    coord_clvl2 = readdlm("data/clvl2.txt")

    I = [i for i in 1:size(coord_terminaux,1)]
    J = [i for i in 1:size(coord_clvl1,1)]
    K = [i for i in 1:size(coord_clvl2,1)]

    Q::Int64 = 7

    a::Float64 = 0.4

    P = 20

    #génération coûts
    distances = generation_matrice_distance(coord_terminaux,coord_clvl1)
    couts_clvl1 = generation_couts_ouverture_clvl(size(coord_clvl1,1))
    couts_clvl2 = generation_couts_ouverture_clvl(size(coord_clvl2,1))

    b = generation_matrice_b(distances,couts_clvl1)
    s = generation_couts_ouverture_clvl(size(coord_clvl2,1))

    @time display(grasp(I, J, K, Q, b, distances, s, distances, a, P))

end

main()