#=
main:
- Julia version: 1.8.0
- Author: sujin
- Date: 2022-11-02
=#
using Random, Profile, Plots, LaTeXStrings
Random.seed!(1235)

include("generation_couts.jl")
include("recup_data.jl")
include("grasp.jl")
include("tabu.jl")
include("refset.jl")

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

    λ::Float64 = 0.5

    P = 60

    β = Int(P/2)

    #génération coûts
    distances = generation_matrice_distance(coord_terminaux,coord_clvl1)
    
    distances_concentrators = generation_matrice_distance(coord_clvl1,coord_clvl2)

    couts_clvl1 = generation_couts_ouverture_clvl(size(coord_clvl1,1))

    c = 1 ./ distances
    b = generation_matrice_b(distances_concentrators,couts_clvl1)

    s = generation_couts_ouverture_clvl(size(coord_clvl2,1))

    solutions = @time grasp(I, J, K, Q, b, c, s, distances, a, λ, P)
    
    
    objective_values = [evaluate_solution(k,distances,c,b,s) for k in solutions]
    z = [[objective_values[k][1] for k in 1:length(objective_values)], [objective_values[k][2] for k in 1:length(objective_values)] ] 
    pop = scatter([z[1][k] for k in 1:Int(P/3)],[z[2][k] for k in 1:Int(P/3)],label="lead obj 1", mc=:blue)
    scatter!([z[1][k] for k in 1+Int(P/3):Int(2*P/3)],[z[2][k] for k in 1+Int(P/3):Int(2*P/3)],label="lead obj 2", mc=:red)
    scatter!([z[1][k] for k in 1+Int(2*P/3):P],[z[2][k] for k in 1+Int(2*P/3):P],label="compromis", mc=:violet)
    xlabel!(L"z_1")
    ylabel!(L"z_2")

    savefig(pop,"out/population.png")

    

    #tabu_test = @time tabu(1, solutions[1],Q,7,0.5, c,b,s, distances)

    first_improvment = vcat([tabu(1,solutions[k],Q,7,0.2,c,b,s,distances) for k in 1:Int(P/2)],
                            [tabu(2,solutions[k],Q,7,0.2,c,b,s,distances) for k in 1+Int(P/2):P])
    #=
    objective_values = [evaluate_solution(k,distances,c,b,s) for k in first_improvment]                       
    z = [[objective_values[k][1] for k in 1:length(objective_values)], [objective_values[k][2] for k in 1:length(objective_values)] ] 
    pop_improv = scatter([z[1][k] for k in 1:Int(P/3)],[z[2][k] for k in 1:Int(P/3)],label="lead obj 1", mc=:blue)
    scatter!([z[1][k] for k in 1+Int(P/3):Int(2*P/3)],[z[2][k] for k in 1+Int(P/3):Int(2*P/3)],label="lead obj 2", mc=:red)
    scatter!([z[1][k] for k in 1+Int(2*P/3):P],[z[2][k] for k in 1+Int(2*P/3):P],label="compromis", mc=:violet)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    savefig(pop_improv,"out/population_improv.png")
    =#

    refSets = update_refset(first_improvment,β, b,c, s, distances)
    println(length(refSets[1]))
    println(length(refSets[2]))
    objective_values1 = [evaluate_solution(k,distances,c,b,s) for k in refSets[1]]                       
    objective_values2 = [evaluate_solution(k,distances,c,b,s) for k in refSets[2]]

    pop_refset = scatter([objective_values1[k][1] for k in 1:β], [objective_values1[k][2] for k in 1:β], label="RefSet1", mc=:blue)
    scatter!([objective_values2[k][1] for k in 1:β], [objective_values2[k][2] for k in 1:β], label="RefSet2", mc=:red)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    savefig(pop_refset,"out/refSets.png")
end

main()