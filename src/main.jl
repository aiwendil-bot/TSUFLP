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
include("plot_test.jl")
include("path_relinking.jl")
include("isFeasible.jl")
#include("SkipList.jl")
include("scatter_search.jl")

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

    P = 150

    tenure = 7

    k = 0.2

    β = 4
    #génération coûts
    distances = generation_matrice_distance(coord_terminaux,coord_clvl1)
    
    distances_concentrators = generation_matrice_distance(coord_clvl1,coord_clvl2)

    couts_clvl1 = generation_couts_ouverture_clvl(size(coord_clvl1,1))

    c = 1 ./ distances
    b = generation_matrice_b(distances_concentrators,couts_clvl1)

    s = generation_couts_ouverture_clvl(size(coord_clvl2,1))

    show(scatter_search(c,b,s,distances,Q,a,P,tenure,k,β))
    #=
    println("generation pop de taille $P w/ grasp")
    solutions = @time grasp(I, J, K, Q, b, c, s, distances, a, λ, P)

    #graph_sol(coord_terminaux, coord_clvl1, coord_clvl2, solutions)
    
    objective_values = [evaluate_solution(3,k,distances,c,b,s) for k in solutions]
    z = [[objective_values[k][1] for k in eachindex(objective_values)], [objective_values[k][2] for k in eachindex(objective_values)] ] 
    pop = scatter([z[1][k] for k in 1:Int(P/5)],[z[2][k] for k in 1:Int(P/5)],label="λ = 1", mc=:blue)
    scatter!([z[1][k] for k in 1+Int(P/5):Int(2*P/5)],[z[2][k] for k in 1+Int(P/3):Int(2*P/3)],label="λ = 0.75", mc=:cyan)
    scatter!([z[1][k] for k in 1+Int(2*P/5):Int(3*P/5)],[z[2][k] for k in 1+Int(2*P/5):Int(3*P/5)],label="λ = 0.5", mc=:green)
    scatter!([z[1][k] for k in 1+Int(3*P/5):Int(4*P/5)],[z[2][k] for k in 1+Int(3*P/5):Int(4*P/5)],label="λ = 0.25", mc=:orange)
    scatter!([z[1][k] for k in 1+Int(4*P/5):P],[z[2][k] for k in 1+Int(4*P/5):P],label="λ = 0", mc=:red)

    xlabel!(L"z_1")
    ylabel!(L"z_2")
    xlims!(150,300)
    ylims!(8,35)
    savefig(pop,"out/swap/population.png")

    

    #tabu_test = @time tabu(1, solutions[1],Q,7,0.5, c,b,s, distances)
    println("tabu 1st population, k = $k, tenure = $tenure")
    first_improvment = @time vcat([ tabu(1,solutions[i],Q,tenure,k,c,b,s,distances) for i in 1:Int(P/2)],
                            [tabu(2,solutions[i],Q,tenure,k,c,b,s,distances) for i in 1+Int(P/2):P])

    
    objective_values = [evaluate_solution(3,k,distances,c,b,s) for k in first_improvment]                       
    z = [[objective_values[k][1] for k in eachindex(objective_values)], [objective_values[k][2] for k in eachindex(objective_values)] ] 
    pop_improv = scatter([z[1][k] for k in 1:Int(P/5)],[z[2][k] for k in 1:Int(P/5)],label="λ = 1", mc=:blue)
    scatter!([z[1][k] for k in 1+Int(P/5):Int(2*P/5)],[z[2][k] for k in 1+Int(P/3):Int(2*P/3)],label="λ = 0.75", mc=:cyan)
    scatter!([z[1][k] for k in 1+Int(2*P/5):Int(3*P/5)],[z[2][k] for k in 1+Int(2*P/5):Int(3*P/5)],label="λ = 0.5", mc=:green)
    scatter!([z[1][k] for k in 1+Int(3*P/5):Int(4*P/5)],[z[2][k] for k in 1+Int(3*P/5):Int(4*P/5)],label="λ = 0.25", mc=:orange)
    scatter!([z[1][k] for k in 1+Int(4*P/5):P],[z[2][k] for k in 1+Int(4*P/5):P],label="λ = 0", mc=:red)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    xlims!(150,300)
    ylims!(8,35)
    savefig(pop_improv,"out/swap/population_improv.png")
    
    println("calcul refsets, β = $β")
    refSets = @time create_refset(first_improvment,β, b,c, s, distances)

    objective_values1 = [evaluate_solution(3,k,distances,c,b,s) for k in refSets[1]]                       
    objective_values2 = [evaluate_solution(3,k,distances,c,b,s) for k in refSets[2]]

    pop_refset = scatter([objective_values1[k][1] for k in 1:β], [objective_values1[k][2] for k in 1:β], label="RefSet1", mc=:blue,legend=false)
    scatter!([objective_values2[k][1] for k in 1:β], [objective_values2[k][2] for k in 1:β], label="RefSet2", mc=:red)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    xlims!(150,300)
    ylims!(8,35)
    savefig(pop_refset,"out/swap/refSets.png")

    for i in refSets[1]
        for j in refSets[2]
            sols_intermediaires = @time path_relinking!(i,j,[],Q,c,b,distances,s)
            #println(length(sols_intermediaires))
            #println(length([k for k in sols_intermediaires if isFeasible(k,Q)]))


            scatter!([evaluate_solution(1,sols_intermediaires[k],distances,c,b,s) for k in 1:(length(sols_intermediaires)-1)],
             [evaluate_solution(2,sols_intermediaires[k],distances,c,b,s) for k in 1:(length(sols_intermediaires)-1)],mc=:orange) 
        end    
    end     
    
    savefig(pop_refset,"out/path_relinking.png")

=#
end

main()