
#= différentes fonctions pour afficher les résultats de 
grasp, 
tabu, 
front de pareto, 
comparaison entre les stratégies
=#
using LaTeXStrings, Plots
function plot_grasp(solutions::Vector{Solution},name_instance::String, c::Array{Int64,2},b::Array{Int64,2},s::Vector{Int64},distances::Array{Int64,2})
    
    P = length(solutions)

    objective_values = [evaluate_solution(3,k,distances,c,b,s) for k in solutions]
    z = [[objective_values[k][1] for k in eachindex(objective_values)], [objective_values[k][2] for k in eachindex(objective_values)] ] 

    pop = scatter([z[1][k] for k in 1:Int(P/5)],[z[2][k] for k in 1:Int(P/5)],label="λ = 1", mc=:blue)
    scatter!([z[1][k] for k in 1+Int(P/5):Int(2*P/5)],[z[2][k] for k in 1+Int(P/3):Int(2*P/3)],label="λ = 0.75", mc=:cyan)
    scatter!([z[1][k] for k in 1+Int(2*P/5):Int(3*P/5)],[z[2][k] for k in 1+Int(2*P/5):Int(3*P/5)],label="λ = 0.5", mc=:green)
    scatter!([z[1][k] for k in 1+Int(3*P/5):Int(4*P/5)],[z[2][k] for k in 1+Int(3*P/5):Int(4*P/5)],label="λ = 0.25", mc=:orange)
    scatter!([z[1][k] for k in 1+Int(4*P/5):P],[z[2][k] for k in 1+Int(4*P/5):P],label="λ = 0", mc=:red)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    title!("population initiale de $name_instance")

    if !ispath("out/$name_instance")
        mkdir("out/$name_instance")
    end    

    savefig(pop,"out/$name_instance/population_initiale.png")

end

function plot_tabu(solutions::Vector{Solution},name_instance::String,c::Array{Int64,2},b::Array{Int64,2},s::Vector{Int64},distances::Array{Int64,2})

    P = length(solutions)

    objective_values = [evaluate_solution(3,k,distances,c,b,s) for k in solutions]
    z = [[objective_values[k][1] for k in eachindex(objective_values)], [objective_values[k][2] for k in eachindex(objective_values)] ] 

    pop = scatter([z[1][k] for k in 1:Int(P/2)],[z[2][k] for k in 1:Int(P/2)],label="lead obj1", mc=:blue)
    scatter!([z[1][k] for k in Int(P/2)+1:P],[z[2][k] for k in Int(P/2)+1:P],label="lead obj2", mc=:red)
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    title!("premier tabu de $name_instance")
    savefig(pop,"out/$name_instance/tabu.png")
end

function plot_pareto(skiplist::SkipList,name_instance::String,cpt_iteration::Int64)

    pareto = scatter([elem.point.x for elem in get_elems(skiplist)],[elem.point.y for elem in get_elems(skiplist)], label=L"Y_N")
    xlabel!(L"z_1")
    ylabel!(L"z_2")

    savefig(pareto,"out/$name_instance/pareto_$cpt_iteration.png")
end

function comparaison_scatter_vOpt(instance_name::String,k::Float64)

    YN_scatter = readdlm("out/$instance_name/YN_scatter_$k.txt")[2:end,:]
    YN_scatter_cross = readdlm("out/$instance_name/YN_scatter_cross_$k.txt")[2:end,:]
    YN_scatter_tabu_cross = readdlm("out/$instance_name/YN_scatter_tabu_cross_$k.txt")[2:end,:]
    YN_vOpt = readdlm("out/$instance_name/YN_vOpt.txt")[2:end,:]

    comparaison = scatter(YN_scatter[:,1],YN_scatter[:,2],color=:blue,label=L"$Y_N$ scatter")
    scatter!(YN_scatter_cross[:,1],YN_scatter_cross[:,2],color=:green,label=L"$Y_N$ scatter_cross")
    scatter!(YN_scatter_tabu_cross[:,1],YN_scatter_tabu_cross[:,2],color=:orange,label=L"$Y_N$ scatter_tabu_cross")
    
    scatter!(YN_vOpt[:,1],YN_vOpt[:,2],color=:red,label=L"$Y_N$ vOpt")
    title!("$instance_name avec k = $k")
    xlabel!(L"z_1")
    ylabel!(L"z_2")
    savefig(comparaison,"out/$instance_name/comparaison_YN_$k.png")
end