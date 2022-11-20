#=
generation_couts:
- Julia version: 1.8.0
- Author: adrien
- Date: 2022-11-20
=#

include("utilities.jl")

using Random
Random.seed!(1234)

function generation_couts_ouverture_clvl(nb_concentrateurs::Int64)::Vector{Float64}
    return rand(Float64,nb_concentrateurs) .* 10
end

function generation_matrice_distance(liste1,liste2)::Array{Float64,2}

    distances = zeros((length(liste1[:,1]),length(liste2[:,1])))
    for i in 1:length(liste1[:,1])
        for j in 1:length(liste2[:,1])
            distances[i,j] = distance(liste1[i,:],liste2[j,:])
        end
    end
    return distances
end

function generation_matrice_b(distances::Array{Float64,2},couts::Vector{Float64})::Array{Float64,2}

    b = distances
    for i in 1:size(distances,1)
        b[i,:] = b[i,:] .+ couts
    end
    return b
end
