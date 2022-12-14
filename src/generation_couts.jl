include("utilities.jl")

using Random
Random.seed!(1234)

#=

fonctions servant à générer les coûts c b z pour une instance type angers (depuis dataset)

=#

function generation_couts_ouverture_clvl(nb_concentrateurs::Int64)::Vector{Int64}
    return rand(Int64,nb_concentrateurs) .* 10
end

function generation_matrice_distance(liste1,liste2)::Array{Int64,2}

    distances = zeros((length(liste1[:,1]),length(liste2[:,1])))
    for i in 1:length(liste1[:,1])
        for j in 1:length(liste2[:,1])
            distances[i,j] = Int(floor(distance(liste1[i,:],liste2[j,:]) * 1000))
        end
    end
    return distances
end

function generation_matrice_b(distances::Array{Int64,2},couts::Vector{Int64})::Array{Int64,2}

    b = distances
    for i in 1:size(distances,2)
        b[:,i] = b[:,i] .+ couts
    end
    return b
end
