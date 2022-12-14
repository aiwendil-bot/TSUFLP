Random.seed!(1235)

#=

fonctions servant à générer les coûts pour les instances de Sanchez 

=#


function distance_euclidienne(coord1::Vector{Float64},coord2::Vector{Float64})::Int64

    return Int(floor((coord1[1]-coord2[1])^2 + (coord1[2]-coord2[2])^2 ))
    
end

function generation_distance_sanchez(vec1::Vector{Vector{Float64}}, vec2::Vector{Vector{Float64}})::Array{Int64,2}

    distances = zeros(Int64,(length(vec1),length(vec2)))

    for i in eachindex(vec1)
        for j in eachindex(vec2)
            distances[i,j] = distance_euclidienne(vec1[i],vec2[j])
        end
    end
    return distances
end

function generation_c_sanchez(size1::Int64,size2::Int64,min::Int64,max::Int64)::Array{Int64,2}

    c = zeros(Int64,(size1,size2))

    for i in 1:size1
        for j in 1:size2
            c[i,j] = rand([i for i in min:max])
        end
    end

    return c
    
end
