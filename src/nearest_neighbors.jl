"""
    nearest_neighbors()

- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-31

# Arguments


"""

#calcule les Q plus proches voisins de origin
#d : lignes = concentrateurs & colonnes = terminaux
#O(n)
function nearest_neighbors(d::Matrix{Float64},origin::Int64,Q::Int64,available_neighbors::Vector{Int64})

    neighbors = [[k,d[k, origin]] for k in available_neighbors]
    nearest = Int64[]
    for n in 1:length(neighbors)
        if length(nearest) < Q
            push!(nearest,available_neighbors[n])
        else
            
            if neighbors[n][2] < maximum([neighbors[i][2] for i in 1:length(neighbors) if (neighbors[i][1] in nearest)])
                nearest[argmax([neighbors[i][2] for i in 1:length(neighbors) if (neighbors[i][1] in nearest)])] = available_neighbors[n]
            end
        end
    end
    return nearest
end