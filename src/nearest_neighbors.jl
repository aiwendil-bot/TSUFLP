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
function nearest_neighbors(d::Vector{Vector{Int64}},origin::Int64,Q::Int64,available_neighbors::Vector{Int64})
    neighbors = [d[origin][k] for k in available_neighbors]
    nearest = Int64[]
    for n in 1:length(neighbors)
        if length(nearest) < Q
            push!(nearest,n)
        else
            if neighbors[n] < maximum([neighbors[k] for k in nearest])
                nearest[argmax([neighbors[k] for k in nearest])] = n
            end
        end
    end
    return nearest
end