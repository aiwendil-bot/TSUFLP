"""
    isFeasible()

- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-30

# Arguments

# Examples

```jldoctest
julia>
```
"""
function isFeasible(solution::Vector{Vector{Int64}}, Q::Int64)::Bool

    #=
    vérifier :
        - chaque terminal est relié à un CLVL1
        - chaque CLVL1 est relié à max Q terminaux
        - chaque CLVL1 est relié à un CLVL2
    =#

    compteur = zeros(Int64,length(solution[2]))

    for k in solution[1]
        if !(k in solution[3])
            return false
        end
        compteur[k] +=1
    end
    
    for k in compteur
        if k > Q
            return false
        end
    end
    
    for k in solution[3]
        if !(solution[2][k] in solution[4])
            return false

        end
    end
    
    return true
    

end