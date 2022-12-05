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
function isFeasible(solution::Solution, Q::Int64)::Bool

    #=
    vérifier :
        - chaque terminal est relié à un CLVL1
        - chaque CLVL1 est relié à max Q terminaux
        - chaque CLVL1 est relié à un CLVL2
    =#

    compteur = zeros(Int64,length(solution.assign_conclvl1))

    for k in solution.assign_term
        if !(k in solution.conclvl1_ouverts)

            return false
        end
        compteur[k] +=1
    end
    
    for k in compteur
        if k > Q


            return false
        end
    end
    
    for k in solution.conclvl1_ouverts
        if !(solution.assign_conclvl1[k] in solution.conclvl2_ouverts)


            return false

        end
    end
    
    return true
    

end