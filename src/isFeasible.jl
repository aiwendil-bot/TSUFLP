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
end