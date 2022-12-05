struct Solution

    assign_term::Vector{Int64}
    assign_conclvl1::Vector{Int64}
    conclvl1_ouverts::Vector{Int64}
    conclvl2_ouverts::Vector{Int64}

end

function Base.:(==)(sol1::Solution, sol2::Solution)::Bool

    return ((sol1.assign_term == sol2.assign_term) && (sol1.assign_conclvl1 == sol2.assign_conclvl1)
            &&(Set(sol1.conclvl1_ouverts) == Set(sol2.conclvl1_ouverts)) && (Set(sol1.conclvl2_ouverts) == Set(sol2.conclvl2_ouverts)))
end