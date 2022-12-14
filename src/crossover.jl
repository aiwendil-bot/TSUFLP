
#=
échange la moitié des terminaux entre deux concentrateurs (ouverts)
renvoie deux nouvelles solutions (selon si on échange la 1ère ou la 2nde moitié)
=#

function crossover(sol::Solution, conc1::Int64, conc2::Int64, Q::Int64)::Vector{Solution}

    res1 = deepcopy(sol)
    res2 = deepcopy(sol)


    terms1 = findall(x -> x == conc1, res1.assign_term)
    terms2 = findall(x -> x == conc2, res1.assign_term)

    if (length(terms1) == Q) && (length(terms2) == Q)
        for k in 1:Int(floor(Q))
            res1.assign_term[terms1[k]] = conc2
            res1.assign_term[terms2[k]] = conc1
        end

        for k in 1+Int(floor(Q)):Q
            res2.assign_term[terms1[k]] = conc2
            res2.assign_term[terms2[k]] = conc1
        end

    end

    return [res1, res2]

end
