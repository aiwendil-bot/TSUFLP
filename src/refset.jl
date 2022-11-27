

function update_refset(solutions::Vector{Vector{Vector{Int64}}},β::Int64,
                            b::Array{Float64,2},c::Array{Float64,2}, s::Vector{Float64}, 
                            d::Array{Float64,2})::Vector{Vector{Vector{Vector{Int64}}}}

    refSet2 = Vector{Vector{Int64}}[]

    meilleures_z1 = sort(solutions,by = sol -> evaluate_solution(1,sol,d,c,b,s))
    meilleures_z2 = sort(solutions,by = sol -> evaluate_solution(2,sol,d,c,b,s))

    refSet1 = [meilleures_z1[k] for k in 1:Int(β/2)]
    i = 1
    while length(refSet2) < Int(β/2)
            
        if !(meilleures_z2[i] in refSet1)
                push!(refSet2,meilleures_z2[i])
            end
        i += 1    
    end


    pull_refset1 = sort(setdiff(solutions,vcat(refSet1,refSet2)),by= sol -> minimum([distance_solutions(sol,k) for k in refSet1]),rev=true)
    pull_refset2 = sort(setdiff(solutions,vcat(refSet1,refSet2)),by= sol -> minimum([distance_solutions(sol,k) for k in refSet2]),rev=true)
    
    k = 1
    while length(refSet1) < β 


        if !(pull_refset1[k] in vcat(refSet1,refSet2))
            push!(refSet1,pull_refset1[k])
        end     
        
        k += 1
    end
    k = 1
    while length(refSet2) < β

        if !(pull_refset2[k] in vcat(refSet1,refSet2))
            push!(refSet2,pull_refset2[k])
        end

        k+= 1
    end
    return [refSet1,refSet2]



end
