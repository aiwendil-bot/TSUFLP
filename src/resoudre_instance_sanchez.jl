

function resoudre_instance_sanchez(filename::String, Q::Int64 , 
                    a::Float64, P::Int64,  k::Float64 , β::Int64, crossover_on::Bool, tabu_cross_on::Bool)

    nameinstance = String(split(Vector(split(filename, '/'))[end],'.')[1])               

    terms, conclvl1, conclvl2 = parsing_instance(filename)

    # génération des matrices de distance, c , b, s
    
    distances = generation_distance_sanchez(terms, conclvl1)

    c = generation_c_sanchez(length(terms),length(conclvl1),minimum(distances),maximum(distances))

    b = generation_distance_sanchez(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in eachindex(conclvl2)]

    starting_time = time()

    YN = scatter_search(nameinstance,c,b,s,distances,Q,a,P,k,β,crossover_on,tabu_cross_on)

    ending_time = round(time()- starting_time, digits=4)

    if crossover_on

        if tabu_cross_on

            writedlm("out/$nameinstance/YN_scatter_tabu_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))


        else

            writedlm("out/$nameinstance/YN_scatter_cross_$k.txt",vcat([[ending_time,length(YN)]],YN))

        end

    else    

        writedlm("out/$nameinstance/YN_scatter_$k.txt",vcat([[ending_time,length(YN)]],YN))

    end

end

