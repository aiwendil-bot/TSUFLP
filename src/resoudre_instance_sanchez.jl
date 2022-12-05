

function resoudre_instance_sanchez(filename::String, Q::Int64 , 
                    a::Float64, P::Int64, tenure::Int64 , k::Float64 , β::Int64)

    nameinstance = String(split(Vector(split(filename, '/'))[end],'.')[1])               

    terms, conclvl1, conclvl2 = parsing_instance(filename)

    # génération des matrices de distance, c , b, s
    
    distances = generation_distance_sanchez(terms, conclvl1)

    c = generation_c_sanchez(length(terms),length(conclvl1),minimum(distances),maximum(distances))

    b = generation_distance_sanchez(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in eachindex(conclvl2)]

    scatter_search(nameinstance,c,b,s,distances,Q,a,P,tenure,k,β)


    
end

