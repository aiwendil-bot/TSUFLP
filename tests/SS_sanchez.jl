include("../src/resoudre_instance_sanchez.jl")


function main(args)
    if length(args) == 0
        resoudre_instance_sanchez("data/files/small2.txt")

    elseif  length(args) ==1
        resoudre_instance_sanchez("data/files/$(args[1]).txt")
    else

        Q = parse(Int64,args[2])
        a = parse(Float64,args[3])
        P = parse(Int64,args[4])
        k = parse(Float64,args[5])
        β = parse(Int64,args[6])
        crossover_on = parse(Bool,args[7])
        tabu_crossover_on = parse(Bool,args[8])

        resoudre_instance_sanchez("data/files/$(args[1]).txt", Q,a,P,k,β,crossover_on,tabu_crossover_on )
    end
end

main(ARGS)

