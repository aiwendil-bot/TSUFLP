include("../src/resoudre_instance_via_dataset.jl")




function main(args)
    if length(args) == 0
        resoudre_instance_via_dataset([100,25,4])
    else

        taille1 = parse(Int64,args[1])
        taille2 = parse(Int64,args[2])
        taille3 = parse(Int64,args[3])
        Q = parse(Int64,args[4])
        a = parse(Float64,args[5])
        P = parse(Int64,args[6])
        k = parse(Float64,args[7])
        β = parse(Int64,args[8])
        crossover_on = parse(Bool,args[9])
        tabu_cross_on = parse(Bool,args[10])

        resoudre_instance_via_dataset([taille1,taille2,taille3],Q,a,P,k,β,crossover_on,tabu_cross_on)
    end

    
end

main(ARGS)