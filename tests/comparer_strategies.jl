include("../src/resoudre_instance_sanchez.jl")

function main(args)

    if length(args) == 0


        resoudre_instance_sanchez("data/files/small2.txt", 7,0.4,300,0.5,10,false,false)
        resoudre_instance_sanchez("data/files/small2.txt", 7,0.4,300,0.5,10,true,false)
        resoudre_instance_sanchez("data/files/small2.txt", 7,0.4,300,0.5,10,true,true)
        comparaison_scatter_vOpt(args[1],0.5)
    else

        if !ispath("/out/$(args[1])/Y_N_vOpt.txt")
            println("résolution exacte avec vOpt Solver ...")
            vopt_resolve("data/files/$(args[1]).txt")
        end



    Q = parse(Int64,args[2])
    a = parse(Float64,args[3])
    P = parse(Int64,args[4])
    k = parse(Float64,args[5])
    β = parse(Int64,args[6])

    resoudre_instance_sanchez("data/files/$(args[1]).txt", Q,a,P,k,β,false,false)
    resoudre_instance_sanchez("data/files/$(args[1]).txt", Q,a,P,k,β,true,false)
    resoudre_instance_sanchez("data/files/$(args[1]).txt", Q,a,P,k,β,true,true)
    comparaison_scatter_vOpt(args[1],k)
    end
    
end

main(ARGS)