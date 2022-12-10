using vOptGeneric, GLPK, DelimitedFiles, JuMP, Gurobi, Random
Random.seed!(1235)

include("../src/Solution.jl")
include("parsing_instance.jl")
include("../src/generation_couts_sanchez.jl")
include("../src/generation_couts.jl")

#=
terminaux = readdlm("data/terminaux.txt")
coord_terminaux = terminaux[rand(1:size(terminaux,1),200),:]
coord_clvl1 = readdlm("data/clvl1.txt")
coord_clvl2 = readdlm("data/clvl2.txt")
=#

function vopt_resolve(filename::String)

    nameinstance = String(split(Vector(split(filename, '/'))[end],'.')[1])               

    terms, conclvl1, conclvl2 = parsing_instance(filename)

    # génération des matrices de distance, c , b, s

    distances = generation_distance_sanchez(terms, conclvl1)

    c = generation_c_sanchez(length(terms),length(conclvl1),minimum(distances),maximum(distances))

    b = generation_distance_sanchez(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in eachindex(conclvl2)]


    nI = size(terms,1)
    nJ = size(conclvl1,1)
    nK = size(conclvl2,1)


    Q = 7

    model = vModel(Gurobi.Optimizer)

    @variable(model,x[1:nI,1:nJ], Bin)
    @variable(model,y[1:nJ,1:nK], Bin)
    @variable(model,z[1:nK], Bin)
    @variable(model, Z >=0)

    @addobjective(model, Min, sum(c[i,j]*x[i,j] for i in 1:nI for j in 1:nJ) + 
                            sum(b[j,k]*y[j,k] for j in 1:nJ for k in 1:nK) +
                            sum(s[k]*z[k]     for k in 1:nK))

    @addobjective(model, Min, Z)

    
    @constraint(model, [i = 1:nI] ,sum(x[i,j] for j in 1:nJ) == 1)
    
    @constraint(model, [i=1:nI, j=1:nJ] ,x[i,j] <= sum(y[j,k] for k in 1:nK))

    @constraint(model, [j=1:nJ, k=1:nK], y[j,k] <= z[k])

    @constraint(model, [j = 1:nJ] ,sum(y[j,k] for k in 1:nK) <= 1)
    
    #contraintes capacité
    @constraint(model, [j = 1:nJ] ,sum(x[i,j] for i in 1:nI) <= Q)
    
    #contraintes pour linéariser obj2 (min max)
    @constraint(model, [i = 1:nI, j=1:nJ] ,Z >= x[i,j]*distances[i,j])

    getTime = time()

    vSolve(model, method=:epsilon, step = 1.0)
    timevOPt = round(time()- getTime, digits=4)
    YN = getY_N( model )

    if !ispath("out/$nameinstance")
        mkdir("out/$nameinstance")
    end    

    writedlm("out/$nameinstance/YN_vOpt.txt",vcat([[timevOPt,length(YN)]],YN))
    #println("tOpt: $(timevOPt) sec     #YN: $(length(YN)) points\n")
    #printX_E( model )


end

function vopt_resolve_angers()

    nameinstance = "angers"          

    terms, conclvl1, conclvl2 = readdlm("data/terminaux.txt")[1:100,:], readdlm("data/clvl1.txt")[25:50,:], readdlm("data/clvl2.txt") 

 
    # génération des matrices de distance, c , b, s

    distances = generation_matrice_distance(terms, conclvl1)

    c = generation_c_sanchez(length(terms),length(conclvl1),minimum(distances),maximum(distances))

    b = generation_matrice_distance(conclvl1, conclvl2)

    s = [rand([i for i in minimum(b):maximum(b)]) for k in eachindex(conclvl2)]


    nI = size(terms,1)
    nJ = size(conclvl1,1)
    nK = size(conclvl2,1)


    Q = 7

    model = vModel(Gurobi.Optimizer)
    set_silent(model)
    @variable(model,x[1:nI,1:nJ], Bin)
    @variable(model,y[1:nJ,1:nK], Bin)
    @variable(model,z[1:nK], Bin)
    @variable(model, Z >=0)

    @addobjective(model, Min, sum(c[i,j]*x[i,j] for i in 1:nI for j in 1:nJ) + 
                            sum(b[j,k]*y[j,k] for j in 1:nJ for k in 1:nK) +
                            sum(s[k]*z[k]     for k in 1:nK))

    @addobjective(model, Min, Z)

    
    @constraint(model, [i = 1:nI] ,sum(x[i,j] for j in 1:nJ) == 1)
    
    @constraint(model, [i=1:nI, j=1:nJ] ,x[i,j] <= sum(y[j,k] for k in 1:nK))

    @constraint(model, [j=1:nJ, k=1:nK], y[j,k] <= z[k])

    @constraint(model, [j = 1:nJ] ,sum(y[j,k] for k in 1:nK) <= 1)
    
    #contraintes capacité
    @constraint(model, [j = 1:nJ] ,sum(x[i,j] for i in 1:nI) <= Q)
    
    #contraintes pour linéariser obj2 (min max)
    @constraint(model, [i = 1:nI, j=1:nJ] ,Z >= x[i,j]*distances[i,j])
    #optimize!(model)
    getTime = time()

    vSolve(model, method=:epsilon, step = 1.0)
    
    timevOPt = round(time()- getTime, digits=4)
    YN = getY_N( model )

    if !ispath("out/$nameinstance")
        mkdir("out/$nameinstance")
    end    

    writedlm("out/$nameinstance/YN_vOpt.txt",vcat([[timevOPt,length(YN)]],YN))
    #println("tOpt: $(timevOPt) sec     #YN: $(length(YN)) points\n")
    #printX_E( model )
    

end