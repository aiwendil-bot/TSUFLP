using vOptGeneric, GLPK, DelimitedFiles, JuMP


include("../src/generation_couts.jl")


terminaux = readdlm("data/terminaux.txt")
coord_terminaux = terminaux[rand(1:size(terminaux,1),200),:]
coord_clvl1 = readdlm("data/clvl1.txt")
coord_clvl2 = readdlm("data/clvl2.txt")

nI = size(coord_terminaux,1)
nJ = size(coord_clvl1,1)
nK = size(coord_clvl2,1)

distances = generation_matrice_distance(coord_terminaux,coord_clvl1)
distances_concentrators = generation_matrice_distance(coord_clvl1,coord_clvl2)
couts_clvl1 = generation_couts_ouverture_clvl(size(coord_clvl1,1))

c = 1 ./ distances
b = generation_matrice_b(distances_concentrators,couts_clvl1)
s = generation_couts_ouverture_clvl(size(coord_clvl2,1))

Q = 7

model = vModel(GLPK.Optimizer)

@variable(model,x[1:nI,1:nJ], Bin)
@variable(model,y[1:nJ,1:nK], Bin)
@variable(model,z[1:nK], Bin)
@variable(model, Z >=0)

@addobjective(model, Min, sum(c[i,j]*x[i,j] for i in 1:nI for j in 1:nJ) + 
                          sum(b[j,k]*y[j,k] for j in 1:nJ for k in 1:nK) +
                          sum(s[k]*z[k]     for k in 1:nK))

@addobjective(model, Min, Z)

for i in 1:nI
    @constraint(model, sum(x[i,j] for j in 1:nJ) == 1)
end

for i in 1:nI
    for j in 1:nJ
        @constraint(model, x[i,j] <= sum(y[j,k] for k in 1:nK))
    end
end

for j in 1:nJ
    for k in 1:nK
        @constraint(model, y[j,k] <= z[k])
    end
end

for j in 1:nJ
    @constraint(model, sum(y[j,k] for k in 1:nK) <= 1)
end

#contraintes capacité

for j in 1:nJ
    @constraint(model, sum(x[i,j] for i in 1:nI) <= Q)
end    

#contraintes pour linéariser obj2 (min max)

for i in 1:nI
    for j in 1:nJ
        @constraint(model, Z >= x[i,j]*distances[i,j])
    end
end
getTime = time()
vSolve(model, method=:epsilon, step = 0.5)
timevOPt = round(time()- getTime, digits=4)
YN = getY_N( model )
println("tOpt: $(timevOPt) sec     #YN: $(length(YN)) points\n")
printX_E( model )


