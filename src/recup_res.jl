

using DelimitedFiles
include("Solution.jl")
include("utilities.jl")


#= fonctions permettant de récupérer les résultats (moyennes, coverage measures, ...) de nos tests sur les
instances (tests/instances_test)
=# 

smalls = ["small1","small2","small7","small10","small11","small14","small16","small17","small18","small20"]
mediums = ["medium5","medium7","medium9","medium11","medium13","medium17","medium20"]

function distance_euclidienne(coord1::Vector{Float64},coord2::Vector{Float64})::Int64

    return Int(floor((coord1[1]-coord2[1])^2 + (coord1[2]-coord2[2])^2 ))
    
end

# distances entre les sol vopt et les solutions de scatter


function distance_a_vopt(instance_name,k,crossover_on,tabu_cross_on)

    vopt = readdlm("out/$instance_name/YN_vOpt.txt")

    if crossover_on
        if tabu_cross_on
            res = readdlm("out/$instance_name/YN_scatter_tabu_cross_$k.txt")
        else
            res = readdlm("out/$instance_name/YN_scatter_cross_$k.txt")
        end
    else
        res  = readdlm("out/$instance_name/YN_scatter_$k.txt")
    end

    return sum(minimum([distance_euclidienne(res[k1,:],vopt[k2,:]) for k2 in 2:size(vopt,2) ]) for k1 in 2:size(res,1))/(size(res,1)-1)
 
    
end

function domine(point,ensemble)
    for i in 1:size(ensemble,1)
        if point[1] <= ensemble[i,1] && point[2] <= ensemble[i,2]
            return 1
        end
    end
            
    return 0
end

function compare_YN(filename1,filename2)

    yn1 = readdlm(filename1)
    yn2 = readdlm(filename2)

    
    return sum([domine(yn1[k,:],yn2[2:end,:]) for k in 2:size(yn1,1)])/yn1[1,2]*100

end

function moyenne_temps(instances)

    return sum([readdlm(instances)[1,1]])/length(instances)
end

function temps_instance(instancename)
    return readdlm(instancename)[1,1]
    
end

function cardYN_instance(instancename)
    return readdlm(instancename)[1,2]
    
end

function compare_temps()

    noms_instances = [String(split(Vector(split(instance, '/'))[end],'.')[1])   for instance in readdir("tests/instances_test")]

    temps_supp_cross =0
    temps_supp_cross_tabu = 0

    

    for instance_name in noms_instances
        temps_supp_cross += sum( [((temps_instance("out/$instance_name/YN_scatter_$k.txt")-temps_instance("out/$instance_name/YN_scatter_cross_$k.txt"))
                                    /(temps_instance("out/$instance_name/YN_scatter_cross_$k.txt"))) for k in [0.1,0.5,0.9]]) 

        temps_supp_cross_tabu += sum( [((temps_instance("out/$instance_name/YN_scatter_$k.txt") - temps_instance("out/$instance_name/YN_scatter_tabu_cross_$k.txt")) 
                                        / temps_instance("out/$instance_name/YN_scatter_tabu_cross_$k.txt")) for k in [0.1,0.5,0.9]]) 

    end

    println("cross : ", temps_supp_cross/(3*length(noms_instances)))
    println("cross tabu : ", temps_supp_cross_tabu/(3*length(noms_instances)))
    
end

function compare_cardYN()

    noms_instances = [String(split(Vector(split(instance, '/'))[end],'.')[1])   for instance in readdir("tests/instances_test")]

    temps_supp_cross =0
    temps_supp_cross_tabu = 0

    

    for instance_name in noms_instances
        temps_supp_cross += sum( [((cardYN_instance("out/$instance_name/YN_scatter_$k.txt")-cardYN_instance("out/$instance_name/YN_scatter_cross_$k.txt"))
                                    /(cardYN_instance("out/$instance_name/YN_scatter_cross_$k.txt"))) for k in [0.1,0.5,0.9]]) 

        temps_supp_cross_tabu += sum( [((cardYN_instance("out/$instance_name/YN_scatter_$k.txt") - cardYN_instance("out/$instance_name/YN_scatter_tabu_cross_$k.txt")) 
                                        / cardYN_instance("out/$instance_name/YN_scatter_tabu_cross_$k.txt")) for k in [0.1,0.5,0.9]]) 

    end

    println("cross : ", temps_supp_cross/(3*length(noms_instances)))
    println("cross tabu : ", temps_supp_cross_tabu/(3*length(noms_instances)))
    
end

function moyenne_temps(k)

    println("temps moyens pour k = $k")
    tsmalls = 0
    tmediums = 0
    for instance in smalls
        tsmalls += readdlm("out/$instance/YN_scatter_$k.txt")[1,1]
        #tsmalls += readdlm("out/$instance/YN_scatter_cross_$k.txt")[1,1]
        #tsmalls += readdlm("out/$instance/YN_scatter_tabu_cross_$k.txt")[1,1]
    end

    for instance in mediums
        tmediums += readdlm("out/$instance/YN_scatter_$k.txt")[1,1]
        #tmediums += readdlm("out/$instance/YN_scatter_cross_$k.txt")[1,1]
        #tmediums += readdlm("out/$instance/YN_scatter_tabu_cross_$k.txt")[1,1]
    end

    println("temps moyens smalls = ", tsmalls/(length(smalls)))
    println("temps moyens mediums = ", tmediums/(length(mediums)))


end

function temps_strat()
    tstratsmall = [0.0,0.0,0.0]
    tstratmedium = [0.0,0.0,0.0]
    strats = ["YN_scatter","YN_scatter_cross","YN_scatter_tabu_cross"]
    for strat in eachindex(strats)

        for k in [0.5]

            for instance in smalls

            tstratsmall[strat] += readdlm("out/$instance/$(strats[strat])_$k.txt")[1,1]
            end

            for instancem in mediums
                tstratmedium[strat] += readdlm("out/$instancem/$(strats[strat])_$k.txt")[1,1]
            end
        end
    end
    for k in 1:3
        println("k = $k")
    println("temps moyens smalls strat $(strats[k]) : ", tstratsmall[k]/(length(smalls)) )

        println("temps moyens mediums strat $(strats[k]) : ", tstratmedium[k]/(length(mediums)) )
    end



            

end

function dominance_k_scatter()
    dom_smalls = [0.0,0.0,0.0] #0.9 dom 0.5 / 0.9 dom 0.1 / 0.5 dom 0.1
    dom_mediums = [0.0,0.0,0.0]
    for instance in smalls
        dom_smalls[1] += compare_YN("out/$instance/YN_scatter_0.9.txt","out/$instance/YN_scatter_0.5.txt")
        dom_smalls[2] += compare_YN("out/$instance/YN_scatter_0.9.txt","out/$instance/YN_scatter_0.1.txt")
        dom_smalls[3] += compare_YN("out/$instance/YN_scatter_0.5.txt","out/$instance/YN_scatter_0.1.txt")
    end

    for instance in mediums
        dom_mediums[1] += compare_YN("out/$instance/YN_scatter_0.9.txt","out/$instance/YN_scatter_0.5.txt") 
        dom_mediums[2] += compare_YN("out/$instance/YN_scatter_0.9.txt","out/$instance/YN_scatter_0.1.txt")
        dom_mediums[3] += compare_YN("out/$instance/YN_scatter_0.5.txt","out/$instance/YN_scatter_0.1.txt")
    end
    println("pour les smalls")
    println("0.9 dom 0.5 : ", dom_smalls[1]/length(smalls))
    println("0.9 dom 0.1 : ", dom_smalls[2]/length(smalls))
    println("0.5 dom 0.1 : ", dom_smalls[3]/length(smalls))
    println("pour les mediums")
    println("0.9 dom 0.5 : ", dom_mediums[1]/length(mediums))
    println("0.9 dom 0.1 : ", dom_mediums[2]/length(mediums))
    println("0.5 dom 0.1 : ", dom_mediums[3]/length(mediums))
end

function dominance_strats()
    dom_smalls = [0.0,0.0,0.0] #tabu cross dom cross / tabu cross dom scatter / cross dom scatter
    dom_mediums = [0.0,0.0,0.0]
    for instance in smalls
        dom_smalls[1] += compare_YN("out/$instance/YN_scatter_tabu_cross_0.1.txt","out/$instance/YN_scatter_cross_0.5.txt")
        dom_smalls[2] += compare_YN("out/$instance/YN_scatter_tabu_cross_0.1.txt","out/$instance/YN_scatter_0.5.txt")
        dom_smalls[3] += compare_YN("out/$instance/YN_scatter_cross_0.1.txt","out/$instance/YN_scatter_0.5.txt")
    end

    for instance in mediums
        dom_mediums[1] += compare_YN("out/$instance/YN_scatter_tabu_cross_0.1.txt","out/$instance/YN_scatter_cross_0.5.txt")
        dom_mediums[2] += compare_YN("out/$instance/YN_scatter_tabu_cross_0.1.txt","out/$instance/YN_scatter_0.5.txt")
        dom_mediums[3] += compare_YN("out/$instance/YN_scatter_cross_0.1.txt","out/$instance/YN_scatter_0.5.txt")
    end
    println("pour les smalls")
    println("tabu cross dom cross : ", dom_smalls[1]/length(smalls))
    println("tabu cross dom scatter : ", dom_smalls[2]/length(smalls))
    println("cross dom scatter : ", dom_smalls[3]/length(smalls))
    println("pour les mediums")
    println("tabu cross dom cross : ", dom_mediums[1]/length(mediums))
    println("tabu cross dom scatter : ", dom_mediums[2]/length(mediums))
    println("cross dom scatter : ", dom_mediums[3]/length(mediums))
end

function moyenne_distance_vopt(k)
    dsmalls = [0.0,0.0,0.0]
    dmediums = [0.0,0.0,0.0]
    strats = ["YN_scatter","YN_scatter_cross","YN_scatter_tabu_cross"]

    for instance in smalls
        dsmalls[1] += distance_a_vopt(instance,k,true,true)
        dsmalls[2] += distance_a_vopt(instance,k,true,false)
        dsmalls[3] += distance_a_vopt(instance,k,false,false)
    end

    for instance in mediums
        dmediums[1] += distance_a_vopt(instance,k,true,true)
        dmediums[2] += distance_a_vopt(instance,k,true,false)
        dmediums[3] += distance_a_vopt(instance,k,false,false)
    end

    println("pour k = $k, distances à vOpt pour les smalls")
    for k in [1,2,3]
    println("pour la strat $(strats[k]) : ", dsmalls[k]/(10^11*length(smalls)))
    end
    println("pour k = $k, distances à vOpt pour les mediums")
    for k in [1,2,3]
    println("pour la strat $(strats[k]) : ", dmediums[k]/(10^11*length(mediums)))
    end

    
end

function moyenne_temps_vopt()

    
    tsmalls = 0
    tmediums = 0
    for instance in smalls
        tsmalls += readdlm("out/$instance/YN_vOpt.txt")[1,1]
        #tsmalls += readdlm("out/$instance/YN_scatter_cross_$k.txt")[1,1]
        #tsmalls += readdlm("out/$instance/YN_scatter_tabu_cross_$k.txt")[1,1]
    end

    for instance in mediums
        tmediums += readdlm("out/$instance/YN_vOpt.txt")[1,1]
        #tmediums += readdlm("out/$instance/YN_scatter_cross_$k.txt")[1,1]
        #tmediums += readdlm("out/$instance/YN_scatter_tabu_cross_$k.txt")[1,1]
    end

    println("temps moyens smalls = ", tsmalls/(length(smalls)))
    println("temps moyens mediums = ", tmediums/(length(mediums)))


end

moyenne_temps(0.1)
println("---------------")
moyenne_temps(0.5)
println("---------------")
moyenne_temps(0.9)
println("---------------")

dominance_k_scatter()
println("---------------")

temps_strat()
println("---------------")

dominance_strats()
println("---------------")

for k in [0.1,0.5,0.9]
    println("---------------")
    moyenne_distance_vopt(k)
end

moyenne_temps_vopt()