

using DelimitedFiles

include("utilities.jl")

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

    return minimum([distance_euclidienne(res[k1,:],vopt[k2,:]) for k1 in 2:size(res,1) for k2 in 2:size(vopt,2)])
    
    
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

    
    println(domine(yn1[2,:],yn2[2:end,:]))
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
compare_temps()
compare_cardYN()