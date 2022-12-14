#=
partir des instances du papier :

clients = terminaux

CLVL1 = 4/5 des dépôts

CLVL2 = 1/5 des dépôts

File format:
    The first line contains four integers, which corresponds to m n p r
    The next m lines contains information about the candidate location to host a facility:
        Each line is conformed with two floats that identifies the x and y coordinates
    The next n lines contains information about the demand points:
    Each line is conformed with two floats that identifies the x and y coordinates
    and an integer that represents the weight of the demand point (1 in all the examples)

=#



function convert(vector::Vector{SubString{String}})::Vector{Float64}
    return [parse(Float64, vector[1]),parse(Float64, vector[1])]

end
#parser
#return liste de liste des coordonnées [terminaux, clvl1, clvl2])
function parsing_instance(filename::String)
    # Instance
#f = open("home/sujin/IdeaProjects/TSUFLP/data/didactique.dat", "r")
open(filename) do f
    lines = readlines(f)
    line1 = Vector(parse.(Int64, split(lines[1], ' ')))
    m = Int(line1[1])
    n = Int(line1[2])

    concentrateurs = map(elem -> convert(split(elem," ")),lines[2:m+1])
    clvl1 = concentrateurs[1:Int(floor(4*m/5))]
    clvl2 = concentrateurs[Int(floor(4*m/5)+1):m]
    terminaux = map(elem -> convert(split(elem," ")),lines[m+2:length(lines)])

    return [terminaux,clvl1,clvl2]


end

end
