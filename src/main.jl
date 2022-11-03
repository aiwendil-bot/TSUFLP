#=
main:
- Julia version: 1.8.0
- Author: sujin
- Date: 2022-11-02
=#

# Instance
#f = open("home/sujin/IdeaProjects/TSUFLP/data/didactique.dat", "r")
open("data/didactique.dat") do f

    # line_number
    line = 0

    # read till end of file
    tab_str = readlines(f)
    ligne1 = Tuple(parse.(Int64, split(tab_str[1], ' ')))
    NK = ligne1[1]
    NJ = ligne1[2]
    NI = ligne1[3]

    println(NK, " ", NJ, " ", NI)
    coord_ctr_2 = [Vector(parse.(Float64, split(tab_str[i], ' '))) for i in 2:(NK+1)]
    coord_ctr_1 = [Vector(parse.(Float64, split(tab_str[i], ' '))) for i in (NK+2):(NK+1+NJ)]
    coord_terminaux = [Vector(parse.(Float64, split(tab_str[i], ' '))) for i in (NJ+NK+2):(NJ+NK+1+NI)]
    println(coord_ctr_2)
    println(coord_ctr_1)
    println(coord_terminaux)

    d = zeros(Float64, (NJ, NI))

    for i in 1:NJ
       for j in 1:NI
           d_ij[i,j] = (coord_ctr_1[i][1] - coord_terminaux[j][1])^2 + (coord_ctr_1[i][2] - coord_terminaux[j][2])^2
       end
    end

    for i in 1:NK
       for j in 1:NJ
           d_ij[i,j] = (coord_ctr_1[i][1] - coord_terminaux[j][1])^2 + (coord_ctr_1[i][2] - coord_terminaux[j][2])^2
       end
    end

    show(d_ij)
end