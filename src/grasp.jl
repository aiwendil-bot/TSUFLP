"""
    grasp()

- Julia version: 1.8.1
- Author: adrien
- Date: 2022-10-30

# Arguments

I : terminaux
J : CLVL1
K : CLVL2
Q : capacité max des concentrateurs
b,c,s : arrays de coûts (vérifier si les couts sont float ou int)
d : arrays de distances entre terminaux et CLVL1
a : alpha de grasp
P : nb de solutions à générer, supposé pair


"""
#=
solution : vecteur composé de :

[1] : vecteur de nI éléments, valeurs = concentrateurs de lvl 1 assignés
[2] : vecteur d'au plus nJ éléments, valeurs = concentrateurs de lvl 2 assignés

optionnels ? éviterait de reparcourir [1] et [2]
[3] concentrateurs de lvl 1 ouverts
[4] concentrateurs de lvl 2 ouverts
=#



function grasp(I::Vector{Int64},J::Vector{Int64},K::Vector{Int64},Q::Int64,b::Array{Int64,2},
    c::Array{Int64,2},s::Array{Int64,2},d::Array{Float64,2},a::Float64,P::Int64)::Vector{Vector{Int64},3}

    res = Vector{Vector{Int64}}

end