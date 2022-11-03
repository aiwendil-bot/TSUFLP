#=
tabu:
- Julia version: 
- Author: sujin
- Date: 2022-11-03
=#

# Paramètres
# taille de la mémoire tabu
LT = 7
# nombre d'itération sans amélioration autorisé (0.5*m avec m le nombre de concentrateurs)
k = 0.5
# budget en ms (ou nombre d'itération max)
budget = 100000

# durée de la mémoire tabu
tenure =

# Mémoire courte: solutions tabu accompagnées de leur tenurs
TM = []

# Voisinage: solutions non tabu
N = []

function tabu(obj::Int8, J, K)
    while time < budget
        # les candidats sont les concentrateurs qui n'ont pas été ouverts et ceux qui n'ont pas été explorer par TABU
        N = obj == 1 ? [k for k in vcat(J,K) if (!(k in sol[3]) and !(k in sol[4]) and !(k in TM))]
        : [k for k in J if (!(k in sol[3]) and !(k in TM))]



