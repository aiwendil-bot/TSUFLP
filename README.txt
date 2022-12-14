Implémentation de Scatter Search multi-objectif (Sanchez, 2021) par BITHO Sullivan et CALLICO Adrien

organisation du dépôt :

data/ 

    -> les instances utilisées par Sanchez (small, medium, large)
    -> les datasets utilisés pour notre instance "angers"
    -> /instances_test regroupe les instances sur lesquelles nous avons testé
        notre implémentation et dont les résultats sont rapportés 

out/ 

    -> chaque sous dossier porte le nom de l'instance concernée
        -> les plots de grasp, 1er passage de tabu et les comparaison entre stratégies (crossover, crossover + tabu)
        -> les fronts de Pareto à chaque itération du test k = 0.9 ; crossover + tabu
        -> les fichiers XE_STRAT_k.txt : chaque ligne = une solution (affectation des terminaux / affectation des conc. de niveau 1)
        -> les fichiers YN_STRAT_k.txt : 1ère ligne : temps d'exec - card(YN) puis les valeurs de (z1,z2) pour chaque solution
        -> de même pour YN_vOpt.txt


src/ 

    -> les fichiers sources de Scatter Search 
    -> /!\ main lance les comparaisons entre les 3 stratégies pour chaque instance de data/instances_test pour k in [0.1,0.5,0.9] 
    

tests/ 

    fichiers pour tester notre implémentation sur différentes instances
    tous les résultats de chaque test sont sauvegardés dans out/nom_instance

    pour le paramètre nom_instance : donner le nom d'une des instances de tests/files 
    par exemple small8, medium4 etc (pas de chemin, ni d'extension à spécifier)

    Les valeurs des paramètres par défaut sont 
        Q = 7
        a = 0.4
        P = 300
        k = 0.5
        beta = 10
        crossover_on = false
        tabu_crossover_on = false

    ci-après comment les utiliser :

    SS_Sanchez

        -> dans un terminal, se placer dans TSUFLP/

        julia tests/SS_Sanchez.jl

        -> lance SS avec les paramètres par défaut sur small2
        
        julia tests/SS_Sanchez.jl noms_instance

        -> lance SS avec les paramètres par défaut sur nom_instance 

        julia tests/SS_Sanchez.jl noms_instance Q a P k beta crossover_on tabu_crossover_on

        -> lance SS avec les paramètres spécifiés sur nom_instance

        Q : capacité (Int)
        a : paramètre glouton (float)
        P : taille population de base (multiple de 5) (int)
        k : paramètre critère arrêt (float)
        beta : taille refSets (int pair)
        crossover_on : activer le crossover (bool)
        tabu_crossover_on : activer le crossover et la recherche tabu en plus (bool)

    SS_angers
        
        -> dans un terminal, se placer dans TSUFLP/

        julia tests/SS_angers.jl

        -> lance l'instance angers avec les paramètres par défaut et taille = [100,25,4]

        julia tests/SS_angers.jl nb_terms nb_clvl1 nb_clvl2 Q a P k beta crossover_on tabu_crossover_on

        -> lance l'instance angers sur les paramètres spécifiés (nb_terms <= 55800 ; nb_clvl1 <= 71 ; nb_clvl2 <= 4 )

    comparer_strategies

        -> dans un terminal, se placer dans TSUFLP/

        julia tests/comparer_strategies.jl
        
        -> lance la comparaison entre SS, SS + crossover, SS + crossover + tabu 
            sur small2.txt avec les paramètres par défaut

        julia tests/comparer_strategies.jl noms_instance Q a P k beta

        -> lance d'abord la résolution avec vOpt si l'instance n'a pas déjà été testée (peut être long sur les mediums)
        -> lance la comparaison entre les trois stratégies avec les paramètres spécifiés



    