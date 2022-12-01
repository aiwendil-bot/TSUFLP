mutable struct Point
    x
    y
    Point() = new()
    Point(a,b) = new(a,b)
end
#=
struct Point{T<:Real}
 x::T
 y::T
 Point() = new()
end
=#
#=
struct Point
 x::Int64
 y::Int64
 Point() = new(1,2)
end
=#

#=
mutable struct SkipList{T}
 #Hauteur max de la skiplist
 etage_max::Int64
 # Contenu des étages
 etages::Vector{Vector{T}}
 SkipList(etag_max) = new(etag_max, [[]])
end
=#

# Structure SkipListe
mutable struct SkipList
 #Hauteur max de la skiplist
 etage_max::Int64
 
 # Contenu des étages
 etages::Vector{Vector{Point}}
 
 # Constructeur : il faut préciser la hauteur maximal (étage max)
 SkipList(etag_max) = new(etag_max, [[Point(-Inf,-Inf),Point(+Inf,+Inf)]])
end

# Défini si elem1 > elem2 selon l'ordre définit sur les elements
function pt_est_superieur(pt1::Point, pt2::Point)
 if (pt1.x == pt2.x)
     return (pt1.y > pt2.y)
 else
     return (pt1.x > pt2.x)
 end
end

# Défini si elem1 dominé elem2 selon la dominance définit sur les elements
function pt_domine(pt1::Point, pt2::Point)
 return (
             ((pt1.x < pt2.x) && (pt1.y <= pt2.y))
         || 	((pt1.x <= pt2.x) && (pt1.y < pt2.y))
        )
end

function creer_point_etage!(SL, etage, new_point, pos)
 SL_ground = SL.etages[1]
 # Si l'étage n'existe pas on le créer
 if length(SL.etages) == etage - 1
     println("L'étage n'existe pas on le créer")
     push!(SL.etages, [Point(-Inf,-Inf),Point(+Inf,+Inf)])
     # on insère un vecteur de point null correspondant à l'étage 1
     v = [Point() for i in 2:(length(SL_ground)-1)]
     SL.etages[etage] = vcat(SL.etages[etage][1], v, SL.etages[etage][end])
 end
 println("on ajoute le point ", new_point, " en position ", pos)
 # on ajoute le point en position pos
 SL.etages[etage][pos] = new_point
end

#=
head = SL_ground[2:pos-1]
tail = SL_ground[pos:end-1]
=#
# insere une colonne de points null en position pos
function inserer_colonne_null!(SL, pos)
 for etage in SL.etages[2:end]
     insert!(etage, pos, Point())
 end
end

# supprime la colonne position pos
function supprimer_colonne!(SL, pos)
 for etage in SL.etages[2:end]
     if isdefined(etage[pos], 1)
         etage[pos] = Point()
     else
         break
     end
 end
end

function inserer_new_elem!(SL, new_point, pos, pt_domine)
 SL_ground = SL.etages[1]
 println("----------------------")
 println("fonction d'insertion du nouvel element ", new_point)
 println("pos = ", pos)
 println("taille ground = ", length(SL_ground))
 domine = false
 
 if pos < length(SL_ground)
     if pt_domine(new_point, SL_ground[pos])
         domine = true
     end
 end
 
 # Si le nouveau point domine son prochain, on remplace ce dernier
 if domine
     #SL.etages[1][pos+1] = new_point
     # On met à jour les étages supérieur (on écrase la colonne précédente avec une colonne "vide")
     println("le point ", new_point, " domine son prochain ", SL_ground[pos], ". On le remplace lui et toute sa colonne")
     SL_ground[pos] = new_point
     supprimer_colonne!(SL, pos)
     for i in 1:length(SL.etages)
         println(SL.etages[i])
     end
 else
     # Sinon on décale les élément de 1 case vers la gauche avant d'insere le nouvel élément à pos+1
     println("on insere ", new_point, " en position ", pos, " dans ", SL_ground)
     insert!(SL_ground, pos, new_point)
     # On met à jour les étages supérieur (insertion d'une colonne "vide")
     for i in 1:length(SL.etages)
         println(SL.etages[i])
     end
     println("le point ", new_point, " a été inseré. Il faut creer une colonne de point vide au dessus de sa tete")
     inserer_colonne_null!(SL, pos)
     println("apres creation de la colonne de points vide")
     for i in 1:length(SL.etages)
         println(SL.etages[i])
     end
 end
 
 # On lance une piece (flip a coin)
 piece = rand((1,0))
 
 # Tant qu'on tombe sur face, on ajoute new_point dans l'étage supérieur
 piece = 1
 etage = 1
 while (piece == 1) && etage <= SL.etage_max
     etage += 1
     println("Face! on ajoute le point ", new_point, " a l'étage superieur (", etage, ")")
     creer_point_etage!(SL, etage, new_point, pos)
     println("après creation du point a l'étage superieur:")
     for i in 1:length(SL.etages)
         println(SL.etages[i])
     end
     piece = rand((1,0))
 end
end

# Calcule de la position d'insertion d'un nouvel élément pour la skip list
#	> s'il est dominé, on retourne 0 
function position_new_elem(SL, new_elem, est_superieur, domine)

 pos_insertion = 1
 
 # On parcours les étages de haut en bas
 for num_etage in length(SL.etages):-1:1
 #for etage in reverse(SL.etages)
     etage = SL.etages[num_etage]
     # On se positionne au bon endroit dans l'étage (selon l'ordre des elements)
     for pos_test in pos_insertion+1:length(etage)-1
     #for elem in etage[pos_insertion:end-1]
         elem = etage[pos_test]
         # Si le point est non défini, on le saute
         if isdefined(elem, 1)
             println("(etage ", num_etage,") ", elem, " domine-t-il le nouveau point ", new_elem, " ?")
             if domine(elem, new_elem)
                     println(elem, " domine ", new_elem)
                     return 0
             else
                 if est_superieur(new_elem, elem)
                     println(new_elem, " est superieur a ", elem)
                     pos_insertion = pos_test+1
                 else
                     println(new_elem, " est inferieur a ", elem)
                     break
                 end
             end
         else
             pos_test += 1
             println("etage :", num_etage, " element en position : ", pos_test, " (pos_test) non definit, pos_test = ", pos_test)
         end
     end
 end
 if pos_insertion == 1
     pos_insertion = 2
 end
 return pos_insertion
end
 
# Insère l'élément élèm dans la SL
function inserer_point_SL(SL, new_point, pt_est_superieur, pt_domine)
 # On calcule la position d'insertion pour le nouvel élément
 pos_inser = position_new_elem(SL, new_point, pt_est_superieur, pt_domine)
 
 # Si le point est non dominé, on l'insere
 if pos_inser > 0
     inserer_new_elem!(SL, new_point, pos_inser, pt_domine)
 end
end

function supprimer_elem_SL()
end

points = [Point(1,2), Point(3,1), Point(5,0), Point(6,3), Point(1,1)]
SL = SkipList(5)
for point in points
 inserer_point_SL(SL, point, pt_est_superieur, pt_domine)
 println("Apres insertion:")
 for i in 1:length(SL.etages)
     println(SL.etages[i])
 end
 print("\n\n")
end