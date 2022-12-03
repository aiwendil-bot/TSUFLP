#----- Définition des structures -------
struct Solution
    f_lvl_concentrators::Vector{Int64}
    s_lvl_concentrators::Vector{Int64}
    locations::Vector{Int64}
end

mutable struct Point
    x::Float64
    y::Float64
    Point() = new()
    Point(a,b) = new(a,b)
end

mutable struct Elem
	point::Point
	sol::Solution
	Elem() = new()
	Elem(pt) = new(pt)
	Elem(pt,s) = new(pt,s)
end

mutable struct Node
	elem::Elem
	etage::Int64
	prev::Node
	next::Node
	up::Node
	down::Node
	
	# Noeud pointant vers le dernier noeud de l'etage
	last::Node
	
	# Noeud correspondant à l'étage 1 (pour y sauter directement)
	ground::Node
	
	Node() = new()
	
	# Constructeur nouvel element (par defaut ground = lui meme)
	Node(el,et) = (x = new(el,et) ; x.prev = Node() ; x.next = Node() ; x.up = Node() ; x.down = Node() ; x.ground = x)
	
	# Constructeur pour insertion
	Node(el,et,p,n,u,d) = new(el,et,p,n,u,d)
	
end

# Structure SkipListe
mutable struct SkipList
 #Hauteur max de la skiplist
 hauteur_max::Int64
 #Hauteur: nombre d'étage de la SL
 hauteur::Int64
 #Largeur: nombre de noeud à l'étage 0
 largeur::Int64
 
 # Noeud d'entrée de la SL: premier noeud situé le plus en haut
 entree::Node
 
 # Noeud curseur: position du curseur de parcours
 curseur::Node
 
 # Constructeur : il faut préciser la hauteur maximal (étage max)
 SkipList(h_max) = new(h_max)
end

#----- Operateurs de comparaisons des points -----
function egaux(pt1::Point, pt2::Point)
	return (pt1.x == pt2.x) && (pt1.y == pt2.y)
end

# Défini si elem1 > elem2 selon l'ordre définit sur les elements
function est_superieur(pt1::Point, pt2::Point)
 if (pt1.x == pt2.x)
     return (pt1.y >= pt2.y)
 else
     return (pt1.x >= pt2.x)
 end
end

# Défini si elem1 dominé elem2 selon la dominance définit sur les elements
function domine(pt1::Point, pt2::Point)
  return (
 	     pt1.x != -Inf && pt2.x != Inf &&
             (
             ((pt1.x < pt2.x) && (pt1.y <= pt2.y))
         ||  ((pt1.x <= pt2.x) && (pt1.y < pt2.y))
         ||  ((pt1.x == pt2.x) && (pt1.y == pt2.y))
             )
        )
end

#----- Méthodes de la SL -------

# Retire un élément de la SL
function remove_elem!(curseur2)
# o1 <-> o2 <-> o3
# se placer sur o1
# o1.next = o3
# o1.next.prev = o1
	curseur2 = curseur2.prev
	
	liste_t = curseur2.next.next
	curseur2.next = liste_t
	curseur2.next.prev = curseur2
end

function filtre_pts_suivant_domines(search_point, SL)
	while !egaux(SL.curseur.elem.point, Point(Inf,Inf))
		if domine(search_point, SL.curseur.elem.point)
			curseur_2 = SL.curseur
			# On retire l'élément à chaque etage (reboutage des pointeurs)
			remove_elem!(curseur2)
			while (curseur2 != curseur2.ground)
				curseur2 = curseur2.down
				remove_elem!(SL)
			end
		end
	end
end

# search_point, curs1.elem.point
# Recherche un élément:
# - positionne le curseur immédiatement avant la position d'insertion
# - retourne true si l'élément est déjà présent ou dominé
function SL_search!(SL, search_point, etage_fin = 1)
	trouve = false
	
	SL.curseur = SL.entree
	continuer = true
	# On parcours les étages de haut en bas
	while continuer
		# On se positionne au bon endroit dans l'étage (selon l'ordre des elements)
		while est_superieur(search_point, SL.curseur.elem.point)
			# Si le point recherché est dominé, on ne l'insere pas
			if domine(SL.curseur.elem.point, search_point)
				if verbose[] println(SL.curseur.elem.point, " etage ", SL.curseur.etage, " domine ", search_point, " etage ", SL.curseur.etage) end
				return true
			else
				SL.curseur = SL.curseur.next
			end
		end

		SL.curseur = SL.curseur.prev
		
		# On descend d'un etage
		if isdefined(SL.curseur.down,1) && (SL.curseur.etage > etage_fin)
			SL.curseur = SL.curseur.down
		else
			continuer = false
		end
	end
	
	return trouve
end

function creer_SL(SL_hauteur::Int64)
	SL = SkipList(5)
	SL.entree = Node(Elem(Point(-Inf,-Inf)),1)
	SL.entree.next = Node(Elem(Point(Inf,Inf)),1)
	SL.entree.next.prev = SL.entree
	
	SL.entree.last = SL.entree.next
	SL.hauteur = 1
	SL.curseur = SL.entree
	return SL
end

function SL_new_floor!(SL)
	num_etage = SL.entree.etage
	
	SL.entree.up = Node(Elem(Point(-Inf,-Inf)), num_etage + 1)
	SL.entree.up.next = Node(Elem(Point(Inf,Inf)), num_etage + 1)
	SL.entree.up.next.prev = SL.entree.up
	SL.entree.up.down = SL.entree
	SL.entree.up.next.down = SL.entree.last
	SL.entree.last.up = SL.entree.up.next
	SL.entree.up.last = SL.entree.up.next
	
	SL.hauteur += 1
	SL.entree = SL.entree.up
end

function afficher_SL(SL)
	curs = SL.entree
	s1 = ""
	continuer = true
	while continuer
		curs = curs.last
		s2 = string(curs.elem.point)
		while isdefined(curs.prev,1)
			curs = curs.prev
			s2 = string(curs.elem.point) * " -> " * s2
		end
		s1 *= s2 * "\n"
		if isdefined(curs.down,1)
			curs = curs.down
		else
			continuer = false
		end
	end
	println(s1)
end

# Weighted coin: p compris entre 0 et 1
function coin(p)
	return rand() < p
end
	
function SL_insert_up!(SL, new_elem, proba)
	# On lance une piece (flip a coin)
	piece = coin(proba)
	
	# Tant qu'on tombe sur face, et qu'on n'a pas atteint le nombre max d'étage, on ajoute new_elem.point dans l'étage supérieur
	#piece = 1
	while (piece) && (SL.curseur.etage + 1 <= SL.hauteur_max)
		etage = SL.curseur.etage + 1
		if verbose[] println("Face! on ajoute le point ", new_elem.point, " a l'étage superieur (", SL.curseur.etage + 1, ")") end
		# Si l'étage n'existe on le créer
		if (SL.curseur.etage + 1 > SL.hauteur)
			if (SL.hauteur + 1 <= SL.hauteur_max)
				if verbose[] println("L'étage n'existe pas on le créer") end
				SL_new_floor!(SL)
				prev = SL.entree
				suiv = SL.entree.last
			else
				return 1 # Si on a atteint le nombre max d'étage, on sort
			end
		else
			# Sinon on recupere les pointeurs prev et suiv du nouveau point
			curs_elem = SL.curseur
			if !(SL_search!(SL, new_elem.point, SL.curseur.etage + 1))
				prev = SL.curseur
				suiv = prev.next
				
				if verbose[] println("Le point ", new_elem.point, " est > a ", prev.elem, " (etage ", prev.etage, ") et < a ", suiv.elem, " (etage ", suiv.etage," )") end
				etage = prev.etage
				# on recupere le curseur qui pointait sur le nouvel element
				SL.curseur = curs_elem
			else
				if verbose[] println("err: point non trouve pour l'insert_up") end
				return 0
			end
		end
		
		if verbose[] println("On insere le point ", new_elem.point, " au dessus de ", SL.curseur.elem, " (etage ", SL.curseur.etage, ")", " a l'etage ", etage, " apres ", prev.elem, " (etage ", prev.etage, ") et avant ", suiv.elem, " (etage ", suiv.etage, ")") end
		SL.curseur.up = Node(Elem(new_elem.point), etage, prev, suiv, Node(), SL.curseur)
		prev.next = SL.curseur.up
		suiv.prev = SL.curseur.up
		
		if verbose[] println("après creation du point a l'étage superieur:") ; afficher_SL(SL) end
		
		SL.curseur = SL.curseur.up
		piece = coin(proba)
	end
	return 1
end

# On ne peut supprimer un élément qu'une fois qu'on s'est assuré qu'il n'est dominé par aucun élément du premier etage
# Une fois qu'on en est assuré, on enleve les noeuds dominé lors du parcours des étages depuis le premier jusqu'au dernier 
function SL_remove!(SL, curs1, curs2)
	curs1.next = curs2
	curs2.prev = curs1
	while curs1.etage < SL.hauteur
		while !(isdefined(curs1.up,1))
			if !(isdefined(curs1.prev,1))
				break # curs1 = -Inf
			else
				curs1 = curs1.prev
			end
		end
		
		while !(isdefined(curs2.up,1))
			if !(isdefined(curs2.next,1))
				break # curs2 = +Inf
			else
				curs2 = curs2.next
			end
		end
		
		curs1 = curs1.up
		curs2 = curs2.up
		
		if curs1.next != curs2
			curs1.next = curs2
		end
		if curs2.prev != curs1
			curs2.prev = curs1
		end
	end
end
	
function SL_insert!(SL, new_elem, proba)
	# Si le nouveau point n'existe pas et qu'il n'est pas dominé (le curseur aura ete positionne a l'élément immédiatement inférieur)
	SL.curseur = SL.entree
	if !(SL_search!(SL, new_elem.point))
		# on commence par retirer de la SL tous les points dominés par new_elem.point
		curs1 = SL.curseur
		curs2 = curs1.next
		
		while domine(new_elem.point, curs2.elem.point)
			curs2 = curs2.next
		end
		
		if verbose[] println("Le point ", new_elem.point, " domine tous les points entre ", curs1.elem.point, " a ", curs2.elem.point) end
		SL_remove!(SL, curs1, curs2)
		
		# puis on insere le nouveau point
		if verbose[] println("Le point ", new_elem.point, " est > a ", SL.curseur.elem.point, " (etage ", SL.curseur.etage, ") et < a ", SL.curseur.next.elem.point, " (etage ", SL.curseur.next.etage," )") end
		liste_t = SL.curseur.next
		# Node(element, prev, next, up, down)
		liste_t.prev = Node(Elem(new_elem.point), SL.curseur.etage, SL.curseur, liste_t, Node(), Node())
		
		SL.curseur.next = liste_t.prev
		
		SL.curseur = SL.curseur.next
		
		if verbose[] println("Apres insertion du point et avant insert up") ; afficher_SL(SL) end
		ins_status = SL_insert_up!(SL, new_elem, proba)
		
		ret = true
	else
		ret = false
	end
	
	return ret
end

function main(args)
	if length(args) >= 1
		verbose[] = args[1] == "-v"
		verbose0[] = args[1] == "-v0"
	end
	
	# -------------------------------
	# - Paramètrage de la Skip List -
	# -------------------------------
	# Hauteur max
	SL_hauteur = 10
	
	# Probabilité de création d'un nouvel étage [0,1] (poids de la piece pour le flip coin)
	proba = 0.5
	
	# -------------------------------
	# - Création de la Skip List -
	# -------------------------------
	SL = creer_SL(SL_hauteur)
	
	# -----------------------------------------
	# - Insertion de points dans la Skip List -
	# -----------------------------------------
	# liste des points à insérer (dans l'exemple ici, points[7] domine {points[2], points[3], points[4]} il seront donc filtrer à l'ajout de points[7] dans la Skip List)
	points = [Point(1,8), Point(3,7), Point(5,6), Point(6,5), Point(7,4), Point(8,3), Point(3,5)]
	sols = [Solution([2], [1,2,3], [1,2,3,4,5,6]) for i in 1:length(points)]
	elements = [Elem(points[i],sols[i]) for i in 1:length(points)]
	
	for elem in elements
		# insertion de chaque elements
		inser_ok = SL_insert!(SL, elem, proba)
	end
	
	#SL_insert!(SL, points[7], proba)
	
	if verbose0[] afficher_SL(SL) end
end

# Niveau de verbosité du script:
# - pas d'option = pas d'affichage
# - "-v" affichage de la SL à la fin des insertions seulement
# - "-v" max
const verbose = Ref(false)
const verbose0 = Ref(false)
main(ARGS)
