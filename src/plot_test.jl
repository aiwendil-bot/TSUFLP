using Plots; gr()
function ajout_term(coord)
	scatter!((coord[1], coord[2]);color=:black, legend = false, markersize=4)
end

function ajout_clvl1(coord)
	scatter!((coord[1], coord[2]);color=:blue, legend = false, markersize=4)
end

function ajout_clvl2(coord)
	scatter!((coord[1], coord[2]);color=:green, legend = false, markersize=4)
end

function connect_points(pt1, pt2)
	plot!([pt1[1], pt2[1]], [pt1[2], pt2[2]])
end

function graph_sol(coord_t, coord_clvl1, coord_clvl2, sol)
	a = plot(1,1)
	#coord = coord_t[1,:]
	#ajout_term(coord)
	#ajout_term(coord_t[1,:])
	#ajout_clvl1(coord_clvl1[1,:])
	for i in eachindex(coord_t[:,1])
		ajout_term(coord_t[i,:])
	end
	for i in eachindex(coord_clvl1[:,1])
		ajout_clvl1(coord_clvl1[i,:])
	end
	for i in eachindex(coord_clvl2[:,1])
		ajout_clvl2(coord_clvl2[i,:])
	end
	
	#pt1 = coord_t[1,:]
	#pt2 = coord_clvl1[1,:]
	#connect_points(pt1,pt2)
	#display(a)
	#println("press")
	clvl1_aff_term = sol[1][1]
	i = 1
	#println("clv1 aff term")
	#println(clvl1_aff_term)
	#readline()
	for j in clvl1_aff_term
		pt1 = coord_t[i,:]
		pt2 = coord_clvl1[j,:]
		#println(pt1)
		#println(pt2)
		connect_points(pt1,pt2)
		#display(a)
		#println("press")
		#readline()
		i += 1
	end
	display(a)
	println("Affectation clvl1 - terminaux (press any key)")
	readline()

	clvl2_aff_clvl1 = sol[1][2]
	println(clvl2_aff_clvl1)
	readline()
	i = 1
	for j in clvl2_aff_clvl1
		if (j != 0)
			pt1 = coord_clvl1[i,:]
			pt2 = coord_clvl2[j,:]
			connect_points(pt1,pt2)
		end
		i += 1
	end
	display(a)
	println("Affectation clvl2 - clvl2 (press any key)")
	readline()

	#println("press any key to continue")
	#readline()
		#coord = coord_t[i,:]
	#for i in 1:size(coord_t,1)
	#	coord = coord_t[i,:]
	#	ajout_term(x,y)
	#end
	#=a = plot(1,1)
	for i in 1:10
		ajout_term(i,1)
	end
	for i in 2:2:10
		ajout_clvl1(i,5)
	end
	for i in 4:4:10
		ajout_clvl1(i,10)
	end
	display(a)
	readline()=#
end
