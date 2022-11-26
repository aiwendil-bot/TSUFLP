using Plots; gr()
function ajout_term(x,y)
	scatter!((x, y);color=:black, legend = false, markersize=4)
end

function ajout_clvl1(x,y)
	scatter!((x, y);color=:black, legend = false, markersize=10)
end

function ajout_clvl1(x,y)
	scatter!((x, y);color=:black, legend = false, markersize=20)
end

a = plot(1,1)
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
readline()
