using Plots	
using DelimitedFiles
using Dates
gr(aspect_ratio=true)
theme(:dark)

posdata = readdlm("./positiondata.csv", ',', Float64);
x, y = posdata[1, 1:2:end], posdata[1, 2:2:end]

#for testing different cell distributions

X = readdlm("GeneratedData/dumbellPosition.csv", ',')
x = X[:, 1]
y = X[:, 2]


struct VEC2
	x::Float64 
	y::Float64 
end

mean(x) = sum(x)/length(x)
dist(a::VEC2, b::VEC2) = √((b.x-a.x)^2+(b.y-a.y)^2)
Base. -(a::VEC2, b::VEC2) = VEC2(a.x-b.x, a.y-b.y)
Base. +(a::VEC2, b::VEC2) = VEC2(a.x+b.x, a.y+b.y)
dot(A::VEC2, B::VEC2) = A.x*B.x + A.y*B.y 
mag(A::VEC2) = √(A.x^2 + A.y^2)
ang(A::VEC2, B::VEC2) = atan(B.y-A.y, B.x-A.x)

function plotcell(A::VEC2, color=:purple)
	scatter!([A.x], [A.y], color=color,
	markerstrokewidth=0,legend=:none,
	markerstrokesize=5)
end
"""Plot an arrow from A to B"""
function plotvec(A::VEC2, B::VEC2, color::Symbol=:yellow)
    plot!([A.x, B.x], [A.y, B.y], arrow=true, color=color, alpha=0.3)
end

cells = VEC2.(x, y)
rcm = VEC2(mean(x), mean(y))

# find the farthest cell 
f(b::VEC2) = dist(rcm, b)
dists = f.(cells)
i = argmax(dists)
rft = cells[i]
# rft = cells[rand(1:length(cells))]

#method 1 try; perpindicular to rcm
begin
	rt = 300.
	rcur = rft 
	boundcells = VEC2[]
	plot(x, y, seriestype=:scatter, markerstrokewidth=0, aspect_ratio=1, legend=:none)
	plotcell(rcur, :orange)
	# begin
	while true
		candidates = VEC2[]
		for j in 1:lastindex(cells)
			if (cells[j]!=rcur)
				if (dist(cells[j], rcur) < rt)
					if !(cells[j] in boundcells) 
						push!(candidates, cells[j])
					end
				end
			end
		end

		# dots = [dot(candi-rcur, rcur-rcm) for candi in candidates]
		dots = [dot(candi-rcur, rcur-rcm)/(mag(candi-rcur)*mag(rcur-rcm)) for candi in candidates]
		# dots = [ang(candi-rcur, rcur-rcm) for candi in candidates]
		
		# dotsmag = abs.(dots)
		dotsmag = dots
		if isempty(dotsmag)
			# error("empty dotsmag")
			break
		end
		# mindotindex = argmin(dotsmag)
		mindotindex = argmax(dotsmag)

		#tempdebug
		plotvec(rcur, candidates[mindotindex], :lightblue)
		plotvec(rcm, rcur, :yellow)

		rcur = candidates[mindotindex]
		
		push!(boundcells, rcur)
		println(rcur)
		plotcell(rcur)
		#condition
		if (rcur == rft)
			break
		end
	end
	plotcell(rcur, :blue)
end

#angle and distance of each cell wrt to the center of mass
begin
	angs = (x->ang(x, rcm)).(cells)
	dists = (x->dist(x, rcm)).(cells)
	#(-π,π)
	maximum(angs), minimum(angs)

	plot()
	plotcell.(cells, :orange)
	plotcell(cells[1], :orange)

	f(x, y) = VEC2(x^2, y^2) 
	f(A::VEC2) = f(A.x, A.y)
	plotcell.(f.(cells), :orange)
	plotcell(cells[1], :orange)

	plot()
	rθ(x, y) = VEC2(atan(y, x), x^2+y^2 )
	rθ(A::VEC2) = rθ(A.x, A.y)
	cycells = rθ.(cells)
	plt = plotcell.(cycells, :blue)
	plotcell(cycells[1], :blue)

	getx(A::VEC2) = A.x
	gety(A::VEC2) = A.y 
	rs = gety.(cycells)
	θs = getx.(cycells)
	argmax(rs)

	πr = LinRange(-π,π,15)
	V = [Vector{Int}(undef, 0) for i in 1:length(πr)-1]
	function whichπr(θ)
		i = 0
		for i in 1:length(πr)-1
			println(πr[i]," ", θ," ", πr[i+1])
			if πr[i] < θ < πr[i+1] 
				return i
			end
		end
		return i
	end	

	for i in 1:length(θs)
		push!(V[whichπr(θs[i])], i)
	end
	#V now contains the indices of cells in each bin 
	maxinVi(Vi) =  Vi[argmax((i->dists[i]).(Vi))]
	#boundarycellindices
	BCI = maxinVi.(V)

	plot()
	plotcell.(cycells, :blue)
	plotcell.(cycells[BCI], :orange)[end]

end

### local center of mass approach
plot() #clear
plotcell.(cells, :orange)[end] #plot the cells

#find mean (center of mass)
rcm = VEC2(mean((i->i.x).(cells)), mean((i->i.y).(cells)))
plt = plotcell(rcm, :blue) 

#distance from center of mass
dists = (x->dist(x, rcm)).(cells)
ift = argmax(dists) #index of furthest cell 
rft = cells[ift]
plotcell(rft, :red) #plot furthest cell

plotvec(rcm, rft)

#set threshold
rt = 30
#function that returns neighbors inside rt 
function getneighbors(cell::VEC2, rt)
	neighbs = VEC2[]
	for i in 1:length(cells)
		if (dist(cell, cells[i])<rt)
			if (cell!=cells[i])
				push!(neighbs, cells[i])
			end 
		end
	end
	return neighbs
end

rftneighbs = getneighbors(rft, rt)
plotcell.(rftneighbs, :purple)[end]

rlcm = VEC2(mean((i->i.x).(rftneighbs)), mean((i->i.y).(rftneighbs)))
plotcell(rlcm, :green)

plotvec(rlcm, rft)

#algo begin

# function angbetween(A::VEC2, B::VEC2)
# 	atan(B.y-A.y,B.x-A.x)
# end

function angbetween(rl::VEC2, rc::VEC2, rn::VEC2)
	θ1 = atan(rc.y-rl.y, rc.x-rl.x)
	θ2 = atan(rn.y-rc.y, rn.x-rc.x)
	return θ1 + π - θ2
end
rc = rft 

for i in 1:40
rcneighbs = getneighbors(rc, rt)
rncm = VEC2(mean((i->i.x).(rcneighbs)), mean((i->i.y).(rcneighbs)))
rcmc = rc - rncm # lcm to rc

scores = [angbetween(rncm, rc, rcni) for rcni in rcneighbs]
inext = argmax(scores)
#update the candidates
rc = rcneighbs[inext]

plotcell(rc, :green)
end
plotcell(rft, :red)


#=
set rc = farthest cell 
boundcells = []
for i in 1:10
	rns = neighbors(rc)
	rl = COM(rns)
	rns = rns not in boundcells
	rn = largestscore(rns)
	rc = rn 
	push!(bcells, rc)
	plot(rc)
end
=#
"""angle from origin to the vector, counterclockwise"""
natan(y, x) = (y>=0) ? (atan(y, x)) : (2π+atan(y, x))
natan(A::VEC2) = natan(A.y, A.x)
function getscore(rnsnb::Vector{VEC2}, rl::VEC2, rc::VEC2)
	function getscore(rn::VEC2, rl::VEC2, rc::VEC2)
		
	end
	getscore.(rnsnb, rl, rc)
end
rc = rft 
boundcells = VEC2[]
rt = 30 
#for 
	rns = getneighbors(rc, rt)
	rl = rncm = VEC2(mean((i->i.x).(rns)),
				     mean((i->i.y).(rns)))
	#neighbors excluding boundary
	rnsnb = setdiff(rns, boundcells)
	rn = getscore(rnsnb, rl, rc) 
#end