using Plots	
using DelimitedFiles
using Dates
gr(aspect_ratio=true)
theme(:dark)

begin
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
end

function plotcell(A::VEC2, color=:purple)
	scatter!([A.x], [A.y], color=color,
	markerstrokewidth=0,legend=:none,
	markerstrokesize=5)
end

function plotvec(A::VEC2, B::VEC2, color::Symbol=:yellow)
    plot!([A.x, B.x], [A.y, B.y], arrow=true, color=color, alpha=0.3)
end

function getneighbors(cell::VEC2, rt)
	neighbs = VEC2[]
	for i in 1:lastindex(cells)
		if (dist(cell, cells[i])<rt)
			if (cell!=cells[i])
				push!(neighbs, cells[i])
			end 
		end
	end
	return neighbs
end

# posdata = readdlm("./positiondata.csv", ',', Float64);
# x, y = posdata[1, 1:2:end], posdata[1, 2:2:end]


X = readdlm("GeneratedData/dumbellPosition.csv", ',')
x = X[:, 1]
y = X[:, 2]



cells = VEC2.(x, y)
rcm = VEC2(mean(x), mean(y))
dists = (b->dist(rcm, b)).(cells)

#farthest cell
i = argmax(dists)
rft = cells[i]

#= 
TERMINOLOGY
    rl <- local center of mass 
    rc <- current candidate cell 
    rn <- potential next candidate cell
=#

"""Returns the anticlockwise angles between vectors 
(rc-rl) and (rn-rc)"""
function angleAC(rl::VEC2, rc::VEC2, rn::VEC2)
    A = rc-rl
    B = rn-rc
    k = B.y * A.x - A.y * B.x
    cosθ = dot(A, B)/(mag(A) * mag(B))
    θ = acos(cosθ)
    if (k>0)
        return π - θ
    else 
        return 2π - θ
    end
end

#=
function angleAC1(rl::VEC2, rc::VEC2, rn::VEC2)
    A = rc-rl
    B = rn-rc
    k = B.y * A.x - A.y * B.x
    θ = atan(B.y-A.y, B.x-A.x)
    if (k<0)
        return π - θ
    else 
        return 2π - θ
    end
end
=#

# rl = VEC2(0, 0)
# rc = VEC2(0, 1)
# rn = VEC2(-1, 1)
# A = rc - rl
# B = rn - rc

# angleAC(rl, rc, rn)



plot() #clear
plotcell.(cells, :orange)[end] #plot the cells

rc = rft 
boundcells = VEC2[]
Rt = 30 

# begin
while true
	rns = getneighbors(rc, Rt)
	rl = VEC2(mean((i->i.x).(rns)),
			  mean((i->i.y).(rns)))
	# plotcell.(rns, :green) #temp
	#neighbors excluding boundary
	rnsnb = setdiff(rns, boundcells)
	# nscores = (A::VEC2->getscore(rl, rc, A)).(rnsnb) 
	nscores = (A::VEC2->angleAC(rl, rc, A)).(rnsnb) 
	
    #remove clockwise paths
	nscores[nscores.>=π] .= 0.
	rnindex = argmax(nscores)
	rn = rnsnb[rnindex]
	plotvec(rc, rn) #temp
	plotvec(rl, rc)
	rc = rn
	push!(boundcells, rc)
	plotcell(rc)[end]
    if (rc == rft)
        break;
    end
end
plotcell(rc)


function getBoundary(cells, Rt)
    #find farthest cell 
    begin
    rcm = VEC2(mean((i->i.x).(cells)),
			  mean((i->i.y).(cells)))
    dists = (b->dist(rcm, b)).(cells)
    i = argmax(dists)
    rft = cells[i]
    end
    
    rc = rft 
    boundcells = VEC2[]

    while true
        rns = getneighbors(rc, Rt)
        rl = VEC2(mean((i->i.x).(rns)),
                  mean((i->i.y).(rns)))
        rnsnb = setdiff(rns, boundcells)
        nscores = (A::VEC2->angleAC(rl, rc, A)).(rnsnb) 
        nscores[nscores.>=π] .= 0.
        rnindex = argmax(nscores)
        rn = rnsnb[rnindex]
        rc = rn
        push!(boundcells, rc)
        if (rc == rft)
            break;
        end
    end
    return boundcells
end

bounds = getBoundary(cells, 30)
plotcell.(bounds, :green)[end]

