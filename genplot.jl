using Plots	
using DelimitedFiles
using Dates
gr()
theme(:juno)

function generateAnimation(W, H)
	data1 = readdlm("./positiondata.csv", ',', Float64);
	date = Dates.format(Dates.now(), "dd-mm-yy-HHMMSS")
	begin
		anim = @animate for i = 1:size(data1)[1]
			x, y = data1[i,1:2:end], data1[i, 2:2:end] 
			# x = x .- W.*round.(x./W);
			# y = y .- H.*round.(y./H);
			plot(x ,y ,seriestype=:scatter, xlims=[-W/2,W/2], ylims=[-H/2,H/2],
			 legend=nothing, axis=nothing, markerstrokewidth=0,c=:royalblue,
			  aspectratio=:equal)

		end
		gif(anim, "Results/Result$(date).gif", fps = 30)
	end
end

function simulateandplot(N=200, n=100, W=1200, H=1200)
	t1 = @elapsed run(`g++ pf.cpp -o pf`);
	t2 = @elapsed run(`./pf $(N) $(n)`)
	println("Time taken to simulate: $(t1+t2) seconds\n");
	t3 = @elapsed plt = generateAnimation(W, H)
	println("Time taken to plot: $(t3) seconds\n");
	plt
end


function tempfunction(W=1200, H=1200)
	data1 = readdlm("./positiondata.csv", ',', Float64);
	x, y = data1[1, 1:2:end], data1[1, 2:2:end]
	plot(x, y, legend=:none, color=:blue,
		 seriestype=:scatter,
		 xlims=[-W/2,W/2], ylims=[-H/2,H/2]) 
end

function tempplot(N=200, n=199, W=1200, H=1200)
	t1 = @elapsed run(`g++ pf.cpp -o pf`);
	t2 = @elapsed run(`./pf $(N) $(n)`)
	# println("Time taken to simulate: $(t1+t2) seconds\n");
	# t3 = @elapsed plt = generateAnimation(W, H)
	# println("Time taken to plot: $(t3) seconds\n");
	# plt
	#temp 
	tempfunction()
end

# edge detection.
posdata = readdlm("./positiondata.csv", ',', Float64);
x, y = posdata[1, 1:2:end], posdata[1, 2:2:end]
plot(x, y, seriestype=:scatter, markerstrokewidth=0)
struct VEC2
	x::Float64 
	y::Float64 
end
mean(x) = sum(x)/length(x)
cells = VEC2.(x, y)
rcm = VEC2(mean(x), mean(y))
plot!([rcm.x], [rcm.y], color=:blue, markersize=5, seriestype=:scatter, markerstrokewidth=0)
rft = cells[1]
dist(a::VEC2, b::VEC2) = √((b.x-a.x)^2+(b.y-a.y)^2)

f(b::VEC2) = dist(rcm, b)
dists = f.(cells)
i = argmax(dists)
rft = cells[i]

plot!([rft.x], [rft.y], color=:orange, markersize=5, seriestype=:scatter, markerstrokewidth=0)
import Base.-
-(a::VEC2, b::VEC2) = VEC2(a.x-b.x, a.y-b.y)
rcmft = rft-rcm


function plotcell(A::VEC2)
	scatter!([A.x], [A.y], color=:purple,
	markerstrokewidth=0,legend=:none,
	markerstrokesize=5)
end

plot(x, y, seriestype=:scatter, markerstrokewidth=0, aspect_ratio=1)
scatter!([rcm.x], [rcm.y], color=:blue, markersize=5, markerstrokewidth=0)
scatter!([cells[i].x], [cells[i].y], color=:orange, markerstrokewidth=0, markersize=5)
plot!([rcm.x, cells[i].x], [rcm.y, cells[i].y], arrow=true)

#test
dot(A::VEC2, B::VEC2) = A.x*B.x + A.y*B.y 
rt = 300.
rcmft = rft - rcm
for i in 1:length(cells)
	if (cells[i]!=rft)
		d = dist(cells[i], rft)
		if (d < rt)
			rftcn = cells[i]-rft #distance from farthest cell to candidate
			dp = dot(rftcn, rcmft) #dot product
			println("($(dp), $i)")
		end
	end
end

rt = 400.
rcur = rft 
boundcells = VEC2[]
for i in 1:10

	candidates = VEC2[]
	for j in 1:lastindex(cells)
		if (cells[j]!=rcur)
			if (dist(cells[j], rft) < rt)
				if !(cells[j] in boundcells) 
					push!(candidates, cells[j])
				end
			end
		end
	end

	dots = [dot(candi-rcur, rcur-rcm) for candi in candidates]
	dotsmag = abs.(dots)
	if isempty(dotsmag)
		break
	end
	mindotindex = argmin(dotsmag)
	rcur = candidates[mindotindex]
	push!(boundcells, rcur)
	println(rcur)

end

plot(x, y, seriestype=:scatter, markerstrokewidth=0, aspect_ratio=1)
plotcell.(boundcells)
scatter!([rcm.x], [rcm.y], color=:orange, markerstrokewidth=0, markersize=5)

