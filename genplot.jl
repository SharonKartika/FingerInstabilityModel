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
			# plot(x ,y ,seriestype=:scatter, xlims=[-W/2,W/2], ylims=[-H/2,H/2],
			#  legend=nothing, axis=nothing, markerstrokewidth=0,c=:royalblue,
			#   aspectratio=:equal)
			plot(x ,y ,seriestype=:scatter,
				legend=nothing, axis=nothing,
				markerstrokewidth=0,c=:royalblue,
				aspectratio=:equal)
		end
		gif(anim, "Results/Result$(date).gif", fps = 30)
	end
end

function simulateandplot(;N=200, n=100, W=1200, H=1200)
	t1 = @elapsed run(`g++ pfbound.cpp -o pf`);
	t2 = @elapsed run(`./pfbound $(N) $(n)`)
	println("Time taken to simulate: $(t1+t2) seconds\n");
	t3 = @elapsed plt = generateAnimation(W, H)
	println("Time taken to plot: $(t3) seconds\n");
	plt
end

#=Haven't figured out what this is
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
=#

