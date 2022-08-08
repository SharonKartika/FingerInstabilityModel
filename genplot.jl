using Plots	
using DelimitedFiles
gr(legend=false, aspectratio=:equal)

function generateAnimation(W, H)
	data1 = readdlm("./positiondata.csv", ',', Float64);
	begin
		anim = @animate for i = 1:size(data1)[1]
			x, y = data1[i,1:2:end], data1[i, 2:2:end] 
			# x = x .- W.*round.(x./W);
			# y = y .- H.*round.(y./H);
			plot(x ,y ,seriestype=:scatter, xlims=[-W/2,W/2], ylims=[-H/2,H/2],
			 legend=nothing, axis=nothing, markerstrokewidth=0,c=:royalblue)
		end
		gif(anim, "Results/Result$(rand()).gif", fps = 30)
	end
end

function simulateandplot(N=100, n=300, W=1200, H=1200)
	t = @elapsed run(`./p $(N) $(n)`)
	println("Time taken: $(t) seconds\n");
	generateAnimation(W, H)
end