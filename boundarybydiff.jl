using Plots
using DelimitedFiles
using Dates
gr(aspect_ratio=:equal)
theme(:dark)

struct VEC2
    x::Float64
    y::Float64
end

function plotcell(A::VEC2, color=:purple)
    scatter!([A.x], [A.y], color=color,
        markerstrokewidth=0, legend=:none,
        markerstrokesize=5)
end

begin
X = readdlm("boundarypositiondata.csv",',')
X = X[1, :]
x = X[1:2:end]
y = X[2:2:end]
boundarycells = VEC2.(x, y)

X = readdlm("positiondata.csv", ',')
x = X[1:2:end]
y = X[2:2:end]
cells = VEC2.(x, y)
end

begin
plot()
plotcell.(cells, :orange)[end] #plot the cells
plotcell.(boundarycells, :cyan)[end] #plot the cells
end

