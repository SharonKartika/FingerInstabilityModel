#manage posdata (WONT RUN)
file = open("GeneratedData/dumbellPositionDataRaw.csv", "r")

lines = readlines(file)
x = parse.(Int, lines[1:6:end])
y = parse.(Int, lines[3:6:end])

writedlm("dumbellPosition.csv",[x y], ',')

#for testing different cell distributions
#(WONT RUN)
X = readdlm("GeneratedData/dumbellPosition.csv", ',')
x = X[:, 1]
y = X[:, 2]
