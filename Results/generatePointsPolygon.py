# %%
import numpy as np
from shapely.geometry import Point, Polygon
from polygenerator import (
    random_polygon,
    random_star_shaped_polygon,
    random_convex_polygon,
)
import matplotlib.pyplot as plt
import random
import sys


# %%
def plot_polygon(polygon, out_file_name):
    plt.figure()
    plt.gca().set_aspect("equal")

    for i, (x, y) in enumerate(polygon):
        plt.text(x, y, str(i), horizontalalignment="center", verticalalignment="center")

    # just so that it is plotted as closed polygon
    polygon.append(polygon[0])

    xs, ys = zip(*polygon)
    plt.plot(xs, ys, "r-", linewidth=0.4)

    plt.savefig(out_file_name, dpi=300)
    plt.close()


def smoothen(polygon):
    # polygon is a list of tuples
    newpoly = []
    for i in range(len(polygon) - 1):
        qx = (3 / 4) * polygon[i][0] + (1 / 4) * polygon[i + 1][0]
        qy = (3 / 4) * polygon[i][1] + (1 / 4) * polygon[i + 1][1]
        rx = (1 / 4) * polygon[i][0] + (3 / 4) * polygon[i + 1][0]
        ry = (1 / 4) * polygon[i][1] + (3 / 4) * polygon[i + 1][1]
        newpoly += [(qx, qy)]
        newpoly += [(rx, ry)]
    qx = (3 / 4) * polygon[-1][0] + (1 / 4) * polygon[0][0]
    qy = (3 / 4) * polygon[-1][1] + (1 / 4) * polygon[0][1]
    rx = (1 / 4) * polygon[-1][0] + (3 / 4) * polygon[0][0]
    ry = (1 / 4) * polygon[-1][1] + (3 / 4) * polygon[0][1]
    newpoly += [(qx, qy)]
    newpoly += [(rx, ry)]
    return newpoly


# %%
def generateReasonablePolygon(n):
    return smoothen(smoothen(random_polygon(n)))


# %%


def generatePointsWithinPolygon(polygon, n):
    polygon = Polygon(polygon)
    minx, miny, maxx, maxy = polygon.bounds
    pointsinpoly = []
    x = []
    y = []
    # for i in range(n):
    while len(x) < n:
        p = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if p.within(polygon):
            pointsinpoly += [(p.x, p.y)]
            x += [p.x]
            y += [p.y]
    return x, y


# %%
def convertrange(ri, x1, x2, y1, y2):
    runit = (ri - x1) / (x2 - x1)
    rf = runit * (y2 - y1) + y1
    return rf

convertrangefast = np.vectorize(convertrange)

#%%

def main():
    nelts = len(sys.argv)
    W, H = float(sys.argv[1]), float(sys.argv[2])
    npoints = int(sys.argv[3])
    nsample = int(sys.argv[4]) 

    x, y = generatePointsWithinPolygon(generateReasonablePolygon(npoints), nsample)
    x, y = np.array(x), np.array(y)
    x = convertrangefast(x, 0, 1, -W, W)
    y = convertrangefast(y, 0, 1, -H, H)
    
    X = np.array([x, y]).T
    np.savetxt("pointsinshape.csv", X, delimiter=" ")
    # plt.scatter(x, y)


# %%
if __name__ == "__main__":
    main()
