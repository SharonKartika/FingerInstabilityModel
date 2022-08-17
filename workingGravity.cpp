#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

#include "utilityfunctions.cpp"

int N{100}; // #agents
int n{500}; // #time steps
constexpr float W{1200.};
constexpr float H{1200.};
constexpr float dt = 0.01;
constexpr float beta = 60;
float w2, h2;

float unitrand();
float map(float, float, float, float, float);
float Hv(float);
float noisemag(float);


class CELL
{
public:
    float x, y,
        vx, vy,
        ax, ay;
    CELL()
    {
        x = randf(-w2, w2);
        y = randf(-h2, h2);
        vx = 0.;
        vy = 0.;
        ax = ay = 0.;
    }
    void update()
    {
        vx += (ax * dt);
        vy += (ay * dt);
        x += (vx * dt);
        y += (vy * dt);
    }
};

void initaccn(CELL M[])
{ /*
    Sets accn of all cells to zero.
  */
    for (int i = 0; i < N; i++)
    {
        M[i].ax = 0.;
        M[i].ay = 0.;
    }
}

void writeposition(CELL M[], std::ofstream &file)
{ /*
  Writes the coordinates of all the cells to file
  */
    for (int i = 0; i < N - 1; i++)
        file << M[i].x << ',' << M[i].y << ',';
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}


float forcemag(float r)
{ /*
    magnitude of the force
  */
    return 1e4 / r * r;
}

void setaccn(CELL &A, CELL &B)
{ /*
    Sets the acceleration of the cells passed
 */
    float dx = B.x - A.x;
    float dy = B.y - A.y;
    float r = sqrt(dx * dx + dy * dy);
    float f = forcemag(r);
    float fx = f * (dx / r);
    float fy = f * (dy / r);
    A.ax += fx;
    A.ay += fy;
}

void looploop(CELL M[])
{
    initaccn(M);
    // find acceleration
    for (int i = 0; i < N; i++)
    { // every particle
        for (int j = 0; j < N; j++)
        { // every particle
            if (i != j)
            {
                setaccn(M[i], M[j]);
            }
        }
    }
    for (int i = 0; i < N; i++)
    { //position update loop
        M[i].update();
    }
}

int main(int argc, char *argv[]){
    N = std::stoi(argv[1]);
    n = std::stoi(argv[2]);
    w2 = W/2;
    h2 = H/2;
    CELL M[N];

    srand(1);
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    for (int t = 0; t < n; t++){
        looploop(M);
        writeposition(M, posfile);
    }
    std::cout << "simulation complete\n";
}