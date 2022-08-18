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
float w2, h2;

// function prototypes; definitions in utilityfunctions.cpp
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
        /*
        Constructor: initialises the cell.
        */
        x = randf(-w2, w2);
        y = randf(-h2, h2);
        vx = randf(-1000, 1000);
        vy = randf(-1000, 1000);
        ax = ay = 0.;
    }
    void update()
    {
        // applies the Euler method to update position
        vx += (ax * dt);
        vy += (ay * dt);
        x += (vx * dt);
        y += (vy * dt);
    }
};

void writeposition(CELL M[], std::ofstream &file)
{ /*
  Writes the coordinates of all the cells to file
  */
    for (int i = 0; i < N - 1; i++)
        file << M[i].x << ',' << M[i].y << ',';
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}
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

float forcemag(float r)
{ /*
    magnitude of the force
  */
    // return 1e4 / r * r;//gravity

    // lennard jones
    float eps = 21.;
    float sig = 78.;
    float force = 0.;
    force += 6 * pow(sig, 6) / pow(r, 7);
    force -= 12 * pow(sig, 12) / pow(r, 13);
    force *= 4 * eps;

    return force;
}

void setaccn(CELL &A, CELL &B)
{ /*
    Increments the acceleration of A due to B,
    and B due to A.
    */
    float dx = B.x - A.x;
    float dy = B.y - A.y;
    float r = sqrt(dx * dx + dy * dy);
    float f = forcemag(r);
    float fx = f * (dx / r);
    float fy = f * (dy / r);
    A.ax += fx;
    A.ay += fy;
    B.ax -= fx;
    B.ay -= fy;
}

void looploop(CELL M[])
{ /*
    Goes through every interacting pair.
    Finds acceleration.
    Also updates position.
   */
    initaccn(M);
    // find acceleration
    for (int i = 0; i < N - 1; i++)
    { // every particle except the last
        for (int j = i + 1; j < N; j++)
        { // every particle from i
            setaccn(M[i], M[j]);
        }
    }
    // position update loop
    for (int i = 0; i < N; i++)
    {
        M[i].update();
    }
}

int main(int argc, char *argv[])
{
    //execution starts here.    
    N = std::stoi(argv[1]);
    n = std::stoi(argv[2]);
    w2 = W / 2;
    h2 = H / 2;
    CELL M[N];

    srand(1);
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    for (int t = 0; t < n; t++)
    {
        looploop(M);
        writeposition(M, posfile);
    }
    std::cout << "simulation complete\n";
}