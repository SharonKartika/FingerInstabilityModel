#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

#include "utilityfunctions.cpp"

float map(float, float, float, float, float);
float randf(float, float);
float Hv(float);


int N{100}; // # agents
int n{500}; // # time steps
constexpr float W{1200.};
constexpr float H{1200.};
constexpr float beta = 60;
constexpr float dt = 0.01;
float w2, h2;

class MOVER
{
public:
    float x, y;
    float vx, vy;
    float ax, ay;
    MOVER()
    {
        x = randf(-w2, w2);
        y = randf(-h2, h2);
        vx = randf(-10, 10) * 0;
        vy = randf(-10, 10) * 0;
        ax = 0.;
        ay = 0.;
    }
    void checkedgesreflect()
    {
        if (x > w2)
        {
            x = w2;
            vx *= -1;
        }
        else if (x < -w2)
        {
            x = -w2;
            vx *= -1;
        }
        if (y > h2)
        {
            y = h2;
            vy *= -1;
        }
        else if (y < -h2)
        {
            y = -h2;
            vy *= -1;
        }
    }
};

void writeposition(MOVER M[], std::ofstream &file)
{
    for (int i = 0; i < N - 1; i++)
        file << M[i].x << ',' << M[i].y << ',';
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}

float forcemag(float r)
{
    return 1E9 / (r * r); // gravity;

    //  lennard jones
    /*float eps = 1e6;
    float sig = 0.1;
    float f = 4 * eps;
    f *= (12 * pow(sig, 12) * pow(r, -13)) -
         6 * pow(sig, 6) * pow(r, -7);
    return f; */
}
void initaccn(MOVER M[])
{   /* 
    Sets accn of all cells to zero.
    Called in every time step */
    for (int i = 0; i < N; i++)
    {
        M[i].ax = 0.;
        M[i].ay = 0.;
    }
}

void getinteractionforce(MOVER &A, MOVER &B)
{
    float dx = B.x - A.x;
    float dy = B.y - A.y;
    float r = sqrt(dx * dx + dy * dy);
    float f = forcemag(r);
    float fx = f * (dx / r);
    float fy = f * (dy / r);
    A.ax += fx;
    A.ay += fy;
    // B.ax -= fx;
    // B.ay -= fy;
}

void looploop(MOVER M[])
{
    /*Two loops over every interacting pair*/
    initaccn(M);
    for (int i = 0; i < N; i++)
    {   //every particle except the last
        for (int j = 0; j < N; j++)
        {   //every particle starting from i,
            // including the last
            getinteractionforce(M[i], M[j]);
        }
    }
}



void updatepos(MOVER M[])
{
    for (int i = 0; i < N; i++)
    {
        M[i].vx += M[i].ax * dt;
        M[i].vy += M[i].ay * dt;
        M[i].x += M[i].vx * dt;
        M[i].y += M[i].vy * dt;
        // M[i].checkedgesreflect();
    }
}

int main(int argc, char *argv[])
{
    N = std::stoi(argv[1]);
    n = std::stoi(argv[2]);
    w2 = W / 2;
    h2 = H / 2;
    MOVER M[N];

    srand(1);
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    for (int t = 0; t < n; t++)
    {
        // time loop
        looploop(M);
        updatepos(M);
        writeposition(M, posfile);
    }
    std::cout << "simulation complete\n";
}