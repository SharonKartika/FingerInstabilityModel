#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

// void calcthetas(MOVER M[]);
// void updateposition(MOVER M[]);
// void writeposition(MOVER M[], std::ofstream&);
float unitrand();
float map(float ri, float x1, float x2, float y1, float y2);

const double PI = 3.14159265358979323846;
float W = 100.; // total W
float H = 100.; // total H
float w2, h2;   // half W and H
float r = 5.0;  // interaction radius
float r2;       // squared iteraction radius
int N = 500;    // number of particles
int n = 500;    // number of time-steps
float eta;      // randomness parameter

class MOVER
{
public:
    float x, y;
    float vx, vy;
    float vabs, theta;
    MOVER()
    {
        x = map(unitrand(), 0., 1., -w2, w2);
        y = map(unitrand(), 0., 1., -h2, h2);
        theta = map(unitrand(), 0., 1., -PI, PI);
        vabs = 1.;
    }
    void update()
    {
        vx = vabs * cos(theta);
        vy = vabs * sin(theta);

        x = x + vx;
        y = y + vy;
    }
    void display()
    {
        std::cout << x << ' ' << y << '\n';
    }
};

void calcthetas(MOVER M[])
/*Calculates the directions for all the movers*/
{
    int ninr;
    float avtheta, ssin, scos, dx, dy, dxe, dye, d2;

    for (int i = 0; i < N; i++)
    {
        // outer object loop
        ninr = 0;
        avtheta = 0.;
        ssin = scos = 0.;
        for (int j = 0; j < N; j++)
        {
            // inner object loop
            dx = M[j].x - M[i].x;
            dy = M[j].y - M[i].y;

            dxe = dx - W * std::nearbyintf(dx / W);
            dye = dy - H * std::nearbyintf(dy / H);

            d2 = pow(dxe, 2) + pow(dye, 2);
            if (d2 < r2)
            {
                // check radius
                ninr++;
                ssin += sin(M[j].theta);
                scos += cos(M[j].theta);
            }
        }
        ssin /= ninr;                // averaging out
        scos /= ninr;                // "
        avtheta = atan2(ssin, scos); // direction
        M[i].theta = avtheta + map(unitrand(), 0, 1, -eta, eta);        // assigning direction, adding randomness
    }
}

void updateposition(MOVER M[])
/*Updates the position of all the movers */
{
    for (int i = 0; i < N; i++)
    {
        M[i].update();
    }
}

void writeposition(MOVER M[], std::ofstream &file)
{
    for (int i = 0; i < N - 1; i++)
    {
        file << M[i].x << ',' << M[i].y << ',';
    }
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}

void writevelocity(MOVER M[], std::ofstream &file)
{
    for (int i = 0; i < N - 1; i++)
    {
        file << M[i].vx << ',' << M[i].vy << ',';
    }
    file << M[N - 1].vx << ',' << M[N - 1].vy << '\n';
}


float unitrand()
/*Returns a float chosen randomly from [0,1]  */
{
    return float(rand()) / INT_MAX;
}

float map(float ri, float x1, float x2, float y1, float y2)
/*Maps one interval to another */
{
    float runit = (ri - x1) / (x2 - x1);
    float rf = runit * (y2 - y1) + y1;
    return rf;
}

int main()
{
    srand(1);              // random seed
    std::ofstream posfile;    // file object
    std::ofstream velfile;    
    posfile.open("positiondata.csv"); // open file
    velfile.open("velocitydata.csv"); 

    eta = PI/8;
    w2 = W / 2;
    h2 = H / 2;
    r2 = r * r;
    MOVER M[N];
    // time-loop
    for (int i = 0; i < n; i++)
    {
        calcthetas(M);
        updateposition(M);
        writeposition(M, posfile);
        writevelocity(M, velfile);
    }
    posfile.close();
    velfile.close();
    return 0;
}