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

class MOVER
{
public:
    float x, y;
    float vx, vy;
    MOVER()
    {
        x = map(unitrand(), 0., 1., -w2, w2);
        y = map(unitrand(), 0., 1., -h2, h2);
        // vx = 0.;
        // vy = 0.;
        vx = (unitrand() - 0.5)*500;
        vy = (unitrand() - 0.5)*500;
    }
    void update()
    {
        x += vx * dt;
        y += vy * dt;
    }
};

void writeposition(MOVER M[], std::ofstream &file)
{
    for (int i = 0; i < N - 1; i++)
        file << M[i].x << ',' << M[i].y << ',';
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}

float getinteractionforce(float r)
{
    float U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
    float A0 = 8, A1 = 2, A2 = 25, A3 = 26;
    float force = 0;
    force += U0 * r * exp(-(pow((r / A0), 2)));
    force += U2 * exp(-r / A2);
    force -= U3 * pow(r - A3, 2) * Hv(r - A3);
    force += U1 * (r - A1) * Hv(r - A1);
    return force;
    // return 10000. / (r * r); // gravity
}
void calcvelocities(MOVER M[])
{
    float f;
    float dx, dy, vx, vy, ax, ay, r;
    float theta;
    for (int i = 0; i < N; i++)
    {
        //new branch
        ax = 0;
        ay = 0;

        for (int j = 0; j < N; j++)
        { // interaction force

            if (i != j)
            {
                // interactionforce(&ax, &ay, M[i], M[j]);

                dx = M[j].x - M[i].x;
                dy = M[j].y - M[i].y;
                r = sqrt(dx * dx + dy * dy);

                if (r < 70)
                {
                    theta = atan2(dy, dx);
                    f = -1*getinteractionforce(r);

                    ax += f * cos(theta);
                    ay += f * sin(theta);
                }
                // interaction force end
            }
        }

        float axt = 0, ayt = 0;
        int Ni = 1; // number of nearest neighbors of M[i]
        for (int j = 0; j < N; j++)
        { // viscek force
            if (i != j)
            {
                dx = M[j].x - M[i].x;
                dy = M[j].y - M[i].y;
                r = sqrt(dx * dx + dy * dy);

                if (r < 70)
                {
                    Ni += 1;
                    axt += (M[j].vx - M[i].vx);
                    ayt += (M[j].vy - M[i].vy);
                }
            }
        }
        axt *= (beta/Ni);
        ayt *= (beta/Ni);

        ax += axt;
        ay += ayt;

        M[i].vx += ax * dt;
        M[i].vy += ay * dt;
        // std::cout << " vx " << vx << " vy " << vy <<  "\n";
    }
}

void updateposition(MOVER M[])
{
    for (int i = 0; i < N; i++)
        M[i].update();
}

int main(int argc, char *argv[])
{
    N = std::stoi(argv[1]);
    n = std::stoi(argv[2]);

    srand(1);
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    w2 = W / 2;
    h2 = H / 2;

    MOVER M[N];
    for (int t = 0; t < n; t++)
    {
        calcvelocities(M);
        updateposition(M);
        writeposition(M, posfile);
    }
    std::cout << "simulation complete\n";
}