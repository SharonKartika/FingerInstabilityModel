#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

constexpr int N{1000}; // #agents
constexpr int n{100};  // #time steps
constexpr float W{1200.};
constexpr float H{1200.};
constexpr float dt = 1;
float w2, h2;

float unitrand();
float map(float, float, float, float, float);

class MOVER
{
public:
    float x, y;
    float vx, vy;
    MOVER()
    {
        x = map(unitrand(), 0., 1., -w2, w2);
        y = map(unitrand(), 0., 1., -h2, h2);
        vx = 0.;
        vy = 0.;
    }
    void update()
    {
        x = x + vx * dt;
        y = y + vy * dt;
    }
};

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
void writeposition(MOVER M[], std::ofstream &file)
{
    for (int i = 0; i < N - 1; i++)
        file << M[i].x << ',' << M[i].y << ',';
    file << M[N - 1].x << ',' << M[N - 1].y << '\n';
}
float Hv(float r)
{ // heaviside function
    return r > 0;
}

float getforce(float r)
{
    float U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
    float A0 = 8, A1 = 2, A2 = 25, A3 = 26;
    float force = 0;
    force += U0 * r * exp(-(pow((r / A0), 2)));
    force += U2 * exp(-r / A2);
    force -= U3 * pow(r - A3, 2) * Hv(r - A3);
    force += U1 * (r - A1) * Hv(r - A1);
    return force;
    // return 100000000. / (r * r);
}
float tempvar;
void calcvelocities(MOVER M[])
{
    float f;
    float dx, dy, vx, vy, r;
    float theta;
    for (int i = 0; i < N; i++)
    {
        vx = 0;
        vy = 0;
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                dx = M[j].x - M[i].x;
                dy = M[j].y - M[i].y;
                r = sqrt(dx * dx + dy * dy);
                if (r < 70)
                {
                    theta = atan2(dy, dx);
                    f = getforce(r);

                    vx += f * cos(theta) * dt;
                    vy += f * sin(theta) * dt;
                }
            }
        }
        M[i].vx = vx;
        M[i].vy = vx;
        // std::cout << " vx " << vx << " vy " << vy << "\n";
    }
}

void updateposition(MOVER M[])
{
    for (int i = 0; i < N; i++)
        M[i].update();
}

int main()
{
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
    std::cout << "finished running";
    std::cout << getforce(70.);
}