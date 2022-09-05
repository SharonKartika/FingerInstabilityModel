#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

#include "utilityfunctions.cpp"

/*
Units
-----
Distance: micrometres
Time    : Hours
*/

int N{100}; // #agents
int n{500}; // #time steps
constexpr float W{1200.};
constexpr float H{1200.};
constexpr float dt = 0.01;
constexpr float beta = 60;
constexpr float PI = 3.1415926535;
constexpr float rt = 70.; // threshold distance
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
        ax, ay,
        etax, etay;
    CELL()
    {
        x = randf(-w2, w2);
        y = randf(-h2, h2);
        vx = randf(-10, 10);
        vy = randf(-10, 10);
        ax = ay = 0.;
        float theta = randf(0, 2*PI);
        etax = cos(theta);
        etay = sin(theta);
    }
    void update()
    {
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

float interactionforcemag(float r)
{
    if (r < rt)
    {
        float U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
        float A0 = 8, A1 = 2, A2 = 25, A3 = 26;
        float force = 0;
        force += U0 * r * exp(-(pow((r / A0), 2)));
        force += U2 * exp(-r / A2);
        force -= U3 * pow(r - A3, 2) * Hv(r - A3);
        force += U1 * (r - A1) * Hv(r - A1);
        return -force;
    }
    else
    {
        return 0.;
    }
}

void setaccn(CELL &A, CELL &B)
{
    float dx = B.x - A.x;
    float dy = B.y - A.y;
    float r = sqrt(dx * dx + dy * dy);
    float f = interactionforcemag(r);
    float fx = f * (dx / r);
    float fy = f * (dy / r);
    A.ax += fx;
    A.ay += fy;
}
void looploop(CELL M[])
{
    // interaction
    for (int i = 0; i < N; i++)
    {
        M[i].ax = M[i].ay = 0.;
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                setaccn(M[i], M[j]);
            }
        }
    }
    // viscek
    for (int i = 0; i < N; i++)
    {
        float ax = 0.;
        float ay = 0.;
        int ninr = 0;
        float dvx = 0.;
        float dvy = 0.;
        for (int j = 0; j < N; j++)
        {

            float dx = M[j].x - M[i].x;
            float dy = M[j].y - M[i].y;
            float r = sqrt(dx * dx + dy * dy);
            if (r < rt)
            {
                dvx = M[j].vx - M[i].vx;
                dvy = M[j].vy - M[i].vy;
                ax += dvx;
                ay += dvy;
                ninr += 1;
            }
        }
        M[i].ax += (beta / ninr) * dvx;
        M[i].ay += (beta / ninr) * dvy;
    }

    // noise
    for (int i = 0; i < N; i++)
    {
        float sig0 = 150.;
        float sig1 = 300.;
        float rho0 = N / (W * H); // reference density
        float rho = 0.;
        float sig = 0.;
        float tau = 1.39;
        float U1 = randf(0,1);
        float U2 = randf(0,1);
        float xix = sqrt(-2*log(U1))*cos(2*PI*U2);
        float xiy = sqrt(-2*log(U1))*sin(2*PI*U2);
        float theta=0.;
        for (int j = 0; j < N; j++)
        {

            float dx = M[j].x - M[i].x;
            float dy = M[j].y - M[i].y;
            float r = sqrt(dx * dx + dy * dy);
            if (r < rt)
            {
                rho += 1;
            }
        }
        rho /= 2 * PI * rt * rt;
        sig = sig0 + (sig1 - sig0) * (1 - rho / rho0);
        M[i].etax += dt*(-M[i].etax+xix)/tau;
        M[i].etay += dt*(-M[i].etay+xiy)/tau;
        //normalize eta
        theta = atan2(M[i].etay, M[i].etax);
        M[i].ax += sig*cos(theta);
        M[i].ax += sig*sin(theta);
    }

    for (int i = 0; i < N; i++)
    {
        M[i].update();
    }
}

int main(int argc, char *argv[])
{
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