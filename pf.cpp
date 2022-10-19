#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

#include "utilityfunctions.cpp"

/*
Units
-----
Distance: micrometres
Time    : hours

Todo: combine the different forces into a single loop
*/

int N{100}; // #agents
int n{500}; // #time steps
constexpr float W{1200.};
constexpr float H{1200.};
constexpr float dt{0.01};
constexpr float beta{60};
constexpr float PI{3.1415926535};
constexpr float rt{70.}; // threshold distance
float w2, h2;

float unitrand();
float map(float, float, float, float, float);
float Hv(float);
float noisemag(float);

class VEC2
{
public:
    float x, y;
    VEC2(float X, float Y)
    {
        x = X;
        y = Y;
    }
    VEC2()
    {
        x = 0.;
        y = 0.;
    }
    VEC2 operator+(VEC2 const &obj)
    {
        VEC2 temp;
        temp.x = x + obj.x;
        temp.y = y + obj.y;
        return temp;
    }
    VEC2 operator-(VEC2 const &obj)
    {
        VEC2 temp;
        temp.x = x - obj.x;
        temp.y = y - obj.y;
        return temp;
    }
    VEC2 operator*(float const &a)
    {
        VEC2 temp;
        temp.x = a * x;
        temp.y = a * y;
        return temp;
    }
    void operator+=(VEC2 const &obj)
    {
        x += obj.x;
        y += obj.y;
    }
    VEC2 operator/(float const &a)
    {
        VEC2 temp;
        temp.x = x / a;
        temp.y = y / a;
        return temp;
    }
    float mag()
    {
        return sqrt(pow(x, 2) + pow(y, 2));
    }
};

class CELL
{
public:
    VEC2 p, v, a, eta;
    CELL()
    {   
        float theta = randf(0, 2 * PI);
        float r = w2*sqrt(randf(0, 1));
        p = VEC2(r*cos(theta), r*sin(theta)); 
        // p = VEC2(randf(-w2, w2), randf(-h2, h2));
        v = VEC2(randf(-10, 10), randf(-10, 10));
        a = VEC2(0., 0.);
        theta = randf(0, 2 * PI);
        eta = VEC2(cos(theta), sin(theta));
    }
    void update()
    {
        v = v + a * dt;
        p = p + v * dt;
    }
};

void writeposition(CELL M[], std::ofstream &file)
{ /*
  Writes the coordinates of all the cells to file
  */
    for (int i = 0; i < N - 1; i++)
        file << M[i].p.x << ',' << M[i].p.y << ',';
    file << M[N - 1].p.x << ',' << M[N - 1].p.y << '\n';
}

float interactionforcemag(float r)
{
    if (r < rt)
    {
        float U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
        float A0 = 8, A1 = 2, A2 = 25, A3 = 26;
        // A0 = 40; //temp
        float force = 0;
        force += U0 * r * exp(-(pow((r / A0), 2)));
        // force += U2 * exp(-r / A2);
        // force -= U3 * pow(r - A3, 2) * Hv(r - A3);
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
    VEC2 dp = B.p - A.p;
    float r = dp.mag(); 
    float f = interactionforcemag(r);
    VEC2 F(f * (dp.x / r), f * (dp.y / r));
    A.a += F;
}

/*Loops through every cell and again through every cell  */
void looploop(CELL M[])
{
    // interaction
    for (int i = 0; i < N; i++)
    {
        M[i].a = VEC2(0., 0.);
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
        VEC2 a;
        int ninr = 0;
        VEC2 dv;
        for (int j = 0; j < N; j++)
        {

            VEC2 dp = M[j].p - M[i].p;
            float r = dp.mag();
            if (r < rt)
            {
                dv = M[j].v - M[i].v;
                a += dv;
                ninr += 1;
            }
        }

        M[i].a += a * (beta / ninr);
    }

    // noise
    for (int i = 0; i < N; i++)
    {
        float sig0 = 150.;
        float sig1 = 300.;
        float rho0 = N / (W * H); // reference density (change to rho1)
        float rho = 0.;
        float sig = 0.;
        float tau = 1.39;
        // float tau = 0.01;

        VEC2 U = VEC2(randf(0, 1), randf(0, 1));
        VEC2 xi(sqrt(-2 * log(U.x)) * cos(2 * PI * U.y),
                sqrt(-2 * log(U.x)) * sin(2 * PI * U.y));
        float std1 = 1;//temp
        xi = xi * std1;
        float theta = 0.;
        for (int j = 0; j < N; j++)
        {

            VEC2 dp = M[j].p - M[i].p;
            float r = dp.mag();
            if (r < rt)
            {
                rho += 1;
            }
        }
        rho /= PI * rt * rt;
        sig = sig0 + (sig1 - sig0) * (1 - rho / rho0);
        // euler maruyama integration
        M[i].eta = M[i].eta - (M[i].eta) * (dt / tau);
        M[i].eta = M[i].eta + (xi) * (sqrt(dt) / tau);

        M[i].eta = M[i].eta / M[i].eta.mag(); //normalize eta
        M[i].a += M[i].eta * sig;
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

    // srand(1);
    // temp
    srand(time(0));
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    for (int t = 0; t < n; t++)
    {
        looploop(M);
        writeposition(M, posfile);
    }
    std::cout << "simulation complete\n";
}