#include <iostream>
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>

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

float randf(float a = 0., float b = 1.)
/*Returns a float chosen randomly
Defaults to [0,1]*/
{
    float r = float(rand()) / INT_MAX;
    return map(r, 0, 1, a, b);
}


float Hv(float r)
{ // heaviside function
    return r > 0;
}

float noisemag(float rho)
{ // noise magnitude
    float s0 = 150;
    float s1 = 300;
    float rho0 = 2.2e-3;
    return s0 + (s1-s0)*(1-(rho/rho0));
}

