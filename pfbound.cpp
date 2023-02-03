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
constexpr float rt{400.}; // threshold distance
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
        float r = w2 * sqrt(randf(0, 1));
        p = VEC2(r * cos(theta), r * sin(theta));
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
        float std1 = 1; // temp
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

        M[i].eta = M[i].eta / M[i].eta.mag(); // normalize eta
        M[i].a += M[i].eta * sig;
    }

    for (int i = 0; i < N; i++)
    {
        M[i].update();
    }
}

VEC2 mean(CELL M[])
{
    VEC2 m(0, 0);
    for (int i = 0; i < N; i++)
    {
        m += M[i].p;
    }
    return m / N;
}

VEC2 mean2(CELL **M)
{
    VEC2 m(0, 0);
    CELL **p = M;
    float i = 0.;

    while (*p != NULL)
    {
        m += (*p)->p;
        i += 1.;
        p++;
    }
    return m / i;
}

/*Checks if A is in B*/
int in(CELL *p, CELL **m)
{
    while (*m)
    {
        if (p == *m)
        {
            return 1;
        }
        m++;
    }
    return 0;
}

/* Takes in an array and its lengths.
Finds index of largest element.*/
int argmax(float M[], int K)
{
    int mi = 0;
    for (int i = 1; i < K; i++)
    {
        if (M[mi] <= M[i])
        {
            mi = i;
        }
    }
    return mi;
}

float dist(CELL a, CELL b)
{
    return (a.p - b.p).mag();
}

/*  Allocate N pointers.
    Return pointer to pointer array
    */
CELL **getcellarray(int N)
{
    CELL **p = (CELL **)malloc(N * sizeof(p));
    for (int i = 0; i < N; i++)
    {
        p[i] = NULL;
    }
    return p;
}

/*  Get set difference A\B
    Elements in A, not in B.
    */
CELL **setdiff(CELL **A, CELL **B)
{
    CELL **p = getcellarray(N / 5);
    CELL **c = p;
    while (*A)
    {
        if (!in(*A, B))
        {
            *c = *A;
            c++;
        }
        A++;
    }
    return p;
}

CELL **getneighbors(CELL M[], CELL *cell, float rt)
{
    // CELL **rns = getcellarray(N / 5);
    CELL **rns = getcellarray(N);

    CELL **p = rns;
    for (int i = 0; i < N; i++)
    {
        if (dist(*cell, M[i]) < rt)
        {
            if (cell != &M[i])
            {
                *p = &M[i];
                p++;
            }
        }
    }
    return rns;
}

float angleAC(VEC2 *rl, VEC2 *rc, VEC2 *rn)
{
    VEC2 A = *rc - *rl;
    VEC2 B = *rn - *rc;
    float k = B.y * A.x - A.y * B.x;
    float costheta = (A.x * B.x + A.y * B.y) /
                     (A.mag() * B.mag());
    float theta = acos(costheta);
    if (k > 0)
        return PI - theta;
    else
        return 2 * PI - theta;
}

int len(CELL **A)
{
    int i = 0;
    while (*A)
    {
        A++;
        i++;
    }
    return i;
}

void writeboundarycells(CELL **B, std::ofstream &file)
{
    int nb = len(B);
    for (int i = 0; i < nb - 1; i++)
        file << (*(B + i))->p.x << ',' << (*(B + i))->p.y << ',';
    file << (*(B + nb - 1))->p.x << ',' << (*(B + nb - 1))->p.y << '\n';
}

void cout(VEC2 *p)
{
    std::cout << p->x << "  " << p->y << "  ";
}

/*Writes to a file the cells that lie on the border*/
CELL **findBorderCells(CELL M[], float rt)
{
    VEC2 rcm = mean(M);
    float dists[N];
    for (int i = 0; i < N; i++)
    {
        dists[i] = (rcm - M[i].p).mag();
    }
    int i = argmax(dists, N);
    CELL *rft = &M[i];
    CELL *rc = rft;

    CELL **boundcells = getcellarray(N);
    while (true)
    {
        CELL **rns = getneighbors(M, rc, rt);
        VEC2 rl = mean2(rns);
        CELL **rnsnb = setdiff(rns, boundcells);
        int lenrnsnb = len(rnsnb);
        float *nscores = (float *)malloc(lenrnsnb * sizeof(*nscores));
        for (int i = 0; i < lenrnsnb; i++)
        {
            *(nscores + i) = angleAC(&rl, &(rc->p), &(*(rnsnb + i))->p);
            if (*(nscores + i) >= PI)
                (*(nscores + i) = 0.);
            // cout(&rl);
            // cout(&(rc->p));
            // cout(&(*(rnsnb + i))->p);
            // std::cout << std::endl;
            // std::cout << *(nscores + i) << std::endl;
        }
        // debug segfault
        std::cout << len(rns) << "\t" << len(rnsnb) << "\t" << len(boundcells) << std::endl;
        // end debug segfault

        int rnindex = argmax(nscores, lenrnsnb);
        CELL *rn = rnsnb[rnindex];
        rc = rn;
        boundcells[len(boundcells)] = rc;

        if (rc == rft)
        {
            std::cout << "DONE ";
            break;
        }
    }
    return boundcells;
}

int *getQuadrantCount(CELL **NC, CELL *C)
{
    static int qc[4];
    for (int i = 0; i < 4; i++)
    {
        qc[i] = 0;
    }
    int Yv, Xv;
    int counttemp = 0;
    while (*NC)
    {
        if ((*NC)->p.x > C->p.x)
            Xv = 1;
        else
            Xv = 0;
        if ((*NC)->p.y > C->p.y)
            Yv = 1;
        else
            Yv = 0;

        if (Xv == 1 && Yv == 1)
            qc[0]++;
        else if (Xv == 1 && Yv == 0)
            qc[1]++;
        else if (Xv == 0 && Yv == 0)
            qc[2]++;
        else if (Xv == 0 && Yv == 1)
            qc[3]++;
        else
            std::cout << "Some kind of error" << std::endl;
        NC++;
    }
    for (int i = 0; i < 4; i++)
    {
        std::cout << "QV: " << qc[i] << "\t";
    }
    std::cout << std::endl;
    return qc;
}

bool isOnBoundary(int *a)
{
    // Difference between number of cells in one region
    //  greater than 5
    //  than the average number of cells in other 3
    for (int i = 0; i < 4; i++)
    {
        float avg = 0;
        for (int j = 0; j < 4; j++)
        {
            if (i != j)
            {
                avg += a[j];
            }
        }
        avg /= 3;
        if (abs(a[i] - avg) > 5)
        {
            std::cout << "1REGIONTRUE" << std::endl;
            return true;
        }
    }
    // Difference between total particles in 2 regions
    //  greater than 8
    //  than total particles in the other 2
    for (int i = 0; i < 4; i++)

    {
        int c1, c2;
        c1 = a[i % 4] + a[(i + 1) % 4];
        c2 = a[(i + 2) % 4] + a[(i + 3) % 4];
        if (abs(c1 - c2) > 8)
        {
            std::cout << "2REGIONTRUE" << std::endl;
            return true;
        }
    }
    // opposite quadrants, count > 8?
    //  int v1 = a[0] + a[2];
    //  int v2 = a[1] + a[3];
    //  if (abs(v1 - v2) > 8)
    //  {
    //      return true;
    //  }

    // else: not on the boundary
    return false;
}

bool isQuadrantEmpty(int *a)
{
    for (int i = 0; i < 4; i++)
    {
        if (a[i] == 0)
        {
            return true;
        }
    }
    return false;
}

void sort(float *q, int nelt)
{
    float key;
    int j;
    for (int i = 1; i < nelt; i++)
    {
        key = q[i];
        j = i - 1;

        while (j >= 0 && q[j] > key)
        {
            q[j + 1] = q[j];
            j = j - 1;
        }
        q[j + 1] = key;
    }
}

float findLargestGap(float *q, int nelt)
{
    float lg = 0.;
    float t;
    for (int i = 0; i < nelt - 1; i++)
    {
        t = q[i + 1] - q[i];
        std::cout << t << " ";
        if (t > lg)
        {
            lg = t;
        }
    }
    //debug
    t = 2 * PI - q[nelt - 1] + q[0];
    // t = 2*PI - (q[0] - q[nelt - 1]);
    std::cout << t << " " << std::endl;
    if (t > lg)
    {
        lg = t;
    }
    return lg;
}

// -[ ] Optimise using short circuiting
bool isOnBoundaryFOV(CELL **rns, CELL *C, float f)
{
    int nelt = len(rns);
    float *q = (float *)malloc(sizeof(float) * nelt);
    for (int i = 0; i < nelt; i++)
    {
        VEC2 s = (*(rns + i))->p - C->p;
        *(q + i) = atan2(s.y, s.x);
    }
    sort(q, nelt);
    float lg = findLargestGap(q, nelt);
    if (lg > f)
    {
        return true;
    }
    return false;
}

CELL **findBorderCellsByDiff(CELL M[], float rt, float f)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        CELL **rns = getneighbors(M, &M[i], rt);
        int *qc = getQuadrantCount(rns, &M[i]);
        // if (isOnBoundary(qc))
        // if (isQuadrantEmpty(qc))
        if (isOnBoundaryFOV(rns, &M[i], f))
        {
            *p = &M[i];
            p++;
        }
    }
    return boundcells;
}

CELL **findBorderCellsByVecSum(CELL M[], float rt)
{
    CELL **boundcells = getcellarray(N);
    CELL **p = boundcells;
    for (int i = 0; i < N; i++)
    {
        VEC2 vecsum(0, 0);
        CELL **rns = getneighbors(M, &M[i], rt);
        CELL **t = rns;
        while (*t)
        {
            vecsum += (*t)->p;
            // vecsum += ((*t)->p)/((*t)->p).mag();
            t++;
        }
        vecsum = vecsum / len(rns);
        if (vecsum.mag() > 400)
        {
            std::cout << vecsum.mag() << std::endl;
            *p = &M[i];
            p++;
        }
    }
    // if the vecsum magnitude is greater than some threshold,
    // the cell is on the boundary
    return boundcells; // WRORNG ANGSER
}

int main(int argc, char *argv[])
{
    N = std::stoi(argv[1]);
    n = std::stoi(argv[2]);
    w2 = W / 2;
    h2 = H / 2;
    CELL M[N];

    // srand(1);
    /*CELL **rns = getneighbors(M, &M[0], 200);
    CELL **p=rns;

    while (*p != NULL)
    {
        std::cout << (*p)->p.y << ",";
        p++;
    }
    VEC2 m = mean2(rns);
    std::cout<<std::endl<<m.y<<std::endl;
    */

    // (*C==NULL)? std::cout <<"NULL" : std::cout<<"Not null";

    // CELL **B = findBorderCells(M, 200);
    // (B == NULL) ? std::cout << "NULL" : std::cout << "Not NULL";
    // std::cout << len(B);
    // std::cout<<p.x<<"\t"<<p.y<<std::endl;

    // int i;
    // while (true){

    // }
    // end test
    // while (true){
    // std::cout << (rns[0])->p.x;

    // }

    // srand(time(0));
    std::ofstream posfile;
    posfile.open("positiondata.csv");

    // begin boundary test
    std::ofstream boundposfile;
    boundposfile.open("boundarypositiondata.csv");
    CELL **B;
    // end boundary test

    // test isonboundaryPOV

    // end test isonboundaryPOV

    // B = findBorderCellsByDiff(M, 300, 2.);
    // B = findBorderCellsByVecSum(M, rt);
    // std::cout << "LENB: " << len(B) << std::endl;
    // writeboundarycells(B, boundposfile);
    // writeposition(M, posfile);
    
    float delta = 0.1;
    // float a = PI - delta;
    // float b = -PI + delta;
    // float q[2] = {a, b};
    // sort(q, 2);
    // for (int i = 0; i < 2; i++){
    //     std::cout<<q[i]<<" ";
    // }
    // findLargestGap(q, 2);

    float q[3] = {-PI+delta, delta, PI-delta};
    sort(q, 3);
    findLargestGap(q, 3);
    // for (int t = 0; t < n; t++)
    // {
    // looploop(M);
    // writeposition(M, posfile);
    // B = findBorderCells(M, 200);
    // writeboundarycells(B, boundposfile);
    // }
    std::cout << "simulation complete\n";
}