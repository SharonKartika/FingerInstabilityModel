#include <iostream>
#include <math.h>
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
    VEC2 operator/(float const &a)
    {
        VEC2 temp;
        temp.x = x / a;
        temp.y = y / a;
        return temp;
    }
    VEC2 operator+=(VEC2 const &obj)
    {
        x += obj.x;
        y += obj.y;
    }
    float mag()
    {
        return sqrt(pow(x, 2) + pow(y, 2));
    }
};

int main()
{
    VEC2 a(3, 4 );
    // VEC2 b = a;
    // b = b / 2;
    // VEC2 c;
//     std::cout << b.x << " " << b.y;
    // std::cout<<a.mag();
    std::cout<< pow((a/a.mag()).x,2) + pow((a/a.mag()).y,2);
}