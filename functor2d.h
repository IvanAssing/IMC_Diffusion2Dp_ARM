#ifndef FUNCTOR2D_H
#define FUNCTOR2D_H


#include "imc_dfm.h"

// Objeto Funcão (Super Classe Abstrata)
class Functor2D
{
    public:
        Functor2D(){}
        virtual tFloat operator()(tFloat x, tFloat y) = 0;
};

class Constant2D : public Functor2D
{
    public:
        tFloat c;

        Constant2D(tFloat constant){ c = constant;}
        tFloat operator()(tFloat, tFloat)
        {
            return c;
        }
};

class Sine : public Functor2D
{
    public:
        tFloat a, b, c;

        Sine(tFloat _a, tFloat _b, tFloat _c):a(_a), b(_b), c(_c){}

        tFloat operator()(tFloat x, tFloat)
        {
            return a*sinq(b*x+c);
        }
};

class SFAS : public Functor2D
{
    public:
        tFloat a, b;

        SFAS(tFloat lx, tFloat ly):a(lx), b(ly){}

        tFloat operator()(tFloat x, tFloat y)
        {
            return sinq(M_PIq *x/a)*sinhq(M_PIq*y/a)/sinhq(M_PIq*b/a);
        }
};


#endif // FUNCTOR2D_H
