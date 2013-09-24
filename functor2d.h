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

class Sine2 : public Functor2D
{
    public:
        tFloat a, b, c;

        Sine2(tFloat _a, tFloat _b, tFloat _c):a(_a), b(_b), c(_c){}

        tFloat operator()(tFloat x, tFloat y)
        {
            return a*sinq(b*x+c)*sinq(b*y+c);
        }
};

class SFAS2 : public Functor2D
{
    public:
        SFAS2(){}

        tFloat operator()(tFloat x, tFloat y)
        {
            //if (x==0.5q && y == 0.5q) return 0.0q;
            return -(1.q/25.q)*(64.q*sinq(8.q*M_PIq*sqrtq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q)))*M_PIq*M_PIq*(x+1.q)*(x+1.q)
                                +8.q*cosq(8.q*M_PIq*sqrtq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q)))*M_PIq*sqrtq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q))
                                -sinq(8.q*M_PIq*sqrtq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q)))+64.q*sinq(8.q*M_PIq*sqrtq((x+1.q)*(x+1.q)+
                                (y+1.q)*(y+1.q)))*M_PIq*M_PIq*(y+1.q)*(y+1.q))/powq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q),1.5q);

            //            return -(1.q/25.q)*(64.q*sinq(8.q*M_PIq*sqrtq(x*x+y*y))*M_PIq*M_PIq*x*x+8.q*cos(8.q*M_PIq*sqrtq(x*x+y*y))*M_PIq*sqrtq(x*x+y*y)
            //                    -sinq(8.q*M_PIq*sqrtq(x*x+y*y))+64.q*sinq(8.q*M_PIq*sqrtq(x*x+y*y))*M_PIq*M_PIq*y*y)/powq(x*x+y*y,1.5q);
        }

};


class SFAS2a : public Functor2D
{
    public:
        SFAS2a(){}

        tFloat operator()(tFloat x, tFloat y)
        {
            //if (x==0.5q && y == 0.5q) return 0.0q;
            return sinq(8.q*M_PIq*sqrtq((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q)))/(25.q*sqrt((x+1.q)*(x+1.q)+(y+1.q)*(y+1.q)));
            //if (x==0 && y == 0) return 0.1q;
            //return sinq(8.q*M_PIq*sqrtq(x*x+y*y))/(25.q*sqrtq(x*x+y*y));
        }

};



#endif // FUNCTOR2D_H
