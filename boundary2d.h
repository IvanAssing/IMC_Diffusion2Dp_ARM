#ifndef BOUNDARY2D_H
#define BOUNDARY2D_H

#include "functor2d.h"
#include "data.h"

// Tipos de condição de contorno
enum BoundaryCondition{
    Dirichlet,
    Neumann,
    Robin
};

// Objeto para definir parametros de contorno
class Boundary2D
{
    public:
        Functor2D *bcValue;
        BoundaryCondition type;
        DirectionType dir;

        Boundary2D();
        Boundary2D(BoundaryCondition type, Functor2D *bcValue, DirectionType dir = South);
};
#endif // BOUNDARY2D_H
