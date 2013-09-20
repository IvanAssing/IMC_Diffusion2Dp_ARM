#ifndef DIFFUSION2DPAR_H
#define DIFFUSION2DPAR_H


#include "imc_dfm.h"
#include "functor2d.h"
#include "node2d.h"
#include "boundary2d.h"
#include "data.h"
#include "element2d.h"


class Diffusion2DpAR
{
    public:
        Diffusion2DData *data; // Dados do problema
        Node2D *nodes; // Lista de nós
        Element2D *elements; // Lista de elementos
        tInteger nx, ny; // Número de divisões em x e y
        tFloat lx, ly, hx, hy; // Tamanho de malha
        Boundary2D *ccS, *ccN, *ccE, *ccW; // Condições de contorno
        tInteger nnodes, nelements;

        tFloat *T, Tm; // Vetor solução

        Diffusion2DpAR(tFloat lengthX, tFloat lengthY, tInteger nx, tInteger ny, Diffusion2DData *data,
                     Boundary2D *ccSouth, Boundary2D *ccNorth, Boundary2D *ccEast, Boundary2D *ccWest);

        void solver(tInteger iterationMax = 1000, tFloat iterationTolerance = 1.0e-28q, bool plotlog = false); // Discretização + solver do sistema linear

        tInteger direction(tInteger position, DirectionType dir);
        Node2D* nodeDirection(tInteger position, DirectionType dir);
        tInteger position(tInteger i, tInteger j);

        void refine(Element2D *element);
        void refine(Node2D *node);
};

#endif // DIFFUSION2DPAR_H
