#ifndef ELEMENT2D_H
#define ELEMENT2D_H

#include "node2d.h"
#include "imc_dfm.h"

class Element2D
{
    public:
        tInteger index;
        Node2D **nodes;

        Element2D(){}
        Element2D(tInteger index, Node2D *node0, Node2D *node1, Node2D *node2, Node2D *node3);
};

#endif // ELEMENT2D_H
