#ifndef NODE2D_H
#define NODE2D_H

#include "imc_dfm.h"
#include "boundary2d.h"
#include "data.h"

class Node2D
{
    public:
        tFloat x, y;
        tInteger index;
        Node2D **nodes;
        Boundary2D *boundary;

        Node2D();
        Node2D(tInteger index, tFloat x, tFloat y);
        Node2D(tInteger index, tFloat x, tFloat y, Node2D *nodeSouth, Node2D *nodeEast, Node2D *nodeNorth, Node2D *nodeWest);
        Node2D(tInteger index, tFloat x, tFloat y, Boundary2D *boundary, Node2D *nodeSouth, Node2D *nodeEast, Node2D *nodeNorth, Node2D *nodeWest);

        void link(Node2D *node, DirectionType type);

        static void referenceNode(Node2D *A, DirectionType dir, Node2D *B, tInteger *refnode);

};

#endif // NODE2D_H
