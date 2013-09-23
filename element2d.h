#ifndef ELEMENT2D_H
#define ELEMENT2D_H

#include "node2d.h"
#include "imc_dfm.h"

#define NULL_INDEX -123456

class Element2D
{
    public:
        tInteger index;
        Node2D **nodes;

        Element2D(){index = NULL_INDEX;}
        Element2D(tInteger index, Node2D *node0, Node2D *node1, Node2D *node2, Node2D *node3);

        ~Element2D()
        {
            if(index != NULL_INDEX)
                delete [] nodes;
        }
};

#endif // ELEMENT2D_H
