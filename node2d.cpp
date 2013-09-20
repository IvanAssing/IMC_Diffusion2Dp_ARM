#include "node2d.h"

#define DIFFTOL 1.0e-6q

Node2D::Node2D(tInteger _index, tFloat _x, tFloat _y, Node2D *node0, Node2D *node1, Node2D *node2, Node2D *node3)
    :index(_index), x(_x), y(_y)
{
    this->nodes = new Node2D*[4];

    nodes[0] = node0;
    nodes[1] = node1;
    nodes[2] = node2;
    nodes[3] = node3;

    boundary = NULL;
}

Node2D::Node2D(tInteger _index, tFloat _x, tFloat _y, Boundary2D *_boundary, Node2D *node0, Node2D *node1, Node2D *node2, Node2D *node3)
    :index(_index), x(_x), y(_y)
{
    this->nodes = new Node2D*[4];

    nodes[0] = node0;
    nodes[1] = node1;
    nodes[2] = node2;
    nodes[3] = node3;

    boundary = _boundary;
}

Node2D::Node2D(tInteger _index, tFloat _x, tFloat _y)
    :index(_index), x(_x), y(_y)
{
    this->nodes = new Node2D*[4];

    nodes[0] = NULL;
    nodes[1] = NULL;
    nodes[2] = NULL;
    nodes[3] = NULL;

    boundary = NULL;
}

Node2D::Node2D()
{
    //    this->nodes = new Node2D*[4];

    //    nodes[0] = NULL;
    //    nodes[1] = NULL;
    //    nodes[2] = NULL;
    //    nodes[3] = NULL;

    //    boundary = NULL;
}

void Node2D::link(Node2D *node, DirectionType type)
{
    switch (type) {
    case North:
        nodes[2] = node;
        node->nodes[0] = this;
        break;
    case South:
        nodes[0] = node;
        node->nodes[2] = this;
        break;
    case West:
        nodes[3] = node;
        node->nodes[1] = this;
        break;
    case East:
        nodes[1] = node;
        node->nodes[3] = this;
        break;
    }
}

void Node2D::referenceNode(Node2D *A, DirectionType type, Node2D *B, tInteger *refnode)
{
    Node2D *NN = A;
    switch (type) {
    case North:
    {
        tFloat pt = 0.5q*(A->y + B->y);
        while(fabsq(NN->nodes[2]->y - pt)>DIFFTOL)
            NN = NN->nodes[2];
        *refnode = NN->nodes[2]->index;
        break;
    }
    case South:
    {
        tFloat pt = 0.5q*(A->y + B->y);
        while(fabsq(NN->nodes[0]->y - pt)>DIFFTOL)
            NN = NN->nodes[0];
        *refnode = NN->nodes[0]->index;
        break;
    }
    case West:
    {
        tFloat pt = 0.5q*(A->x + B->x);
        while(fabsq(NN->nodes[3]->x - pt)>DIFFTOL)
            NN = NN->nodes[3];
        *refnode = NN->nodes[3]->index;
        break;
    }
    case East:
    {
        tFloat pt = 0.5q*(A->x + B->x);
        while(fabsq(NN->nodes[1]->x - pt)>DIFFTOL)
            NN = NN->nodes[1];
        *refnode = NN->nodes[1]->index;
        break;
    }
    }
}
