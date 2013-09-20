#include "element2d.h"

Element2D::Element2D(tInteger _index, Node2D *node0, Node2D *node1, Node2D *node2, Node2D *node3)
    :index(_index)
{
    this->nodes = new Node2D*[4];

    nodes[0] = node0;
    nodes[1] = node1;
    nodes[2] = node2;
    nodes[3] = node3;
}
