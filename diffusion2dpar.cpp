#include "diffusion2dpar.h"

#include "gaussseidel.h"
#include "solverlu.h"
#include "stdlib.h"

#include "graphics.h"

#define MAXNODES(nx, ny) 50000*nx*ny
#define MAXELEMENTS(nx, ny) 50000*nx*ny

Diffusion2DpAR::Diffusion2DpAR(tFloat _lx, tFloat _ly, tInteger _nx, tInteger _ny, Diffusion2DData *_data,
                               Boundary2D *_bS, Boundary2D *_bN, Boundary2D *_bE, Boundary2D *_bW)
    :lx(_lx), ly(_ly), nx(_nx), ny(_ny), data(_data), ccS(_bS), ccN(_bN), ccE(_bE), ccW(_bW)
{
    hx = lx/(nx - 1.0q);
    hy = ly/(ny - 1.0q);

    // Nodes
    nodes = new Node2D*[MAXNODES(nx,ny)];
    tInteger p = 0;


    for(tInteger j=0; j<ny; j++)
        for(tInteger i=0; i<nx; i++)
            nodes[p++] = new Node2D(p, i*hx, j*hy,
                                    nodeDirection(p, South),
                                    nodeDirection(p, East),
                                    nodeDirection(p, North),
                                    nodeDirection(p, West));

    nnodes = p;

    for(tInteger i=0; i<nnodes; i++){
        nodes[i]->nodes[0] = nodeDirection(i, South);
        nodes[i]->nodes[1] = nodeDirection(i, East);
        nodes[i]->nodes[2] = nodeDirection(i, North);
        nodes[i]->nodes[3] = nodeDirection(i, West);
    }


    for(tInteger i=0; i<nnodes; i++)
        std::cout<<"\n NODE2D "<<i<<"\t"<<nodes[i]<<"\t"<<QtoD(nodes[i]->x)<<"\t"<<
                   QtoD(nodes[i]->y)<<"\t"<<nodes[i]->nodes[0]<<"\t"<<nodes[i]->nodes[1]<<"\t"<<
                   nodes[i]->nodes[2]<<"\t"<<nodes[i]->nodes[3];


    // Elements
    elements = new Element2D*[MAXELEMENTS(nx,ny)];
    p = 0;

    for(tInteger j=0; j<ny-1; j++)
        for(tInteger i=0; i<nx-1; i++){
            tInteger nn = position(i, j);
            elements[p++] = new Element2D(p, nodes[nn],
                                          nodeDirection(nn, East),
                                          nodeDirection(direction(nn, East), North),
                                          nodeDirection(nn, North));
        }

    nelements = p;

    for(tInteger i=0; i<nelements; i++)
        std::cout<<"\n ELEMENT2D "<<i<<"\t"<<elements[i]<<"\t"<<
                   elements[i]->nodes[0]<<"\t"<<elements[i]->nodes[1]<<"\t"<<
                   elements[i]->nodes[2]<<"\t"<<elements[i]->nodes[3];

    T = new tFloat[MAXNODES(nx,ny)];

    for(tInteger i=0; i<MAXNODES(nx,ny); i++)
        T[i] = 0.q;


    // Condições de Contorno
    // Sul
    for(tInteger i=0; i<nx; i++)
        nodes[i]->boundary = ccS;

    // East
    for(tInteger i=1; i<ny; i++)
        nodes[i*nx-1]->boundary = ccE;

    // West
    for(tInteger i=0; i<ny; i++)
        nodes[i*nx]->boundary = ccW;

    // Norte
    for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
        nodes[i]->boundary = ccN;


}


tInteger Diffusion2DpAR::position(tInteger i, tInteger j)
{
    return i+j*nx;
}


void Diffusion2DpAR::updateTm(void)
{
    Tm = 0.0q;

    for(tInteger i=0; i<nelements; i++)
    {
        if(elements[i]->index == -1) continue;
        Tm += (T[elements[i]->nodes[0]->index] + T[elements[i]->nodes[1]->index] +
                T[elements[i]->nodes[2]->index] + T[elements[i]->nodes[3]->index])
                *
                (fabsq(elements[i]->nodes[0]->x - elements[i]->nodes[1]->x)*
                fabsq(elements[i]->nodes[0]->y - elements[i]->nodes[3]->y));
    }

    Tm /= 4.0q*lx*ly;

}

tInteger Diffusion2DpAR::direction(tInteger position, DirectionType dir)
{
    switch (dir) {
    case North:
        if(position+nx > nx*ny)
            throw std::string("error: invalid index");
        else
            return position+nx;
        break;
    case South:
        if(position-nx < 0)
            throw std::string("error: invalid index");
        else
            return position-nx;
        break;
    case West:
        if(position-1 < 0)
            throw std::string("error: invalid index");
        else
            return position-1;
        break;
    case East:
        if(position+1 > nx*ny)
            throw std::string("error: invalid index");
        else
            return position+1;
        break;
    default:
        throw std::string("error: invalid argument");
        break;
    }
}

Node2D* Diffusion2DpAR::nodeDirection(tInteger position, DirectionType dir)
{
    switch (dir) {
    case North:
        if(position+nx > nx*ny)
            return NULL;
        else
            return nodes[position+nx];
        break;
    case South:
        if(position-nx < 0)
            return NULL;
        else
            return nodes[position-nx];
        break;
    case West:
        if(!(position%nx))
            return NULL;
        else
            return nodes[position-1];
        break;
    case East:
        if(!((position+1)%nx))
            return NULL;
        else
            return nodes[position+1];
        break;
    default:
        return NULL;
        break;
    }
}

void Diffusion2DpAR::solver(tInteger iterationMax, tFloat iterationTolerance, bool plotlog)
{
    //    refine(nodes[82]);

    //    refine(elements[63]);
    //    refine(elements[73]);
    //    refine(elements[66]);
    //    refine(elements[76]);

    //    refine(elements[72]);
    //    refine(elements[82]);
    //    refine(elements[77]);
    //    refine(elements[87]);

//    tInteger ne = nelements;

//    for(tInteger i=0; i<ne; i++)
//            refine(elements[i]);


    // Condições de Contorno
    // Sul
    Node2D *ptr;
    for(tInteger i=0; i<nx-1; i++){
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }
    // East
    for(tInteger i=1; i<ny; i++)
    {
        ptr = nodes[i*nx-1];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // West
    for(tInteger i=0; i<ny-1; i++)
    {
        ptr = nodes[i*nx];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // Norte
    for(tInteger i=(ny-1)*nx+nx/2; i>(ny-1)*nx; i--)
    {
        ptr = nodes[i];
        while(ptr->nodes[3]->boundary == NULL){
            ptr->nodes[3]->boundary = ptr->boundary;
            ptr = ptr->nodes[3];
        }
    }

    for(tInteger i=(ny-1)*nx+nx/2; i<nx*ny-1; i++)
    {
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }


    GaussSeidel sys2(T, nnodes, iterationMax, iterationTolerance, 10);

    tFloat ap, he, hw, hn, hs, chx, chy;

    for(tInteger i =0; i<nnodes; i++)
    {
        if(nodes[i]->boundary == NULL) // não está no contorno
        {
            ap = 0.0q;

            // Diff2x
            if(nodes[i]->nodes[3] != NULL && nodes[i]->nodes[1] != NULL){ //CDS

                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);

                if(fabsq(he-hw)<1.0e-10q){ //he == hw
                    chx = 1.0q/(he*he);

                    sys2(i, nodes[i]->nodes[1]->index, chx);
                    sys2(i, nodes[i]->nodes[3]->index, chx);
                    ap += -2.0q*chx;
                }
                else{
                    chx = 2.0q/(he*he*hw + hw*hw*he);
                    sys2(i, nodes[i]->nodes[1]->index, hw*chx);
                    sys2(i, nodes[i]->nodes[3]->index, he*chx);
                    ap += -(he+hw)*chx;
                }
            }
            else if(nodes[i]->nodes[1] != NULL){ //UDS
                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[2]->nodes[3]->x);
                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                chx = 2.0q/(he*he*hw + hw*hw*he);

                sys2(i, nodes[i]->nodes[1]->index, hw*chx);
                sys2(i, nodes[i]->nodes[2]->nodes[3]->index, hs*he*chx/(hn+hs));
                sys2(i, nodes[i]->nodes[0]->nodes[3]->index, hn*he*chx/(hn+hs));
                ap += -(he+hw)*chx;
            }
            else if(nodes[i]->nodes[3] != NULL){ //DDS
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);
                he = fabsq(nodes[i]->x - nodes[i]->nodes[2]->nodes[1]->x);
                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                chx = 2.0q/(he*he*hw + hw*hw*he);

                sys2(i, nodes[i]->nodes[3]->index, he*chx);
                sys2(i, nodes[i]->nodes[2]->nodes[1]->index, hs*hw*chx/(hn+hs));
                sys2(i, nodes[i]->nodes[0]->nodes[1]->index, hn*hw*chx/(hn+hs));
                ap += -(he+hw)*chx;
            }

            // Diff2y
            if(nodes[i]->nodes[0] != NULL && nodes[i]->nodes[2] != NULL){

                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                if(fabsq(hn-hs)<1.0e-10q){ //hn == hs
                    chy = 1.0q/(hn*hn);

                    sys2(i, nodes[i]->nodes[0]->index, chy);
                    sys2(i, nodes[i]->nodes[2]->index, chy);
                    ap += -2.0q*chy;
                }
                else{
                    chy = 2.0q/(hn*hn*hs + hs*hs*hn);
                    sys2(i, nodes[i]->nodes[0]->index, hn*chy);
                    sys2(i, nodes[i]->nodes[2]->index, hs*chy);
                    ap += -(hn+hs)*chy;
                }

            }
            else if(nodes[i]->nodes[0] != NULL){ //UDS
                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);
                hn = fabsq(nodes[i]->y - nodes[i]->nodes[1]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                chy = 2.0q/(hn*hn*hs + hs*hs*hn);

                sys2(i, nodes[i]->nodes[0]->index, hn*chy);
                sys2(i, nodes[i]->nodes[1]->nodes[2]->index, hw*hs*chy/(he+hw));
                sys2(i, nodes[i]->nodes[3]->nodes[2]->index, he*hs*chy/(he+hw));
                ap += -(hs+hn)*chy;
            }
            else if(nodes[i]->nodes[2] != NULL){ //DDS
                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);
                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[1]->nodes[0]->y);

                chy = 2.0q/(hn*hn*hs + hs*hs*hn);

                sys2(i, nodes[i]->nodes[2]->index, hs*chy);
                sys2(i, nodes[i]->nodes[1]->nodes[0]->index, hw*hn*chy/(he+hw));
                sys2(i, nodes[i]->nodes[3]->nodes[0]->index, he*hn*chy/(he+hw));
                ap += -(hs+hn)*chy;
            }

            //std::cout<<"\n"<<i<<"\t"<<QtoD(chx)<<"\t"<<QtoD(chy)<<"\t"<<QtoD(ap);

            sys2(i, i, ap);
            sys2(i, data->heatSource->operator ()(nodes[i]->x, nodes[i]->y));
        }
        else
        {
            sys2(i, nodes[i]->boundary->bcValue->operator ()(nodes[i]->x, nodes[i]->y));
        }
    }


    sys2.solver();

}


void Diffusion2DpAR::refine(Element2D *element)
{
    tInteger id[5];
    tInteger refnode;

    id[0] = nnodes++;
    nodes[id[0]] = new Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[1]->x), 0.5q*(element->nodes[0]->y + element->nodes[3]->y));

    // South
    if(element->nodes[0]->nodes[1] == element->nodes[1])
    {
        id[1] = nnodes++;
        nodes[id[1]] = new Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[1]->x), 0.5q*(element->nodes[0]->y + element->nodes[1]->y));
        nodes[id[1]]->link(nodes[id[0]], North);
        nodes[id[1]]->link(element->nodes[0], West);
        nodes[id[1]]->link(element->nodes[1], East);
    }
    else
    {
        Node2D::referenceNode(element->nodes[0], East, element->nodes[1], &refnode);
        id[1] = refnode;
        nodes[id[1]]->link(nodes[id[0]], North);
    }

    // East
    if(element->nodes[1]->nodes[2] == element->nodes[2])
    {
        id[2] = nnodes++;
        nodes[id[2]] = new Node2D(nnodes-1, 0.5q*(element->nodes[1]->x + element->nodes[2]->x), 0.5q*(element->nodes[1]->y + element->nodes[2]->y));
        nodes[id[2]]->link(nodes[id[0]], West);
        nodes[id[2]]->link(element->nodes[1], South);
        nodes[id[2]]->link(element->nodes[2], North);
    }
    else
    {
        Node2D::referenceNode(element->nodes[1], North, element->nodes[2], &refnode);
        id[2] = refnode;
        nodes[id[2]]->link(nodes[id[0]], West);
    }


    //North
    if(element->nodes[3]->nodes[1] == element->nodes[2])
    {
        id[3] = nnodes++;
        nodes[id[3]] = new Node2D(nnodes-1, 0.5q*(element->nodes[3]->x + element->nodes[2]->x), 0.5q*(element->nodes[3]->y + element->nodes[2]->y));
        nodes[id[3]]->link(nodes[id[0]], South);
        nodes[id[3]]->link(element->nodes[3], West);
        nodes[id[3]]->link(element->nodes[2], East);
    }
    else
    {
        Node2D::referenceNode(element->nodes[3], East, element->nodes[2], &refnode);
        id[3] = refnode;
        nodes[id[3]]->link(nodes[id[0]], South);
    }


    // West
    if(element->nodes[0]->nodes[2] == element->nodes[3])
    {
        id[4] = nnodes++;
        nodes[id[4]] = new Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[3]->x), 0.5q*(element->nodes[0]->y + element->nodes[3]->y));
        nodes[id[4]]->link(nodes[id[0]], East);
        nodes[id[4]]->link(element->nodes[0], South);
        nodes[id[4]]->link(element->nodes[3], North);
    }
    else
    {
        Node2D::referenceNode(element->nodes[0], North, element->nodes[3], &refnode);
        id[4] = refnode;
        nodes[id[4]]->link(nodes[id[0]], East);
    }

    elements[nelements++] = new Element2D(nelements, element->nodes[0], nodes[id[1]], nodes[id[0]], nodes[id[4]]);
    elements[nelements++] = new Element2D(nelements, nodes[id[1]], element->nodes[1], nodes[id[2]], nodes[id[0]]);
    elements[nelements++] = new Element2D(nelements, nodes[id[0]], nodes[id[2]], element->nodes[2], nodes[id[3]]);
    elements[nelements++] = new Element2D(nelements, nodes[id[4]], nodes[id[0]], nodes[id[3]], element->nodes[3]);


    element->index = -1;
}

void Diffusion2DpAR::refine(Node2D *node)
{
    tInteger ne = nelements;

    for(tInteger i=0; i<ne; i++)
        if(elements[i]->nodes[0] == node ||
                elements[i]->nodes[1] == node ||
                elements[i]->nodes[2] == node ||
                elements[i]->nodes[3] == node)
            if(elements[i]->index!= -1)
                refine(elements[i]);
}

void Diffusion2DpAR::solver2(tInteger iterationMax, tFloat iterationTolerance, bool plotlog)
{
    Node2D *ptr;

    refine(nodes[82]);

    // Condições de Contorno
    // Sul
    for(tInteger i=0; i<nx-1; i++){
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }
    // East
    for(tInteger i=1; i<ny; i++)
    {
        ptr = nodes[i*nx-1];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // West
    for(tInteger i=0; i<ny-1; i++)
    {
        ptr = nodes[i*nx];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // Norte
    for(tInteger i=(ny-1)*nx+nx/2; i>(ny-1)*nx; i--)
    {
        ptr = nodes[i];
        while(ptr->nodes[3]->boundary == NULL){
            ptr->nodes[3]->boundary = ptr->boundary;
            ptr = ptr->nodes[3];
        }
    }

    for(tInteger i=(ny-1)*nx+nx/2; i<nx*ny-1; i++)
    {
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }


    //GaussSeidel sys2(T, nnodes, iterationMax, iterationTolerance, 15);

    SolverLU sys2(T, nnodes);

    tFloat ap, he, hw, hn, hs, chx, chy;

    for(tInteger i =0; i<nnodes; i++)
    {
        if(nodes[i]->boundary == NULL) // não está no contorno
        {
            ap = 0.0q;

            // Diff2x
            if(nodes[i]->nodes[3] != NULL && nodes[i]->nodes[1] != NULL){ //CDS

                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);

                if(fabsq(he-hw)<1.0e-10q){ //he == hw
                    chx = 1.0q/(he*he);

                    sys2(i, nodes[i]->nodes[1]->index, chx);
                    sys2(i, nodes[i]->nodes[3]->index, chx);
                    ap += -2.0q*chx;
                }
                else{
                    chx = 2.0q/(he*he*hw + hw*hw*he);
                    sys2(i, nodes[i]->nodes[1]->index, hw*chx);
                    sys2(i, nodes[i]->nodes[3]->index, he*chx);
                    ap += -(he+hw)*chx;
                }
            }
            else if(nodes[i]->nodes[1] != NULL){ //UDS
                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);

                chx = 1.0q/(he*he);

                sys2(i, nodes[i]->nodes[1]->index, -26.q/3.q*chx);
                sys2(i, nodes[i]->nodes[1]->nodes[1]->index, 19.q/2.q*chx);
                sys2(i, nodes[i]->nodes[1]->nodes[1]->nodes[1]->index, -14.q/3.q*chx);
                sys2(i, nodes[i]->nodes[1]->nodes[1]->nodes[1]->nodes[1]->index, 11.q/12.q*chx);
                ap += 35.q/12.q*chx;
            }
            else if(nodes[i]->nodes[3] != NULL){ //DDS
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);

                chx = 1.0q/(hw*hw);

                sys2(i, nodes[i]->nodes[3]->index, -26.q/3.q*chx);
                sys2(i, nodes[i]->nodes[3]->nodes[3]->index, 19.q/2.q*chx);
                sys2(i, nodes[i]->nodes[3]->nodes[3]->nodes[3]->index, -14.q/3.q*chx);
                sys2(i, nodes[i]->nodes[3]->nodes[3]->nodes[3]->nodes[3]->index, 11.q/12.q*chx);
                ap += 35.q/12.q*chx;
            }

            // Diff2y
            if(nodes[i]->nodes[0] != NULL && nodes[i]->nodes[2] != NULL){

                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                if(fabsq(hn-hs)<1.0e-10q){ //hn == hs
                    chy = 1.0q/(hn*hn);

                    sys2(i, nodes[i]->nodes[0]->index, chy);
                    sys2(i, nodes[i]->nodes[2]->index, chy);
                    ap += -2.0q*chy;
                }
                else{
                    chy = 2.0q/(hn*hn*hs + hs*hs*hn);
                    sys2(i, nodes[i]->nodes[0]->index, hn*chy);
                    sys2(i, nodes[i]->nodes[2]->index, hs*chy);
                    ap += -(hn+hs)*chy;
                }

            }
            else if(nodes[i]->nodes[0] != NULL){ //UDS

                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                chy = 1.0q/(hs*hs);

                sys2(i, nodes[i]->nodes[0]->index, -26.q/3.q*chy);
                sys2(i, nodes[i]->nodes[0]->nodes[0]->index, 19.q/2.q*chy);
                sys2(i, nodes[i]->nodes[0]->nodes[0]->nodes[0]->index, -14.q/3.q*chy);
                sys2(i, nodes[i]->nodes[0]->nodes[0]->nodes[0]->nodes[0]->index, 11.q/12.q*chy);
                ap += 35.q/12.q*chy;
            }
            else if(nodes[i]->nodes[2] != NULL){ //DDS

                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);

                chy = 1.0q/(hn*hn);

                sys2(i, nodes[i]->nodes[2]->index, -26.q/3.q*chy);
                sys2(i, nodes[i]->nodes[2]->nodes[2]->index, 19.q/2.q*chy);
                sys2(i, nodes[i]->nodes[2]->nodes[2]->nodes[2]->index, -14.q/3.q*chy);
                sys2(i, nodes[i]->nodes[2]->nodes[2]->nodes[2]->nodes[2]->index, 11.q/12.q*chy);
                ap += 35.q/12.q*chy;
            }

            std::cout<<"\n"<<i<<"\t"<<QtoD(chx)<<"\t"<<QtoD(chy)<<"\t"<<QtoD(ap);

            sys2(i, i, ap);
            sys2(i, data->heatSource->operator ()(nodes[i]->x, nodes[i]->y));
        }
        else
        {
            sys2(i, nodes[i]->boundary->bcValue->operator ()(nodes[i]->x, nodes[i]->y));
        }
    }


    sys2.solver();

}


void Diffusion2DpAR::solver3(tInteger iterationMax, tFloat iterationTolerance, bool plotlog)
{

    // Condições de Contorno

    Node2D *ptr;
    // East
    for(tInteger i=1; i<ny; i++)
    {
        ptr = nodes[i*nx-1];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // West
    for(tInteger i=0; i<ny-1; i++)
    {
        ptr = nodes[i*nx];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // Norte
    for(tInteger i=(ny-1)*nx+nx/2; i>(ny-1)*nx; i--)
    {
        ptr = nodes[i];
        while(ptr->nodes[3]->boundary == NULL){
            ptr->nodes[3]->boundary = ptr->boundary;
            ptr = ptr->nodes[3];
        }
    }

    for(tInteger i=(ny-1)*nx+nx/2; i<nx*ny-1; i++)
    {
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }

    // Sul
    for(tInteger i=nx/2; i<nx-1; i++){
        ptr = nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }

    for(tInteger i=nx/2; i>0; i--){
        ptr = nodes[i];
        while(ptr->nodes[3]->boundary == NULL){
            ptr->nodes[3]->boundary = ptr->boundary;
            ptr = ptr->nodes[3];
        }
    }


    //GaussSeidel sys2(T, nnodes, iterationMax, iterationTolerance, 15);

    SolverLU sys2(T, nnodes);

    tFloat ap, he, hw, hn, hs, chx, chy;

    for(tInteger i =0; i<nnodes; i++)
    {
        if(nodes[i]->boundary == NULL) // não está no contorno
        {
            ap = 0.0q;

            // Diff2x
            if(nodes[i]->nodes[3] != NULL && nodes[i]->nodes[1] != NULL){ //CDS

                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);

                if(fabsq(he-hw)<1.0e-10q){ //he == hw
                    chx = 1.0q/(he*he);

                    sys2(i, nodes[i]->nodes[1]->index, chx);
                    sys2(i, nodes[i]->nodes[3]->index, chx);
                    ap += -2.0q*chx;
                }
                else{
                    chx = 2.0q/(he*he*hw + hw*hw*he);
                    sys2(i, nodes[i]->nodes[1]->index, hw*chx);
                    sys2(i, nodes[i]->nodes[3]->index, he*chx);
                    ap += -(he+hw)*chx;
                }
            }
            else if(nodes[i]->nodes[1] != NULL){ //UDS
                he = fabsq(nodes[i]->x - nodes[i]->nodes[1]->x);

                chx = 1.0q/(he*he);

                sys2(i, nodes[i]->nodes[1]->index, -2.q*chx);
                sys2(i, nodes[i]->nodes[1]->nodes[1]->index, 1.q*chx);
                ap += chx;
            }
            else if(nodes[i]->nodes[3] != NULL){ //DDS
                hw = fabsq(nodes[i]->x - nodes[i]->nodes[3]->x);

                chx = 1.0q/(hw*hw);

                sys2(i, nodes[i]->nodes[3]->index, -2.q*chx);
                sys2(i, nodes[i]->nodes[3]->nodes[3]->index, 1.q*chx);
                ap += chx;
            }

            // Diff2y
            if(nodes[i]->nodes[0] != NULL && nodes[i]->nodes[2] != NULL){

                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);
                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                if(fabsq(hn-hs)<1.0e-10q){ //hn == hs
                    chy = 1.0q/(hn*hn);

                    sys2(i, nodes[i]->nodes[0]->index, chy);
                    sys2(i, nodes[i]->nodes[2]->index, chy);
                    ap += -2.0q*chy;
                }
                else{
                    chy = 2.0q/(hn*hn*hs + hs*hs*hn);
                    sys2(i, nodes[i]->nodes[0]->index, hn*chy);
                    sys2(i, nodes[i]->nodes[2]->index, hs*chy);
                    ap += -(hn+hs)*chy;
                }

            }
            else if(nodes[i]->nodes[0] != NULL){ //UDS

                hs = fabsq(nodes[i]->y - nodes[i]->nodes[0]->y);

                chy = 1.0q/(hs*hs);

                sys2(i, nodes[i]->nodes[0]->index, -2.q*chy);
                sys2(i, nodes[i]->nodes[0]->nodes[0]->index, 1.q*chy);
                ap += chy;
            }
            else if(nodes[i]->nodes[2] != NULL){ //DDS

                hn = fabsq(nodes[i]->y - nodes[i]->nodes[2]->y);

                chy = 1.0q/(hn*hn);


                sys2(i, nodes[i]->nodes[2]->index, -2.q*chy);
                sys2(i, nodes[i]->nodes[2]->nodes[2]->index, 1.q*chy);
                ap += chy;
            }

            //std::cout<<"\n"<<i<<"\t"<<QtoD(chx)<<"\t"<<QtoD(chy)<<"\t"<<QtoD(ap);

            sys2(i, i, ap);
            sys2(i, data->heatSource->operator ()(nodes[i]->x, nodes[i]->y));
        }
        else
        {
            sys2(i, nodes[i]->boundary->bcValue->operator ()(nodes[i]->x, nodes[i]->y));
        }
    }


    sys2.solver();

}

