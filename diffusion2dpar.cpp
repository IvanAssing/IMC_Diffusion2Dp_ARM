#include "diffusion2dpar.h"

#include "gaussseidel.h"
#include "solverlu.h"
#include "stdlib.h"

#define MAXNODES(nx, ny) 100*nx*ny
#define MAXELEMENTS(nx, ny) 50*nx*ny

Diffusion2DpAR::Diffusion2DpAR(tFloat _lx, tFloat _ly, tInteger _nx, tInteger _ny, Diffusion2DData *_data,
                               Boundary2D *_bS, Boundary2D *_bN, Boundary2D *_bE, Boundary2D *_bW)
    :lx(_lx), ly(_ly), nx(_nx), ny(_ny), data(_data), ccS(_bS), ccN(_bN), ccE(_bE), ccW(_bW)
{
    hx = lx/(nx - 1.0q);
    hy = ly/(ny - 1.0q);

    // Nodes
    nodes = new Node2D[MAXNODES(nx,ny)]();
    tInteger p = 0;


    for(tInteger j=0; j<ny; j++)
        for(tInteger i=0; i<nx; i++)
            nodes[p] = Node2D(p++, i*hx, j*hy,
                              nodeDirection(p, South),
                              nodeDirection(p, East),
                              nodeDirection(p, North),
                              nodeDirection(p, West));

    nnodes = p;

    // Elements
    elements = new Element2D[MAXELEMENTS(nx,ny)]();
    p = 0;

    for(tInteger j=0; j<ny-1; j++)
        for(tInteger i=0; i<nx-1; i++){
            tInteger nn = position(i, j);
            elements[p] = Element2D(p++, &nodes[nn],
                                    nodeDirection(nn, East),
                                    nodeDirection(direction(nn, East), North),
                                    nodeDirection(nn, North));
        }

    nelements = p;

    T = new tFloat[nx*ny];

    for(tInteger i=0; i<nx*ny; i++)
        T[i] = 0.q;


    // Condições de Contorno
    // Sul
    for(tInteger i=0; i<nx; i++)
        nodes[i].boundary = ccS;

    // Norte
    for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
        nodes[i].boundary = ccN;

    // East
    for(tInteger i=1; i<ny; i++)
        nodes[i*nx-1].boundary = ccE;

    // West
    for(tInteger i=0; i<ny; i++)
        nodes[i*nx].boundary = ccW;

}


tInteger Diffusion2DpAR::position(tInteger i, tInteger j)
{
    return i+j*nx;
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
            return &nodes[position+nx];
        break;
    case South:
        if(position-nx < 0)
            return NULL;
        else
            return &nodes[position-nx];
        break;
    case West:
        if(!(position%nx))
            return NULL;
        else
            return &nodes[position-1];
        break;
    case East:
        if(!((position+1)%nx))
            return NULL;
        else
            return &nodes[position+1];
        break;
    default:
        return NULL;
        break;
    }
}

void Diffusion2DpAR::solver(tInteger iterationMax, tFloat iterationTolerance, bool plotlog)
{

    //GaussSeidel sys(T, nnodes, iterationMax, iterationTolerance);
    SolverLU sys(T, nnodes);

    tFloat chx, chy;

    for(tInteger i =0; i<nnodes; i++)
        if(nodes[i].boundary == NULL) // não está no contorno
        {

            // Diff2x
            if(nodes[i].nodes[3] != NULL && nodes[i].nodes[1] != NULL){ //CDS
                chx = 2.0q/(nodes[i].nodes[1]->x - nodes[i].nodes[3]->x);
                sys(i, nodes[i].nodes[1]->index, chx*chx);
                sys(i, nodes[i].nodes[3]->index, chx*chx);
            }
            else if(nodes[i].nodes[1] != NULL){ //UDS
                chx = 1.0q/(-nodes[i].x + nodes[i].nodes[1]->x);
                sys(i, nodes[i].index, 2.0q*chx); //Tp
                sys(i, nodes[i].nodes[1]->index, -5.0q*chx); //Tw
                sys(i, nodes[i].nodes[1]->nodes[1]->index, 4.0q*chx); //Tww
                sys(i, nodes[i].nodes[1]->nodes[1]->nodes[1]->index, -1.0q*chx); //Twww
            }
            else if(nodes[i].nodes[3] != NULL){ //DDS
                chx = 1.0q/(+nodes[i].x - nodes[i].nodes[3]->x);
                sys(i, nodes[i].index, 2.0q*chx); //Tp
                sys(i, nodes[i].nodes[3]->index, -5.0q*chx); //Te
                sys(i, nodes[i].nodes[3]->nodes[3]->index, 4.0q*chx); //Tee
                sys(i, nodes[i].nodes[3]->nodes[3]->nodes[3]->index, -1.0q*chx); //Teee
            }

            // Diff2y
            if(nodes[i].nodes[0] != NULL && nodes[i].nodes[2] != NULL){
                chy = 2.0q/(nodes[i].nodes[2]->y - nodes[i].nodes[0]->y);
                sys(i, nodes[i].nodes[0]->index, chy*chy);
                sys(i, nodes[i].nodes[2]->index, chy*chy);
            }
            else if(nodes[i].nodes[0] != NULL){ //UDS
                chy = 1.0q/(+nodes[i].y - nodes[i].nodes[0]->y);
                sys(i, nodes[i].index, 2.0q*chy); //Tp
                sys(i, nodes[i].nodes[0]->index, -5.0q*chy); //Tw
                sys(i, nodes[i].nodes[0]->nodes[0]->index, 4.0q*chy); //Tww
                sys(i, nodes[i].nodes[0]->nodes[0]->nodes[0]->index, -1.0q*chy); //Twww
            }
            else if(nodes[i].nodes[2] != NULL){ //DDS
                chy = 1.0q/(-nodes[i].y + nodes[i].nodes[2]->y);
                sys(i, nodes[i].index, 2.0q*chy); //Tp
                sys(i, nodes[i].nodes[2]->index, -5.0q*chy); //Te
                sys(i, nodes[i].nodes[2]->nodes[2]->index, 4.0q*chy); //Tee
                sys(i, nodes[i].nodes[2]->nodes[2]->nodes[2]->index, -1.0q*chy); //Teee
            }

            std::cout<<"\n"<<i<<"\t"<<QtoD(1.0q/chx)<<"\t"<<QtoD(1.0q/chy);
            sys(i, i, -2.0q*(chx*chx+chy*chy));
            sys(i, data->heatSource->operator ()(nodes[i].x, nodes[i].y));
        }
        else
        {
            sys(i, nodes[i].boundary->bcValue->operator ()(nodes[i].x, nodes[i].y));
        }


    sys.solver();

    //

    //    Tm = 0.0q;
    //    int p;
    //    for(tInteger j=0; j<ny-1; j++)
    //        for(tInteger i=0; i<nx-1; i++){
    //            p = position(i,j);
    //            Tm += T[p] + T[direction(p, East)] + T[direction(p, North)] + T[direction(direction(p, North), East)];
    //        }
    //    Tm *= hx*hy/(4.0q*lx*ly);

    //    if(plotlog)
    //        sys.plotIterationLog();
    tFloat Tref = T[82];

    tInteger ne = nelements;

    //    for(int i = 0; i<ne; i++)
    //        refine(&elements[i]);

    ne = nelements;

    //        for(int i = 0; i<4*ne; i++)
    //            refine(&elements[i]);

    //    for(tInteger j=2; j<ny-2; j++)
    //        for(tInteger i=2; i<nx-2; i++)
    //            refine(&elements[position(i,j)]);

    //refine(&nodes[104]);

    //refine(&nodes[15]);

    //refine(&elements[62]);
    //    refine(&nodes[84]);
    //    refine(&nodes[104]);
    //    refine(&nodes[106]);
    refine(&nodes[82]);

    //        refine(&nodes[82]);
    //        refine(&nodes[104]);
    //    refine(&nodes[59]);

    //    refine(&nodes[10])*/


    // Condições de Contorno
    // Sul
    Node2D *ptr;
    for(tInteger i=0; i<nx-1; i++){
        ptr = &nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }

    // Norte
    for(tInteger i=(ny-1)*nx; i<nx*ny-1; i++)
    {
        ptr = &nodes[i];
        while(ptr->nodes[1]->boundary == NULL){
            ptr->nodes[1]->boundary = ptr->boundary;
            ptr = ptr->nodes[1];
        }
    }

    // East
    for(tInteger i=1; i<ny; i++)
    {
        ptr = &nodes[i*nx-1];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    // West
    for(tInteger i=0; i<ny-1; i++)
    {
        ptr = &nodes[i*nx];
        while(ptr->nodes[2]->boundary == NULL){
            ptr->nodes[2]->boundary = ptr->boundary;
            ptr = ptr->nodes[2];
        }
    }

    for(tInteger i=0; i<nnodes; i++)
        if(nodes[i].nodes[2] == NULL)
            nodes[i].boundary = ccN;

    /////////////////////////////////////////////////

    nodes[132].boundary = NULL;
    nodes[136].boundary = NULL;

    SolverLU sys2(T, nnodes);

    tFloat ap;
    tFloat he, hw, hn, hs;

    for(tInteger i =0; i<nnodes; i++)
    {
        if(i==110 || i==121)
            ap = 0.0q;

        if(nodes[i].boundary == NULL) // não está no contorno
        {
            ap = 0.0q;



            // Diff2x
            if(nodes[i].nodes[3] != NULL && nodes[i].nodes[1] != NULL){ //CDS

                he = fabsq(nodes[i].x - nodes[i].nodes[1]->x);
                hw = fabsq(nodes[i].x - nodes[i].nodes[3]->x);

                if(fabsq(he-hw)<1.0e-10q){ //he == hw
                    chx = 1.0q/(he*he);

                    sys2(i, nodes[i].nodes[1]->index, chx);
                    sys2(i, nodes[i].nodes[3]->index, chx);
                    ap += -2.0q*chx;
                }
                else{
                    chx = 2.0q/(he*he*hw + hw*hw*he);
                    sys2(i, nodes[i].nodes[1]->index, hw*chx);
                    sys2(i, nodes[i].nodes[3]->index, he*chx);
                    ap += -(he+hw)*chx;
                }
            }
            else if(nodes[i].nodes[1] != NULL){ //UDS
                he = fabsq(nodes[i].x - nodes[i].nodes[1]->x);
                chx = 1.0q/(he*he);

                sys2(i, nodes[i].nodes[1]->index, -5.0q*chx); //Tw
                sys2(i, nodes[i].nodes[1]->nodes[1]->index, +4.0q*chx); //Tww
                sys2(i, nodes[i].nodes[1]->nodes[1]->nodes[1]->index, -1.0q*chx); //Twww
                ap += 2.0q*chx;
            }
            else if(nodes[i].nodes[3] != NULL){ //DDS
                hw = fabsq(nodes[i].x - nodes[i].nodes[3]->x);
                chx = 1.0q/(hw*hw);

                sys2(i, nodes[i].nodes[3]->index, -5.0q*chx); //Te
                sys2(i, nodes[i].nodes[3]->nodes[3]->index, +4.0q*chx); //Tee
                sys2(i, nodes[i].nodes[3]->nodes[3]->nodes[3]->index, -1.0q*chx); //Teee
                ap += 2.0q*chx;
            }

            // Diff2y
            if(nodes[i].nodes[0] != NULL && nodes[i].nodes[2] != NULL){

                hn = fabsq(nodes[i].y - nodes[i].nodes[2]->y);
                hs = fabsq(nodes[i].y - nodes[i].nodes[0]->y);

                if(fabsq(hn-hs)<1.0e-10q){ //hn == hs
                    chy = 1.0q/(hn*hn);

                    sys2(i, nodes[i].nodes[0]->index, chy);
                    sys2(i, nodes[i].nodes[2]->index, chy);
                    ap += -2.0q*chy;
                }
                else{
                    chy = 2.0q/(hn*hn*hs + hs*hs*hn);
                    sys2(i, nodes[i].nodes[0]->index, hn*chy);
                    sys2(i, nodes[i].nodes[2]->index, hs*chy);
                    ap += -(hn+hs)*chy;
                }

            }
            else if(nodes[i].nodes[0] != NULL){ //UDS
                hs = fabsq(nodes[i].y - nodes[i].nodes[0]->y);
                chy = 1.0q/(hs*hs);

                sys2(i, nodes[i].nodes[0]->index, -5.0q*chy); //Tw
                sys2(i, nodes[i].nodes[0]->nodes[0]->index, +4.0q*chy); //Tww
                sys2(i, nodes[i].nodes[0]->nodes[0]->nodes[0]->index, -1.0q*chy); //Twww
                ap += 2.0q*chy;
            }
            else if(nodes[i].nodes[2] != NULL){ //DDS
                hn = fabsq(nodes[i].y - nodes[i].nodes[2]->y);
                chy = 1.0q/(hn*hn);

                sys2(i, nodes[i].nodes[2]->index, -5.0q*chy); //Te
                sys2(i, nodes[i].nodes[2]->nodes[2]->index, +4.0q*chy); //Tee
                sys2(i, nodes[i].nodes[2]->nodes[2]->nodes[2]->index, -1.0q*chy); //Teee
                ap += 2.0q*chy;
            }

            std::cout<<"\n"<<i<<"\t"<<QtoD(chx)<<"\t"<<QtoD(chy)<<"\t"<<QtoD(ap);

            sys2(i, i, ap);
            sys2(i, data->heatSource->operator ()(nodes[i].x, nodes[i].y));
        }
        else
        {
            sys2(i, nodes[i].boundary->bcValue->operator ()(nodes[i].x, nodes[i].y));
        }
    }


    //    for(tInteger i =0; i<nnodes; i++)
    //        if(nodes[i].boundary != NULL)
    //            sys2(i, nodes[i].boundary->bcValue->operator ()(nodes[i].x, nodes[i].y));

    //sys2.printindex();

    //sys2.diagonalcheck();


    //sys2.printindex();

    //return;
    tInteger neqmax = nnodes;

    tFloat *b = new tFloat[neqmax];

    tFloat **A = new tFloat*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        A[i] = new tFloat[neqmax];

    for(tInteger i = 0; i<neqmax; i++){
        b[i] = sys2.b[i];
        for(tInteger j = 0; j<neqmax; j++)
            A[i][j] = sys2.A[i][j];
    }

    sys2.solver();

    tFloat sum, sumtotal = 0.0q;
    for(tInteger k = 0; k<neqmax; k++)
    {
        sum = -b[k];
        for(tInteger i = 0; i<neqmax; i++){
                sum += T[i]* A[k][i];
        }
        sumtotal += sum;
        std::cout<<"\n  "<<k<<"\t\t"<<print(sum)<<"\t\t"<<print(sumtotal);
    }



    //    if(plotlog)
    //        sys2.plotIterationLog();

    std::cout<<"\n****"<<print(Tref)<<"\t"<<print(T[82]);


}


//// Diff2x
//if(nodes[i].nodes[3] != NULL && nodes[i].nodes[1] != NULL){ //CDS
//    chx = 2.0q/(nodes[i].nodes[1]->x - nodes[i].nodes[3]->x);
//    sys2(i, nodes[i].nodes[1]->index, chx);
//    sys2(i, nodes[i].nodes[3]->index, chx);
//}
//else if(nodes[i].nodes[1] != NULL){ //UDS
//    chx = 1.0q/(-nodes[i].x + nodes[i].nodes[1]->x);
//    sys2(i, nodes[i].nodes[1]->index, +5.0q*chx); //Tw
//    sys2(i, nodes[i].nodes[1]->nodes[1]->index, -4.0q*chx); //Tww
//    sys2(i, nodes[i].nodes[1]->nodes[1]->nodes[1]->index, 1.0q*chx); //Twww
//}
//else if(nodes[i].nodes[3] != NULL){ //DDS
//    chx = -1.0q/(+nodes[i].x - nodes[i].nodes[3]->x);
//    sys2(i, nodes[i].nodes[3]->index, 5.0q*chx); //Te
//    sys2(i, nodes[i].nodes[3]->nodes[3]->index, -4.0q*chx); //Tee
//    sys2(i, nodes[i].nodes[3]->nodes[3]->nodes[3]->index, 1.0q*chx); //Teee
//}

//// Diff2y
//if(nodes[i].nodes[0] != NULL && nodes[i].nodes[2] != NULL){
//    chy = 2.0q/(nodes[i].nodes[2]->y - nodes[i].nodes[0]->y);
//    sys2(i, nodes[i].nodes[0]->index, chy);
//    sys2(i, nodes[i].nodes[2]->index, chy);
//}
//else if(nodes[i].nodes[0] != NULL){ //UDS
//    chy = -1.0q/(+nodes[i].y - nodes[i].nodes[0]->y);
//    sys2(i, nodes[i].nodes[0]->index, +5.0q*chy); //Tw
//    sys2(i, nodes[i].nodes[0]->nodes[0]->index, -4.0q*chy); //Tww
//    sys2(i, nodes[i].nodes[0]->nodes[0]->nodes[0]->index, +1.0q*chy); //Twww
//}
//else if(nodes[i].nodes[2] != NULL){ //DDS
//    chy = -1.0q/(-nodes[i].y + nodes[i].nodes[2]->y);
//    sys2(i, nodes[i].nodes[2]->index, +5.0q*chy); //Te
//    sys2(i, nodes[i].nodes[2]->nodes[2]->index, -4.0q*chy); //Tee
//    sys2(i, nodes[i].nodes[2]->nodes[2]->nodes[2]->index, +1.0q*chy); //Teee
//}



//// Diff2x
//if(nodes[i].nodes[3] != NULL && nodes[i].nodes[1] != NULL){ //CDS
//    chx = 2.0q/(nodes[i].nodes[1]->x - nodes[i].nodes[3]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[1]->index, chx);
//    sys2(i, nodes[i].nodes[3]->index, chx);
//    ap += -2.0q*chx;
//}
//else if(nodes[i].nodes[1] != NULL){ //UDS
//    chx = 1.0q/(-nodes[i].x + nodes[i].nodes[1]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[1]->index, -26.0q/3.0q*chx); //Tw
//    sys2(i, nodes[i].nodes[1]->nodes[1]->index, +19.0q/2.0q*chx); //Tww
//    sys2(i, nodes[i].nodes[1]->nodes[1]->nodes[1]->index, -14.0q/3.0q*chx); //Twww
//    sys2(i, nodes[i].nodes[1]->nodes[1]->nodes[1]->nodes[1]->index, +11.0q/12.0q*chx); //Twwww
//    ap += 35.0q/12.0q*chx;
//}
//else if(nodes[i].nodes[3] != NULL){ //DDS
//    chx = 1.0q/(+nodes[i].x - nodes[i].nodes[3]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[3]->index, -26.0q/3.0q*chx); //Te
//    sys2(i, nodes[i].nodes[3]->nodes[3]->index, +19.0q/2.0q*chx); //Tee
//    sys2(i, nodes[i].nodes[3]->nodes[3]->nodes[3]->index, -14.0q/3.0q*chx); //Teee
//    sys2(i, nodes[i].nodes[3]->nodes[3]->nodes[3]->nodes[3]->index, +11.0q/12.0q*chx); //Teeee
//    ap += 35.0q/12.0q*chx;
//}

//// Diff2y
//if(nodes[i].nodes[0] != NULL && nodes[i].nodes[2] != NULL){
//    chy = 2.0q/(nodes[i].nodes[2]->y - nodes[i].nodes[0]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[0]->index, chy);
//    sys2(i, nodes[i].nodes[2]->index, chy);
//    ap += -2.0q*chy;
//}
//else if(nodes[i].nodes[0] != NULL){ //UDS
//    chy = 1.0q/(+nodes[i].y - nodes[i].nodes[0]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[0]->index, -26.0q/3.0q*chy); //Tw
//    sys2(i, nodes[i].nodes[0]->nodes[0]->index, +19.0q/2.0q*chy); //Tww
//    sys2(i, nodes[i].nodes[0]->nodes[0]->nodes[0]->index, -14.0q/3.0q*chy); //Twww
//    sys2(i, nodes[i].nodes[0]->nodes[0]->nodes[0]->nodes[0]->index, +11.0q/12.0q*chy); //Twwww
//    ap += 35.0q/12.0q*chy;
//}
//else if(nodes[i].nodes[2] != NULL){ //DDS
//    chy = 1.0q/(-nodes[i].y + nodes[i].nodes[2]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[2]->index, -26.0q/3.0q*chy); //Te
//    sys2(i, nodes[i].nodes[2]->nodes[2]->index, +19.0q/2.0q*chy); //Tee
//    sys2(i, nodes[i].nodes[2]->nodes[2]->nodes[2]->index, -14.0q/3.0q*chy); //Teee
//    sys2(i, nodes[i].nodes[2]->nodes[2]->nodes[2]->nodes[2]->index, +11.0q/12.0q*chy); //Teeee
//    ap += 35.0q/12.0q*chy;
//}


//// Diff2x
//if(nodes[i].nodes[3] != NULL && nodes[i].nodes[1] != NULL){ //CDS
//    chx = 2.0q/(nodes[i].nodes[1]->x - nodes[i].nodes[3]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[1]->index, chx);
//    sys2(i, nodes[i].nodes[3]->index, chx);
//    ap += -2.0q*chx;
//}
//else if(nodes[i].nodes[1] != NULL){ //UDS
//    chx = 1.0q/(-nodes[i].x + nodes[i].nodes[1]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[1]->index, -2.0q*chx); //Tw
//    sys2(i, nodes[i].nodes[1]->nodes[1]->index, chx); //Tww
//    ap += chx;
//}
//else if(nodes[i].nodes[3] != NULL){ //DDS
//    chx = 1.0q/(+nodes[i].x - nodes[i].nodes[3]->x);
//    chx *= chx;
//    sys2(i, nodes[i].nodes[3]->index, -2.0q*chx); //Te
//    sys2(i, nodes[i].nodes[3]->nodes[3]->index, chx); //Tee
//    ap += chx;
//}

//// Diff2y
//if(nodes[i].nodes[0] != NULL && nodes[i].nodes[2] != NULL){
//    chy = 2.0q/(nodes[i].nodes[2]->y - nodes[i].nodes[0]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[0]->index, chy);
//    sys2(i, nodes[i].nodes[2]->index, chy);
//    ap += -2.0q*chy;
//}
//else if(nodes[i].nodes[0] != NULL){ //UDS
//    chy = 1.0q/(+nodes[i].y - nodes[i].nodes[0]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[0]->index, -2.0q*chy); //Tw
//    sys2(i, nodes[i].nodes[0]->nodes[0]->index, chy); //Tww
//    ap += chy;
//}
//else if(nodes[i].nodes[2] != NULL){ //DDS
//    chy = 1.0q/(-nodes[i].y + nodes[i].nodes[2]->y);
//    chy *= chy;
//    sys2(i, nodes[i].nodes[2]->index, -2.0q*chy); //Te
//    sys2(i, nodes[i].nodes[2]->nodes[2]->index, chy); //Tee
//    ap += chy;
//}


void Diffusion2DpAR::refine(Element2D *element)
{
    tInteger id[5];
    tInteger refnode;

    id[0] = nnodes++;
    nodes[id[0]] = Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[1]->x), 0.5q*(element->nodes[0]->y + element->nodes[3]->y));

    // South
    if(element->nodes[0]->nodes[1] == element->nodes[1])
    {
        id[1] = nnodes++;
        nodes[id[1]] = Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[1]->x), 0.5q*(element->nodes[0]->y + element->nodes[1]->y));
        nodes[id[1]].link(&nodes[id[0]], North);
        nodes[id[1]].link(element->nodes[0], West);
        nodes[id[1]].link(element->nodes[1], East);
    }
    else
    {
        Node2D::referenceNode(element->nodes[0], East, element->nodes[1], &refnode);
        id[1] = refnode;
        nodes[id[1]].link(&nodes[id[0]], North);
    }

    // East
    if(element->nodes[1]->nodes[2] == element->nodes[2])
    {
        id[2] = nnodes++;
        nodes[id[2]] = Node2D(nnodes-1, 0.5q*(element->nodes[1]->x + element->nodes[2]->x), 0.5q*(element->nodes[1]->y + element->nodes[2]->y));
        nodes[id[2]].link(&nodes[id[0]], West);
        nodes[id[2]].link(element->nodes[1], South);
        nodes[id[2]].link(element->nodes[2], North);
    }
    else
    {
        Node2D::referenceNode(element->nodes[1], North, element->nodes[2], &refnode);
        id[2] = refnode;
        nodes[id[2]].link(&nodes[id[0]], West);
    }


    //North
    if(element->nodes[3]->nodes[1] == element->nodes[2])
    {
        id[3] = nnodes++;
        nodes[id[3]] = Node2D(nnodes-1, 0.5q*(element->nodes[3]->x + element->nodes[2]->x), 0.5q*(element->nodes[3]->y + element->nodes[2]->y));
        nodes[id[3]].link(&nodes[id[0]], South);
        nodes[id[3]].link(element->nodes[3], West);
        nodes[id[3]].link(element->nodes[2], East);
    }
    else
    {
        Node2D::referenceNode(element->nodes[3], East, element->nodes[2], &refnode);
        id[3] = refnode;
        nodes[id[3]].link(&nodes[id[0]], South);
    }


    // West
    if(element->nodes[0]->nodes[2] == element->nodes[3])
    {
        id[4] = nnodes++;
        nodes[id[4]] = Node2D(nnodes-1, 0.5q*(element->nodes[0]->x + element->nodes[3]->x), 0.5q*(element->nodes[0]->y + element->nodes[3]->y));
        nodes[id[4]].link(&nodes[id[0]], East);
        nodes[id[4]].link(element->nodes[0], South);
        nodes[id[4]].link(element->nodes[3], North);
    }
    else
    {
        Node2D::referenceNode(element->nodes[0], North, element->nodes[3], &refnode);
        id[4] = refnode;
        nodes[id[4]].link(&nodes[id[0]], East);
    }

    elements[nelements] = Element2D(nelements++, element->nodes[0], &nodes[id[1]], &nodes[id[0]], &nodes[id[4]]);
    elements[nelements] = Element2D(nelements++, &nodes[id[1]], element->nodes[1], &nodes[id[2]], &nodes[id[0]]);
    elements[nelements] = Element2D(nelements++, &nodes[id[0]], &nodes[id[2]], element->nodes[2], &nodes[id[3]]);
    elements[nelements] = Element2D(nelements++, &nodes[id[4]], &nodes[id[0]], &nodes[id[3]], element->nodes[3]);


    element->index = -1;
}

void Diffusion2DpAR::refine(Node2D *node)
{
    tInteger ne = nelements;
    for(int i=0; i<ne; i++)
        if(elements[i].nodes[0] == node ||
                elements[i].nodes[1] == node ||
                elements[i].nodes[2] == node ||
                elements[i].nodes[3] == node)
            if(elements[i].index!= -1)
                refine(&elements[i]);
}

