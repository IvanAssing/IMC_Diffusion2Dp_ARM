#include <iostream>

#include <QApplication>

#include "imc_dfm.h"
#include "diffusion2dpar.h"
#include "gaussseidel.h"

#include "graphics.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    // MANAGE_EXCEPTIONS

    Diffusion2DData data;

    //data.heatSource = new Constant2D(0.0q); // Sem geração de calor
    //data.heatSource = new Sine2(-4.q*M_PIq*M_PIq, 2.q*M_PIq, 0.0q);
    data.heatSource = new SFAS2();
    data.k = 1.0q;



    Constant2D cc0(0.0);

    Sine ccS(1.0q, M_PIq, 0.0q);

    Boundary2D south(Dirichlet, &cc0);
    //Boundary2D north(Dirichlet, &ccS);
    Boundary2D north(Dirichlet, &cc0);
    Boundary2D east(Dirichlet, &cc0);
    Boundary2D west(Dirichlet, &cc0);


    tFloat lx = 1.0q;
    tFloat ly = 1.0q;
    tInteger nx = 11;
    tInteger ny = 11;

    tFloat TmAS = 2.0q*lx * (coshq(M_PIq*ly/lx) - 1.0q) / (M_PIq*M_PIq*ly*sinhq(M_PIq*ly/lx));
    //tFloat TmAS = 0.0q;
    Diffusion2DpAR mesh(lx, ly, nx, ny, &data, &south, &north, &east, &west);


    //SFAS as(1.0q, 1.0q);
    //Sine2 as(1.0q, 2.q*M_PIq, 0.0q);
    SFAS2a as;

    // 1
    mesh.solver3(100000, 1.0e-28q, true); // itmax, itol, plotlog?
    mesh.updateTm();
    std::cout<<"\n *** TM: "<<mesh.nnodes<<"\t"<<print(mesh.Tm)<<"\t"<<print(TmAS)<<"\t"<<print(TmAS-mesh.Tm);


    tInteger nne = mesh.nelements;


    for(int z=1; z<4; z++)
    {
        nne = mesh.nelements;
        for(int k=0; k<nne; k++)
        {
            if(mesh.elements[k]->index == -1) continue;

            tFloat Ta = 0.0q, Tn = 0.0q;

            for(tInteger i=0; i<4; i++)
            {
                Tn += mesh.T[mesh.elements[k]->nodes[i]->index]/4.0q;
                Ta += as(mesh.elements[k]->nodes[i]->x, mesh.elements[k]->nodes[i]->y)/4.0q;
            }

            //if(fabsq(Ta - Tn)>(1.0e-4q/((z+1)*1.0e-2q)))
                if(fabsq(Ta - Tn)>4.0e-3q)
            //if(fabsq(Ta - Tn)>fabsq(1.0e-2q*(z*1.0e-1q)))
                mesh.refine(mesh.elements[k]);
        }


        std::cout<<"\n *** TM: "<<mesh.nnodes<<std::flush;
        mesh.solver3(100000, 1.0e-28q, true); // itmax, itol, plotlog?
        mesh.updateTm();
        std::cout<<"\t"<<print(mesh.Tm)<<"\t"<<print(TmAS)<<"\t"<<print(TmAS-mesh.Tm);
    }



    tFloat *erro = new tFloat[mesh.nnodes];

    // Erro númerico
    for(int i=0; i<mesh.nnodes; i++)
    {
        erro[i] = fabsq(as(mesh.nodes[i]->x, mesh.nodes[i]->y)-mesh.T[i]);
        //erro[i] = as(mesh.nodes[i]->x, mesh.nodes[i]->y);
        //std::cout<<"\n"<<i<<"\t"<<print(as(mesh.nodes[i]->x, mesh.nodes[i]->y))<<"\t"<<print(mesh.T[i])<<"\t"<<print(erro[i]);
    }


    Graphics w11(&mesh, mesh.T, QString("Temperatura"));
    w11.show();

    Graphics w2(&mesh, erro, QString("Erro"));
    w2.show();

    int nn = 10;

    // Resultados númericos
    std::cout<<"\n\nT("<<QtoD(mesh.nodes[nn]->x)<<", "<<QtoD(mesh.nodes[nn]->y)<<"):";
    std::cout<<std::setw(15)<<std::right<<"\nNumérica: "<<print(mesh.T[nn]);
    std::cout<<std::setw(15)<<std::right<<"\nAnalítica: "<<print(as(mesh.nodes[nn]->x, mesh.nodes[nn]->y));
    std::cout<<std::setw(15)<<std::right<<"\nErro: "<<print(as(mesh.nodes[nn]->x, mesh.nodes[nn]->y) - mesh.T[nn])<<std::endl;

    //END_EXCEPTIONS

    return a.exec();
}
