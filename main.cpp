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

    data.heatSource = new Constant2D(0.0q); // Sem geração de calor
    data.k = 1.0q;



    Constant2D cc0(0.0);
    Constant2D cc1(1.0);
    Constant2D cc2(-1.0);

    Sine ccS(1.0q, M_PIq, 0.0q);

    // Condições de contorno do problema 7.1
    Boundary2D south(Dirichlet, &cc0);
    Boundary2D north(Dirichlet, &ccS);
    Boundary2D east(Dirichlet, &cc0);
    Boundary2D west(Dirichlet, &cc0);

    //Condições de contorno de Neumann
    //Boundary2D east(Neumann, &cc0);
    //Boundary2D west(Neumann, &cc0);
    //Boundary2D south(Dirichlet, &cc0);
    //Boundary2D north(Dirichlet, &cc1);

    tFloat lx = 1.0q;
    tFloat ly = 1.0q;
    tInteger nx = 11;
    tInteger ny = 11;

    tFloat TmAS = 2.0q*lx * (coshq(M_PIq*ly/lx) - 1.0q) / (M_PIq*M_PIq*ly*sinhq(M_PIq*ly/lx));
    Diffusion2DpAR mesh(lx, ly, nx, ny, &data, &south, &north, &east, &west);

    mesh.solver(100000, 1.0e-28q, true); // itmax, itol, plotlog?

    SFAS as(1.0q, 1.0q);

//    mesh.plotX(ny/2, as);
//    mesh.plotY(nx/2, as);
//    mesh.plot(as);

    tFloat *erro = new tFloat[mesh.nnodes];

    // Erro númerico
    for(int i=0; i<mesh.nnodes; i++)
    {
        erro[i] = fabsq(as(mesh.nodes[i].x, mesh.nodes[i].y)-mesh.T[i]);
        std::cout<<"\n"<<i<<"\t"<<print(as(mesh.nodes[i].x, mesh.nodes[i].y))<<"\t"<<print(mesh.T[i])<<"\t"<<print(erro[i]);
    }

//    // Gerar mapa de cores
//    Graphics w1, w2;
//    w1.mesh = &mesh;
//    w1.X = mesh.T;
//    w2.mesh = &mesh;
//    w2.X = erro;
//    w1.show();
//    w2.show();

    Graphics w1(&mesh, mesh.T, QString("Temperatura"));
    w1.show();

    Graphics w2(&mesh, erro, QString("Erro"));
    w2.show();

    int nn = 82;

    // Resultados númericos
    std::cout<<"\n\nT("<<QtoD(mesh.nodes[nn].x)<<", "<<QtoD(mesh.nodes[nn].y)<<"):";
    std::cout<<std::setw(15)<<std::right<<"\nNumérica: "<<print(mesh.T[nn]);
    std::cout<<std::setw(15)<<std::right<<"\nAnalítica: "<<print(as(mesh.nodes[nn].x, mesh.nodes[nn].y));
    std::cout<<std::setw(15)<<std::right<<"\nErro: "<<print(as(mesh.nodes[nn].x, mesh.nodes[nn].y) - mesh.T[nn])<<std::endl;

//    std::cout<<"\n\nTemperatura Média: "<<std::setfill(' ');
//    std::cout<<std::setw(15)<<std::right<<"\nNumérica: "<<print(mesh.Tm);
//    std::cout<<std::setw(15)<<std::right<<"\nAnalítica: "<<print(TmAS);
//    std::cout<<std::setw(15)<<std::right<<"\nErro: "<<print(TmAS - mesh.Tm)<<std::endl;

    //END_EXCEPTIONS

    return a.exec();
}
