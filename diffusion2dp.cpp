#include "diffusion2dp.h"
#include "gaussseidel.h"

Diffusion2Dp::Diffusion2Dp(tFloat _lx, tFloat _ly, tInteger _nx, tInteger _ny, Diffusion2DData *_data,
                           Boundary2D *_bS, Boundary2D *_bN, Boundary2D *_bE, Boundary2D *_bW)
    :lx(_lx), ly(_ly), nx(_nx), ny(_ny), data(_data), ccS(_bS), ccN(_bN), ccE(_bE), ccW(_bW)
{
    hx = lx/(nx - 1.0q);
    hy = ly/(ny - 1.0q);

    nodes = new Node2D[nx*ny]();
    tInteger p = 0;

//    for(tInteger j=0; j<ny; j++)
//        for(tInteger i=0; i<nx; i++)
//            nodes[p] = Node2D(i*hx, j*hy, i, j, p++);

    T = new tFloat[nx*ny];

    for(tInteger i=0; i<nx*ny; i++)
        T[i] = 10.q;

}


void Diffusion2Dp::solver(tInteger iterationMax, tFloat iterationTolerance, bool plotlog)
{


    GaussSeidel sys(T, nx*ny, iterationMax, iterationTolerance);
    tFloat a[4];

    for(tInteger j=1; j<ny-1; j++)
        for(tInteger i=1; i<nx-1; i++)
        {
            tInteger p = position(i,j);
            a[0] = a[1] = 1.0q/(hy*hy);
            a[2] = a[3] = 1.0q/(hx*hx);

            sys(p, p, -(a[0]+a[1]+a[2]+a[3]));

            sys(p, direction(p, South), a[0]);
            sys(p, direction(p, North), a[1]);
            sys(p, direction(p, East), a[2]);
            sys(p, direction(p, West), a[3]);

            sys(p, data->heatSource->operator ()(nodes[p].x, nodes[p].y));
        }

    // Condições de Contorno
    // Sul
    for(tInteger i=0; i<nx; i++)
        if(ccS->type == Dirichlet)
            sys(i, ccS->bcValue->operator ()(nodes[i].x, nodes[i].y));
        else{
            sys(i, i, -3.0q);
            sys(i, direction(i, North), 4.0q);
            sys(i, direction(direction(i, North), North), -1.0q);
            sys(i, -2.0q*hy*ccS->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
        }

    // Norte
    for(tInteger i=(ny-1)*nx; i<nx*ny; i++)
        if(ccN->type == Dirichlet)
            sys(i, ccN->bcValue->operator ()(nodes[i].x, nodes[i].y));
        else{
            sys(i, i, 3.0q);
            sys(i, direction(i, South), -4.0q);
            sys(i, direction(direction(i, South), South), 1.0q);
            sys(i, -2.0q*hy*ccN->bcValue->operator ()(nodes[i].x, nodes[i].y)/data->k);
        }
    // East
    for(tInteger i=1; i<ny; i++)
        if(ccE->type == Dirichlet)
            sys(i*nx-1, ccE->bcValue->operator ()(nodes[i*nx-1].x, nodes[i*nx-1].y));
        else{
            int p = i*nx-1;
            sys(p, p, 3.0q);
            sys(p, direction(p, West), -4.0q);
            sys(p, direction(direction(p, West), West), 1.0q);
            sys(p, -2.0q*hx*ccE->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
        }

    // West
    for(tInteger i=0; i<ny; i++)
        if(ccW->type == Dirichlet)
            sys(i*nx, ccW->bcValue->operator ()(nodes[i*nx].x, nodes[i*nx].y));
        else{
            int p = i*nx;
            sys(p, p, -3.0q);
            sys(p, direction(p, East), 4.0q);
            sys(p, direction(direction(p, East), East), -1.0q);
            sys(p, -2.0q*hx*ccW->bcValue->operator ()(nodes[p].x, nodes[p].y)/data->k);
        }

    sys.solver();

    Tm = 0.0q;
    int p;
    for(tInteger j=0; j<ny-1; j++)
        for(tInteger i=0; i<nx-1; i++){
            p = position(i,j);
            Tm += T[p] + T[direction(p, East)] + T[direction(p, North)] + T[direction(direction(p, North), East)];
        }
    Tm *= hx*hy/(4.0q*lx*ly);

    if(plotlog)
        sys.plotIterationLog();
}


tInteger Diffusion2Dp::position(tInteger i, tInteger j)
{
    return i+j*nx;
}


tInteger Diffusion2Dp::direction(tInteger position, DirectionType dir)
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


tFloat Diffusion2Dp::operator()(tInteger i, tInteger j)
{
    return T[position(i, j)];
}


void Diffusion2Dp::plotX(tInteger p, Functor2D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "plotx.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nx; i++)
        file1<<QtoD(nodes[p*nx+i].x)<<"\t"<<QtoD(T[p*nx+i])<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat h_10 = lx/(10*nx-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*nx; i++){
        tFloat xp = i*h_10;
        file2<<QtoD(xp)<<"\t"<<QtoD(analyticalSolution(xp,p*hy))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \"DIFUSÃO DE CALOR BIDIMENSIONAL EM REGIME PERMANENTE\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'x'\n"
             "set ylabel 'T( x, "<<QtoD(p*hy)<<")'\n"
             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}


void Diffusion2Dp::plotY(tInteger p, Functor2D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "ploty.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<ny; i++)
        file1<<QtoD(nodes[p+nx*i].y)<<"\t"<<QtoD(T[p+nx*i])<<std::endl;
    file1.close();

    // Solução Analítica
    std::ofstream file2(dat2_filename.c_str());
    tFloat h_10 = ly/(10*ny-1); // Aumenta o número de pontos em 10X
    for(tInteger i=0; i<10*ny; i++){
        tFloat yp = i*h_10;
        file2<<QtoD(yp)<<"\t"<<QtoD(analyticalSolution(p*hx, yp))<<std::endl;
    }
    file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \"DIFUSÃO DE CALOR BIDIMENSIONAL EM REGIME PERMANENTE\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set xlabel 'y'\n"
             "set ylabel 'T( "<<QtoD(p*hx)<<", y)'\n"
             "plot '" <<dat2_filename<<"' t\"Solução Analítica\" with lines lt 2 lc 2 lw 2, "
             "'" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

void Diffusion2Dp::plot(Functor2D &analyticalSolution)
{
    const std::string cmd_filename = "plotconfig.gnu";
    const std::string pic_filename = "plot.png";
    const std::string dat1_filename = "data1.txt";
    const std::string dat2_filename = "data2.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nx*ny; i++)
        file1<<QtoD(nodes[i].x)<<"\t"<<QtoD(nodes[i].y)<<"\t"<<QtoD(T[i])<<std::endl;
    file1.close();

        // Solução Analítica
        std::ofstream file2(dat2_filename.c_str());
        int f = 10;
        tFloat hy_10 = ly/(f*ny-1); // Aumenta o número de pontos em f X
        tFloat hx_10 = lx/(f*nx-1);
        for(tInteger j=0; j<f*nx; j++)
        for(tInteger i=0; i<f*ny; i++)
        {
            tFloat yp = i*hy_10;
            tFloat xp = j*hx_10;
            file2<<QtoD(xp)<<"\t"<<QtoD(yp)<<"\t"<<QtoD(analyticalSolution(xp, yp))<<std::endl;
        }
        file2.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "#set terminal wxt font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set title \"DIFUSÃO DE CALOR BIDIMENSIONAL EM REGIME PERMANENTE\\nResolução com MDF / Aproximação com CDS-2\"\n"
             "set lmargin 8\n"
             "set xlabel 'x'\n"
             "set ylabel 'y'\n"
             "set zlabel 'T( x, y)'\n"
             "splot '" <<dat2_filename<<"' t\"Solução Analítica\" with points lt 2 lc 0 lw 2 pt 0, "
             " '" <<dat1_filename<<"' t\"Solução Numérica\" with points lt 2 lc 1 pt 13 lw 5\n"
             "#pause -1";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

