#include "gaussseidel.h"

#include <stdlib.h>


GaussSeidel::GaussSeidel(tFloat *results, tInteger equationMax, tInteger iterationMax, tFloat iterationTolerance, tInteger bandWidth)
    :x(results), neqmax(equationMax), itol(iterationTolerance), imax(iterationMax), bwidth(bandWidth)
{
    neIndex = new tInteger[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        neIndex[i] = 0;

    neq = 0;
    nit = 0;


    AIndex = new tInteger*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        AIndex[i] = new tInteger[bwidth-1];


    ax = new tFloat[neqmax];
    b = new tFloat[neqmax];
    L = new tFloat[imax];

    for(tInteger i = 0; i<neqmax; i++)
        ax[i] = 1.0q, b[i] = 0.0q;


    for(tInteger i = 0; i<imax; i++)
        L[i] = 0.0q;

    A = new tFloat*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        A[i] = new tFloat[bwidth-1];


}

void GaussSeidel::operator()(tInteger i, tInteger j, tFloat value)
{
    if(i<0 || i>neqmax || j<0 || j>neqmax)
        throw std::string("error: invalid index");
    else if(i==j)
        ax[i] = value;
    else{
        A[i][neIndex[i]] = value;
        AIndex[i][neIndex[i]] = j;
        neIndex[i]++;
    }
}

void GaussSeidel::operator()(tInteger i, tFloat _b)
{
    if(i<0 || i>neqmax)
        throw std::string("error: invalid index");
    else
        b[i] = _b;
}

void GaussSeidel::solver()
{
    tFloat sum, residual;

    do{
        // Solver
        for(tInteger i = 0; i<neqmax; i++){
            sum = b[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum -= A[i][j]*x[AIndex[i][j]];
            x[i] = sum/ax[i];
        }

        // Residuo
        residual = 0.0q;
        for(tInteger i = 0; i<neqmax; i++){
            sum = b[i]-ax[i]*x[i];
            for(tInteger j=0; j<neIndex[i]; j++)
                sum -= A[i][j]*x[AIndex[i][j]];
            residual += sum*sum;
        }
        L[nit] = sqrtq(residual);

        std::cout<<"\n"<<nit<<"\t"<<print(L[nit]);

    }while(L[nit] > itol && nit++ < imax);

}

#define DIV0_TOL 1.0e-25



void GaussSeidel::diagonalcheck(void)
{
    nswap = 0;
    for(tInteger i=0; i<neqmax; i++)
        if(fabsq(ax[i]) < DIV0_TOL)
        {
            for(tInteger k=0; k<neIndex[i]; k++)
                for(tInteger q=0; q<neIndex[AIndex[i][k]]; q++)
                    if(i == AIndex[AIndex[i][k]][q])
                    {
                        tInteger t = AIndex[i][k];


                        vswap[nswap][0] = i;
                        vswap[nswap++][1] = t;

                        tInteger ni = neIndex[i];
                        tInteger nt = neIndex[t];

                        tInteger AIi[100], AIt[100];
                        tFloat Ai[100], At[100];

                        for(tInteger kk=0; kk<neIndex[i]; kk++)
                        {
                            AIi[kk] = AIndex[i][kk];
                            Ai[kk] = A[i][kk];
                        }

                        for(tInteger kk=0; kk<neIndex[t]; kk++)
                        {
                            AIt[kk] = AIndex[t][kk];
                            At[kk] = A[t][kk];
                        }

                        neIndex[i] = 0;

                        this->operator ()(i, t, ax[t]);
                        for(tInteger kk=0; kk<nt; kk++)
                            if(AIt[kk] != i)
                                this->operator ()(i, AIt[kk], At[kk]);
                            else
                                this->operator ()(i, i, At[kk]);

                        neIndex[t] = 0;

                        this->operator ()(t, i, ax[i]);
                        for(tInteger kk=0; kk<ni; kk++)
                            if(AIi[kk] != t)
                                this->operator ()(t, AIi[kk], Ai[kk]);
                            else
                                this->operator ()(t, t, Ai[kk]);




                        tFloat swap = b[t];
                        b[t] = b[i];
                        b[i] = swap;


                        swap = x[t];
                        x[t] = x[i];
                        x[i] = swap;






                    }
                    else
                        std::cout<<"ERRO";
        }
}

//tInteger GaussSeidel::findSwap(tInteger pivot, tInteger column)
//{
//    for(tInteger i=0; i<neqmax; i++)
//        for(tInteger k=0; k<neIndex[i]; k++)
//            if(column == AIndex[i][k])
//                for(tInteger kk=0; kk<neIndex[i]; kk++)
//                    if(pivot == AIndex[i][kk])
//                        return i;
//    return -1;
//}

void GaussSeidel::plotIterationLog()
{
    const std::string cmd_filename = "plotconfig_it.gnu";
    const std::string pic_filename = "iteractive_it.png";
    const std::string dat1_filename = "data_it.txt";

    // Solução numérica
    std::ofstream file1(dat1_filename.c_str());
    for(tInteger i=0; i<nit; i++)
        file1<<i<<"\t"<<QtoD(L[i]/L[0])<<std::endl;
    file1.close();

    std::ofstream file3(cmd_filename.c_str());
    file3 <<
             "set terminal pngcairo enhanced font \"arial,12\" size 1600, 1000 \n"
             "set output '" << pic_filename <<"'\n"
             "#set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 1.000\n"
             "set grid\n"
             "set logscale y\n"
             "set format y \"10^{%L}\" \n"
             "set lmargin 10 \n"
             "set title \"DESEMPENHO DE ITERAÇÃO \\n itmax = "<<imax<<"  itol = "<<QtoD(itol)<<"  nit = "<<nit-1<<"\"\n"
             "set ylabel \"L^{n}/L^{0}\" \n"
             "set xlabel 'Número de iterações'\n"

             "plot '" <<dat1_filename<<"' t\"\" with lines lt 2 lc 1 lw 2";
    file3.close();

    const std::string cmd1 = "gnuplot " + cmd_filename; // Gráfico com GNUPLOT
    const std::string cmd2 = "eog " + pic_filename; // Visualizador de imagem

    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
}

void GaussSeidel::printindex()
{
    for(tInteger i = 0; i<neqmax; i++)
    {
        std::cout<<"\n"<<i<<"\t"<<print(ax[i]);
        for(tInteger j=0; j<neIndex[i]; j++)
            std::cout<<"\t"<<AIndex[i][j];
    }

}
