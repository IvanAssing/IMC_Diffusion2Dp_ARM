#ifndef GAUSSSEIDEL_H
#define GAUSSSEIDEL_H

#include "imc_dfm.h"

class GaussSeidel
{
    public:

        tInteger neq,; // Contador de equações

        tInteger **AIndex; // Matriz de indices
        tInteger *neIndex; // Número de termos por linha


        tFloat **A; // Coeficientes
        tFloat *x; // Resultados
        tFloat *ax; // Coeficientes da diagonal
        tFloat *b; // Vetor constante
        tFloat *L; // Resíduos
        tFloat itol; // Tolerância mímima sobre resíduo
        tInteger nit;
        tInteger vswap[100][2];
        tInteger nswap;

        tInteger bwidth; // Largura da matriz de banda
        tInteger neqmax; // Número máximo de equações
        tInteger imax; // Número máximo de iterações


        GaussSeidel(tFloat *results, tInteger equationMax, tInteger iterationMax = 100, tFloat iterationTolerance = 1.0e-28q, tInteger bandWidth = 5);

        void operator()(tInteger equation, tInteger index, tFloat value);
        void operator()(tInteger equation, tFloat b);

        void diagonalcheck(void);

        void solver();

        void printindex();

        void plotIterationLog();
};

#endif // GAUSSSEIDEL_H
