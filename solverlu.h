#ifndef SOLVERLU_H
#define SOLVERLU_H

#include "imc_dfm.h"

class SolverLU
{
    public:

        tInteger neq; // Contador de equações
        tInteger neqmax;
        tFloat **A; // Coeficientes
        tFloat *x; // Resultados
        tFloat *b; // Vetor constante


        SolverLU(tFloat *results, tInteger equationMax);

        void operator()(tInteger equation, tInteger index, tFloat value);
        void operator()(tInteger equation, tFloat b);


        void solver();

        void printindex();

};

#endif // SOLVERLU_H
