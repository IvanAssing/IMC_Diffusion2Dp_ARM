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

//        ~SolverLU()
//        {

//            for(tInteger i = 0; i<neqmax; i++)
//                delete [] A[i];

//            delete [] A;
//            delete [] b;

//        }

};

#endif // SOLVERLU_H
