#include "solverlu.h"

#include <stdlib.h>


SolverLU::SolverLU(tFloat *results, tInteger equationMax)
    :x(results), neqmax(equationMax)
{

    neq = 0;

    b = new tFloat[neqmax];

    A = new tFloat*[neqmax];
    for(tInteger i = 0; i<neqmax; i++)
        A[i] = new tFloat[neqmax];

    for(tInteger i = 0; i<neqmax; i++){
        b[i] = 0.0q;
        for(tInteger j = 0; j<neqmax; j++)
            A[i][j] = 0.0q;
        A[i][i] = 1.0q;
    }
}

void SolverLU::operator()(tInteger i, tInteger j, tFloat value)
{
    if(i<0 || i>neqmax || j<0 || j>neqmax)
        throw std::string("error: invalid index");
    else{
        A[i][j] = value;

    }

}

void SolverLU::operator()(tInteger i, tFloat _b)
{
    if(i<0 || i>neqmax)
        throw std::string("error: invalid index");
    else
        b[i] = _b, neq++;
}

void SolverLU::solver()
{
    //void gaussp(int n, double** m, double* v, double** l)
    //Definir vector  permutacion
    tInteger* per;
    per = new tInteger[neqmax];
    for (tInteger i=1;i<neqmax;i++) per[i]=i;

//    for (tInteger i=0;i<neqmax; i++){
//        l[i][i]=1;
//        if (i<n-1) for(tInteger j=i+1; j<n; j++) l[i][j]=0;
//    }


    //ciclo eliminacion
    for (tInteger i=0;i<neqmax;i++){
        //search for the pivot
        tFloat piv = fabsq(A[i][i]);
        tInteger ipiv = i;
        for (tInteger j=i;j<neqmax;j++){
            if( fabsq(A[j][i])>piv){
                ipiv=j;
                piv= fabsq(A[j][i]);
            }
        }
        //Permutar filas i y  ipiv
        for (tInteger j=i;j<neqmax;j++) {
            tFloat temp = A[i][j];
            A[i][j] = A[ipiv][j];
            A[ipiv][j] = temp;
        }
        tFloat temp = b[i];
        b[i]=b[ipiv];
        b[ipiv]=temp;

        //Permutar vector permutacion
        tInteger itemp=per[i];
        per[i]=per[ipiv];
        per[ipiv]=itemp;

        //Eliminacion
        for (tInteger j=i+1; j<neqmax; j++){
            tFloat   piv=A[j][i]/A[i][i]; //m[j][i] se anula despues de la primera elim.
            //l[j][i]=piv;
            for (tInteger k=i; k<neqmax; k++){
                A[j][k] = A[j][k] - piv*A[i][k];
            }
            b[j] = b[j]- piv*b[i];
        }
    }


    //solucion de un sistema triangular superior por sustitucion hacia atras
    for (int i=neqmax-1;i>=0;i--){
       x[i]=b[i];
      for (int j= neqmax-1;j>i;j--)
        x[i]=x[i]-A[i][j]*x[j];
        x[i]=x[i]/A[i][i];
      }

}



void SolverLU::printindex()
{
    //    for(tInteger i = 0; i<neqmax; i++)
    //    {
    //        std::cout<<"\n"<<i<<"\t"<<print(ax[i]);
    //        for(tInteger j=0; j<neIndex[i]; j++)
    //            std::cout<<"\t"<<AIndex[i][j];
    //    }

}
