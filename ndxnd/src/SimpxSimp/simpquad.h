#ifndef SIMPQUAD_H
#define SIMPQUAD_H

#include "QuadRule.h"
#include "AffineTrafo.h"

void PermRefl(int k, int d, QuadRule* QP, AffineTrafo A, int* P, double* Q, int flag, int whichF);
void Quadrature(int d, int k,int whichF, double* ptrQ, Quad2 q);
int NextPerm(int j, int k, int d, int* U);
void Quad2RefS(int k, int d,  Quad2* qp);
Quad2 Quad2PhyS(int k, int d, Quad2 qp, AffineTrafo A );
double squadtrafoj_jacobi(int whichF, AffineTrafo A, QuadRule *QP, int d, int k, int j, int* ndof);
double squadtrafoj(int whichF, AffineTrafo A, QuadRule* QP, int d, int k, int j, int* ndof);

#endif /*SIMPQUAD_H*/
