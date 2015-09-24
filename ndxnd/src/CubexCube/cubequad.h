#ifndef CUBEQUAD_H
#define CUBEQUAD_H
#include "QuadRule.h"
#include "AffineTrafo.h"

void Quadrature(int d, int k,int whichF, double* ptrQ, QuadRule q, int whichoutput);
int NextSubset(int k, int numN, int* U);

double cubequad(int k, int d, AffineTrafo A, int whichF, int whichoutput,int *ndof, QuadRule *QP);
void cubetransform(int k, int d, FILE *f, QuadRule *QP);
void cubeaffine(int k, int d, AffineTrafo A, int whichoutput, QuadRule *QP, FILE *f);

QuadRule Quad2PhyP(int k, int d, QuadRule *QP, AffineTrafo A , int whichoutput);
QuadRule BasicSubQuad(int j, int k, int d, int numN, QuadRule* QP);

void Step5(int j, int k, int d, QuadRule* QP, QuadRule* qp);
void Step4(int numN, int k, QuadRule* QP);
void Step3(int numN, int k, int d, QuadRule* QP);
void Step2(int k, int d, QuadRule* QP);

#endif /*CUBEQUAD_H*/
