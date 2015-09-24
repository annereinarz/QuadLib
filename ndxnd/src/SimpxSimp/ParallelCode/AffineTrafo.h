#ifndef AFFTR_H
#define AFFTR_H

#include "quadpoints.h"

typedef struct {
  double** A1;
  double** A2;
  double* v10;
  double* v20;
  int lvtx;
  int d;
} AffineTrafo;

typedef struct {
    int s1;            //size of vertexlist in first dimension
    int s2;            //size of vertexlist in second dimension
    double** vtxlist;     //vertex list, has the size s1 x s2
} Vertexlist;  

void PermRefl(int k, int d, QuadraturePoints* QP, AffineTrafo A, int* P, double* Q, int flag, int whichF);
double determinant(double**A, int m);
void Quadrature(int d, int k,int whichF, double* ptrQ, Quad2 q);
AffineTrafo determineAffineTrafo(int d,int k, Vertexlist vtxlist);
void free_afftrafo(AffineTrafo A);
int NextPerm(int j, int k, int d, int* U);
void Quad2RefS(int k, int d,  Quad2* qp);
Quad2 Quad2PhyS(int k, int d, Quad2 qp, AffineTrafo A );

#endif /*AFFTR_H*/
