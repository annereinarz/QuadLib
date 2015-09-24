#ifndef SQUAD_H
#define SQUAD_H
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
    int s1;            /*size of vertexlist in first dimension*/
    int s2;            /*size of vertexlist in second dimension*/
    double** vtxlist;     /*vertex list, has the size s1 x s2*/
} Vertexlist;

Vertexlist get_vertexlist(int sim);
void free_vertexlist(Vertexlist vtxlist);
double determinant(double**A, int m);
void Quadrature(int d, int k,int whichF, double* ptrQ, QuadraturePoints q);
AffineTrafo determineAffineTrafo(int d,int k, Vertexlist vtxlist);
QuadraturePoints QuadratureRule(int k, int d, int nr, int ns);
int NextSubset(int k, int numN, int* U);
QuadraturePoints Quad2PhyP(int k, int d, QuadraturePoints *QP, AffineTrafo A );
double squad(int nr, int ns, int whichvtxlist, int whichF);
/*QuadraturePoints squad(int nr, int ns, int whichvtxlist, int whichF);*/

QuadraturePoints Step5(int j, int k, int d, QuadraturePoints* QP);
void Step4(int numN, int k, QuadraturePoints* QP);
void Step3(int numN, int k, int d, QuadraturePoints* QP);
void Step2(int k, int d, QuadraturePoints* QP);

#endif /*SQUAD_H*/
