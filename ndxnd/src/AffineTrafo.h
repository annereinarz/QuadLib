#ifndef AFFINETRAFO_H
#define AFFINETRAFO_H

#include "vertexlist.h"

typedef struct {
  double** A1;
  double** A2;
  double* v10;
  double* v20;
  int lvtx;
  int d;
} AffineTrafo;

double determinant(double**A, int m);

AffineTrafo determineAffineTrafo(int d,int k, Vertexlist vtxlist);
void freeAffineTrafo(AffineTrafo AfTr);

#endif /*AFFINETRAFO_H*/
