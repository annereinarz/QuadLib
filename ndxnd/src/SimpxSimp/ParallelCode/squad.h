#ifndef SQUAD_H
#define SQUAD_H
#include "AffineTrafo.h"
#include "quadpoints.h"

double squad_nonperm(int whichF, AffineTrafo A,int d, int k, int nr, int ns);
double squad_perm(int whichF, AffineTrafo A,int d, int k, int nr, int ns, int j, int perm, int my_rank);
QuadraturePoints QuadratureRule(int k, int d, int nr, int ns);

#endif /*SQUAD_H*/
