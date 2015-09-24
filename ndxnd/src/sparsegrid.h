#ifndef SPARSEGRID_H
#define SPARSEGRID_H

#include "QuadRule.h"

typedef int bool;
#define true  (1)
#define false (0)

bool elemIl(int d, int* K, int l, int p, int T);
QuadRule Sparse(int d, int l, int p, int whichquadrule, int T);
void set_up_1dquadrules(int l, QuadRule *QP1D, int regular);
int next_sparse(QuadRule *QPtemp,int l,QuadRule* QP1D, int* K, int d, int p, int T);

#endif
