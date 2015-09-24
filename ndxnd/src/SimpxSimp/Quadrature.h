#ifndef QUAD_H
#define QUAD_H
#include "QuadRule.h"

QuadRule QuadratureRule(int k, int d, int nr, int ns);
void SJacobi(QuadRule* SP, int n, int d);
QuadRule QuadratureRule_gaujac(int k, int d, int nr, int ns, int j, double alpha);
QuadRule QuadratureRule_gaujacsing(int k, int d,  int ns, int j, double alpha);
QuadRule QuadratureRule_gaujacreg(int k, int d, int nr, int j, int opt);
void QuadratureRule_gaujacreg2(int k, int d, int j, int nr, QuadRule *qp);
QuadRule QuadratureRule_gaulegsing(int ns);
QuadRule QuadratureRule_cgl_gl(int k, int d, int nr, int ns);
QuadRule gauleg_reg(int d, int nr);

#endif /*QUAD_H*/
