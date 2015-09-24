#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "QuadRule.h"

#define GL  (2)
#define CC  (3)
#define KP  (1)
#define CGL (3)
#define GJ  (1)

QuadRule QuadRule_sobol(int d, int l);
QuadRule QuadRule_reg_direct(int d, int l, QuadRule qp);
QuadRule QuadRule_reg(int d, int l, int regular);

void KPquad(int level, QuadRule*qp);
void GLquad(QuadRule* qp, int n, double a, double b);
void CGLquad(QuadRule* qp, int n);
void gaujac(double* x, double* w, int n, double alf, double bet);
double gammln(double xx);
double *vector(long nl, long nh);
void GJquad01(QuadRule* qp,int n, double alpha, double beta);

#endif /*QUADRATURE_H*/
