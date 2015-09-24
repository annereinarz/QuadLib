#ifndef TENSORQUAD_H
#define TENSORQUAD_H

void GLquad(QuadraturePoints* qp, int n, double a, double b);
void CGLquad(QuadraturePoints* qp, int n);
void TensorQuad(int it, int lx, QuadraturePoints* QP, QuadraturePoints qp);

#endif /*TENSORQUAD_H*/