#ifndef TENSORQUAD_H
#define TENSORQUAD_H

void GLquad(QuadraturePoints* qp, int n, double a, double b);
void CGLquad(QuadraturePoints* qp, int n);
void TensorQuad(int it, int ly, QuadraturePoints*QP, QuadraturePoints qp);
void GJquad(QuadraturePoints* qp,int n, double alpha, double beta);
void GJquad01(QuadraturePoints* qp,int n, double alpha, double beta);
void TensorQuad_nonit(QuadraturePoints* SP, QuadraturePoints QP, QuadraturePoints qp);
void TensorQuad_nonit_ip(QuadraturePoints* QP, QuadraturePoints qp);

#endif



