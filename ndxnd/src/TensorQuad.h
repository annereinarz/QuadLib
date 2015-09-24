#ifndef TENSORQUAD_H
#define TENSORQUAD_H
#include "QuadRule.h"

void TensorQuad(int it, int lx, QuadRule* QP, QuadRule qp);
void TensorQuad_sparse(QuadRule *QP, QuadRule qp);
void TensorQuad_new(QuadRule *SP,QuadRule QP, QuadRule qp);
void TensorQuad_nonit(QuadRule* SP, QuadRule QP, QuadRule qp);
void TensorQuad_nonit_ip(QuadRule* QP, QuadRule qp);

#endif /*TENSORQUAD_H*/
