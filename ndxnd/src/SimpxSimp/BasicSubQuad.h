#ifndef BASICSUBQUAD_H
#define BASICSUBQUAD_H

void BasicSubQuad_jacobi( int j, int k, int d, QuadRule* quadpoints);
void OriginQuad(int* active, int size, QuadRule* qp);
void C2S(QuadRule* qp,int startidx, int endidx);
void BasicSubQuad( int j, int k, int d, QuadRule* quadpoints);

#endif /*BASICSUBQUAD_H*/
