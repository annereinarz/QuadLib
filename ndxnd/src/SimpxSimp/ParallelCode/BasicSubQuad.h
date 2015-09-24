#ifndef BASICSUBQUAD_H
#define BASICSUBQUAD_H

void OriginQuad(int* active,int size, QuadraturePoints* quadpoints);
void C2S(QuadraturePoints* qp,int startidx, int endidx);
void BasicSubQuad( int j, int k, int d, QuadraturePoints* quadpoints);

#endif /*BASICSUBQUAD_H*/
