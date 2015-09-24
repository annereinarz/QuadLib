#ifndef QUADPOINTS_H
#define QUADPOINTS_H

typedef struct {
    int n;            /*Number of Quadrature Points*/
    int d;            /*Spatial Dimension*/
    double** t;       /*Quadrature Points, has the size nxd*/
    double* wt;       /*Quadrature Weights, has length n*/
} QuadraturePoints;


typedef struct {
  double** A1;
  double** A2;
  double* v10;
  double* v20;
  int lvtx;
  int d;
} AffineTrafo;


typedef struct {
    int s1;            /*size of vertexlist in first dimension*/
    int s2;            /*size of vertexlist in second dimension*/
    double** vtxlist;     /*vertex list, has the size s1 x s2*/
} Vertexlist;  

AffineTrafo determineAffineTrafo(int d, int k, Vertexlist vtxlist);
void resizequadpointsn(QuadraturePoints* quadpoints, int size);
void resizequadpointsd(QuadraturePoints* quadpoints, int size);
void init_quadpoints(QuadraturePoints* qp, int n, int d); 
void free_qp(QuadraturePoints qp);

void GLquad(QuadraturePoints* qp, int n, double a, double b);
void CGLquad(QuadraturePoints* qp, int n);
void TensorQuad(QuadraturePoints* SP, QuadraturePoints QP, QuadraturePoints qp);
void init_vtxlist(Vertexlist *vtx, int s1, int s2);
void free_vtxlist(Vertexlist *vtx);
void QuadRule_cgl_gl(QuadraturePoints *SP, int nr, int ns, int s2);

#endif /*QUADPOINTS_H*/
