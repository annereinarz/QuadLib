#ifndef QUADPOINTS_H
#define QUADPOINTS_H

typedef struct {
    int n;            /*Number of Quadrature Points*/
    int d;            /*Spatial Dimension*/
    double** t;       /*Quadrature Points, has the size nxd*/
    double* wt;       /*Quadrature Weights, has length n*/
} QuadraturePoints;

QuadraturePoints resizequadpointsn(QuadraturePoints quadpoints, int size);
QuadraturePoints resizequadpointsd(QuadraturePoints quadpoints, int size);
void init_quadpoints(QuadraturePoints* qp, int n, int d); 
void import_quadppoints(char* fname, QuadraturePoints* qp); 
void free_qp(QuadraturePoints qp);

#endif /*QUADPOINTS_H*/
