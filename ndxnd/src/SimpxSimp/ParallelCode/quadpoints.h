#ifndef QUADPOINTS_H
#define QUADPOINTS_H

typedef struct {
    int n;            //Number of Quadrature Points
    int d;            //Spatial Dimension
    double** t;       //Quadrature Points, has the size nxd
    double* wt;       //Quadrature Weights, has length n
} QuadraturePoints;

typedef struct {
    int luv;
    int ruv;       
    int lzh;
    int rzh;          
    double **u;          // vector of size luvxruv
    double **v;           // vector of size luvxruv
    double **zh;        // vector of size lzhxrzh
    double *w;        //Quadrature weights, ruv in length
} Quad2;

QuadraturePoints resizequadpointsn(QuadraturePoints quadpoints, int size);
QuadraturePoints resizequadpointsd(QuadraturePoints quadpoints, int size);
QuadraturePoints resizequadpoints(QuadraturePoints quadpoints, int size1, int size2);
void init_quadpoints(QuadraturePoints* qp, int n, int d); 
void init_Quad2(Quad2* qp, int luv, int ruv, int lzh, int rzh);
void import_quadpoints(char* fname, QuadraturePoints* qp); 
void export_quadpoints(char* fname, QuadraturePoints* qp); 
void export_quad2(char* fname, Quad2* qp); 
void free_qp(QuadraturePoints qp);
void free_Quad2(Quad2 qp);

#endif /*QUADPOINTS_H*/
