#ifndef QUADRULE_H
#define QUADRULE_H

typedef struct {
    int n;            /*Number of Quadrature Points*/
    int d;            /*Spatial Dimension*/
    double** t;       /*Quadrature Points, has the size nxd*/
    double* wt;       /*Quadrature Weights, has length n*/
} QuadRule;


void write_quadrule(FILE* f, QuadRule *qp);
QuadRule resize_quadrule_n(QuadRule quadpoints, int size);
QuadRule resize_quadrule_d(QuadRule quadpoints, int size);
void init_quadrule(QuadRule* qp, int n, int d);
void import_quadrule(char* fname, QuadRule* qp);
void free_quadrule(QuadRule qp);
QuadRule resize_quadrule(QuadRule quadpoints, int size1, int size2);
void import_1d_quadrule(char* fname, QuadRule* qp);
void write_1d_quadrule(FILE* f, QuadRule *qp);
void print_quadrule(QuadRule qr);

#endif /*QuadRule_H*/
