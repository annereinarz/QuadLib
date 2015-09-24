#ifndef QUAD2_H
#define QUAD2_H

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

void init_Quad2(Quad2* qp, int luv, int ruv, int lzh, int rzh);
void export_quad2(char* fname, Quad2* qp); 
void free_Quad2(Quad2 qp);

#endif /*QUAD2_H*/
