#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "Ref2PhyJ2d.h"


/*
 * Transfers quadrature rules calculated by the sing_*2dtria rules onto physical Cubes
 * INPUT:
 * n_gl  : the number of gauss legendre points
 * n_cgl : number for the composite gauss legendre quadrature
 * number of quadrature points is 0.5*n_cgl*(n_cgl+1)
 * k      : dimension of the intersection
 * vtx    : list of vertices of the two quadrilateral.  All vertices should be given in a clockwise order
 * starting with the lower left
 * vertex if there are no common vertices. If there are common vertices they are given first in clockwise
 * order, then the remaining vertices of the first cube are given in clockwise order starting
 * next to a common vertex and then the vertices of the second cube are given, also in
 * clockwise order starting next to a common vertex.
 * SP     : Quadrature points on the reference cube
 *
 * OUTPUT:
 * [t,wt]  : quadrature rule on the two input cubes, with singularities as described in the sing_*2dquad routines.
 *
 */

void Ref2PhyJ2d(QuadraturePoints* SP, Vertexlist vtx, int k){
    
    int i;
    double** A1, ** A2;
    double *v1, *v2;
    
    v1=(double*)malloc(2*sizeof(double));
    v2=(double*)malloc(2*sizeof(double));
    
    A1=(double**)malloc(2*sizeof(double*));
    for(i=0; i<2; i++){
        A1[i]=(double*)malloc(2*sizeof(double));
    }
    A2=(double**)malloc(2*sizeof(double*));
    for(i=0; i<2; i++){
        A2[i]=(double*)malloc(2*sizeof(double));
    }
    
    A1[0][0]=vtx.vtxlist[0][1]-vtx.vtxlist[0][0];  A1[0][1]=vtx.vtxlist[0][2]-vtx.vtxlist[0][0];
    A1[1][0]=vtx.vtxlist[1][1]-vtx.vtxlist[1][0];  A1[1][1]=vtx.vtxlist[1][2]-vtx.vtxlist[1][0];

    v1[0]=vtx.vtxlist[0][0];  v1[1]=vtx.vtxlist[1][0];   
            
    switch(k){
        case -1:
            A2[0][0]=vtx.vtxlist[0][4]-vtx.vtxlist[0][3];  A2[0][1]=vtx.vtxlist[0][5]-vtx.vtxlist[0][3];
            A2[1][0]=vtx.vtxlist[1][4]-vtx.vtxlist[1][3];  A2[1][1]=vtx.vtxlist[1][5]-vtx.vtxlist[1][3];
            
            v2[0]=vtx.vtxlist[0][3];  v2[1]=vtx.vtxlist[1][3];
            break;
            
        case 0:            
            A2[0][0]=vtx.vtxlist[0][3]-vtx.vtxlist[0][0];  A2[0][1]=vtx.vtxlist[0][4]-vtx.vtxlist[0][0];
            A2[1][0]=vtx.vtxlist[1][3]-vtx.vtxlist[1][0];  A2[1][1]=vtx.vtxlist[1][4]-vtx.vtxlist[1][0];
            
            v2[0]=vtx.vtxlist[0][0];  v2[1]=vtx.vtxlist[1][0];
            break;
            
        case 1:
            A2[0][0]=vtx.vtxlist[0][1]-vtx.vtxlist[0][0];  A2[0][1]=vtx.vtxlist[0][3]-vtx.vtxlist[0][0];
            A2[1][0]=vtx.vtxlist[1][1]-vtx.vtxlist[1][0];  A2[1][1]=vtx.vtxlist[1][3]-vtx.vtxlist[1][0];
            
            v2[0]=vtx.vtxlist[0][0];  v2[1]=vtx.vtxlist[1][0];
            break;
            
        case 2:
            A2[0][0]=vtx.vtxlist[0][1]-vtx.vtxlist[0][0];  A2[0][1]=vtx.vtxlist[0][2]-vtx.vtxlist[0][0];
            A2[1][0]=vtx.vtxlist[1][1]-vtx.vtxlist[1][0];  A2[1][1]=vtx.vtxlist[1][2]-vtx.vtxlist[1][0];
            
            v2[0]=vtx.vtxlist[0][0];  v2[1]=vtx.vtxlist[1][0];
            break;
    }
    
    double th;    
    /*z=*/
    switch(k){ 
        case 2:
            for(i=0; i<SP->n; i++){
                th=SP->t[4][i];
                SP->t[4][i]=A1[0][0]*th+A1[0][1]*SP->t[5][i];
                SP->t[5][i]=A1[1][0]*th+A1[1][1]*SP->t[5][i];
            }
            break;
        case 1:
            for(i=0; i<SP->n; i++){
                SP->t[4][i] = SP->t[4][i]*A1[0][0] + SP->t[3][i]*A2[0][1] - SP->t[1][i]*A1[0][1] ;
            }            
            break;
        /*zhat variable does not exist in the cases k=-1,0*/
    }
    
    /*x=A1u+v1, u=SP->t[0:1]*/
    for(i=0; i<SP->n; i++){
        th=SP->t[0][i];
        SP->t[0][i]=A1[0][0]*th+A1[0][1]*SP->t[1][i]+v1[0];
        SP->t[1][i]=A1[1][0]*th+A1[1][1]*SP->t[1][i]+v1[1];
    }
    
    /*y=A2v+v2, v=SP->t[2:3]*/
    for(i=0; i<SP->n; i++){
        th=SP->t[2][i];
        SP->t[2][i]=A2[0][0]*th+A2[0][1]*SP->t[3][i]+v2[0];
        SP->t[3][i]=A2[1][0]*th+A2[1][1]*SP->t[3][i]+v2[1];       
    }  
   
    /*wt=|detA1||detA2|wt*/
    for(i=0; i<SP->n; i++){
        SP->wt[i]*=fabs(A1[1][1]*A1[0][0]-A1[1][0]*A1[0][1])*fabs(A2[1][1]*A2[0][0]-A2[1][0]*A2[0][1]);
    }
    
    free(v1);
    free(v2);
    free(A1[0]); free(A1[1]); free(A2[0]); free(A2[1]);
    free(A1); free(A2);
    
    return;


}


