#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../quadpoints.h"
#include "sing_identical2dtria.h"


void sing_identical2dtria(QuadraturePoints* SP){

     int i;
     double *pth;
     pth=SP->t[0];
     SP->t[0]=SP->t[3];
     SP->t[3]=pth;
     
    for(i=0; i<SP->n; i++){
        SP->wt[i]=SP->wt[i]*(SP->t[3][i]);
        SP->t[2][i]=SP->t[2][i]*SP->t[3][i];
        SP->t[3][i]=SP->t[3][i]-SP->t[2][i];
     }  
          
     for(i=0; i<SP->n; i++){
        SP->wt[i]=SP->wt[i]*(1-SP->t[1][i]);
        SP->t[0][i]=SP->t[0][i]*(1-SP->t[1][i]);
     }
    
     double *sp3w;
     double **sp3t;
     sp3w=(double*)malloc(SP->n*sizeof(double));
     sp3t=(double**)malloc(SP->d*sizeof(double*));
     for(i=0; i<SP->d; i++){
        sp3t[i]=(double*)malloc(SP->n*sizeof(double));
     }
     
     for(i=0; i<SP->n; i++){
        sp3w[i]=SP->wt[i]*pow(1-SP->t[2][i]-SP->t[3][i],2);
        sp3t[0][i]=SP->t[0][i]*(1-SP->t[2][i]-SP->t[3][i]);
        sp3t[1][i]=SP->t[1][i]*(1-SP->t[2][i]-SP->t[3][i]);
        sp3t[2][i]=SP->t[2][i];
        sp3t[3][i]=SP->t[3][i];
     }
     
     for(i=0; i<SP->n; i++){
         SP->t[3][i]= SP->t[3][i]+SP->t[2][i];
         SP->t[2][i]=-SP->t[2][i];
     } 
     
     for(i=0; i<SP->n; i++){
         SP->wt[i]=SP->wt[i]*(1-SP->t[3][i])*(1-SP->t[3][i]);
         SP->t[0][i]=SP->t[0][i]*(1-SP->t[3][i])-SP->t[2][i];
         SP->t[1][i]=SP->t[1][i]*(1-SP->t[3][i]);
     }

     resizequadpointsn(SP, 3);
     
     for(i=0; i<SP->n/3; i++){
         SP->t[0][i+SP->n/3]=SP->t[0][i];
         SP->t[1][i+SP->n/3]=SP->t[1][i];
         SP->t[2][i+SP->n/3]=SP->t[2][i];
         SP->t[3][i+SP->n/3]=SP->t[3][i];
         SP->t[0][i]=sp3t[0][i];
         SP->t[1][i]=sp3t[1][i];
         SP->t[2][i]=sp3t[2][i];
         SP->t[3][i]=sp3t[3][i];
         SP->t[0][i+2*SP->n/3]=SP->t[1][i+SP->n/3];
         SP->t[1][i+2*SP->n/3]=SP->t[0][i+SP->n/3];
         SP->t[2][i+2*SP->n/3]=SP->t[3][i+SP->n/3];
         SP->t[3][i+2*SP->n/3]=SP->t[2][i+SP->n/3];
         SP->wt[i+SP->n/3]=SP->wt[i];
         SP->wt[i+2*SP->n/3]=SP->wt[i];
         SP->wt[i]=sp3w[i];
     }
       
     resizequadpointsn(SP,2);
     
     for(i=0; i<SP->n/2; i++){
         SP->wt[i+SP->n/2]=SP->wt[i];
         SP->t[0][i+SP->n/2]=SP->t[0][i]+SP->t[2][i];
         SP->t[1][i+SP->n/2]=SP->t[1][i]+SP->t[3][i];
         SP->t[2][i+SP->n/2]=-SP->t[2][i];
         SP->t[3][i+SP->n/2]=-SP->t[3][i];
     }
     
    resizequadpointsd(SP, 2);
    for(i=0; i<SP->n; i++){
         SP->t[4][i]=SP->t[2][i];
         SP->t[5][i]=SP->t[3][i];
     }
          
     for(i=0; i<SP->n; i++){
         SP->t[2][i]=SP->t[0][i]+SP->t[2][i];
         SP->t[3][i]=SP->t[1][i]+SP->t[3][i];
     }

}

     