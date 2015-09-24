#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "quadpoints.h"
#include "Ref2PhyS1d_circle.h"

void Ref2PhyS1d_circle(QuadraturePoints *SP, Vertexlist vtx, int k){

     int i;
     double * th;
     double thelp;
 
     resizequadpointsd(SP, 2);   /* Add second coordinate*/   
     th=SP->t[1];
     SP->t[1]=SP->t[2];
     SP->t[2]=th;
     
     switch(k){ /*map to planar element*/
         case 1:/*identical elements*/             
               for(i=0; i<SP->n; i++){
                  SP->t[1][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][1]-vtx.vtxlist[1][0])*SP->t[0][i];
                  SP->t[0][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][1]-vtx.vtxlist[0][0])*SP->t[0][i];
                  SP->t[3][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][1]-vtx.vtxlist[1][0])*SP->t[2][i];
                  SP->t[2][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][1]-vtx.vtxlist[0][0])*SP->t[2][i];
                  
                  SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][0]-vtx.vtxlist[0][0],2)+pow(vtx.vtxlist[1][1]-vtx.vtxlist[0][1],2)))*SP->wt[i];
               }
               break;
               
         case 0:/*common vertex*/             
               for(i=0; i<SP->n; i++){
                  SP->t[1][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][1]-vtx.vtxlist[1][0])*SP->t[0][i];
                  SP->t[0][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][1]-vtx.vtxlist[0][0])*SP->t[0][i];
                  
                  SP->t[3][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][2]-vtx.vtxlist[1][0])*SP->t[2][i];
                  SP->t[2][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][2]-vtx.vtxlist[0][0])*SP->t[2][i];                  
               }
               for(i=0; i<SP->n/2; i++){
                  SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][1]-vtx.vtxlist[1][0],2)+
                                       pow(vtx.vtxlist[0][1]-vtx.vtxlist[0][0],2)))*SP->wt[i];
                  SP->wt[i+SP->n/2]=(sqrt(pow(vtx.vtxlist[1][2]-vtx.vtxlist[1][0],2)+
                                       pow(vtx.vtxlist[0][2]-vtx.vtxlist[0][0],2)))*SP->wt[i];
               }
               break;
               
        case -1:/*disant elements*/             
               for(i=0; i<SP->n; i++){
                  SP->t[1][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][1]-vtx.vtxlist[1][0])*SP->t[0][i];
                  SP->t[0][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][1]-vtx.vtxlist[0][0])*SP->t[0][i];
                  
                  SP->t[3][i]=vtx.vtxlist[1][2]+(vtx.vtxlist[1][3]-vtx.vtxlist[1][2])*SP->t[2][i];
                  SP->t[2][i]=vtx.vtxlist[0][2]+(vtx.vtxlist[0][3]-vtx.vtxlist[0][2])*SP->t[2][i];  
               }
               for(i=0; i<SP->n; i++){
                  SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][1]-vtx.vtxlist[1][0],2)+
                                       pow(vtx.vtxlist[0][1]-vtx.vtxlist[0][0],2)))*SP->wt[i];
                  SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][3]-vtx.vtxlist[1][2],2)+
                                       pow(vtx.vtxlist[0][3]-vtx.vtxlist[0][2],2)))*SP->wt[i];
               }
               break;
     }
     /* raise to the circle */
     for(i=0; i<SP->n; i++){
         thelp=SP->t[0][i];
         SP->wt[i]=SP->wt[i]*sqrt(pow(SP->t[1][i],2)+pow(SP->t[0][i],2));
         SP->t[0][i]=SP->t[0][i]/(sqrt(pow(SP->t[0][i],2)+pow(SP->t[1][i],2)));
         SP->t[1][i]=SP->t[1][i]/(sqrt(pow(thelp,2)+pow(SP->t[1][i],2)));

         thelp=SP->t[2][i];
         SP->wt[i]=SP->wt[i]*sqrt(pow(SP->t[2][i],2)+pow(SP->t[3][i],2));
         SP->t[2][i]=SP->t[2][i]/(sqrt(pow(SP->t[2][i],2)+pow(SP->t[3][i],2)));
         SP->t[3][i]=SP->t[3][i]/(sqrt(pow(thelp,2)+pow(SP->t[3][i],2))); 
     }
}
