#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "quadpoints.h"
#include "Ref2PhyS1d_affine.h"

void Ref2PhyS1d_affine(QuadraturePoints *SP, Vertexlist vtx, int k){
    
    int i;
    if(vtx.s1==1){
        AffineTrafo A;
        
        A=determineAffineTrafo(1, k, vtx);
        
        for(i=0; i<SP->n; i++){
            SP->t[0][i]=A.A1[0][0]*SP->t[0][i] +A.v10[0];
            SP->t[1][i]=A.A2[0][0]*SP->t[1][i] +A.v20[0];
            SP->wt[i]=fabs(A.A2[0][0]*A.A1[0][0])*SP->wt[i];
        }
    }
    else if(vtx.s1==2){
        double th1, th2;
        resizequadpointsd(SP, 2);   /* Add second coordinate*/
        
        switch(k){ /*map to planar element*/
            case 1:/*identical elements*/
                for(i=0; i<SP->n; i++){
                    th1=SP->t[0][i];
                    th2=SP->t[1][i];
                    
                    SP->t[0][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[1][0]-vtx.vtxlist[0][0])*th1;
                    SP->t[1][i]=vtx.vtxlist[0][1]+(vtx.vtxlist[1][1]-vtx.vtxlist[0][1])*th1;
                    
                    SP->t[2][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[1][0]-vtx.vtxlist[0][0])*th2;
                    SP->t[3][i]=vtx.vtxlist[0][1]+(vtx.vtxlist[1][1]-vtx.vtxlist[0][1])*th2;
                    
                    SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][0]-vtx.vtxlist[0][0], 2)+pow(vtx.vtxlist[1][1]-vtx.vtxlist[0][1], 2)))*SP->wt[i];
                }
                break;
                
            case 0:/*common vertex*/
                
                for(i=0; i<SP->n; i++){                                    
                    th1=SP->t[0][i];
                    th2=SP->t[1][i];
                    
                    SP->t[0][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][1]-vtx.vtxlist[0][0])*th1;
                    SP->t[1][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][1]-vtx.vtxlist[1][0])*th1;
                    
                    SP->t[2][i]=vtx.vtxlist[0][0]+(vtx.vtxlist[0][2]-vtx.vtxlist[0][0])*th2;
                    SP->t[3][i]=vtx.vtxlist[1][0]+(vtx.vtxlist[1][2]-vtx.vtxlist[1][0])*th2;
                }
                for(i=0; i<SP->n/2; i++){
                    SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][1]-vtx.vtxlist[1][0], 2)+
                            pow(vtx.vtxlist[0][1]-vtx.vtxlist[0][0], 2)))*SP->wt[i];
                    SP->wt[i+SP->n/2]=(sqrt(pow(vtx.vtxlist[1][2]-vtx.vtxlist[1][0], 2)+
                            pow(vtx.vtxlist[0][2]-vtx.vtxlist[0][0], 2)))*SP->wt[i];
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
                    SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][1]-vtx.vtxlist[1][0], 2)+
                            pow(vtx.vtxlist[0][1]-vtx.vtxlist[0][0], 2)))*SP->wt[i];
                    SP->wt[i]=(sqrt(pow(vtx.vtxlist[1][3]-vtx.vtxlist[1][2], 2)+
                            pow(vtx.vtxlist[0][3]-vtx.vtxlist[0][2], 2)))*SP->wt[i];
                }
                break;
        }
        
    }
    
}
