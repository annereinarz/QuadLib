#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "quadpoints.h"
#include "BasicSubQuad.h"
#include "Tensorquad.h"
#include "squad.h"

#include "mpi.h"

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
}

void assert_correct_number_of_inputs (int actual, int expected);
int  get_integer_argument(char* argumentString, char* errorMsg);

Vertexlist get_vertexlist(int sim);
void free_vertexlist(Vertexlist vtxlist);
Vertexlist validate_results(int i, int sim);   //Input: sim=d

int factorial( int n )
{
    if ( n <= 1 )
        return 1;
    else
        return  n * factorial( n-1 );
}

int main(int argc, char** argv){
   int i;
   // Set up Mpi
   int my_rank;      // Rank of the process
   int p;            // Number of processes
   int tag=50;        // tag for messages
   MPI_Status status;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &p);

   int whichF=1;        // select which function to integrate, 1,2,3 or 4

  /*Get Input data */
   assert_correct_number_of_inputs (argc, 3);
   // number of regular quadrature points1
   int nrbg = get_integer_argument(argv[1], "first input argument must be integer!");
   int nred = get_integer_argument(argv[2], "second input argument must be integer!");
   // select dimensions: first number gives dimension in which to integrate and the
   // second gives the dimension of the overlap
   int sim = get_integer_argument(argv[3], "third input argument must be integer!");

   double Q;
   Vertexlist vtxlist;
   vtxlist=get_vertexlist(sim);

   // compute space dimension d and dimension of the intersection k
   int d= min(vtxlist.s1, vtxlist.s2);
   int k= 2*d+1 - max(vtxlist.s1, vtxlist.s2);

   clock_t begin=clock();
   //Calculate Affine Transformation
   AffineTrafo A;
   A=determineAffineTrafo(d,k,vtxlist);

   if(k<=0){
      if(my_rank==0){
        for(i=nrbg; i<=nred; i++){
            Q=squad_nonperm(whichF, A, d, k, i, 2*i);
            printf("sim=%d\nQ=%16.15lf\n",sim,Q);
        }
      }
      else{
         printf("Didn't need this process \n");
      }
   }
   else if(k>0){
      for(i=nrbg; i<=nred; i++){
         if(p!=(int)pow(2,k)-1){
            printf("Wanted %d processes \n Have %d processes\n",(int)pow(2,k)-1, p);
            return 1;
         }
         //printf("preparing for quadrature\n");
         int j, perm, N;
         if(my_rank==0){
            int cnt=0;
            for(j=0; j<=k-1; j++){
               N= factorial(k)/(factorial(j)*factorial(k-j));
               for(perm=0; perm<N; perm++){
                  MPI_Send(&j,    1, MPI_INT, cnt, 1, MPI_COMM_WORLD);
                  MPI_Send(&perm, 1, MPI_INT, cnt, 2, MPI_COMM_WORLD);
                  cnt++;
               }
            }
         }
         MPI_Recv(&j,    1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
         MPI_Recv(&perm, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
         Q=squad_perm(whichF, A, d, k, i, 2*i, j, perm,my_rank);
         //printf("Process %d has calculated quadrature value %lf\n", my_rank, Q);
         //Communication
         if(my_rank!=0){
            MPI_Send(&Q, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            //printf("Process %d has sent quadrature value\n",my_rank);
         }
         else{
            int j;
            double Qhelp=0;
            for(j=1; j<(int)pow(2,k)-1; j++){
               printf("Waiting for message from %d \n", j);
               MPI_Recv(&Qhelp, 1, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, &status);
               Q=Q+Qhelp;
            }
            printf("sim=%d\nQ=%16.15lf\n",sim,Q);
         }
      }
   }
   clock_t end=clock();
   if(my_rank==0){
      double diff=diffclock(end,begin);
      printf("time taken: %e\n",diff);
   }
   free_vertexlist(vtxlist);
   MPI_Finalize();
   return 0;
}


Vertexlist get_vertexlist(int sim){
 int i;
 Vertexlist vtxlist;

 switch(sim){
    case 1:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=-1;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=sqrt(3);

         break;
    case 2:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0.5;       vtxlist.vtxlist[1][2]=0.5*sqrt(3);

         break;

    case 4:
         vtxlist.s1=2;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=0.5;       vtxlist.vtxlist[1][1]=0.5*sqrt(3);
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-0.5;      vtxlist.vtxlist[1][3]=0.5*sqrt(3);

         break;

    case 5:
         vtxlist.s1=2;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0.5;       vtxlist.vtxlist[1][0]=0.5*sqrt(3);
         vtxlist.vtxlist[0][1]=0;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=1;         vtxlist.vtxlist[1][2]=0;
         vtxlist.vtxlist[0][3]=-0.5;      vtxlist.vtxlist[1][3]=0.5*sqrt(3);
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=sqrt(3);

         break;

      case 19:
         vtxlist.s1=2;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=-1;        vtxlist.vtxlist[1][3]=-1;
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=-1;
         vtxlist.vtxlist[0][5]=-1;        vtxlist.vtxlist[1][5]=0;

         break;

      case 20:
         vtxlist.s1=2;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=-1;        vtxlist.vtxlist[1][3]=0;
         vtxlist.vtxlist[0][4]=0;         vtxlist.vtxlist[1][4]=-1;
         break;

      case 21:
         vtxlist.s1=2;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         vtxlist.vtxlist[0][3]=0;         vtxlist.vtxlist[1][3]=-1;
         break;

      case 22:
         vtxlist.s1=2;
         vtxlist.s2=3;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;         vtxlist.vtxlist[1][0]=0;
         vtxlist.vtxlist[0][1]=1;         vtxlist.vtxlist[1][1]=0;
         vtxlist.vtxlist[0][2]=0;         vtxlist.vtxlist[1][2]=1;
         break;

      case 29:
         vtxlist.s1=3;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=-1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;
         vtxlist.vtxlist[0][6]=-1;   vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=-1;     vtxlist.vtxlist[2][7]=0;
         break;

       case 30:
         vtxlist.s1=3;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=-1;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;
         break;

       case 31:
         vtxlist.s1=3;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=-1;     vtxlist.vtxlist[2][4]=0;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=-1;
         break;


      case 32:
         vtxlist.s1=3;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         vtxlist.vtxlist[0][4]=0;   vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=-1;
         break;

      case 33:
         vtxlist.s1=3;
         vtxlist.s2=4;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;
         break;

      case 39:
         vtxlist.s1=4;
         vtxlist.s2=10;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=-1;   vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=-1;    vtxlist.vtxlist[3][5]=-1;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;     vtxlist.vtxlist[2][6]=-1;    vtxlist.vtxlist[3][6]=-1;
         vtxlist.vtxlist[0][7]=-1;   vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=-1;     vtxlist.vtxlist[3][7]=-1;
         vtxlist.vtxlist[0][8]=-1;   vtxlist.vtxlist[1][8]=-1;     vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=-1;
         vtxlist.vtxlist[0][9]=-1;   vtxlist.vtxlist[1][9]=-1;     vtxlist.vtxlist[2][9]=-1;     vtxlist.vtxlist[3][9]=0;

          break;

      case 40:
         vtxlist.s1=4;
         vtxlist.s2=9;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=-1;   vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=-1;     vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=-1;     vtxlist.vtxlist[3][7]=0;
         vtxlist.vtxlist[0][8]=0;    vtxlist.vtxlist[1][8]=0;      vtxlist.vtxlist[2][8]=0;     vtxlist.vtxlist[3][8]=-1;
         break;

      case 41:
         vtxlist.s1=4;
         vtxlist.s2=8;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=-1;     vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=-1;    vtxlist.vtxlist[3][6]=0;
         vtxlist.vtxlist[0][7]=0;    vtxlist.vtxlist[1][7]=0;      vtxlist.vtxlist[2][7]=0;     vtxlist.vtxlist[3][7]=-1;
         break;

      case 42:
         vtxlist.s1=4;
         vtxlist.s2=7;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=-1;    vtxlist.vtxlist[3][5]=0;
         vtxlist.vtxlist[0][6]=0;    vtxlist.vtxlist[1][6]=0;      vtxlist.vtxlist[2][6]=0;     vtxlist.vtxlist[3][6]=-1;
         break;

      case 43:
         vtxlist.s1=4;
         vtxlist.s2=6;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=0;     vtxlist.vtxlist[3][5]=-1;
         break;

      case 44:
         vtxlist.s1=4;
         vtxlist.s2=5;

         vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
         for(i=0; i<vtxlist.s1; i++){
            vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
         }

         vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;     vtxlist.vtxlist[3][0]=0;
         vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;     vtxlist.vtxlist[3][1]=0;
         vtxlist.vtxlist[0][2]=0;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;     vtxlist.vtxlist[3][2]=0;
         vtxlist.vtxlist[0][3]=0;    vtxlist.vtxlist[1][3]=0;      vtxlist.vtxlist[2][3]=1;     vtxlist.vtxlist[3][3]=0;
         vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=0;     vtxlist.vtxlist[3][4]=1;
         break;

   }
 return vtxlist;

}

 void free_vertexlist(Vertexlist vtxlist){
   int i;
   for(i=0; i<vtxlist.s1; i++){
       free(vtxlist.vtxlist[i]);
   }
   free(vtxlist.vtxlist);
}

void assert_correct_number_of_inputs (int actual, int expected) {
   if (actual-1 != expected){
      printf("Requires %d input arguments, got %d.\n", expected, actual-1);
      exit(1);
   }
}

int get_integer_argument(char* argumentString, char* errorMsg) {
   int argument;
   if (1 != sscanf(argumentString, "%d", &argument)) {
      printf("%s\n", errorMsg);
      exit(1);
   }
   return argument;
}


Vertexlist validate_results(int i, int sim){
   Vertexlist vtxlist;

   switch(i){
       case 1:   // identical elements
       vtxlist.s1=3;
       vtxlist.s2=4;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=1;
       break;

       case 2:   // elements with a common face
       vtxlist.s1=3;
       vtxlist.s2=5;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=0;      vtxlist.vtxlist[2][1]=0;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=1;      vtxlist.vtxlist[2][2]=1;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=1;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       break;

       case 3:   // elements with a common edge 1
       vtxlist.s1=3;
       vtxlist.s2=6;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=1;      vtxlist.vtxlist[2][1]=1;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       vtxlist.vtxlist[0][5]=1;    vtxlist.vtxlist[1][5]=0;      vtxlist.vtxlist[2][5]=1;
       break;

       case 4:   // elements with a common edge 2
       vtxlist.s1=3;
       vtxlist.s2=6;

       vtxlist.vtxlist = (double**)malloc(vtxlist.s1*sizeof(double*));
       for(i=0; i<vtxlist.s1; i++){
          vtxlist.vtxlist[i] = (double*)malloc(vtxlist.s2*sizeof(double));
       }
       vtxlist.vtxlist[0][0]=0;    vtxlist.vtxlist[1][0]=0;      vtxlist.vtxlist[2][0]=0;
       vtxlist.vtxlist[0][1]=1;    vtxlist.vtxlist[1][1]=1;      vtxlist.vtxlist[2][1]=1;
       vtxlist.vtxlist[0][2]=1;    vtxlist.vtxlist[1][2]=0;      vtxlist.vtxlist[2][2]=0;
       vtxlist.vtxlist[0][3]=1;    vtxlist.vtxlist[1][3]=1;      vtxlist.vtxlist[2][3]=0;
       vtxlist.vtxlist[0][4]=0;    vtxlist.vtxlist[1][4]=0;      vtxlist.vtxlist[2][4]=1;
       vtxlist.vtxlist[0][5]=0;    vtxlist.vtxlist[1][5]=1;      vtxlist.vtxlist[2][5]=1;
       break;
   }

   return vtxlist;
}

