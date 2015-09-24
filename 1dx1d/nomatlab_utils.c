#include <stdio.h>
#include <stdlib.h>

#include "quadpoints.h"
#include "nomatlab_utils.h"

#define BOOL    int
#define FALSE   0
#define TRUE    1

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

double get_double_argument(char* argumentString, char* errorMsg) { 
   double argument;
   if (1 != sscanf(argumentString, "%lf", &argument)) {
      printf("%s\n", errorMsg);
      exit(1);
   }
   return argument;
}

void read_vertexlist(Vertexlist *vtx){
   int m = 1;
   int n = 1;
   BOOL firstRow = TRUE;
   double* values;
   int count = 0;
   char chr;
   
   values = malloc (n*m*sizeof(double));
  
   printf("Enter a vertexlist");
 
   while (1) {
        if (1 != scanf("%lf", &values[count])) {
            printf("double expected on input (%d values read successfully into %dx%d matrix)!\n", count, n, m);
            exit(1);
        }
        count += 1;
        if (feof(stdin)) {
            break;
        }
        if (' ' == getchar()) { /* same row */
            if (firstRow) {
                n += 1;
            }
        } else { /* new row */
            m += 1;
            firstRow = FALSE;
        
        }
        chr = getchar(); /* FIXME: for some reason, feof(stdin) does not work here. */
        if (EOF == chr) { /* newline at end of input */
            m -= 1;
            break;
        }
        ungetc(chr, stdin);
        /* TOD0: Kill whitespace? */
        values = realloc (values, n*m*sizeof(double));

   }
      if (count != n*m) {
        printf("input must be a matrix!\n");
        exit(1);
   }
   
   /* Copy values into vertex list */
   init_vtxlist(vtx, m, n);
   
   {
   int i;
   for(i=0; i<vtx->s1*vtx->s2; i++){
       vtx->vtxlist[i/vtx->s2][i%vtx->s2]=values[i];
   }
   }
   
   free(values);
}

void print_quadpoints(QuadraturePoints SP) {
  int i,j;
   printf("\n");
  for (i=0; i<SP.n; i++) {
    for (j=0; j<SP.d; j++) {
      printf("%lf ", SP.t[j][i]);
    }
    printf("\n");
  }
  printf("\n");
  
  for (i=0; i<SP.n; i++) {
      printf("%lf\n", SP.wt[i]);
  }
}
