#include <stdio.h>
#include <stdlib.h>
#include "time.h"
#include "math.h"
#include <string.h>
#include <argp.h>

#include "QuadRule.h"
#include "TensorQuad.h"
#include "cubequad.h"
#include "NextPoint.h"
#include "unique.h"
#include "sparsegrid.h"
#include "quadrature.h"
#include "vertexlist.h"
#include <unistd.h>

#define eps 0.3183099

struct arguments {
    int whichF;
    int d;
    int k;
    char* vtxlist;
    int regular;
    char* regular_file;
    int p;
    int T;
    int nr;
    int singular;
    int method;
    char *qrregular;
    char *qrsingular;
    int whichoutput;
    int portionwise;
    char *output;
};

#define QR_SOBOL (-2)
#define QR_FULL (-3)
#define QR_SPARSE (-1)

#define qr     (256+0)
#define qs     (256+1)
#define no_x   (256+2)
#define no_y   (256+3)
#define no_z   (256+4)
#define pw     (256+5)

static struct argp_option options[] = {
       {"function",                'f',   "NUMBER", 0, "Use a predefined function." },
       {"dimension",               'd',   "DIM",    0, "Spatial dimension (Don't specify with --vertexlist)." },
       {"intersection",            'k',   "DIM",    0, "Dimension of the intersection (Don't specify with --vertexlist)." },
       {"vertexlist",              'l',   "FILE",   0, "Read vertexlist from FILE to transform to physical parallelotope." },
       {"regular",                 'r',   "NAME",   0, "One dimensional quadrature rule for regular coordinates (Don't specify with --combination-method=sobol).\nNAME must be gauss-legendre, kronrod-patterson, or clenshaw-curtis or the path to a file containing a one-dimensional quadrature rule." },
       {"regular-sparse-p",        'p',   "NUMBER", 0, "Parameter used for sparse grid quadrature. (Needed iff --regular=sparse)" },
       {"regular-sparse-T",        'T',   "NUMBER", 0, "Parameter used for sparse grid quadrature. (Needed iff --regular=sparse and --regular-sparse-p=-4)" },
       {"n-regular",               'n',   "NUMBER", 0, "Number of quadrature points in each regular coordinate or in highest level of sparse-grid quadrature (Needs --regular)."},
       {"singular",                's',   "NAME",   0, "One dimensional quadrature rule for singular coordinate. NAME must be composite-gauss-legendre or gauss-jacobi" },
       {"combination-method",      'm',   "METHOD", 0, "Method for producing high-dimensional quadrature in regular coordinates. METHOD must be sparse, full-tensor, or sobol." },
       {"quadrature-rule-regular",  qr,   "FILE",   0, "Read regular quadrature rule from FILE (Don't specify with --regular or --combination-method)." },
       {"quadrature-rule-singular", qs,   "FILE",   0, "Read singular quadrature rule from FILE (Don't specify with --singular)." },
       {"no-x",                     no_x, 0,        0, "Don't calculate x-part of the quadrature." },
       {"no-y",                     no_y, 0,        0, "Don't calculate y-part of the quadrature." },
       {"no-z",                     no_z, 0,        0, "Don't calculate z-part of the quadrature." },
       {"portionwise",              pw,   0,        0, "Optimise memory usage by calculating only portions of fixed size of the quadrature rule before evaluating. If not set, the whole quadrature rule is held in memory, which is usually faster. Only use with --combination-method=sparse and when --function is set. " },
       {"output",                   'o',  "FILE",   0, "Output quadrature rules into FILE (Don't specify with --function). Note that the quadrature rules will be on the unit cube unless --vertexlist is set." },
       { 0 }
};

static error_t  parse_opt (int key, char *arg, struct argp_state *state){
   struct arguments *arguments = state->input;
   switch(key){
        case 'f' :
           if(1!=sscanf(arg,"%d",&arguments->whichF)){
              argp_error(state, "Argument to function must be number.\n");
           }
           break;
        case 'd':
           if(1!=sscanf(arg,"%d",&arguments->d)){
              argp_error(state, "Argument to dimension must be number.\n");
           }
           break;
        case 'k':
           if(1!=sscanf(arg,"%d",&arguments->k)){
              argp_error(state, "Argument to intersection must be number.\n");
           }
           break;
        case 'l':
           arguments->vtxlist = arg;
           break;
        case 'r':
           if(0==strncmp(arg,"gauss-legendre",sizeof("gauss-legendre"))){
                arguments->regular = GL;
           }
           else if(0==strncmp(arg,"clenshaw-curtis",sizeof("clenshaw-curtis"))){
                arguments->regular = CC;
           }
           else if(0==strncmp(arg,"kronrod-patterson",sizeof("kronrod-patterson"))){
                arguments->regular = KP;
           }
           else {
                arguments->regular_file = arg;
                if (!access(arg, F_OK) == 0){
                  argp_error(state, "The file \"%s\" cannot be accessed. Did you try to use one of the predefined quadrature rules clenshaw-curtis, kronrod-patterson, or gauss-legendre? Make sure the spelling is correct.\n");
                }
           }
           break;
        case 'p':
           if(1!=sscanf(arg,"%d",&arguments->p)){
              argp_error(state, "Argument to regular-sparse-p must be a number.\n");
           }
           break;
        case 'T':
           if(1!=sscanf(arg,"%d",&arguments->T)){
              argp_error(state, "Argument to regular-sparse-T must be a number.\n");
           }
           break;
        case 'n':
           if(1!=sscanf(arg,"%d",&arguments->nr)){
              argp_error(state, "Argument to nr must be number.\n");
           }
           break;
       case 's':
           if(0==strncmp(arg,"composite-gauss-legendre",sizeof("composite-gauss-legendre"))){
                arguments->singular = CGL;
           }
           else if(0==strncmp(arg,"gauss-jacobi",sizeof("gauss-jacobi"))){
                arguments->singular = GJ;
           }
           else {
            argp_error(state, "Argument singular must be gauss-jacobi, or composite-gauss-legendre\n");
           }
           break;
       case 'm'  :
           if(0==strncmp(arg,"sobol",sizeof("sobol"))){
                arguments->method = QR_SOBOL;
           }
           else if(0==strncmp(arg,"sparse",sizeof("sparse"))){
                arguments->method = QR_SPARSE;
           }
           else if(0==strncmp(arg,"full-tensor",sizeof("full-tensor"))){
                arguments->method = QR_FULL;
           }
           else {
            argp_error(state, "Argument combination-method must be sobol, sparse, or full-tensor\n");
           }
           break;
       case qr   :
           arguments->qrregular = arg;
           break;
       case qs   :
           arguments->qrsingular = arg;
           break;
       case no_x :
           arguments->whichoutput &= ~1;
           break;
       case no_y :
           arguments->whichoutput &= ~2;
           break;
       case no_z :
           arguments->whichoutput &= ~4;
           break;
       case pw :
           arguments->portionwise = 1;
       case 'o'  :
           arguments->output = arg;
           break;
        default:
           return ARGP_ERR_UNKNOWN;
     }
   return 0;
}

void set_nr(struct arguments arguments, int* nr){
   if(arguments.nr) {
      *nr = arguments.nr;
   }

   if(!*nr){
      fprintf(stderr,"Error: Need number of quadrature points. You have to set them using --n-regular=NUMBER or enter a file using --quadrature-rule-regular=FILE or --regular=FILE.\n");
      exit(1);
   }
}

void sparse_check(struct arguments *arguments){
   if(! arguments->p) {
      fprintf(stderr,"Warning: No value given for --regular-sparse-p. Using standard value p=1.\n");
      arguments->p = 1;
   }
   if(arguments->T==1 && arguments->p==-4) {
      fprintf(stderr,"Warning: No value given for --regular-sparse-T. Using standard value T=0.\n");
      arguments->T = 0;
   }
   if(! arguments->regular) {
      fprintf(stderr,"error: when using --combination-method=sparse, --regular can only be set to gauss-legendre, clenshaw-curtis or kronrod-patterson. importing custom rules from files is not implemented for sparse yet.\n");
   }
}

int get_quadrature_rule_regular(QuadRule* QR, struct arguments arguments, int* nr){
   if(arguments.nr) {
      *nr = arguments.nr;
   }
   if(!arguments.method && !arguments.regular && !arguments.qrregular && !arguments.regular_file){
      arguments.method  = QR_FULL;
      arguments.regular = GL;
      fprintf(stderr, "Warning: No rules for the regular coordinates were set. Assuming full tensor product Gauss-Legendre is wanted, i.e. --combination-method=full-tensor --regular=gauss-legendre.\n");
   }
   if(arguments.qrregular){
       //Einlesen;
       import_quadrule(arguments.qrregular, QR);
       if(! arguments.nr) {
         *nr = floor(pow(QR->n,1/(2*arguments.d-1)));
       }
       if(! arguments.qrsingular){
          fprintf(stderr,"Warning: Number of points for singular quadrature might not be optimal. You may want to choose your own singular quadrature rule with --quadrature-rule-singular.\n");
       }
       return 0;
   }

   QuadRule QR_1d;
   if(arguments.regular_file){
       import_quadrule(arguments.qrregular, &QR_1d);
       if(arguments.nr) {
          fprintf(stderr,"Warning: Setting --n-regular when using --regular=FILENAME has no effect.\n");
       }
       *nr = QR_1d.n;

       if(QR_1d.d != 1) {
          fprintf(stderr,"Error: Rule given by --regular=FILENAME is not 1-dimensional. If you want to provide the entire 2d-1-dimensional regular part of the quadrature, use --quadrature-rule-regular=FILENAME.\n");
       }
   }

   if(!*nr){
      fprintf(stderr,"Error: Need number of quadrature points. You have to set them using --n-regular=NUMBER or enter a file using --quadrature-rule-regular=FILE or --regular=FILE.\n");
      exit(1);
   }

   if(arguments.method==QR_SPARSE){ //sparse
      sparse_check(&arguments);
      *QR = Sparse( arguments.d,  arguments.nr,  arguments.p,  arguments.regular, arguments.T);
      return 0;
   }
   if(arguments.p){
      printf("Warning: Setting --regular-sparse-p without --combination-method=sparse has no effect\n");
   }
   if(arguments.method==QR_FULL){  //full TP
      if(arguments.regular){
         *QR = QuadRule_reg(arguments.d, arguments.nr, arguments.regular);
         return 0;
      }
      if(arguments.regular_file){
         *QR = QuadRule_reg_direct(arguments.d, arguments.nr, QR_1d);
         return 0;
      }
      fprintf(stderr,"Error: Either --quadrature-rule-regular or --regular must be set.\n");
      exit(1);
   }

   if(arguments.method==QR_SOBOL){
      if(arguments.regular || arguments.regular_file){
         fprintf(stderr,"Warning: Setting --regular with --combination-method=sobol has no effect.\n");
      }
      int N=(int)ceil(0.125*pow(arguments.nr,2*arguments.d-2));
      *QR = QuadRule_sobol(arguments.d,N);
      return 0;
   }

   fprintf(stderr,"Internal Error: Something went wrong when setting arguments.method!\n");
   exit(1);
}

QuadRule get_quadrature_rule_singular(struct arguments arguments, int nr, double alpha){
   QuadRule QR;
   //printf("nr=%d\n",nr);
   int ns;
   if(!arguments.singular && !arguments.qrsingular){
      arguments.singular = CGL;
      fprintf(stderr, "Warning: No rule for the singular coordinate was set. Assuming composite Gauss-Legendre is wanted, i.e. --singular=composite-gauss-legendre.\n");
   }

   if(arguments.qrsingular){
       //Einlesen
       import_quadrule(arguments.qrsingular, &QR);
       return QR;
   }

   if(arguments.k==-1){
      init_quadrule(&QR, nr,1);
      GLquad(&QR, nr,0,1);
      return QR;
   }

   if(arguments.singular == CGL){
      if(arguments.method==QR_SOBOL){
         ns=nr+1;
      }
      else if(arguments.method==QR_SPARSE){
         ns=2*((int)pow(2,nr)-1);
      }
      else{
         ns=2*nr;
      }
      init_quadrule(&QR, (int)(0.5*ns*(ns+1)),1);
      CGLquad(&QR, ns);
      return QR;
   }

   if (arguments.singular == GJ){
      if(arguments.method==QR_SOBOL){
         ns=nr+1;
      }
      else if(arguments.method==QR_SPARSE){
         ns=(int)pow(2,nr)-1;
      }
      else{
         ns=nr;
      }
      init_quadrule(&QR, ns, 1);
      GJquad01(&QR, ns, 0, 2*arguments.d-arguments.k-1+alpha);
      int i;
      for(i=0; i<QR.n; i++){
         QR.wt[i]=QR.wt[i]/pow(QR.t[0][i],alpha+2*arguments.d-arguments.k-1);
      }
      return QR;
   }
   fprintf(stderr,"Internal Error: Something went wrong when setting arguments.singular!\n");
   exit(1);
}

static struct argp argp = { options, parse_opt, "", "" };

char* generate_filename(struct arguments arguments){
   //char fname[FILENAME_MAX+1];
   if (arguments.vtxlist) {
      return "quadrature-rule-parallelotope.dat";
   }
   return "quadrature-rule-cube.dat";
}

int main(int argc, char **argv){
   int i;
   /* int a[22] = {19,20,21,22,29,30,31,32,33,39,40,41,42,43,44,49,50,51,52,53,54,55};
    Vertexlist vtxlist;
    char fname[] = "vertexlists/__.dat";
    for (int i=0; i<22; i++) {
      vtxlist=get_vertexlist(a[i]);
      sprintf(fname, "vertexlists/%d.dat", a[i]);
      export_vertexlist(vtxlist, fname);
    }*/
   clock_t start, end;
   double elapsed;
   int ndof=0;

   /*Set default values*/
   struct arguments arguments;
    arguments.T=1;
    arguments.whichF = 0;
    arguments.d = 0;
    arguments.k = 0;
    arguments.vtxlist = 0;
    arguments.regular = 0;
    arguments.regular_file = 0;
    arguments.p = 0;
    arguments.nr = 0;
    arguments.singular = 0;
    arguments.method = 0;
    arguments.qrregular = 0;
    arguments.qrsingular = 0;
    arguments.whichoutput = 7;
    arguments.portionwise = 0;
    arguments.output = 0;

    /*Get Input Data*/
    argp_parse(&argp, argc, argv, 0, 0, &arguments);


    if(argc==1){  //no input arguments
      fprintf(stderr, "Error: Need input arguments.\n");
      argp_help(&argp, stderr,ARGP_HELP_USAGE,argv[0]);
      exit(1);
    }


    /*Process Input Data*/
    if(arguments.portionwise && (! arguments.whichF || arguments.method != QR_SPARSE)){
       fprintf(stderr,"Error: You can only use --portionwise with --combination-method=sparse --function=NUMBER.\n");
    }

    if(arguments.whichF){
       if(arguments.whichoutput!=7){
          fprintf(stderr,"Warning: Using --no-x, --no-y or --no-z with --function has no effect.\n");
       }
       if(arguments.whichF==2 || arguments.whichF==4){
          arguments.whichoutput=6;
       }
       else if(arguments.whichF==11){
          arguments.whichoutput=7;
       }
       else{
          arguments.whichoutput=1;
       }

       if(! arguments.vtxlist)  {
          fprintf(stderr,"Error: --vertexlist must be set when using --function. If --function is set, the quadrature value on parallelotopes is calculated, so the position of parallelotopes is needed.\n");
          exit(1);
       }
    }
    AffineTrafo A;
    if(arguments.vtxlist) {
      Vertexlist vtxlist = import_vertexlist(arguments.vtxlist);
      if(arguments.d || arguments.k){
         fprintf(stderr,"Warning: Using --dimension or --intersection with --vertexlist has no effect.\n");
      }

      //compute space dimension d and dimension of the intersection k
      int d= vtxlist.s1;
      int k= 2*(d+1)-vtxlist.s2-1;

      if(arguments.d && arguments.d!=d){
         fprintf(stderr,"Warning: Spatial dimension set with --dimension is not identical with spatial dimension given by --vertexlist!\n");
      }
      arguments.d = d;
      if(arguments.k && arguments.k!=k){
         fprintf(stderr,"Warning: Dimension of intersection set with --intersection is not identical with dimension given by --vertexlist!\n");
      }
      arguments.k = k;

      //Calculate Affine Transformation
      A=determineAffineTrafo(d, k, vtxlist);
      free_vertexlist(vtxlist);
    }
    double alpha=-2*arguments.d+arguments.k+eps;

    FILE *f;
    if(!arguments.whichF){
       if(!arguments.output){
          arguments.output = generate_filename(arguments);
          fprintf(stderr, "Warning: No output file name given, writing in \"%s\". Use --output=FILENAME to set file name.\n",arguments.output);
       }
       f=fopen(arguments.output, "wb");
       if(f==NULL){
          fprintf(stderr,"Error: Could not open file \"%s\" for writing.\n", arguments.output);
          exit(1);
       }
    }

    int K[2*arguments.d-1];
    QuadRule QP1D[arguments.nr];
    if(arguments.portionwise){
         set_up_1dquadrules(arguments.nr,QP1D,arguments.regular);
         memset(K,0,(2*arguments.d-1)*sizeof(int));
    }

    QuadRule QR_reg, QR_sing, QR;
    int nr=0;
    if(arguments.portionwise){
       sparse_check(&arguments);
       set_nr(arguments, &nr);
    }
    double Q=0;

    int done = 0;
    while (! done) {
       if(arguments.portionwise){
           if (! next_sparse(&QR_reg, arguments.nr, QP1D, K, arguments.d, arguments.p, arguments.T)) {
              break;
           }
       }
       else{
          get_quadrature_rule_regular (&QR_reg, arguments, &nr);
          done = 1;
       }
       QR_sing = get_quadrature_rule_singular(arguments,  nr, alpha);
       //print_quadpoints(QR_sing);

       for(i=0; i<QR_sing.n; i++){
          init_quadrule(&QR, QR_reg.n, QR_reg.d+1);
          int m,l;
          for(m=0; m<QR_reg.n; m++){
             for(l=0; l<QR_reg.d; l++){
                QR.t[l+1][m]=QR_reg.t[l][m];
                QR.wt[m]=QR_sing.wt[i]*QR_reg.wt[m];
             }
             QR.t[0][m]=QR_sing.t[0][i];
          }
          if(arguments.whichF){
             start=clock();
             Q=Q+cubequad(arguments.k,arguments.d, A, arguments.whichF, arguments.whichoutput, &ndof, &QR);
             end=clock();
             elapsed = ((double)(end-start))/CLOCKS_PER_SEC;
          }
          else{
             if(!arguments.vtxlist){
                cubetransform(arguments.k,arguments.d,f,&QR);
             }
             else{
                cubeaffine(arguments.k,arguments.d,A,arguments.whichoutput, &QR, f);
             }
          }
          free_quadrule(QR);
       }
       free_quadrule(QR_sing);
       free_quadrule(QR_reg);
    }

    if(arguments.whichF){
      printf("%3.16lf, %d, %lf; \n", Q, ndof, elapsed);
    }
    if(!arguments.whichF){
       fclose(f);
    }
    return 0;
}
