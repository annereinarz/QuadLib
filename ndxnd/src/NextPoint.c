#include "NextPoint.h"

int NextPoint(int* I, int d, int l){
  int k=0;
  while(k<d){
     I[k] = (I[k] + 1) % l;
     if(I[k]!=0){
        return 1;
     }
     k++;
  }
  I[0]=-1;
  return 0;
}

int morePoints(int* I) {
  return I[0]!=-1;
}
