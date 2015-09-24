#ifndef VERTEXLIST_H
#define VERTEXLIST_H


typedef struct {
    int s1;            /*size of vertexlist in first dimension*/
    int s2;            /*size of vertexlist in second dimension*/
    double** vtxlist;     /*vertex list, has the size s1 x s2*/
} Vertexlist;

//Vertexlist get_vertexlist(int sim);
void free_vertexlist(Vertexlist vtxlist);
Vertexlist import_vertexlist(char* filename);
void export_vertexlist(Vertexlist vtxlist, char *fname);
void free_vertexlist(Vertexlist vtxlist);

#endif
