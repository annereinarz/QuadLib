#ifndef VTX_H
#define VTX_H

#include "AffineTrafo.h"

Vertexlist validate_results(int i, int sim);
Vertexlist get_vertexlist(int sim);
void free_vertexlist(Vertexlist vtxlist);  
void import_vtx(Vertexlist *vtx, FILE* f);

#endif /*VTX_H*/
