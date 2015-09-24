#ifndef MXMEMORY_H
#define MXMEMORY_H

void *mxMalloc(int n);
void *mxCalloc(int n, int size);
void *mxRealloc(void *ptr, int size);
void mxFree(void *ptr);

int mexPrintf(const char *message, ...);

#endif /* MXMEMORY_H */