#include <stdlib.h>
#include "mxMemory.h"

void *mxMalloc(int n) {
    return malloc(n);
}

void *mxCalloc(int n, int size) {
    return calloc(n, size);
}

void *mxRealloc(void *ptr, int size) {
    return realloc(ptr, size);
}

void mxFree(void *ptr) {
    free(ptr);
}

/*
int mexPrintf(const char *message, ...){
    return printf(message, ...);
}
*/