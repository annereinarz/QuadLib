#ifndef NOMATLAB_UTILS_H
#define NOMATLAB_UTILS_H

#include "quadpoints.h"

void assert_correct_number_of_inputs (int actual, int expected);
int    get_integer_argument(char* argumentString, char* errorMsg);
double get_double_argument (char* argumentString, char* errorMsg);
void read_vertexlist(Vertexlist *vtx);
void print_quadpoints(QuadraturePoints SP);

#endif /*NOMATLAB_UTILS_H*/
