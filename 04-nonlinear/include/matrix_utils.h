#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include "HYPRE_utilities.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

void init_hypre();
void finalize_hypre();
void solve_with_hypre(HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x);

#endif