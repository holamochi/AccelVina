#ifndef VINA_PACKED_BFGS_H
#define VINA_PACKED_BFGS_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
double vina_bfgs_dot(size_t n, const double *a, const double *b);
void vina_bfgs_product(size_t n, const double *matrix, const double *x, double *y);
void vina_bfgs_update(size_t n, double *matrix, const double *p, const double *y,
                      double *minus_hy, double alpha, double yp);
#ifdef __cplusplus
}
#endif
#endif
