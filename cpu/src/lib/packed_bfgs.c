#include "packed_bfgs.h"

double vina_bfgs_dot(size_t n, const double *a, const double *b) {
    double sum = 0;
    for (size_t i = 0; i < n; ++i)
        sum += a[i] * b[i];
    return sum;
}

void vina_bfgs_product(size_t n, const double *restrict matrix, const double *restrict x,
                       double *restrict y) {
    /* Column j holds entries (0,j)..(j,j). Each y[i] receives its j terms
       from left to right, with the same initial +0 as the upstream dot. */
    const double *column = matrix;
    for (size_t j = 0; j < n; ++j) {
        double diagonal_row = 0;
        for (size_t i = 0; i < j; ++i) {
            diagonal_row += column[i] * x[i];
            y[i] += column[i] * x[j];
        }
        y[j] = diagonal_row + column[j] * x[j];
        column += j + 1;
    }
    for (size_t i = 0; i < n; ++i)
        y[i] = -y[i];
}

void vina_bfgs_update(size_t n, double *restrict matrix, const double *restrict p,
                      const double *restrict y, double *restrict minus_hy,
                      double alpha, double yp) {
    vina_bfgs_product(n, matrix, y, minus_hy);
    const double yhy = -vina_bfgs_dot(n, y, minus_hy);
    const double r = 1 / (alpha * yp);
    const double linear = alpha * r;
    const double quadratic = alpha * alpha * (r * r * yhy + r);
    double *column = matrix;
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i <= j; ++i)
            column[i] += linear * (minus_hy[i] * p[j] + minus_hy[j] * p[i])
                         + quadratic * p[i] * p[j];
        column += j + 1;
    }
}
