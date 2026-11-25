/* SPDX-License-Identifier: Apache-2.0 */
#include "c_kernel.h"

#include <float.h>
#include <immintrin.h>

static inline __m256d load_xyz(const unsigned char *address)
{
    return _mm256_maskz_loadu_pd(0x7, address);
}

static inline void store_xyz(unsigned char *address, __m256d xyz)
{
    _mm256_mask_storeu_pd(address, 0x7, xyz);
}

static inline __m256d add_force(__m256d force_a, unsigned char *force_b, __m256d force)
{
    force_a = _mm256_sub_pd(force_a, force);
    store_xyz(force_b, _mm256_add_pd(load_xyz(force_b), force));
    return force_a;
}

double vina_ck_pairs_deriv(const vina_ck_pair_row *rows, size_t row_count,
                          const vina_ck_pair_chunk *chunks, const double *smooth,
                          const double *x, const double *y, const double *z,
                          void *forces, double factor, double cutoff_sqr, double cap)
{
    unsigned char *const force_bytes = forces;
    const size_t stride = 3 * sizeof(double);
    const __m256d zero = _mm256_setzero_pd();
    const __m256d cutoff = _mm256_set1_pd(cutoff_sqr);
    const __m256d scale = _mm256_set1_pd(factor);
    const __m256d vc = _mm256_set1_pd(cap);
    const int curl_enabled = cap < 0.1 * DBL_MAX;
    double energy = 0;

    for (size_t i = 0; i < row_count; ++i) {
        const vina_ck_pair_row *row = &rows[i];
        const __m256d ax = _mm256_set1_pd(x[row->a]);
        const __m256d ay = _mm256_set1_pd(y[row->a]);
        const __m256d az = _mm256_set1_pd(z[row->a]);
        unsigned char *const force_a = forces ? force_bytes + stride * row->a : NULL;
        __m256d fa = forces ? load_xyz(force_a) : zero;

        for (size_t j = row->begin; j < row->end; ++j) {
            const vina_ck_pair_chunk *chunk = &chunks[j];
            const size_t b = chunk->b;
            const __mmask8 slots = (__mmask8)chunk->mask;
            const __m256d dx = _mm256_maskz_sub_pd(slots, _mm256_loadu_pd(x + b), ax);
            const __m256d dy = _mm256_maskz_sub_pd(slots, _mm256_loadu_pd(y + b), ay);
            const __m256d dz = _mm256_maskz_sub_pd(slots, _mm256_loadu_pd(z + b), az);
            const __m256d dx2 = _mm256_mul_pd(dx, dx);
            const __m256d dy2 = _mm256_mul_pd(dy, dy);
            const __m256d dz2 = _mm256_mul_pd(dz, dz);
            const __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
            const __mmask8 active = _mm256_mask_cmp_pd_mask(slots, r2, cutoff, _CMP_LT_OQ);
            if (!active)
                continue;

            const __m256d scaled = _mm256_maskz_mul_pd(active, scale, r2);
            const __m256i bin = _mm256_cvttpd_epi64(scaled);
            const __m256d rem = _mm256_sub_pd(scaled, _mm256_cvtepi64_pd(bin));
            const __m256i table = _mm256_loadu_si256((const __m256i_u *)chunk->table);
            const __m256i sample = _mm256_add_epi64(table, _mm256_slli_epi64(bin, 1));
            const __m256d e0 = _mm256_mmask_i64gather_pd(zero, active, sample, smooth, 8);
            const __m256d e1 = _mm256_mmask_i64gather_pd(zero, active, sample, smooth + 2, 8);
            __m256d e = _mm256_add_pd(e0, _mm256_mul_pd(rem, _mm256_sub_pd(e1, e0)));
            const __mmask8 positive = curl_enabled ? _mm256_mask_cmp_pd_mask(active, e, zero, _CMP_GT_OQ) : 0;
            __m256d t2 = zero;
            if (positive) {
                const __m256d t = cap < DBL_EPSILON ? zero :
                    _mm256_maskz_div_pd(positive, vc, _mm256_maskz_add_pd(positive, vc, e));
                e = _mm256_mask_mul_pd(e, positive, e, t);
                t2 = _mm256_mul_pd(t, t);
            }
            if (!forces) {
                double lane[4];
                _mm256_storeu_pd(lane, e);
                for (unsigned k = 0; k < 4; ++k)
                    if (active & (1u << k)) energy += lane[k];
                continue;
            }
            const __m256d d0 = _mm256_mmask_i64gather_pd(zero, active, sample, smooth + 1, 8);
            const __m256d d1 = _mm256_mmask_i64gather_pd(zero, active, sample, smooth + 3, 8);
            const __m256d dor = _mm256_add_pd(d0, _mm256_mul_pd(rem, _mm256_sub_pd(d1, d0)));
            __m256d fx = _mm256_maskz_mul_pd(active, dor, dx);
            __m256d fy = _mm256_maskz_mul_pd(active, dor, dy);
            __m256d fz = _mm256_maskz_mul_pd(active, dor, dz);
            if (positive) {
                fx = _mm256_mask_mul_pd(fx, positive, fx, t2);
                fy = _mm256_mask_mul_pd(fy, positive, fy, t2);
                fz = _mm256_mask_mul_pd(fz, positive, fz, t2);
            }

            /* SIMD lanes never participate in a horizontal energy or force reduction. */
            const __m256d xy_even = _mm256_unpacklo_pd(fx, fy);
            const __m256d xy_odd = _mm256_unpackhi_pd(fx, fy);
            const __m256d z_even = _mm256_unpacklo_pd(fz, zero);
            const __m256d z_odd = _mm256_unpackhi_pd(fz, zero);
            const __m128d e_lo = _mm256_castpd256_pd128(e);
            const __m128d e_hi = _mm256_extractf128_pd(e, 1);
            unsigned char *const fb = force_bytes + stride * b;
            if (active & 0x1) {
                energy += _mm_cvtsd_f64(e_lo);
                fa = add_force(fa, fb, _mm256_permute2f128_pd(xy_even, z_even, 0x20));
            }
            if (active & 0x2) {
                energy += _mm_cvtsd_f64(_mm_unpackhi_pd(e_lo, e_lo));
                fa = add_force(fa, fb + stride, _mm256_permute2f128_pd(xy_odd, z_odd, 0x20));
            }
            if (active & 0x4) {
                energy += _mm_cvtsd_f64(e_hi);
                fa = add_force(fa, fb + 2 * stride, _mm256_permute2f128_pd(xy_even, z_even, 0x31));
            }
            if (active & 0x8) {
                energy += _mm_cvtsd_f64(_mm_unpackhi_pd(e_hi, e_hi));
                fa = add_force(fa, fb + 3 * stride, _mm256_permute2f128_pd(xy_odd, z_odd, 0x31));
            }
        }

        if (forces) store_xyz(force_a, fa);
    }

    return energy;
}

typedef struct {
    __m512d maps01;
    __m256d map2;
} grid_value;

static inline __m512d join_maps(__m256d a, __m256d b)
{
    return _mm512_insertf64x4(_mm512_castpd256_pd512(a), b, 1);
}

static inline grid_value grid_term(const double *const maps[3], size_t corner,
                                   __m256d wx, __m256d wy, __m256d wz)
{
    const __m512d fs = join_maps(_mm256_broadcast_sd(maps[0] + corner),
                                _mm256_broadcast_sd(maps[1] + corner));
    const __m256d f2 = _mm256_broadcast_sd(maps[2] + corner);
    const __m512d vx = _mm512_mul_pd(fs, join_maps(wx, wx));
    const __m512d vxy = _mm512_mul_pd(vx, join_maps(wy, wy));
    const __m256d v2x = _mm256_mul_pd(f2, wx);
    const __m256d v2xy = _mm256_mul_pd(v2x, wy);
    return (grid_value){_mm512_mul_pd(vxy, join_maps(wz, wz)), _mm256_mul_pd(v2xy, wz)};
}

static inline grid_value grid_add(grid_value a, grid_value b)
{
    return (grid_value){_mm512_add_pd(a.maps01, b.maps01), _mm256_add_pd(a.map2, b.map2)};
}

static inline double finish_grid(__m256d value, double penalty, double slope,
                                 double cap, const double factor[3],
                                 const int region[3], double gradient[3])
{
    double v[4];
    _mm256_storeu_pd(v, value);
    if (v[0] > 0 && cap < 0.1 * DBL_MAX) {
        const double t = cap < DBL_EPSILON ? 0 : cap / (cap + v[0]);
        v[0] *= t;
        const double t2 = t * t;
        v[1] *= t2;
        v[2] *= t2;
        v[3] *= t2;
    }
    for (size_t i = 0; i < 3; ++i)
        gradient[i] = factor[i] * (region[i] == 0 ? v[i + 1] : 0) + slope * region[i];
    return v[0] + penalty;
}

void vina_ck_ad4_deriv(const double *const maps[3], const size_t dims[3],
                      const double origin[3], const double factor[3],
                      const double factor_inv[3], const double last[3],
                      const double location[3], double slope, double cap,
                      double values[3], double *const gradients[3])
{
    double s[3];
    double miss[3] = {0, 0, 0};
    int region[3] = {0, 0, 0};
    size_t a[3];
    for (size_t i = 0; i < 3; ++i) {
        s[i] = (location[i] - origin[i]) * factor[i];
        if (s[i] < 0) {
            miss[i] = -s[i];
            region[i] = -1;
            a[i] = 0;
            s[i] = 0;
        } else if (s[i] >= last[i]) {
            miss[i] = s[i] - last[i];
            region[i] = 1;
            a[i] = dims[i] - 2;
            s[i] = 1;
        } else {
            a[i] = (size_t)s[i];
            s[i] -= (double)a[i];
        }
    }
    const double penalty = slope * ((miss[0] * factor_inv[0] + miss[1] * factor_inv[1])
                                    + miss[2] * factor_inv[2]);
    const size_t row = dims[0];
    const size_t plane = row * dims[1];
    const size_t base = a[0] + row * (a[1] + dims[1] * a[2]);
    const size_t p01 = base + row;
    const size_t p10 = base + plane;
    const size_t p11 = p10 + row;

    /* Lanes within each map are energy, dx, dy, dz, in that order. */
    const double mx = 1 - s[0], my = 1 - s[1], mz = 1 - s[2];
    const __m256d x0 = _mm256_setr_pd(mx, -1, mx, mx);
    const __m256d x1 = _mm256_setr_pd(s[0], 1, s[0], s[0]);
    const __m256d y0 = _mm256_setr_pd(my, my, -1, my);
    const __m256d y1 = _mm256_setr_pd(s[1], s[1], 1, s[1]);
    const __m256d z0 = _mm256_setr_pd(mz, mz, mz, -1);
    const __m256d z1 = _mm256_setr_pd(s[2], s[2], s[2], 1);

    grid_value sum =    grid_term(maps, base,     x0, y0, z0);
    sum = grid_add(sum, grid_term(maps, base + 1, x1, y0, z0));
    sum = grid_add(sum, grid_term(maps, p01,      x0, y1, z0));
    sum = grid_add(sum, grid_term(maps, p01 + 1,  x1, y1, z0));
    sum = grid_add(sum, grid_term(maps, p10,      x0, y0, z1));
    sum = grid_add(sum, grid_term(maps, p10 + 1,  x1, y0, z1));
    sum = grid_add(sum, grid_term(maps, p11,      x0, y1, z1));
    sum = grid_add(sum, grid_term(maps, p11 + 1,  x1, y1, z1));

    values[0] = finish_grid(_mm512_castpd512_pd256(sum.maps01), penalty, slope, cap, factor,
                            region, gradients[0]);
    values[1] = finish_grid(_mm512_extractf64x4_pd(sum.maps01, 1), penalty, slope, cap, factor,
                            region, gradients[1]);
    values[2] = finish_grid(sum.map2, penalty, slope, cap, factor, region, gradients[2]);
}
