/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_BFGS_H
#define VINA_BFGS_H

#include "matrix.h"
#include <cstring>
#include "packed_bfgs.h"

typedef triangular_matrix<fl> flmat;

template<typename Change>
void minus_mat_vec_product(const flmat& m, const Change& in, Change& out, flv& in_flat) {
	sz n = m.dim();
	VINA_FOR(j, n)
		in_flat[j] = in(j);
	VINA_FOR(i, n) {
		fl sum = 0;
		VINA_FOR(j, n)
			sum += m(m.index_permissive(i, j)) * in_flat[j];
		out(i) = -sum;
	}
}

template<typename Change>
void minus_mat_vec_product(const flmat& m, const Change& in, Change& out) {
	if (&in == &out) {
		sz n = m.dim();
		VINA_FOR(i, n) {
			fl sum = 0;
			VINA_FOR(j, n)
				sum += m(m.index_permissive(i, j)) * in(j);
			out(i) = -sum;
		}
		return;
	}
	flv in_flat(m.dim());
	minus_mat_vec_product(m, in, out, in_flat);
}

template<typename Change>
inline fl scalar_product(const Change& a, const Change& b, sz n) {
	fl tmp = 0;
	VINA_FOR(i, n)
		tmp += a(i) * b(i);
	return tmp;
}

template<typename Change>
inline bool bfgs_update(flmat& h, const Change& p, const Change& y, const fl alpha,
                        flv& in_flat, flv& p_flat, flv& minus_hy_flat) {
	const fl yp  = scalar_product(y, p, h.dim());
	if(alpha * yp < epsilon_fl) return false; // FIXME?
	Change minus_hy(y); minus_mat_vec_product(h, y, minus_hy, in_flat);
	const fl yhy = - scalar_product(y, minus_hy, h.dim());
	const fl r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
	const sz n = p.num_floats();
	VINA_FOR(i, n) {
		p_flat[i] = p(i);
		minus_hy_flat[i] = minus_hy(i);
	}
	VINA_FOR(i, n)
		VINA_RANGE(j, i, n) // includes i
			h(i, j) +=   alpha * r * (minus_hy_flat[i] * p_flat[j]
	                                + minus_hy_flat[j] * p_flat[i]) +
			           + alpha * alpha * (r*r * yhy  + r) * p_flat[i] * p_flat[j]; // s * s == alpha * alpha * p * p
	return true;
}

template<typename Change>
inline bool bfgs_update(flmat& h, const Change& p, const Change& y, const fl alpha) {
	flv in_flat(h.dim()), p_flat(p.num_floats()), minus_hy_flat(p.num_floats());
	return bfgs_update(h, p, y, alpha, in_flat, p_flat, minus_hy_flat);
}

#ifndef FUSION_PACKED_BFGS
// Keep the public triangular matrix API; only the optimizer owns mirrored storage.
class bfgs_dense_matrix {
	flv m_data;
	sz m_dim;
public:
	bfgs_dense_matrix(sz n, fl filler) : m_data(n * n, filler), m_dim(n) {}
	sz dim() const { return m_dim; }
	const fl& operator()(sz i, sz j) const { return m_data[i + j * m_dim]; }
	fl& operator()(sz i, sz j) { return m_data[i + j * m_dim]; }
};

typedef fl bfgs_lane2 __attribute__((vector_size(16)));

template<typename Change>
void minus_mat_vec_product(const bfgs_dense_matrix& m, const Change& in,
                           Change& out, flv& in_flat) {
	const sz n = m.dim();
	VINA_FOR(j, n) in_flat[j] = in(j);
	sz i = 0;
	for (; i + 1 < n; i += 2) {
		bfgs_lane2 sums = {0, 0};
		VINA_FOR(j, n) {
			bfgs_lane2 coefficients;
			std::memcpy(&coefficients, &m(i, j), sizeof(coefficients));
			const bfgs_lane2 input = {in_flat[j], in_flat[j]};
			sums += coefficients * input;
		}
		fl lanes[2];
		std::memcpy(lanes, &sums, sizeof(sums));
		out(i) = -lanes[0];
		out(i + 1) = -lanes[1];
	}
	for (; i < n; ++i) {
		fl sum = 0;
		VINA_FOR(j, n) sum += m(i, j) * in_flat[j];
		out(i) = -sum;
	}
}

template<typename Change>
inline bool bfgs_update(bfgs_dense_matrix& h, const Change& p, const Change& y,
                        const fl alpha, flv& in_flat, flv& p_flat,
                        flv& minus_hy_flat) {
	const fl yp = scalar_product(y, p, h.dim());
	if (alpha * yp < epsilon_fl) return false;
	Change minus_hy(y);
	minus_mat_vec_product(h, y, minus_hy, in_flat);
	const fl yhy = -scalar_product(y, minus_hy, h.dim());
	const fl r = 1 / (alpha * yp);
	const sz n = p.num_floats();
	VINA_FOR(i, n) {
		p_flat[i] = p(i);
		minus_hy_flat[i] = minus_hy(i);
	}
	VINA_FOR(j, n)
		for (sz i = 0; i <= j; ++i) {
			h(i, j) += alpha * r * (minus_hy_flat[i] * p_flat[j]
			                         + minus_hy_flat[j] * p_flat[i]) +
			             + alpha * alpha * (r*r * yhy + r) * p_flat[i] * p_flat[j];
			h(j, i) = h(i, j);
		}
	return true;
}
#endif

template<typename F, typename Conf, typename Change>
fl line_search_trial(F& f, Conf& x, Change& g) {
	return f(x, g);
}

template<typename F, typename Conf, typename Change>
void line_search_finish(F&, Conf&, Change&) {}

template<typename F, typename Conf, typename Change>
fl line_search(F& f, sz n, const Conf& x, const Change& g, const fl f0, const Change& p, Conf& x_new, Change& g_new, fl& f1, int& evalcount) { // returns alpha
	const fl c0 = 0.0001;
	const unsigned max_trials = 10;
	const fl multiplier = 0.5;
	fl alpha = 1;

	const fl pg = scalar_product(p, g, n);

	VINA_U_FOR(trial, max_trials) {
		x_new = x; x_new.increment(p, alpha);
		f1 = line_search_trial(f, x_new, g_new);
		evalcount++;
		if(f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
			break;
		alpha *= multiplier;
	}
	line_search_finish(f, x_new, g_new);
	return alpha;
}

template<typename Matrix>
inline void set_diagonal(Matrix& m, fl x) {
	VINA_FOR(i, m.dim())
		m(i, i) = x;
}

template<typename Change>
void subtract_change(Change& b, const Change& a, sz n) { // b -= a
	VINA_FOR(i, n)
		b(i) -= a(i);
}

template<typename F, typename Conf, typename Change>
fl bfgs(F& f, Conf& x, Change& g, const unsigned max_steps, const fl average_required_improvement, const sz over,
        int& evalcount) {
    const sz n = g.num_floats();
    flmat h(n, 0);
    set_diagonal(h, 1);
    std::vector<fl> storage(5*n+1);
    fl* gf = storage.data();
    fl* gn = gf + n;
    fl* pf = gn + n;
    fl* yf = pf + n;
    fl* hy = yf + n;
    Change g_new(g), p(g);
    Conf x_new(x);
    fl f0 = f(x, g);
    ++evalcount;
    const fl f_orig = f0;
    const Change g_orig(g);
    const Conf x_orig(x);
    g.to_flat(gf);
    for (unsigned step = 0; step < max_steps; ++step) {
        const fl* triangle = n ? &h(0,0) : NULL;
        vina_bfgs_product(n, triangle, gf, pf);
        p.from_flat(pf);
        fl f1 = 0;
        const fl alpha = line_search(f, n, x, g, f0, p, x_new, g_new, f1, evalcount);
        g_new.to_flat(gn);
        for (sz i = 0; i < n; ++i) yf[i] = gn[i] - gf[i];
        f0 = f1;
        x = x_new;
        if (!(std::sqrt(vina_bfgs_dot(n, gf, gf)) >= 1e-5)) break;
        g = g_new;
        std::swap(gn, gf);
        const fl yp = vina_bfgs_dot(n, yf, pf);
        if (step == 0) {
            const fl yy = vina_bfgs_dot(n, yf, yf);
            if (std::abs(yy) > epsilon_fl) set_diagonal(h, alpha * yp / yy);
        }
        if (!(alpha * yp < epsilon_fl) && n)
            vina_bfgs_update(n, &h(0,0), pf, yf, hy, alpha, yp);
    }
    if (!(f0 <= f_orig)) {
        f0 = f_orig;
        x = x_orig;
        g = g_orig;
    }
    return f0;
}

#endif
