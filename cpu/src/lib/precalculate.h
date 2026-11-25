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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include <map>
#include <atomic>
#include <cstring>
#include "scoring_function.h"
#include "matrix.h"


//Forward declaration
struct model;

struct precalculate_element
{
public:
    precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) { }
    fl eval_fast(fl r2) const{
        assert(r2 * factor < fast.size());
        sz i = sz(factor * r2); // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        assert(i < fast.size());
        return fast[i];
    };
    pr eval_deriv(fl r2) const{
        fl r2_factored = factor * r2;
        assert(r2_factored + 1 < smooth.size());
        sz i1 = sz(r2_factored);
        sz i2 = i1 + 1; // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
        assert(i1 < smooth.size());
        assert(i2 < smooth.size());
        fl rem = r2_factored - i1;
        assert(rem >= -epsilon_fl);
        assert(rem < 1 + epsilon_fl);
        const pr &p1 = smooth[i1];
        const pr &p2 = smooth[i2];
        fl e = p1.first + rem * (p2.first - p1.first);
        fl dor = p1.second + rem * (p2.second - p1.second);
        return pr(e, dor);
    };
    void init_from_smooth_fst(const flv &rs){
        sz n = smooth.size();
        VINA_CHECK(rs.size() == n);
        VINA_CHECK(fast.size() == n);
        VINA_FOR(i, n) {
            // calculate dor's
            fl &dor = smooth[i].second;
            if (i == 0 || i == n - 1)
                dor = 0;
            else {
                fl delta = rs[i + 1] - rs[i - 1];
                fl r = rs[i];
                dor = (smooth[i + 1].first - smooth[i - 1].first) / (delta * r);
            }
            // calculate fast's from smooth.first's
            fl f1 = smooth[i].first;
            fl f2 = (i + 1 >= n) ? 0 : smooth[i + 1].first;
            fast[i] = (f2 + f1) / 2;
        }
    };
    sz min_smooth_fst() const{
        sz tmp = 0; // returned if smooth.empty()
        VINA_FOR_IN(i_inv, smooth) {
            sz i = smooth.size() - i_inv - 1; // i_inv < smooth.size()  => i_inv + 1 <= smooth.size()
            if (i_inv == 0 || smooth[i].first < smooth[tmp].first)
                tmp = i;
        }
        return tmp;
    };
    void widen_smooth_fst(const flv &rs, fl left, fl right){
        flv tmp(smooth.size(), 0); // the new smooth[].first's
        sz min_index = min_smooth_fst();
        VINA_CHECK(min_index < rs.size()); // won't hold for n == 0
        VINA_CHECK(rs.size() == smooth.size());
        fl optimal_r = rs[min_index];
        VINA_FOR_IN(i, smooth)
        {
            fl r = rs[i];
            if (r < optimal_r - left)
                r += left;
            else if (r > optimal_r + right)
                r -= right;
            else
                r = optimal_r;

            if (r < 0)
                r = 0;
            if (r > rs.back())
                r = rs.back();

            tmp[i] = eval_deriv(sqr(r)).first;
        }
        VINA_FOR_IN(i, smooth)
        smooth[i].first = tmp[i];
    };
    void widen(const flv &rs, fl left, fl right){
        widen_smooth_fst(rs, left, right);
        init_from_smooth_fst(rs);
    };
    prv smooth; // [(e, dor)]
    fl get_factor() const { return factor; }
    const flv& fast_samples() const { return fast; }
private:
    flv fast;
    fl factor;
};

struct precalculate 
{
public:
    precalculate() { }
    precalculate(const ScoringFunction& sf, fl v=max_fl, fl factor=32){
        m_factor = factor;
        m_cutoff_sqr = sqr(sf.get_cutoff());
        m_max_cutoff_sqr = sqr(sf.get_max_cutoff());
        m_n = sz(m_factor * m_max_cutoff_sqr) + 3; // sz(factor * r^2) + 1 <= sz(factor * max_cutoff_sqr) + 2 <= n-1 < n  // see assert below
        //std::cout << "-- DEBUG -- sf.cutoff^2 in precalculate = " << m_max_cutoff_sqr << "\n";
        triangular_matrix<precalculate_element> data(num_atom_types(sf.get_atom_typing()), precalculate_element(m_n, m_factor));

        VINA_CHECK(m_factor > epsilon_fl);
        VINA_CHECK(sz(m_max_cutoff_sqr * m_factor) + 1 < m_n); // max_cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
        VINA_CHECK(m_max_cutoff_sqr * m_factor + 1 < m_n);

        flv rs = calculate_rs();

        VINA_FOR(t1, data.dim())
        {
            VINA_RANGE(t2, t1, data.dim())
            {
                precalculate_element &p = data(t1, t2);
                // init smooth[].first
                VINA_FOR_IN(i, p.smooth)
                {
                    p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));
                }

                // init the rest
                p.init_from_smooth_fst(rs);
            }
        }

        m_data = data;
    };
    fl eval_fast(sz type_pair_index, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(type_pair_index).eval_fast(r2);
    };
    pr eval_deriv(sz type_pair_index, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        return m_data(type_pair_index).eval_deriv(r2);
    };
    sz index_permissive(sz t1, sz t2) const { return m_data.index_permissive(t1, t2); }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    void widen(fl left, fl right){
        flv rs = calculate_rs();
        VINA_FOR(t1, m_data.dim())
        VINA_RANGE(t2, t1, m_data.dim())
        m_data(t1, t2).widen(rs, left, right);
    };
private:
    flv calculate_rs() const{
        flv tmp(m_n, 0);
        VINA_FOR(i, m_n)
            tmp[i] = std::sqrt(i / m_factor);
        return tmp;
    };

    fl m_cutoff_sqr;
    fl m_max_cutoff_sqr;
    sz m_n;
    fl m_factor;

    triangular_matrix<precalculate_element> m_data;
};

struct precalculate_byatom
{
public:
    precalculate_byatom() : m_cutoff_sqr(0), m_max_cutoff_sqr(0), m_n(0), m_factor(0), m_serial(0), m_dim(0) { }
    precalculate_byatom(const ScoringFunction &sf, const model &model, fl v=max_fl, fl factor=32){
        // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
        m_factor = factor;
        m_cutoff_sqr = sqr(sf.get_cutoff());
        m_max_cutoff_sqr = sqr(sf.get_max_cutoff());
        m_n = sz(m_factor * m_max_cutoff_sqr) + 3; // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
        // std::cout << "-- DEBUG -- sf.cutoff^2 in precalculate = " << m_cutoff_sqr << "\n";
        sz n_atoms = model.num_atoms();
        atomv atoms = model.get_atoms();

        VINA_CHECK(m_factor > epsilon_fl);
        VINA_CHECK(sz(m_max_cutoff_sqr * m_factor) + 1 < m_n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
        VINA_CHECK(m_max_cutoff_sqr * m_factor + 1 < m_n);

        flv rs = calculate_rs();

        // ACCEL: the potentials only read the atom typing (el, ad, xs, sy) and the
        // partial charge of the two atoms. Atoms sharing all of these produce
        // bit-identical tables, so one table is computed and stored per ordered
        // pair of atom classes instead of one per atom pair. This keeps the
        // lookup tables small enough to stay in cache during the search.
        std::vector<sz> atom_class(n_atoms);
        {
            std::map<atom_class_key, sz> classes;
            VINA_FOR(i, n_atoms) {
                // The Vina and Vinardo atom-based potentials, including glue, read only xs.
                // AD4 also reads charge and other atom attributes.
                atom_class_key k = sf.get_atom_typing() == atom_type::XS
                    ? atom_class_key{0, 0, atoms[i].xs, 0, 0} : make_atom_class_key(atoms[i]);
                std::map<atom_class_key, sz>::const_iterator it = classes.find(k);
                if (it == classes.end()) {
                    sz id = classes.size();
                    classes[k] = id;
                    atom_class[i] = id;
                } else
                    atom_class[i] = it->second;
            }
        }

        m_serial = next_serial();
        m_dim = n_atoms;
        m_index.assign(n_atoms * (n_atoms + 1) / 2, 0);
        std::map<std::pair<sz, sz>, unsigned> table_of_classes;
        std::vector<std::pair<sz, sz> > representatives;

        VINA_FOR(i, n_atoms)
        {
            VINA_RANGE(j, i, n_atoms)
            {
                // Preserve the order of the two atoms passed to sf.eval.
                const std::pair<sz, sz> ck(atom_class[i], atom_class[j]);
                auto entry = table_of_classes.lower_bound(ck);
                if (entry == table_of_classes.end() || entry->first != ck) {
                    VINA_CHECK(representatives.size() <= max_unsigned);
                    entry = table_of_classes.emplace_hint(entry, ck, unsigned(representatives.size()));
                    representatives.emplace_back(i, j);
                }
                m_index[triangular_matrix_index(m_dim, i, j)] = entry->second;
            }
        }

        // A single double array gives masked gathers one base address for all tables.
        VINA_CHECK(representatives.size() <= m_smooth.max_size() / m_n / 2);
        m_smooth.resize(2 * m_n * representatives.size());
        m_fast.resize(m_n * representatives.size());
        precalculate_element table(m_n, m_factor);
        VINA_FOR_IN(t, representatives) {
            const std::pair<sz, sz>& atoms_for_table = representatives[t];
            VINA_FOR_IN(k, table.smooth)
                table.smooth[k].first = (std::min)(v, sf.eval(atoms[atoms_for_table.first],
                                                            atoms[atoms_for_table.second], rs[k]));
            table.init_from_smooth_fst(rs);
            store_table(t, table);
        }
    };
    fl eval_fast(sz i, sz j, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        const sz bin = sz(m_factor * r2);
        assert(bin < m_n);
        return m_fast[m_n * table_index(i, j) + bin];
    };
    pr eval_deriv(sz i, sz j, fl r2) const{
        assert(r2 <= m_max_cutoff_sqr);
        const fl r2_factored = m_factor * r2;
        const sz bin = sz(r2_factored);
        const fl rem = r2_factored - bin;
        assert(bin + 1 < m_n);
        const fl* sample = smooth_data() + smooth_offset(i, j) + 2 * bin;
        const fl e = sample[0] + rem * (sample[2] - sample[0]);
        const fl dor = sample[1] + rem * (sample[3] - sample[1]);
        return pr(e, dor);
    };
    const fl* smooth_data() const { return m_smooth.data(); }
    sz smooth_offset(sz i, sz j) const { return 2 * m_n * table_index(i, j); }
    sz index_permissive(sz t1, sz t2) const { return triangular_matrix_index_permissive(m_dim, t1, t2); }
    fl cutoff_sqr() const { return m_cutoff_sqr; }
    fl max_cutoff_sqr() const { return m_max_cutoff_sqr; }
    sz get_factor() const { return m_factor; }
    void widen(fl left, fl right){
        if (m_n == 0) return;
        flv rs = calculate_rs();
        // identical tables stay identical under the same transformation
        precalculate_element table(m_n, m_factor);
        VINA_FOR(t, m_fast.size() / m_n) {
            const fl* samples = m_smooth.data() + 2 * m_n * t;
            VINA_FOR(k, m_n)
                table.smooth[k] = pr(samples[2 * k], samples[2 * k + 1]);
            table.widen(rs, left, right);
            store_table(t, table);
        }
        m_serial = next_serial();
    };
    // Rebuilt or modified tables invalidate the model's compiled offsets.
    sz serial() const { return m_serial; }
    fl factor() const { return m_factor; }
private:
    sz table_index(sz i, sz j) const { return m_index[triangular_matrix_index(m_dim, i, j)]; }

    void store_table(sz t, const precalculate_element& table) {
        fl* samples = m_smooth.data() + 2 * m_n * t;
        VINA_FOR(k, m_n) {
            samples[2 * k] = table.smooth[k].first;
            samples[2 * k + 1] = table.smooth[k].second;
        }
        std::copy(table.fast_samples().begin(), table.fast_samples().end(), m_fast.data() + m_n * t);
    }

    flv calculate_rs() const{
        flv tmp(m_n, 0);
        VINA_FOR(i, m_n)
            tmp[i] = std::sqrt(i / m_factor);
        return tmp;
    };

    fl m_cutoff_sqr;
    fl m_max_cutoff_sqr;
    sz m_n;
    fl m_factor;

    struct atom_class_key {
        sz el, ad, xs, sy;
        unsigned long long charge_bits; // exact bit pattern of the charge
        bool operator<(const atom_class_key& o) const {
            if (el != o.el) return el < o.el;
            if (ad != o.ad) return ad < o.ad;
            if (xs != o.xs) return xs < o.xs;
            if (sy != o.sy) return sy < o.sy;
            return charge_bits < o.charge_bits;
        }
    };
    static atom_class_key make_atom_class_key(const atom& a) {
        atom_class_key k;
        k.el = a.el; k.ad = a.ad; k.xs = a.xs; k.sy = a.sy;
        const fl c = a.charge;
        std::memcpy(&k.charge_bits, &c, sizeof(k.charge_bits));
        return k;
    }

    static sz next_serial() {
        static std::atomic<sz> counter(0);
        return ++counter;
    }
    sz m_serial;
    sz m_dim;                                  // number of atoms
    std::vector<unsigned> m_index;             // triangular atom pairs -> unique table
    flv m_smooth;                              // interleaved energy and derivative samples
    flv m_fast;
};

#endif
