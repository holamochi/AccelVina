/* SPDX-License-Identifier: Apache-2.0 */
#ifndef VINA_PAIRS_ACCEL_H
#define VINA_PAIRS_ACCEL_H

#include "common.h"
#include "c_kernel.h"

struct precalculate_byatom;

struct compiled_pairs {
	bool usable = false;
	std::vector<vina_ck_pair_row> rows;
	std::vector<vina_ck_pair_chunk> chunks;
};

struct model_accel {
	bool ready = false;
	const precalculate_byatom* p = NULL;
	sz p_serial = 0;
	std::vector<compiled_pairs> ligands;
	compiled_pairs inter, other, glue;
	std::vector<fl> xs, ys, zs;
	std::vector<unsigned char> grid_type;
	bool grid_type_ready = false;
	enum { grid_type_none = 255 };
};

bool accel_have_avx();
bool accel_have_c_kernel();

#endif
