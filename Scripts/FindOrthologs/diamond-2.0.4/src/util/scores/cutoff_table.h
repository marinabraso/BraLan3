/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once
#include "../../basic/score_matrix.h"
#include "../intrin.h"

namespace Util { namespace Scores {

struct CutoffTable {
	
	CutoffTable(double evalue) {
		for (int b = 1; b <= MAX_BITS; ++b) {
			data_[b] = score_matrix.rawscore(score_matrix.bitscore_norm(evalue, 1 << (b - 1)));
		}
	}

	int operator()(int query_len) const {
		const int b = 32 - clz((uint32_t)query_len);
		return data_[b];
	}

private:

	enum { MAX_BITS = 31 };

	int data_[MAX_BITS + 1];

};

}}