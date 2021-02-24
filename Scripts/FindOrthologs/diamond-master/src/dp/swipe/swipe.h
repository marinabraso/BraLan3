/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include <assert.h>
#include <limits>
#include <vector>
#include "../score_vector.h"
#include "../score_vector_int8.h"
#include "../score_vector_int16.h"
#include "../../basic/value.h"
#include "../../util/simd/vector.h"
#include "../../util/system.h"
#include "../../util/memory/alignment.h"
#include "../../util/simd/transpose.h"

static inline uint8_t cmp_mask(int x, int y) {
	return x == y;
}

static inline int blend(int v, int w, int mask) {
	return mask ? w : v;
}

template<typename _sv>
struct TraceStat {
	_sv length;
	_sv gapopen;
	_sv qstart;
	_sv sstart;
	_sv ident;
	_sv mismatch;
};

struct DummyRowCounter {
	DummyRowCounter() {}
	DummyRowCounter(int) {}
	void store(void*) {}
	template<typename _sv>
	FORCE_INLINE void inc(const _sv&, const _sv&) {}
	enum { MAX_LEN = INT_MAX };
};

template<typename _sv>
struct RowCounter {
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
	RowCounter(int i):
		i(::DISPATCH_ARCH::ScoreTraits<_sv>::zero_score() + Score(i)),
		i_max()
	{
	}
	FORCE_INLINE void inc(const _sv& best, const _sv& current_cell) {
		i_max = blend(i_max, i, best == current_cell);
		i += _sv(Score(1));
	}
	void store(Score *ptr) {
		store_sv(i_max, ptr);
	}
	static constexpr int MAX_LEN = ::DISPATCH_ARCH::ScoreTraits<_sv>::max_int_score();
	_sv i;
	_sv i_max;
};

template<typename _sv>
MSC_INLINE _sv add_cbs(const _sv &v, void*) {
	return v;
}

template<typename _sv>
MSC_INLINE _sv add_cbs(const _sv& v, const _sv& query_bias) {
	return v + query_bias;
}

template<typename _score>
_score add_cbs_scalar(_score x, int8_t b) {
	return x + _score(b);
}

template<typename _score>
_score add_cbs_scalar(_score x, void *b) {
	return x;
}

template<typename _sv>
MSC_INLINE void make_gap_mask(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask *trace_mask, const _sv& current_cell, const _sv& vertical_gap, const _sv& horizontal_gap) {
	trace_mask->gap = ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask::make(cmp_mask(current_cell, vertical_gap), cmp_mask(current_cell, horizontal_gap));
}

template<typename _sv>
MSC_INLINE void make_gap_mask(std::nullptr_t, const _sv&, const _sv&, const _sv&) {
}

template<typename _sv>
MSC_INLINE void make_open_mask(typename ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask *trace_mask, const _sv& open, const _sv& vertical_gap, const _sv& horizontal_gap) {
	trace_mask->open = ::DISPATCH_ARCH::ScoreTraits<_sv>::TraceMask::make(cmp_mask(vertical_gap, open), cmp_mask(horizontal_gap, open));
}

template<typename _sv>
MSC_INLINE void make_open_mask(std::nullptr_t, const _sv&, const _sv&, const _sv&) {
}

namespace DP {
	struct NoCBS;
}

template<typename _sv, typename _cbs>
struct CBSBuffer {
	CBSBuffer(const DP::NoCBS&, int, uint32_t) {}
	void* operator()(int i) const {
		return nullptr;
	}
};

template<typename _sv>
struct CBSBuffer<_sv, const int8_t*> {
	CBSBuffer(const int8_t* v, int l, uint32_t channel_mask) {
		typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
		data.reserve(l);
		for (int i = 0; i < l; ++i)
			data.push_back(load_sv(Score(v[i]), (Score)0, channel_mask));
	}
	_sv operator()(int i) const {
		return data[i];
	}
	std::vector<_sv, Util::Memory::AlignmentAllocator<_sv, 32>> data;
};


template<typename _sv, typename _cbs, typename _trace_mask, typename _row_counter>
MSC_INLINE _sv swipe_cell_update(const _sv &diagonal_cell,
	const _sv &scores,
	_cbs query_bias,
	const _sv &gap_extension,
	const _sv &gap_open,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best,
	void*,
	void*,
	void*,
	_trace_mask trace_mask,
	_row_counter& row_counter)
{
	using std::max;

	_sv current_cell = diagonal_cell + add_cbs(scores, query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	make_gap_mask(trace_mask, current_cell, vertical_gap, horizontal_gap);

	best = max(best, current_cell);

	row_counter.inc(best, current_cell);

	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	make_open_mask(trace_mask, open, vertical_gap, horizontal_gap);

	return current_cell;
}

/*template<typename _sv>
static inline _sv swipe_cell_update(const _sv& diagonal_cell,
	const _sv& scores,
	const _sv& query_bias,
	const _sv& gap_extension,
	const _sv& gap_open,
	_sv& horizontal_gap,
	_sv& vertical_gap,
	_sv& best,
	TraceStat<_sv> &trace_stat_diag,
	TraceStat<_sv> &trace_stat_vertical,
	TraceStat<_sv> &trace_stat_horizontal,
	void*,
	const RowCounter<_sv>&)
{
	typedef typename ::DISPATCH_ARCH::ScoreTraits<_sv>::Score Score;
	using std::max;
	_sv current_cell = diagonal_cell + (scores + query_bias);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);

	const _sv one = _sv(Score(1)), zero = _sv(), zero2 = _sv(Score(0));
	const _sv vgap_mask = current_cell == vertical_gap, hgap_mask = current_cell == horizontal_gap, zero_mask = current_cell == zero;

	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);

	trace_stat_vertical.length += one;
	trace_stat_horizontal.length += one;
	trace_stat_diag.length += one;
	trace_stat_diag.length = blend(trace_stat_diag.length, trace_stat_vertical.length, vgap_mask);
	trace_stat_diag.length = blend(trace_stat_diag.length, trace_stat_horizontal.length, hgap_mask);
	trace_stat_diag.length = blend(trace_stat_diag.length, zero2, zero_mask);
	
	return current_cell;
}*/

template<typename _sv>
static inline _sv cell_update(const _sv &diagonal_cell,
	const _sv &shift_cell0,
	const _sv &shift_cell1,
	const _sv &scores,
	const _sv &gap_extension,
	const _sv &gap_open,
	const _sv &frame_shift,
	_sv &horizontal_gap,
	_sv &vertical_gap,
	_sv &best)
{
	using std::max;
	_sv current_cell = diagonal_cell + scores;
	const _sv f = scores - frame_shift;
	current_cell = max(current_cell, shift_cell0 + f);
	current_cell = max(current_cell, shift_cell1 + f);
	current_cell = max(max(current_cell, vertical_gap), horizontal_gap);
	::DISPATCH_ARCH::ScoreTraits<_sv>::saturate(current_cell);
	best = max(best, current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const _sv open = current_cell - gap_open;
	vertical_gap = max(vertical_gap, open);
	horizontal_gap = max(horizontal_gap, open);
	return current_cell;
}

namespace DISPATCH_ARCH {

template<typename _sv>
struct SwipeProfile
{

#ifdef __SSSE3__
	inline void set(typename ScoreTraits<_sv>::Vector seq)
	{
		assert(sizeof(data_) / sizeof(_sv) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = _sv(j, seq);
	}
#else
	inline void set(uint64_t seq)
	{
		assert(sizeof(data_) / sizeof(_sv) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < AMINO_ACID_COUNT; ++j)
			data_[j] = _sv(j, seq);
	}
#endif

	inline const _sv& get(Letter i) const
	{
		return data_[(int)i];
	}

	void set(const int8_t** target_scores) {
#if ARCH_ID == 2
		transpose(target_scores, 32, (int8_t*)data_, __m256i());
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i].expand_from_8bit();
#elif defined(__SSE2__)
		transpose(target_scores, 16, (int8_t*)data_, __m128i());
		for (int i = 0; i < 16; ++i)
			target_scores[i] += 16;
		transpose(target_scores, 16, (int8_t*)(&data_[16]), __m128i());
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i].expand_from_8bit();
#else
		for (int i = 0; i < AMINO_ACID_COUNT; ++i)
			data_[i] = target_scores[0][i];
#endif
	}

	void set(const int32_t** target_scores) {
		typename ScoreTraits<_sv>::Score s[ScoreTraits<_sv>::CHANNELS];
		for (size_t i = 0; i < AMINO_ACID_COUNT; ++i) {
			for (size_t j = 0; j < ScoreTraits<_sv>::CHANNELS; ++j)
				s[j] = target_scores[j][i];
			data_[i] = load_sv(s);
		}
	}

	//_sv data_[AMINO_ACID_COUNT];
	_sv data_[32];

};

template<>
struct SwipeProfile<int32_t>
{
#ifdef __AVX2__
	void set(const __m256i& seq)
	{
		int16_t s[32];
		_mm256_storeu_si256((__m256i*)s, seq);
		const int* row = score_matrix.row((char)s[0]);
		std::copy(row, row + 32, this->row);
	}
#endif
#ifdef __SSE2__
	void set(const __m128i& seq)
	{
		int16_t s[8];
		_mm_storeu_si128((__m128i*)s, seq);
		const int* row = score_matrix.row((char)s[0]);
		std::copy(row, row + 32, this->row);
	}
#endif
	void set(uint64_t seq)
	{
		const int* row = score_matrix.row((char)seq);
		std::copy(row, row + 32, this->row);
	}
	void set(const int8_t** target_scores) {
		for (int i = 0; i < 32; ++i)
			row[i] = target_scores[0][i];
	}
	void set(const int32_t * *target_scores) {
		for (int i = 0; i < 32; ++i)
			row[i] = target_scores[0][i];
	}
	int32_t get(char i) const
	{
		return row[(int)i];
	}
	int32_t row[32];
};

}
