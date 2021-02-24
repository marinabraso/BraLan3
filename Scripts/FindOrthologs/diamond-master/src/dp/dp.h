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
#include <list>
#include <vector>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../stats/hauser_correction.h"
#include "../basic/statistics.h"
#include "../basic/config.h"
#include "../util/dynamic_iterator.h"
#include "../stats/cbs.h"

int smith_waterman(const Sequence&query, const Sequence&subject, unsigned band, unsigned padding, int op, int ep);

struct Local {};
struct Global {};

template<typename _t>
struct Fixed_score_buffer
{

	inline void init(size_t col_size, size_t cols, _t init)
	{
		col_size_ = col_size;
		data_.clear();
		data_.reserve(col_size * cols);
		data_.resize(col_size);
		for (size_t i = 0; i<col_size; ++i)
			data_[i] = init;
	}
	
	std::pair<int, int> find(_t s) const
	{
		const int i = int(std::find(data_.begin(), data_.end(), s) - data_.begin());
		return std::pair<int, int>(int(i%col_size_), int(i / col_size_));
	}

	inline std::pair<_t*, _t*> get()
	{
		data_.resize(data_.size() + col_size_);
		_t* ptr = last();
		return std::pair<_t*, _t*>(ptr - col_size_, ptr);
	}

	inline _t* last()
	{
		return &*(data_.end() - col_size_);
	}

	const _t* column(int col) const
	{
		return &data_[col_size_*col];
	}

	_t operator()(int i, int j) const
	{
		return data_[j*col_size_ + i];
	}

	friend std::ostream& operator<<(std::ostream &s, const Fixed_score_buffer &buf)
	{
		s << '\t';
		for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
			s << j << '\t';
		s << std::endl;
		for (int i = 0; i < int(buf.col_size_); ++i) {
			s << i << '\t';
			for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
				s << buf(i, j) << '\t';
			s << std::endl;
		}
		return s;
	}

private:
	std::vector<_t> data_;
	size_t col_size_;

};

template<typename _score, typename _mode>
const Fixed_score_buffer<_score>& needleman_wunsch(Sequence query, Sequence subject, int &max_score, const _mode&, const _score&);

struct Band
{
	void init(int diags, int cols)
	{
		diags_ = diags;
		cols_ = cols;
		data_.clear();
		data_.resize((size_t)diags*cols);
	}
	struct Iterator {
		Iterator(uint8_t *p, int diags) :
			diags_(diags),
			p_(p)			
		{}
		uint8_t& operator[](int i)
		{
			return p_[i*diags_];
		}
	private:
		const int diags_;
		uint8_t *p_;
	};
	Iterator diag(int o)
	{
		return Iterator(&data_[o], diags_);
	}
	int cols() const
	{
		return cols_;
	}
	int diags() const
	{
		return diags_;
	}
	uint8_t* data()
	{
		return data_.data();
	}
	bool check(uint8_t *ptr) const
	{
		return ptr >= data_.data() && ptr <= data_.data() + data_.size();
	}
private:
	int diags_, cols_;
	vector<uint8_t> data_;
};

extern size_t cells;

struct DpTarget
{
	DpTarget():
		target_idx(-1),
		matrix(nullptr),
		previous_i1(0),
		previous_j1(0)
	{}
	DpTarget(const Sequence &seq, int d_begin, int d_end, int j_begin, int j_end, int target_idx = 0, int qlen = 0, const Stats::TargetMatrix* matrix = nullptr) :
		seq(seq),
		d_begin(d_begin),
		d_end(d_end),
		j_begin(j_begin),
		j_end(j_end),
		target_idx(target_idx),
		previous_i1(0),
		previous_j1(0),
		matrix(matrix)
	{
		int pos = std::max(d_end - 1, 0) - (d_end - 1);
		const int d0 = d_begin;
		const int j1 = std::min(qlen - 1 - d0, (int)(seq.length() - 1)) + 1;
		cols = j1 - pos;
	}
	DpTarget(const Sequence& seq, int target_idx, int previous_i1 = 0, int previous_j1 = 0):
		seq(seq),
		target_idx(target_idx),
		previous_i1(previous_i1),
		previous_j1(previous_j1)
	{}
	int left_i1() const
	{
		return std::max(d_end - 1, 0);
	}
	int band() const {
		return d_end - d_begin;
	}
	bool operator<(const DpTarget &x) const
	{
		const int i = left_i1(), j = x.left_i1(), b1 = band(), b2 = x.band(), bin_b1 = b1 / config.band_bin, bin_b2 = b2 / config.band_bin,
			t1 = cols, t2 = x.cols, bin_t1 = t1 / config.col_bin, bin_t2 = t2 / config.col_bin;
		return bin_b1 < bin_b2 || (bin_b1 == bin_b2 && (bin_t1 < bin_t2 || (bin_t1 == bin_t2 && i < j)));
		//return i < j || (i == j && (target_idx < x.target_idx || (target_idx == x.target_idx && d_begin < x.d_begin)));
	}
	bool blank() const {
		return target_idx == -1;
	}
	bool adjusted_matrix() const {
		return matrix != nullptr;
	}
	int matrix_scale() const {
		return adjusted_matrix() ? config.cbs_matrix_scale : 1;
	}
	Sequence seq;
	int d_begin, d_end, j_begin, j_end, target_idx, cols, previous_i1, previous_j1;
	const Stats::TargetMatrix* matrix;
};

struct DpStat
{
	DpStat():
		gross_cells(0),
		net_cells(0)
	{}
	DpStat& operator+=(DpStat &x)
	{
		mtx_.lock();
		gross_cells += x.gross_cells;
		net_cells += x.net_cells;
		mtx_.unlock();
		return *this;
	}
	size_t gross_cells, net_cells;
private:
	std::mutex mtx_;
};

extern DpStat dp_stat;

void smith_waterman(Sequence q, Sequence s, Hsp &out);
int score_range(Sequence query, Sequence subject, int i, int j, int j_end);

namespace DP {

struct Traceback {};
struct StatTraceback {};
struct VectorTraceback {};
struct ScoreOnly {};
struct ScoreWithCoords {};

enum { TRACEBACK = 1, PARALLEL = 2, FULL_MATRIX = 4, WITH_COORDINATES = 8 };

struct NoCBS {
	constexpr void* operator[](int i) const { return nullptr; }
};
	
namespace Swipe {

//DECL_DISPATCH(std::list<Hsp>, swipe, (const sequence &query, const sequence *subject_begin, const sequence *subject_end, int score_cutoff))

}

namespace BandedSwipe {

DECL_DISPATCH(std::list<Hsp>, swipe, (const Sequence&query, std::vector<DpTarget> &targets8, const std::vector<DpTarget> &targets16, const std::vector<DpTarget>& targets32, DynamicIterator<DpTarget>* targets, Frame frame, const Bias_correction *composition_bias, int flags, Statistics &stat))

}

}

void banded_sw(const Sequence&query, const Sequence&subject, int d_begin, int d_end, int j_begin, int j_end, Hsp &out);

void anchored_3frame_dp(const TranslatedSequence &query, Sequence&subject, const DiagonalSegment &anchor, Hsp &out, int gap_open, int gap_extend, int frame_shift);
int sw_3frame(const TranslatedSequence &query, Strand strand, const Sequence&subject, int gap_open, int gap_extend, int frame_shift, Hsp &out);

DECL_DISPATCH(std::list<Hsp>, banded_3frame_swipe, (const TranslatedSequence &query, Strand strand, vector<DpTarget>::iterator target_begin, vector<DpTarget>::iterator target_end, DpStat &stat, bool score_only, bool parallel))
