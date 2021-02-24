/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef SHAPE_CONFIG_H_
#define SHAPE_CONFIG_H_

#include <vector>
#include <algorithm>
#include <stdint.h>
#include <string>
#include "shape.h"

extern std::vector<std::vector<std::string>> shape_codes;

class shape_config
{

public:
	
	shape_config():
		n_ (0),
		mode_ (0)
	{ }

	shape_config(unsigned mode, unsigned count, const vector<string> &shape_mask):
		n_ (0),
		mode_ (mode)
	{
		if (shape_mask.size() == 0) {
			unsigned maxShapes = count == 0 ? (unsigned)shape_codes[mode_].size() : count;
			for (unsigned i = 0; i < maxShapes; ++i)
				shapes_[n_++] = Shape(shape_codes[mode_][i].c_str(), i);
		}
		else {
			for (unsigned i = 0; i < (count == 0 ? shape_mask.size() : std::min((unsigned)shape_mask.size(), count)); ++i)
				shapes_[n_++] = Shape(shape_mask[i].c_str(), i);
		}
	}

	unsigned count() const
	{ return n_; }

	const Shape& operator[](size_t i) const
	{ return shapes_[i]; }

	unsigned mode() const
	{ return mode_; }

	friend std::ostream& operator<<(std::ostream&s, const shape_config &cfg)
	{
		for (unsigned i = 0; i < cfg.n_; ++i)
			s << cfg.shapes_[i] << (i<cfg.n_-1?",":"");
		return s;
	}

	std::vector<uint32_t> patterns(unsigned begin, unsigned end) const {
		std::vector<uint32_t> v;
		for (unsigned i = begin; i < end; ++i)
			v.push_back(shapes_[i].mask_);
		return v;
	}

private:

	Shape shapes_[Const::max_shapes];
	unsigned n_, mode_;

};

extern shape_config shapes;
extern unsigned shape_from, shape_to;

#endif /* SHAPE_CONFIG_H_ */
