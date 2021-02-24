/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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

#include <string.h>
#include <assert.h>
#include "value.h"
#include "../util/util.h"

const Letter Char_representation::invalid = '\xff';

invalid_sequence_char_exception::invalid_sequence_char_exception(char ch) :
	msg(std::string("Invalid character (") + print_char(ch) + ") in sequence")
{ }

Char_representation::Char_representation(unsigned size, const char* chars, char mask, const char* mask_chars)
{
	memset(data_, invalid, sizeof(data_));
	for (unsigned i = 0; i < size; ++i) {
		assert(chars[i] != (char)invalid);
		data_[(long)chars[i]] = i;
		data_[(long)tolower(chars[i])] = i;
	}
	while (*mask_chars != 0) {
		const char ch = *mask_chars;
		data_[(long)ch] = mask;
		data_[(long)tolower(ch)] = mask;
		++mask_chars;
	}
}

Value_traits::Value_traits(const char* alphabet, Letter mask_char, const char* ignore, Sequence_type seq_type) :
	alphabet(alphabet),
	alphabet_size((unsigned)strlen(alphabet)),
	mask_char(mask_char),
	from_char(Char_representation((unsigned)alphabet_size, alphabet, mask_char, ignore)),
	seq_type(seq_type)
{}

const Value_traits amino_acid_traits(AMINO_ACID_ALPHABET, 23, "UO-", Sequence_type::amino_acid);
const Value_traits nucleotide_traits("ACGTN", 4, "MRWSYKVHDBX", Sequence_type::nucleotide);
Value_traits value_traits(amino_acid_traits);
Value_traits input_value_traits(amino_acid_traits);