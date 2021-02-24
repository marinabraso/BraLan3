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

#include <set>
#include "output_format.h"
#include "../data/taxonomy.h"

using std::set;

void Taxon_format::print_match(const Hsp_context &r, const Metadata &metadata, TextBuffer &out)
{
	const vector<unsigned> taxons((*metadata.taxon_list)[r.orig_subject_id]);
	if (taxons.empty())
		return;
	evalue = std::min(evalue, r.evalue());
	try {
		for (vector<unsigned>::const_iterator i = taxons.begin(); i != taxons.end(); ++i)
			taxid = metadata.taxon_nodes->get_lca(taxid, *i);
	}
	catch (std::exception &) {
		std::cerr << "Query=" << r.query_name << endl << "Subject=" << r.subject_name << endl;
		throw;
	}
}

void Taxon_format::print_query_epilog(TextBuffer &out, const char *query_title, bool unaligned, const Parameters &params) const
{
	out.write_until(query_title, Const::id_delimiters);
	out << '\t' << taxid << '\t';
	if (taxid != 0)
		out.print_e(evalue);
	else
		out << '0';
	out << '\n';
}