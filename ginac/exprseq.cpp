/** @file exprseq.cpp
 *
 *  Implementation of GiNaC's exprseq. */

/*
 *  GiNaC Copyright (C) 1999-2005 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "exprseq.h"

namespace GiNaC {

template <> GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(exprseq, basic,
  print_func<print_context>(&exprseq::do_print).
  print_func<print_tree>(&exprseq::do_print_tree))

/** Specialization of container::info() for exprseq. */
bool exprseq::info(unsigned inf) const
{
	if (inf == info_flags::exprseq)
		return true;
	else
		return inherited::info(inf);
}

} // namespace GiNaC
