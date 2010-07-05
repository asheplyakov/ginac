/** @file parse_context.h
 *
 *  Interface to parser context. */

/*
 *  GiNaC Copyright (C) 1999-2010 Johannes Gutenberg University Mainz, Germany
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
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef GINAC_PARSE_CONTEXT_H
#define GINAC_PARSE_CONTEXT_H

#include "ex.h"
#include "symbol.h"

#include <cstddef> // for size_t
#include <map>
#include <string>
#include <utility>

namespace GiNaC {

/**
 * Establishes correspondence between the strings and expressions.
 * The parser will create missing symbols (if not instructed otherwise,
 * in which case it fails if the expression contains unknown symbols).
 */
typedef std::map<std::string, ex> symtab;

/**
 * Find the symbol (or abbreviation) with the @a name in the symbol table @a syms.
 *
 * If symbol is missing and @a strict = false, insert it, otherwise
 * throw an exception.
 */
extern ex 
find_or_insert_symbol(const std::string& name, symtab& syms,
	              const bool strict);

/**
 * Function (or class ctor) prototype
 * .first is  the name of function(or ctor),
 * .second is the number of arguments (each of type ex)
 */
typedef std::pair<std::string, std::size_t> prototype;

/**
 * A (C++) function for reading functions and classes from the stream.
 *
 * The parser uses (an associative array of) such functions to construct
 * (GiNaC) classes and functions from a sequence of characters.
 */
typedef ex (*reader_func)(const exvector& args);

/**
 * Prototype table.
 *
 * If parser sees an expression which looks like a function call (e.g.
 * foo(x+y, z^2, t)), it looks up such a table to find out which
 * function (or class) corresponds to the given name and has the given
 * number of the arguments.
 */
typedef std::map<prototype, unsigned> prototype_table;

/**
 * Creates a default prototype table containing all defined GiNaC functions.
 *
 * The data referenced by the return value is only created once when this
 * function is called for the first time. This might cause problems in very
 * rare stituations (i.e. if functions are added after this first call). In
 * that case, a new initialization can be forced with an "true" argument.
 */
extern const prototype_table& get_default_reader(bool force_init = false);

} // namespace GiNaC

#endif // GINAC_PARSE_CONTEXT_H
