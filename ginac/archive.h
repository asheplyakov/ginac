/** @file archive.h
 *
 *  Archiving of GiNaC expressions. */

/*
 *  GiNaC Copyright (C) 1999 Johannes Gutenberg University Mainz, Germany
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

#ifndef __GINAC_ARCHIVE_H__
#define __GINAC_ARCHIVE_H__

#include <ginac/ex.h>

#include <string>
#include <vector>

class ostream;
class istream;

#ifndef NO_GINAC_NAMESPACE
namespace GiNaC {
#endif // ndef NO_GINAC_NAMESPACE

class lst;
class archive;


/** Numerical ID value to refer to an archive_node. */
typedef unsigned int archive_node_id;

/** Numerical ID value to refer to a string. */
typedef unsigned int archive_atom;


/** This class stores all properties needed to record/retrieve the state
 *  of one object of class basic (or a derived class). Each property is
 *  addressed by its name and data type. */
class archive_node
{
	friend ostream &operator<<(ostream &os, const archive_node &ar);
	friend istream &operator>>(istream &is, archive_node &ar);

public:
	archive_node(archive &ar) : a(ar), has_expression(false) {}
	archive_node(archive &ar, const ex &expr);
	~archive_node() {}

	archive_node(const archive_node &other);
	const archive_node &operator=(const archive_node &other);

	bool has_same_ex_as(const archive_node &other) const;

	void add_bool(const string &name, bool value);
	void add_unsigned(const string &name, unsigned int value);
	void add_string(const string &name, const string &value);
	void add_ex(const string &name, const ex &value);

	bool find_bool(const string &name, bool &ret) const;
	bool find_unsigned(const string &name, unsigned int &ret) const;
	bool find_string(const string &name, string &ret) const;
	bool find_ex(const string &name, ex &ret, const lst &sym_lst, unsigned int index = 0) const;

	ex unarchive(const lst &sym_lst) const;

	void forget(void);
	void printraw(ostream &os) const;

private:
	/** Property data types */
	enum property_type {
		PTYPE_BOOL,
		PTYPE_UNSIGNED,
		PTYPE_STRING,
		PTYPE_NODE
	};

	/** Archived property (data type, name and associated data) */
	struct property {
		property() {}
		property(archive_atom n, property_type t, unsigned int v) : type(t), name(n), value(v) {}
		~property() {}

		property(const property &other) : type(other.type), name(other.name), value(other.value) {}
		const property &operator=(const property &other);

		property_type type;	/**< Data type of property. */
		archive_atom name;	/**< Name of property. */
		unsigned int value;	/**< Stored value. */
	};

	/** Reference to the archive to which this node belongs. */
	archive &a;

	/** Vector of stored properties. */
	vector<property> props;

	/** Flag indicating whether a cached unarchived representation of this node exists. */
	mutable bool has_expression;

	/** The cached unarchived representation of this node (if any). */
	mutable ex e;
};


/** This class holds archived versions of GiNaC expressions (class ex).
 *  An archive can be constructed from an expression and then written to
 *  a stream; or it can be read from a stream and then unarchived, yielding
 *  back the expression. Archives can hold multiple expressions which can
 *  be referred to by name or index number. The main component of the
 *  archive class is a vector of archive_nodes which each store one object
 *  of class basic (or a derived class). */
class archive
{
	friend ostream &operator<<(ostream &os, const archive &ar);
	friend istream &operator>>(istream &is, archive &ar);

public:
	archive() {}
	~archive() {}

	/** Construct archive from expression using the default name "ex". */
	archive(const ex &e) {archive_ex(e, "ex");}

	/** Construct archive from expression using the specified name. */
	archive(const ex &e, const char *n) {archive_ex(e, n);}

	archive_node_id add_node(const archive_node &n);
	archive_node &get_node(archive_node_id id);

	void archive_ex(const ex &e, const char *name);
	ex unarchive_ex(const lst &sym_lst, const char *name) const;
	ex unarchive_ex(const lst &sym_lst, unsigned int index = 0) const;
	ex unarchive_ex(const lst &sym_lst, string &name, unsigned int index = 0) const;
	unsigned int num_expressions(void) const;

	void clear(void);

	void forget(void);
	void printraw(ostream &os) const;

private:
	/** Vector of archived nodes. */
	vector<archive_node> nodes;

	/** Archived expression descriptor. */
	struct archived_ex {
		archived_ex() {}
		archived_ex(archive_atom n, archive_node_id node) : name(n), root(node) {}

		archive_atom name;		/**< Name of expression. */
		archive_node_id root;	/**< ID of root node. */
	};

	/** Vector of archived expression descriptors. */
	vector<archived_ex> exprs;

public:
	archive_atom atomize(const string &s) const;
	const string &unatomize(archive_atom id) const;

private:
	/** Vector of atomized strings (using a vector allows faster unarchiving). */
	mutable vector<string> atoms;
};


ostream &operator<<(ostream &os, const archive &ar);
istream &operator>>(istream &is, archive &ar);


#ifndef NO_GINAC_NAMESPACE
} // namespace GiNaC
#endif // ndef NO_GINAC_NAMESPACE

#endif // ndef __GINAC_ARCHIVE_H__