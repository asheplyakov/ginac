/** @file tensor.h
 *
 *  Interface to GiNaC's special tensors. */

/*
 *  GiNaC Copyright (C) 1999-2001 Johannes Gutenberg University Mainz, Germany
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

#ifndef __GINAC_TENSOR_H__
#define __GINAC_TENSOR_H__

#include "ex.h"

namespace GiNaC {


/** This class holds one of GiNaC's predefined special tensors such as the
 *  delta and the metric tensors. They are represented without indices.
 *  To attach indices to them, wrap them in an object of class indexed. */
class tensor : public basic
{
	GINAC_DECLARE_REGISTERED_CLASS(tensor, basic)

	// other constructors
protected:
	tensor(unsigned ti);

	// functions overriding virtual functions from bases classes
protected:
	unsigned return_type(void) const { return return_types::noncommutative_composite; }
};


/** This class represents the delta tensor. If indexed, it must have exactly
 *  two indices of the same type. */
class tensdelta : public tensor
{
	GINAC_DECLARE_REGISTERED_CLASS(tensdelta, tensor)

	// functions overriding virtual functions from bases classes
public:
	void print(const print_context & c, unsigned level = 0) const;
	ex eval_indexed(const basic & i) const;
	bool contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const;
};


/** This class represents a general metric tensor which can be used to
 *  raise/lower indices. If indexed, it must have exactly two indices of the
 *  same type which must be of class varidx or a subclass. */
class tensmetric : public tensor
{
	GINAC_DECLARE_REGISTERED_CLASS(tensmetric, tensor)

	// functions overriding virtual functions from bases classes
public:
	void print(const print_context & c, unsigned level = 0) const;
	ex eval_indexed(const basic & i) const;
	bool contract_with(exvector::iterator self, exvector::iterator other, exvector & v) const;
};


/** This class represents a Minkowski metric tensor. It has all the
 *  properties of a metric tensor and is (as a matrix) equal to
 *  diag(1,-1,-1,...) or diag(-1,1,1,...). */
class minkmetric : public tensmetric
{
	GINAC_DECLARE_REGISTERED_CLASS(minkmetric, tensmetric)

	// other constructors
public:
	/** Construct Lorentz metric tensor with given signature. */
	minkmetric(bool pos_sig);

	// functions overriding virtual functions from bases classes
public:
	void print(const print_context & c, unsigned level = 0) const;
	ex eval_indexed(const basic & i) const;

	// member variables
private:
	bool pos_sig; /**< If true, the metric is diag(-1,1,1...). Otherwise it is diag(1,-1,-1,...). */
};


/** This class represents the totally antisymmetric epsilon tensor. If
 *  indexed, all indices must be of the same type and their number must
 *  be equal to the dimension of the index space. */
class tensepsilon : public tensor
{
	GINAC_DECLARE_REGISTERED_CLASS(tensepsilon, tensor)

	// other constructors
public:
	tensepsilon(bool minkowski, bool pos_sig);

	// functions overriding virtual functions from bases classes
public:
	void print(const print_context & c, unsigned level = 0) const;
	ex eval_indexed(const basic & i) const;

	// member variables
private:
	bool minkowski; /**< If true, tensor is in Minkowski-type space. Otherwise it is in a Euclidean space. */
	bool pos_sig;  /**< If true, the metric is assumed to be diag(-1,1,1...). Otherwise it is diag(1,-1,-1,...). This is only relevant if minkowski = true. */
};


// utility functions
inline const tensor &ex_to_tensor(const ex &e)
{
	return static_cast<const tensor &>(*e.bp);
}

/** Create a delta tensor with specified indices. The indices must be of class
 *  idx or a subclass. The delta tensor is always symmetric and its trace is
 *  the dimension of the index space.
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @return newly constructed delta tensor */
ex delta_tensor(const ex & i1, const ex & i2);

/** Create a symmetric metric tensor with specified indices. The indices
 *  must be of class varidx or a subclass. A metric tensor with one
 *  covariant and one contravariant index is equivalent to the delta tensor.
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @return newly constructed metric tensor */
ex metric_tensor(const ex & i1, const ex & i2);

/** Create a Minkowski metric tensor with specified indices. The indices
 *  must be of class varidx or a subclass. The Lorentz metric is a symmetric
 *  tensor with a matrix representation of diag(1,-1,-1,...) (negative
 *  signature, the default) or diag(-1,1,1,...) (positive signature).
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @param pos_sig Whether the signature is positive
 *  @return newly constructed Lorentz metric tensor */
ex lorentz_g(const ex & i1, const ex & i2, bool pos_sig = false);

/** Create an epsilon tensor in a Euclidean space with two indices. The
 *  indices must be of class idx or a subclass, and have a dimension of 2.
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @return newly constructed epsilon tensor */
ex epsilon_tensor(const ex & i1, const ex & i2);

/** Create an epsilon tensor in a Euclidean space with three indices. The
 *  indices must be of class idx or a subclass, and have a dimension of 3.
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @param i3 Third index
 *  @return newly constructed epsilon tensor */
ex epsilon_tensor(const ex & i1, const ex & i2, const ex & i3);

/** Create an epsilon tensor in a Minkowski space with four indices. The
 *  indices must be of class varidx or a subclass, and have a dimension of 4.
 *
 *  @param i1 First index
 *  @param i2 Second index
 *  @param i3 Third index
 *  @param i4 Fourth index
 *  @param pos_sig Whether the signature of the metric is positive
 *  @return newly constructed epsilon tensor */
ex lorentz_eps(const ex & i1, const ex & i2, const ex & i3, const ex & i4, bool pos_sig = false);

} // namespace GiNaC

#endif // ndef __GINAC_TENSOR_H__