/** @file function.cpp
 *
 *  Implementation of class of symbolic functions. */

/*
 *  GiNaC Copyright (C) 1999-2021 Johannes Gutenberg University Mainz, Germany
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

#include "function.h"
#include "operators.h"
#include "fderivative.h"
#include "ex.h"
#include "lst.h"
#include "symmetry.h"
#include "print.h"
#include "power.h"
#include "archive.h"
#include "inifcns.h"
#include "utils.h"
#include "hash_seed.h"
#include "remember.h"
#include "relational.h"

#include <iostream>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>

namespace GiNaC {
template<typename F, typename V>
static inline auto visit(F&& f, V&& v)
#if __cplusplus < 201403L
-> decltype(::mapbox::util::apply_visitor(std::forward<F>(f), std::forward<V>(v)))
#endif
{
	return ::mapbox::util::apply_visitor(std::forward<F>(f), std::forward<V>(v));
}

template<int N, typename A>
static inline typename A::const_reference get(const A& array) {
	return array[N];
}

template<typename F, typename... Rest> struct adjusted_arity {
	static constexpr unsigned value = get_arity_s<F>::value - sizeof...(Rest);
};

template<typename F, typename... Rest> struct genseq2 : genseq<adjusted_arity<F, Rest...>::value> { };
template<typename F> struct return_type;
template<typename R, typename... Args> struct return_type<R(Args...)> {
	using type = R;
};
template<typename R, typename... Args> struct return_type<R(*)(Args...)> {
	using type = R;
};

/* Unpacks arguments and calls a function
 * f, {1, 2} -> f(1, 2)
 * g, {1, 2, 3}, 'a', "b" -> g(1, 2, 3, 'a', "b")
 */
template<bool PassExvector> struct exvector_args_dispatcher {
	template<typename F, typename C, int... S, typename... Rest>
	static inline auto do_call(const F& func, const C& args, seq<S...>, Rest&&... rest)
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
	{
		return func(get<S>(args)..., std::forward<Rest>(rest)...);
	}

	template<typename F, typename ArgT, typename... Rest>
	static inline auto call(F fptr, const ArgT& args, Rest&&... rest)
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
	{
		return do_call(fptr, args, typename genseq2<F, Rest...>::type{}, std::forward<Rest>(rest)...);
	}
};

/* Calls function taking exvector arguments
 * f, {1, 2} -> f({1, 2})
 * g, {1, 2, 3}, 'a', "b" -> g({1, 2, 3}, 'a', "b")
 */
template<> struct exvector_args_dispatcher<true> {
	/** Calls a function taking const exvector&, const Rest&... */
	template<typename F, typename ArgT, typename... Rest>
	static inline auto call(F fptr, const ArgT& args, Rest&&... rest)
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
	{
		return fptr(args, std::forward<Rest>(rest)...);
	}
};

template<typename F, typename Args, typename... Rest>
static inline auto funcp_call(F fptr, const Args& args, Rest&&... rest)
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
{
	return exvector_args_dispatcher<accepts_exvector_args_s<F>::value>::call(fptr, args, std::forward<Rest>(rest)...);
}

template<typename... Rest>
struct eval_visitor {
	const exvector& args;
	std::tuple<Rest...> rest;

	explicit eval_visitor(const exvector& args_, Rest&&... rest_) : args(args_),
		rest{std::forward<Rest>(rest_)...} { }

	template<typename F, int... S>
	inline auto do_call(F fptr, seq<S...>) 
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
	{
		return funcp_call(fptr, this->args, std::get<S>(rest)...);
	}

	template<typename F>
	inline auto operator()(F fptr)
#if __cplusplus < 201402L
	-> typename return_type<F>::type
#endif
	{
		return this->do_call(fptr, typename genseq<sizeof...(Rest)>::type());
	}
};

template<typename Variant, typename... Rest>
static inline auto call_funcp(const Variant& f, const exvector& args, Rest&&... rest)
#if __cplusplus < 201402L
-> decltype(visit(eval_visitor<Rest...>{args, std::forward<Rest>(rest)...}, f))
#endif
{
	return visit(eval_visitor<Rest...>{args, std::forward<Rest>(rest)...}, f);
}

struct is_null_visitor {
	template<typename F>
	inline bool operator()(F fptr) {
		return fptr == nullptr;
	}
};

template<typename Variant>
inline bool is_null(const Variant& v) {
	return visit(is_null_visitor{}, v);
}

//////////
// helper class function_options
//////////

function_options::function_options()
{
	initialize();
}

function_options::function_options(std::string const & n, std::string const & tn)
{
	initialize();
	set_name(n, tn);
}

function_options::function_options(std::string const & n, unsigned np)
{
	initialize();
	set_name(n, std::string());
	nparams = np;
}

function_options::~function_options()
{
	// nothing to clean up at the moment
}

void function_options::initialize()
{
	set_name("unnamed_function", "\\mbox{unnamed}");
	nparams = 0;
	evalf_params_first = true;
	use_return_type = false;
	use_remember = false;
	functions_with_same_name = 1;
	symtree = 0;
}

function_options & function_options::set_name(std::string const & n,
                                              std::string const & tn)
{
	name = n;
	if (tn==std::string())
		TeX_name = "\\mbox{"+name+"}";
	else
		TeX_name = tn;
	return *this;
}

function_options & function_options::latex_name(std::string const & tn)
{
	TeX_name = tn;
	return *this;
}

function_options & function_options::set_return_type(unsigned rt, const return_type_t* rtt)
{
	use_return_type = true;
	return_type = rt;
	if (rtt != nullptr)
		return_type_tinfo = *rtt;
	else
		return_type_tinfo = make_return_type_t<function>();
	return *this;
}

function_options & function_options::do_not_evalf_params()
{
	evalf_params_first = false;
	return *this;
}

function_options & function_options::remember(unsigned size,
                                              unsigned assoc_size,
                                              unsigned strategy)
{
	use_remember = true;
	remember_size = size;
	remember_assoc_size = assoc_size;
	remember_strategy = strategy;
	return *this;
}

function_options & function_options::overloaded(unsigned o)
{
	functions_with_same_name = o;
	return *this;
}

function_options & function_options::set_symmetry(const symmetry & s)
{
	symtree = s;
	return *this;
}
	
void function_options::test_and_set_nparams(unsigned n)
{
	if (nparams==0) {
		nparams = n;
	} else if (nparams!=n) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// be caught anyhow
		std::cerr << "WARNING: " << name << "(): number of parameters ("
		          << n << ") differs from number set before (" 
		          << nparams << ")" << std::endl;
	}
}

void function_options::set_print_func(unsigned id, print_funcp_variant&& f)
{
	if (id >= print_dispatch_table.size())
		print_dispatch_table.resize(id + 1);
	print_dispatch_table[id] = f;
}


GINAC_IMPLEMENT_REGISTERED_CLASS(function, exprseq)

//////////
// default constructor
//////////

// public

function::function() : serial(0)
{
}

//////////
// other constructors
//////////

// public

function::function(unsigned ser) : serial(ser)
{
}

function::function(unsigned ser, const exprseq & es) : exprseq(es), serial(ser)
{

	// Force re-evaluation even if the exprseq was already evaluated
	// (the exprseq copy constructor copies the flags)
	clearflag(status_flags::evaluated);
}

function::function(unsigned ser, const exvector & v)
  : exprseq(v), serial(ser)
{
}

function::function(unsigned ser, exvector && v)
  : exprseq(std::move(v)), serial(ser)
{
}

//////////
// archiving
//////////

/** Construct object from archive_node. */
void function::read_archive(const archive_node& n, lst& sym_lst)
{
	inherited::read_archive(n, sym_lst);
	// Find serial number by function name and number of parameters
	unsigned np = seq.size();
	std::string s;
	if (n.find_string("name", s)) {
		unsigned int ser = 0;
		for (auto & it : registered_functions()) {
			if (s == it.name && np == registered_functions()[ser].nparams) {
				serial = ser;
				return;
			}
			++ser;
		}
		throw (std::runtime_error("unknown function '" + s +
		                          "' with " + std::to_string(np) + " parameters in archive"));
	} else
		throw (std::runtime_error("unnamed function in archive"));
}

/** Archive the object. */
void function::archive(archive_node &n) const
{
	inherited::archive(n);
	GINAC_ASSERT(serial < registered_functions().size());
	n.add_string("name", registered_functions()[serial].name);
}

GINAC_BIND_UNARCHIVER(function);

//////////
// functions overriding virtual functions from base classes
//////////

// public

void function::print(const print_context & c, unsigned level) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	const auto &pdt = opt.print_dispatch_table;

	// Dynamically dispatch on print_context type
	const print_context_class_info *pc_info = &c.get_class_info();

next_context:
	unsigned id = pc_info->options.get_id();
	if (id >= pdt.size() || is_null(pdt[id])) {

		// Method not found, try parent print_context class
		const print_context_class_info *parent_pc_info = pc_info->get_parent();
		if (parent_pc_info) {
			pc_info = parent_pc_info;
			goto next_context;
		}

		// Method still not found, use default output
		if (is_a<print_tree>(c)) {

			c.s << std::string(level, ' ') << class_name() << " "
			    << opt.name << " @" << this
			    << std::hex << ", hash=0x" << hashvalue << ", flags=0x" << flags << std::dec
			    << ", nops=" << nops()
			    << std::endl;
			unsigned delta_indent = static_cast<const print_tree &>(c).delta_indent;
			for (size_t i=0; i<seq.size(); ++i)
				seq[i].print(c, level + delta_indent);
			c.s << std::string(level + delta_indent, ' ') << "=====" << std::endl;

		} else if (is_a<print_csrc>(c)) {

			// Print function name in lowercase
			std::string lname = opt.name;
			size_t num = lname.size();
			for (size_t i=0; i<num; i++)
				lname[i] = tolower(lname[i]);
			c.s << lname;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());

		} else if (is_a<print_latex>(c)) {
			c.s << opt.TeX_name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		} else {
			c.s << opt.name;
			printseq(c, '(', ',', ')', exprseq::precedence(), function::precedence());
		}

	} else {
		// Method found, call it
		call_funcp(pdt[id], seq, std::cref(c));
	}
}

ex function::eval() const
{
	if (flags & status_flags::evaluated) {
		return *this;
	}

	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Canonicalize argument order according to the symmetry properties
	if (seq.size() > 1 && !(opt.symtree.is_zero())) {
		exvector v = seq;
		GINAC_ASSERT(is_a<symmetry>(opt.symtree));
		int sig = canonicalize(v.begin(), ex_to<symmetry>(opt.symtree));
		if (sig != std::numeric_limits<int>::max()) {
			// Something has changed while sorting arguments, more evaluations later
			if (sig == 0)
				return _ex0;
			return ex(sig) * thiscontainer(std::move(v));
		}
	}

	if (is_null(opt.eval_f)) {
		return this->hold();
	}

	bool use_remember = opt.use_remember;
	ex eval_result;
	if (use_remember && lookup_remember_table(eval_result)) {
		return eval_result;
	}
	eval_result = call_funcp(opt.eval_f, seq);
	if (use_remember) {
		store_remember_table(eval_result);
	}
	return eval_result;
}

ex function::evalf() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	// Evaluate children first
	exvector eseq;
	if (!opt.evalf_params_first)
		eseq = seq;
	else {
		eseq.reserve(seq.size());
		for (auto & it : seq) {
			eseq.push_back(it.evalf());
		}
	}

	if (is_null(opt.evalf_f)) {
		return function(serial,eseq).hold();
	}
	return call_funcp(opt.evalf_f, eseq);
}

/**
 *  This method is defined to be in line with behavior of function::return_type()
 */
ex function::eval_ncmul(const exvector & v) const
{
	// If this function is called then the list of arguments is non-empty
	// and the first argument is non-commutative, see  function::return_type()
	return seq.begin()->eval_ncmul(v);
}

unsigned function::calchash() const
{
	unsigned v = golden_ratio_hash(make_hash_seed(typeid(*this)) ^ serial);
	for (size_t i=0; i<nops(); i++) {
		v = rotate_left(v);
		v ^= this->op(i).gethash();
	}

	if (flags & status_flags::evaluated) {
		setflag(status_flags::hash_calculated);
		hashvalue = v;
	}
	return v;
}

ex function::thiscontainer(const exvector & v) const
{
	return function(serial, v);
}

ex function::thiscontainer(exvector && v) const
{
	return function(serial, std::move(v));
}

/** Implementation of ex::series for functions.
 *  \@see ex::series */
ex function::series(const relational & r, int order, unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (is_null(opt.series_f)) {
		return basic::series(r, order);
	}
	ex res;
	try {
		res = call_funcp(opt.series_f, seq, std::cref(r), order, options);
	} catch (do_taylor) {
		res = basic::series(r, order, options);
	}
	return res;
}

/** Implementation of ex::conjugate for functions. */
ex function::conjugate() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (is_null(opt.conjugate_f))
		return conjugate_function(*this).hold();

	return call_funcp(opt.conjugate_f, seq);
}

/** Implementation of ex::real_part for functions. */
ex function::real_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (is_null(opt.real_part_f))
		return basic::real_part();

	return call_funcp(opt.real_part_f, seq);
}

/** Implementation of ex::imag_part for functions. */
ex function::imag_part() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (is_null(opt.imag_part_f))
		return basic::imag_part();

	return call_funcp(opt.imag_part_f, seq);
}

/** Implementation of ex::info for functions. */
bool function::info(unsigned inf) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options & opt = registered_functions()[serial];

	if (is_null(opt.info_f)) {
		return basic::info(inf);
	}
	return call_funcp(opt.info_f, seq, inf);
}

// protected

/** Implementation of ex::diff() for functions. It applies the chain rule,
 *  except for the Order term function.
 *  \@see ex::diff */
ex function::derivative(const symbol & s) const
{
	ex result;

	try {
		// Explicit derivation
		result = expl_derivative(s);
	} catch (...) {
		// Chain rule
		ex arg_diff;
		size_t num = seq.size();
		for (size_t i=0; i<num; i++) {
			arg_diff = seq[i].diff(s);
			// We apply the chain rule only when it makes sense.  This is not
			// just for performance reasons but also to allow functions to
			// throw when differentiated with respect to one of its arguments
			// without running into trouble with our automatic full
			// differentiation:
			if (!arg_diff.is_zero())
				result += pderivative(i)*arg_diff;
		}
	}
	return result;
}

int function::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return serial < o.serial ? -1 : 1;
	else
		return exprseq::compare_same_type(o);
}

bool function::is_equal_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	if (serial != o.serial)
		return false;
	else
		return exprseq::is_equal_same_type(o);
}

bool function::match_same_type(const basic & other) const
{
	GINAC_ASSERT(is_a<function>(other));
	const function & o = static_cast<const function &>(other);

	return serial == o.serial;
}

unsigned function::return_type() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return return_types::commutative;
		else
			return seq.begin()->return_type();
	}
}

return_type_t function::return_type_tinfo() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (opt.use_return_type) {
		// Return type was explicitly specified
		return opt.return_type_tinfo;
	} else {
		// Default behavior is to use the return type of the first
		// argument. Thus, exp() of a matrix behaves like a matrix, etc.
		if (seq.empty())
			return make_return_type_t<function>();
		else
			return seq.begin()->return_type_tinfo();
	}
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

// protected

ex function::pderivative(unsigned diff_param) const // partial differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	if (!is_null(opt.derivative_f)) {
		// Invoke the defined derivative function.
		return call_funcp(opt.derivative_f, seq, diff_param);
	}
	// No derivative defined? Fall back to abstract derivative object.
	return fderivative(serial, diff_param, seq);
}

ex function::expl_derivative(const symbol & s) const // explicit differentiation
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (!is_null(opt.expl_derivative_f)) {
		// Invoke the defined explicit derivative function.
		return call_funcp(opt.expl_derivative_f, seq, std::cref(s));
	}
	// There is no fallback for explicit derivative.
	throw(std::logic_error("function::expl_derivative(): explicit derivation is called, but no such function defined"));
}

ex function::power(const ex & power_param) const // power of function
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];
	
	if (!is_null(opt.power_f)) {
		// Invoke the defined power function.
		return call_funcp(opt.power_f, seq, power_param);
	}
	// No power function defined? Fall back to returning a power object.
	return dynallocate<GiNaC::power>(*this, power_param).setflag(status_flags::evaluated);
}

ex function::expand(unsigned options) const
{
	GINAC_ASSERT(serial<registered_functions().size());
	const function_options &opt = registered_functions()[serial];

	if (!is_null(opt.expand_f)) {
		// Invoke the defined expand function.
		return call_funcp(opt.expand_f, seq, options);

	}
	// No expand function defined? Return the same function with expanded arguments (if required)
	if (options & expand_options::expand_function_args)
		return inherited::expand(options);
	else
		return (options == 0) ? setflag(status_flags::expanded) : *this;
}

std::vector<function_options> & function::registered_functions()
{
	static std::vector<function_options> rf = std::vector<function_options>();
	return rf;
}

bool function::lookup_remember_table(ex & result) const
{
	return remember_table::remember_tables()[this->serial].lookup_entry(*this,result);
}

void function::store_remember_table(ex const & result) const
{
	remember_table::remember_tables()[this->serial].add_entry(*this,result);
}

// public

unsigned function::register_new(function_options const & opt)
{
	size_t same_name = 0;
	for (auto & i : registered_functions()) {
		if (i.name==opt.name) {
			++same_name;
		}
	}
	if (same_name>=opt.functions_with_same_name) {
		// we do not throw an exception here because this code is
		// usually executed before main(), so the exception could not
		// caught anyhow
		std::cerr << "WARNING: function name " << opt.name
		          << " already in use!" << std::endl;
	}
	registered_functions().push_back(opt);
	if (opt.use_remember) {
		remember_table::remember_tables().
			push_back(remember_table(opt.remember_size,
			                         opt.remember_assoc_size,
			                         opt.remember_strategy));
	} else {
		remember_table::remember_tables().push_back(remember_table());
	}
	return registered_functions().size()-1;
}

/** Find serial number of function by name and number of parameters.
 *  Throws exception if function was not found. */
unsigned function::find_function(const std::string &name, unsigned nparams)
{
	unsigned serial = 0;
	for (auto & it : function::registered_functions()) {
		if (it.get_name() == name && it.get_nparams() == nparams)
			return serial;
		++serial;
	}
	throw (std::runtime_error("no function '" + name + "' with " + std::to_string(nparams) + " parameters defined"));
}

/** Return the print name of the function. */
std::string function::get_name() const
{
	GINAC_ASSERT(serial<registered_functions().size());
	return registered_functions()[serial].name;
}

bool function_options::has_derivative() const {
	return !is_null(derivative_f);
}

bool function_options::has_power() const {
	return !is_null(power_f);
}

} // namespace GiNaC

