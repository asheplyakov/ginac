/** @file function.h
 *
 *  Interface to class of symbolic functions. */

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

#ifndef GINAC_FUNCTION_H
#define GINAC_FUNCTION_H

#include "exprseq.h"
#include "relational.h"
#include "symbol.h"

#include <string>
#include <vector>
#include <type_traits>
#include "mapbox/variant.hpp"

#define DECLARE_FUNCTION(NAME, ARITY) \
class NAME##_SERIAL { public: static unsigned serial; }; \
const unsigned NAME##_NPARAMS = ARITY; \
template<typename... Args> \
typename std::enable_if<sizeof...(Args) == ARITY, const GiNaC::function>::type \
NAME(Args const&... args) { \
	return GiNaC::function(NAME##_SERIAL::serial, {GiNaC::ex{args}...}); \
}
 
// compatibility
#define DECLARE_FUNCTION_1P(NAME) DECLARE_FUNCTION(NAME, 1)
#define DECLARE_FUNCTION_2P(NAME) DECLARE_FUNCTION(NAME, 2)
#define DECLARE_FUNCTION_3P(NAME) DECLARE_FUNCTION(NAME, 3)
#define DECLARE_FUNCTION_4P(NAME) DECLARE_FUNCTION(NAME, 4)
#define DECLARE_FUNCTION_5P(NAME) DECLARE_FUNCTION(NAME, 5)
#define DECLARE_FUNCTION_6P(NAME) DECLARE_FUNCTION(NAME, 6)
#define DECLARE_FUNCTION_7P(NAME) DECLARE_FUNCTION(NAME, 7)
#define DECLARE_FUNCTION_8P(NAME) DECLARE_FUNCTION(NAME, 8)
#define DECLARE_FUNCTION_9P(NAME) DECLARE_FUNCTION(NAME, 9)
#define DECLARE_FUNCTION_10P(NAME) DECLARE_FUNCTION(NAME, 10)
#define DECLARE_FUNCTION_11P(NAME) DECLARE_FUNCTION(NAME, 11)
#define DECLARE_FUNCTION_12P(NAME) DECLARE_FUNCTION(NAME, 12)
#define DECLARE_FUNCTION_13P(NAME) DECLARE_FUNCTION(NAME, 13)
#define DECLARE_FUNCTION_14P(NAME) DECLARE_FUNCTION(NAME, 14)

#define REGISTER_FUNCTION(NAME,OPT) \
unsigned NAME##_SERIAL::serial = \
	GiNaC::function::register_new(GiNaC::function_options(#NAME, NAME##_NPARAMS).OPT);

namespace GiNaC {
constexpr unsigned FUNCTION_MAX_ARGS = 14;

class function;
class symmetry;

template<int ...> struct seq {};

template<int N, int ...S> struct genseq : genseq<N-1, N-1, S...> { };

template<int ...S> struct genseq<0, S...> {
	typedef seq<S...> type;
};

namespace variant_ns = ::mapbox::util;

/* Declare a function ValT (*)(ArgsT, ArgsT, ... (N times), Rest...) */
template<typename ValT, typename ArgsT, int N, typename... Rest> struct nary_func_declarator {
	template<int x> struct dummy {
		typedef ArgsT type;
	};
	template<typename X> struct helper;
	template<int... S> struct helper<seq<S...>> {
		typedef ValT (*type)(typename dummy<S>::type..., Rest...);
	};
	typedef typename helper<typename genseq<N>::type>::type type;
};


template<typename ValT, typename... Rest> struct variant_funcp {
	// using vector_funcp = ValT (*)(const exvector&, Rest...);
	typedef ValT (*vector_funcp)(const exvector&, Rest...);

	template<typename X> struct helper;
	template<int... S> struct helper<seq<S...>> {
		typedef variant_ns::variant<
			vector_funcp,
			typename nary_func_declarator<ValT, const ex&, S, Rest...>::type...
			> type;
	};
	typedef typename helper<typename genseq<FUNCTION_MAX_ARGS>::type>::type type;
};

template<typename T> struct accepts_exvector_args_s {
	static constexpr bool value = false;
};

template<typename ValT, typename... Rest>
struct accepts_exvector_args_s<ValT(*)(const exvector&, Rest...)> {
	static constexpr bool value = true;
};

template<typename ValT, typename... Rest>
struct accepts_exvector_args_s<ValT(const exvector&, Rest...)> {
	static constexpr bool value = true;
};

template<typename F>
inline constexpr bool accepts_exvector_args(const F& func) {
	return accepts_exvector_args_s<F>::value;
}

template<typename T> struct get_arity_s;

template<typename R, typename ...Args>
struct get_arity_s<R(*)(Args...)> {
	static constexpr unsigned value = sizeof...(Args);
};

template<typename R, typename ...Args>
struct get_arity_s<R(Args...)> {
	static constexpr unsigned value = sizeof...(Args);
};

template<typename F>
inline constexpr unsigned get_arity(const F& func) {
	return get_arity_s<F>::value;
}

typedef variant_funcp<ex>::type eval_funcp_variant;
typedef variant_funcp<ex>::type evalf_funcp_variant;
typedef variant_funcp<ex>::type conjugate_funcp_variant;
typedef variant_funcp<ex>::type real_part_funcp_variant;
typedef variant_funcp<ex>::type imag_part_funcp_variant;
typedef variant_funcp<ex, unsigned>::type expand_funcp_variant;
typedef variant_funcp<ex, unsigned>::type derivative_funcp_variant;
typedef variant_funcp<ex, const symbol&>::type expl_derivative_funcp_variant;
typedef variant_funcp<ex, const ex&>::type power_funcp_variant;
typedef variant_funcp<ex, const relational&, int, unsigned>::type series_funcp_variant;
typedef variant_funcp<void, const print_context&>::type print_funcp_variant;
typedef variant_funcp<bool, unsigned>::type info_funcp_variant;

class function_options
{
	friend class function;
	friend class fderivative;
public:
	function_options();
	function_options(std::string const & n, std::string const & tn=std::string());
	function_options(std::string const & n, unsigned np);
	~function_options();
	void initialize();

	function_options & dummy() { return *this; }
	function_options & set_name(std::string const & n, std::string const & tn=std::string());
	function_options & latex_name(std::string const & tn);

	template<typename F>
	function_options& eval_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value);
		}
		eval_f = e;
		return *this;
	}
	template<typename F>
	function_options& evalf_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value);
		}
		evalf_f = e;
		return *this;
	}
	template<typename F>
	function_options& conjugate_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value);
		}
		conjugate_f = e;
		return *this;
	}
	template<typename F>
	function_options& real_part_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value);
		}
		real_part_f = e;
		return *this;
	}
	template<typename F>
	function_options& imag_part_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value);
		}
		imag_part_f = e;
		return *this;
	}
	template<typename F>
	function_options& expand_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 1);
		}
		expand_f = e;
		return *this;
	}
	template<typename F>
	function_options& derivative_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 1);
		}
		derivative_f = e;
		return *this;
	}
	template<typename F>
	function_options& expl_derivative_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 1);
		}
		expl_derivative_f = e;
		return *this;
	}
	template<typename F>
	function_options& power_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 1);
		}
		power_f = e;
		return *this;
	}
	template<typename F>
	function_options& series_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 3);
		}
		series_f = e;
		return *this;
	}
	template<typename F>
	function_options& info_func(F e) {
		if (!accepts_exvector_args_s<F>::value) {
			test_and_set_nparams(get_arity_s<F>::value - 1);
		}
		info_f = e;
		return *this;
	}

	template <class Ctx, typename F> function_options & print_func(F p)
	{
		if (!accepts_exvector_args(p)) {
			test_and_set_nparams(get_arity(p) - 1);
		}
		set_print_func(Ctx::get_class_info_static().options.get_id(), print_funcp_variant{p});
		return *this;
	}

	function_options & set_return_type(unsigned rt, const return_type_t* rtt = nullptr);
	function_options & do_not_evalf_params();
	function_options & remember(unsigned size, unsigned assoc_size=0,
	                            unsigned strategy=remember_strategies::delete_never);
	function_options & overloaded(unsigned o);
	function_options & set_symmetry(const symmetry & s);

	std::string get_name() const { return name; }
	unsigned get_nparams() const { return nparams; }

protected:
	bool has_derivative() const;
	bool has_power() const;
	void test_and_set_nparams(unsigned n);
	void set_print_func(unsigned id, print_funcp_variant&& f);

	std::string name;
	std::string TeX_name;

	unsigned nparams;

	eval_funcp_variant eval_f;
	evalf_funcp_variant evalf_f;
	conjugate_funcp_variant conjugate_f;
	real_part_funcp_variant real_part_f;
	imag_part_funcp_variant imag_part_f;
	expand_funcp_variant expand_f;
	derivative_funcp_variant derivative_f;
	expl_derivative_funcp_variant expl_derivative_f;
	power_funcp_variant power_f;
	series_funcp_variant series_f;
	std::vector<print_funcp_variant> print_dispatch_table;
	info_funcp_variant info_f;

	bool evalf_params_first;

	bool use_return_type;
	unsigned return_type;
	return_type_t return_type_tinfo;

	bool use_remember;
	unsigned remember_size;
	unsigned remember_assoc_size;
	unsigned remember_strategy;

	unsigned functions_with_same_name;

	ex symtree;
};


/** Exception class thrown by classes which provide their own series expansion
 *  to signal that ordinary Taylor expansion is safe. */
class do_taylor {};


/** The class function is used to implement builtin functions like sin, cos...
	and user defined functions */
class function : public exprseq
{
	GINAC_DECLARE_REGISTERED_CLASS(function, exprseq)

	friend class remember_table_entry;

// member functions

	// other constructors
public:
	function(unsigned ser);
	function(unsigned ser, std::initializer_list<ex> il) : exprseq{il}, serial(ser) { }
	template<typename... T> function(unsigned ser, const T&... args) : function{ser, {args...}} { }

	function(unsigned ser, const exprseq & es);
	function(unsigned ser, const exvector & v);
	function(unsigned ser, exvector && v);
	
	// functions overriding virtual functions from base classes
public:
	void print(const print_context & c, unsigned level = 0) const override;
	unsigned precedence() const override {return 70;}
	ex expand(unsigned options=0) const override;
	ex eval() const override;
	ex evalf() const override;
	ex eval_ncmul(const exvector & v) const override;
	unsigned calchash() const override;
	ex series(const relational & r, int order, unsigned options = 0) const override;
	ex thiscontainer(const exvector & v) const override;
	ex thiscontainer(exvector && v) const override;
	ex conjugate() const override;
	ex real_part() const override;
	ex imag_part() const override;
	void archive(archive_node& n) const override;
	void read_archive(const archive_node& n, lst& syms) override;
	bool info(unsigned inf) const override;
protected:
	ex derivative(const symbol & s) const override;
	bool is_equal_same_type(const basic & other) const override;
	bool match_same_type(const basic & other) const override;
	unsigned return_type() const override;
	return_type_t return_type_tinfo() const override;
	
	// new virtual functions which can be overridden by derived classes
	// none
	
	// non-virtual functions in this class
protected:
	ex pderivative(unsigned diff_param) const; // partial differentiation
	ex expl_derivative(const symbol & s) const; // partial differentiation
	static std::vector<function_options> & registered_functions();
	bool lookup_remember_table(ex & result) const;
	void store_remember_table(ex const & result) const;
public:
	ex power(const ex & exp) const;
	static unsigned register_new(function_options const & opt);
	static unsigned find_function(const std::string &name, unsigned nparams);
	static std::vector<function_options> get_registered_functions() { return registered_functions(); };
	unsigned get_serial() const {return serial;}
	std::string get_name() const;

// member variables

protected:
	unsigned serial;
};
GINAC_DECLARE_UNARCHIVER(function);

// utility functions/macros

template <typename T>
inline bool is_the_function(const ex & x)
{
	return is_exactly_a<function>(x)
	    && ex_to<function>(x).get_serial() == T::serial;
}

// Check whether OBJ is the specified symbolic function.
#define is_ex_the_function(OBJ, FUNCNAME) (GiNaC::is_the_function<FUNCNAME##_SERIAL>(OBJ))

} // namespace GiNaC

#endif // ndef GINAC_FUNCTION_H

