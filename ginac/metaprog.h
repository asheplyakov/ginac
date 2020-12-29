#pragma once

template<int ...> struct seq {};

template<int N, int ...S> struct genseq : genseq<N-1, N-1, S...> { };

template<int ...S> struct genseq<0, S...> {
	typedef seq<S...> type;
};

template<int N, typename A>
static inline auto get(const A& array) /* XXX: C++14 does not need this: */ ->  decltype(array[N]) {
	return array[N];
}

template<int N>
struct array_unpacker {
	using seq_type = typename genseq<N>::type;

	template<typename F, typename C, int... S>
	inline auto do_call(const F& func, const C& params, seq<S...>) /* XXX: C++14 does not need this: */ -> decltype(func(get<S>(params)...)) {
		return func(get<S>(params)...);
	}

	template<typename F, typename C>
	inline auto call(const F& func, const C& params) /* XXX: C++14 does not need this: */ -> decltype(this->do_call(func, params, seq_type{})) {
		return do_call(func, params, seq_type{});
	}

};

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
constexpr unsigned get_arity(const F& func) {
	return get_arity_s<F>::value;
}

template<typename F, typename C>
static inline auto call_with_array(const F& f, const C& packed_args) -> decltype(array_unpacker<get_arity_s<F>::value>{}.call(f, packed_args)) {
	return array_unpacker<get_arity_s<F>::value>{}.call(f, packed_args);
}

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

