#include <thread>
#include <vector>
#include <array>
#include <typeinfo>
#include <iostream>
#include <algorithm> 
#include <string>
#include <utility>
#include <functional>
#include <future>
#include <cassert>
#include <type_traits>
#include <list>
#include <ranges>
#include <stdexcept>


// THE MAIN OBJECTIF OF THIS PACKAGE 
// - to provide a range of META FUNCTIONS
// -



//================================================================================================================================
// META-FUNCTIONS FOR EXTRACTING THE n-th TYPE OF A PARAMETER PACK

// Declare primary template
template<int I, typename... Ts>
struct nth_type_of
{
};

// Base step
template<typename T, typename... Ts>
struct nth_type_of<0, T, Ts...>
{
    using type = T;
};

// Induction step
template<int I, typename T, typename... Ts>
struct nth_type_of<I, T, Ts...>
{
    using type = typename nth_type_of<I - 1, Ts...>::type;
};

// Helper meta-function for retrieving the first type in a parameter pack
template<typename... Ts>
struct first_type_of
{
    using type = typename nth_type_of<0, Ts...>::type;
};

// Helper meta-function for retrieving the last type in a parameter pack
template<typename... Ts>
struct last_type_of
{
    using type = typename nth_type_of<sizeof...(Ts) - 1, Ts...>::type;
};

//================================================================================================================================
// FUNCTIONS FOR EXTRACTING THE n-th VALUE OF AN ARGUMENT PACK

// Base step
template<int I, typename T, typename... Ts>
auto nth_value_of(T&& t, Ts&&... args) ->
    typename std::enable_if<(I == 0), decltype(std::forward<T>(t))>::type
{
    return std::forward<T>(t);
}

// Induction step
template<int I, typename T, typename... Ts>
auto nth_value_of(T&& t, Ts&&... args) ->
    typename std::enable_if<(I > 0), decltype(
        std::forward<typename nth_type_of<I, T, Ts...>::type>(
            std::declval<typename nth_type_of<I, T, Ts...>::type>()
            )
        )>::type
{
    using return_type = typename nth_type_of<I, T, Ts...>::type;
    return std::forward<return_type>(nth_value_of<I - 1>((std::forward<Ts>(args))...));
}

// Helper function for retrieving the first value of an argument pack
template<typename... Ts>
auto first_value_of(Ts&&... args) ->
    decltype(
        std::forward<typename first_type_of<Ts...>::type>(
            std::declval<typename first_type_of<Ts...>::type>()
            )
        )
{
    using return_type = typename first_type_of<Ts...>::type;
    return std::forward<return_type>(nth_value_of<0>((std::forward<Ts>(args))...));
}

// Helper function for retrieving the last value of an argument pack
template<typename... Ts>
auto last_value_of(Ts&&... args) ->
    decltype(
        std::forward<typename last_type_of<Ts...>::type>(
            std::declval<typename last_type_of<Ts...>::type>()
            )
        )
{
    using return_type = typename last_type_of<Ts...>::type;
    return std::forward<return_type>(nth_value_of<sizeof...(Ts) - 1>((std::forward<Ts>(args))...));
}

//================================================================================================================================
// METAFUNCTION FOR COMPUTING THE UNDERLYING TYPE OF HOMOGENEOUS PARAMETER PACKS

// Used as the underlying type of non-homogeneous parameter packs
struct null_type
{
};

// Declare primary template
template<typename... Ts>
struct homogeneous_type;

// Base step
template<typename T>
struct homogeneous_type<T>
{
    using type = T;
    static const bool isHomogeneous = true;
};

// Induction step
template<typename T, typename... Ts>
struct homogeneous_type<T, Ts...>
{
    // The underlying type of the tail of the parameter pack
    using type_of_remaining_parameters = typename homogeneous_type<Ts...>::type;

    // True if each parameter in the pack has the same type
    static const bool isHomogeneous = std::is_same<T, type_of_remaining_parameters>::value;

    // If isHomogeneous is "false", the underlying type is the ficTaskflow_HPCious null_type
    using type = typename std::conditional<isHomogeneous, T, null_type>::type;
};

// Meta-function to determine if a parameter pack is homogeneous
template<typename... Ts>
struct is_homogeneous_pack
{
    static const bool value = homogeneous_type<Ts...>::isHomogeneous;
};

//================================================================================================================================
// META-FUNCTIONS FOR CREATING INDEX LISTS

// The structure that encapsulates index lists
template <unsigned... Is>
struct index_list
{
};

// Collects internal details for generating index ranges [MIN, MAX)
namespace detail
{
    // Declare primary template for index range builder
    template <unsigned MIN, unsigned N, unsigned... Is>
    struct range_builder;

    // Base step
    template <unsigned MIN, unsigned... Is>
    struct range_builder<MIN, MIN, Is...>
    {
        typedef index_list<Is...> type;
    };

    // Induction step
    template <unsigned MIN, unsigned N, unsigned... Is>
    struct range_builder : public range_builder<MIN, N - 1, N - 1, Is...>
    {
    };
}

// Meta-function that returns a [MIN, MAX) index range
template<unsigned MIN, unsigned MAX>
using index_range = typename detail::range_builder<MIN, MAX>::type;



//================================================================================================================================
// CLASSES AND FUNCTIONS FOR REALIZING LOOPS ON ARGUMENT PACKS

// Collects internal details for implementing functor invocation
namespace detail
{
    // Functor invocation is realized through variadic inheritance.
    // The constructor of each base class invokes an input functor.
    // An functor invoker for an argument pack has one base class
    // for each argument in the pack

    // Realizes the invocation of the functor for one parameter
    template<unsigned I, typename T>
    struct invoker_base
    {
        template<typename F, typename U>
        invoker_base(F&& f, U&& u) { f(u); }
    };

    // Necessary because a class cannot inherit the same class twice
    template<unsigned I, typename T>
    struct indexed_type
    {
        static const unsigned int index = I;
        using type = T;
    };

    // The functor invoker: inherits from a list of base classes.
    // The constructor of each of these classes invokes the input
    // functor with one of the arguments in the pack.
    template<typename... Ts>
    struct invoker : public invoker_base<Ts::index, typename Ts::type>...
    {
        template<typename F, typename... Us>
        invoker(F&& f, Us&&... args)
            :
            invoker_base<Ts::index, typename Ts::type>(std::forward<F>(f), std::forward<Us>(args))...
        {
        }
    };
}

// The functor provided in the first argument is invoked for each
// argument in the pack whose index is contained in the index list
// specified in the second argument
template<typename F, unsigned... Is, typename... Ts>
void for_each_in_arg_pack_subset(F&& f, index_list<Is...> const& i, Ts&&... args)
{
    // Constructors of invoker's sub-objects will invoke the functor.
    // Note that argument types must be paired with numbers because the
    // implementation is based on inheritance, and one class cannot
    // inherit the same base class twice.
    detail::invoker<detail::indexed_type<Is, typename nth_type_of<Is, Ts...>::type>...> invoker(
        f,
        (nth_value_of<Is>(std::forward<Ts>(args)...))...
        );
}

// The functor provided in the first argument is invokedcd for each
// argument in the pack
template<typename F, typename... Ts>
void for_each_in_arg_pack(F&& f, Ts&&... args)
{
    for_each_in_arg_pack_subset(f, index_range<0, sizeof...(Ts)>(), std::forward<Ts>(args)...);
}

// The functor provided in the first argument is given in input the
// arguments in whose index is contained in the index list specified
// as the second argument.
template<typename F, unsigned... Is, typename... Ts>
void forward_subpack(F&& f, index_list<Is...> const& i, Ts&&... args)
{
    f((nth_value_of<Is>(std::forward<Ts>(args)...))...);
}

// The functor provided in the first argument is given in input all the
// arguments in the pack.
template<typename F, typename... Ts>
void forward_pack(F&& f, Ts&&... args)
{
    f(std::forward<Ts>(args)...);
}


//================================================================================================================================
// Returns a reference tuple... like this we can modify directly on the variable...

template<typename ...T, size_t... I>
auto makeTupleReferencesSub(std::tuple<T...>& t ,  std::index_sequence<I...>) { return std::tie(*std::get<I>(t)...) ;}

template<typename ...T>
auto makeTupleReferences( std::tuple<T...>& t ){ return makeTupleReferencesSub<T...>(t, std::make_index_sequence<sizeof...(T)>{}); }

//================================================================================================================================
// Writes a list of arguments of a tuple according to a sequence...

template< class T, T... Ints > 
class integer_sequence;

template<std::size_t... Ints>
using index_sequence = std::integer_sequence<std::size_t, Ints...>;

template <typename T>
void writeElem(const T& x) { std::cout << x << ','; };

template <typename TupleT, std::size_t... Is>
void writeTupleSequence(const TupleT& tp, std::index_sequence<Is...>) { (writeElem(std::get<Is>(tp)), ...); }


//================================================================================================================================
// Returns the value of a tuple without having the problems of the template td::get<...>(...) function during precompilation

template<
    typename Tuple,
    typename Indices=std::make_index_sequence<std::tuple_size<Tuple>::value>>
struct runtime_get_func_table;

template<typename Tuple,size_t ... Indices>
struct runtime_get_func_table<Tuple,std::index_sequence<Indices...>>{
    using return_type=typename std::tuple_element<0,Tuple>::type&;
    using get_func_ptr=return_type (*)(Tuple&) noexcept;
    static constexpr get_func_ptr table[std::tuple_size<Tuple>::value]={
        &std::get<Indices>...
    };
};

template<typename Tuple,size_t ... Indices>
constexpr typename
runtime_get_func_table<Tuple,std::index_sequence<Indices...>>::get_func_ptr
runtime_get_func_table<Tuple,std::index_sequence<Indices...>>::table[std::tuple_size<Tuple>::value];

template<typename Tuple>
constexpr
typename std::tuple_element<0,typename std::remove_reference<Tuple>::type>::type&
runtime_get_value_tuple(Tuple&& t,size_t index){
    using tuple_type=typename std::remove_reference<Tuple>::type;
    if(index>=std::tuple_size<tuple_type>::value)
        throw std::runtime_error("Out of range");
    return runtime_get_func_table<tuple_type>::table[index](t);
}		

//================================================================================================================================
// Constructs an tuple of index

template <int... Indices> struct indices;
template <> struct indices<-1> { typedef indices<> type; };
template <int... Indices>
struct indices<0, Indices...>
{
    typedef indices<0, Indices...> type;
};
template <int Index, int... Indices>
struct indices<Index, Indices...>
{
    typedef typename indices<Index - 1, Index, Indices...>::type type;
};

template <typename T>
typename indices<std::tuple_size<T>::value - 1>::type const* make_indices() { return 0; }

//================================================================================================================================
// Implementation of different methods to access the values of a tuple.
// A short tutorial with different methods to access it.

template <typename Function, typename... T, int... N>
void call_FwT1_impl(Function && fun, std::tuple<T...>&& t) {
    fun(std::get<N>(t)...);
}

template <typename Function, typename Tuple>
void call_FwT1(Function && fun, Tuple&& t)
{
    call_FwT1_impl(std::forward<Function>(fun), std::forward<Tuple>(t), make_indices<Tuple>());
}

template<typename Function, typename Tuple, size_t ... I>
auto call_FwT2(Function f, Tuple t, std::index_sequence<I ...>)
{
     return f(std::get<I>(t) ...);
}

template<typename Function, typename Tuple>
auto call_FwT2(Function f, Tuple t)
{
    static constexpr auto size = std::tuple_size<Tuple>::value;
    return call_FwT2(f, t, std::make_index_sequence<size>{});
}

 template<typename Function, typename Tuple, size_t ...S > 
 decltype(auto) apply_tuple_impl(Function && fn, Tuple&& t, std::index_sequence<S...>) 
 {
      return std::forward<Function>(fn)(std::get<S>(std::forward<Tuple>(t))...);
 }

 template<typename Function, typename Tuple>
 decltype(auto) apply_from_tuple(Function && fn, Tuple&& t)
 {
    std::size_t constexpr tSize=std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
    return apply_tuple_impl(std::forward<Function>(fn),std::forward<Tuple>(t),std::make_index_sequence<tSize>());
 }

 template<class Tuple, std::size_t N>
struct TuplePrinter
{
    static void write_tuple(const Tuple& t)
    {
        TuplePrinter<Tuple, N - 1>::write_tuple(t);
        std::cout << ", " << std::get<N-1>(t);
    }
};
 
template<class Tuple>
struct TuplePrinter<Tuple, 1>
{
    static void write_tuple(const Tuple& t)
    {
        std::cout << std::get<0>(t);
    }
};
 
template<typename... Args, std::enable_if_t<sizeof...(Args) == 0, int> = 0>
void write_tuple(const std::tuple<Args...>& t)
{
    std::cout << "()\n";
}
 
template<typename... Args, std::enable_if_t<sizeof...(Args) != 0, int> = 0>
void write_tuple(const std::tuple<Args...>& t)
{
    std::cout << "(";
    TuplePrinter<decltype(t), sizeof...(Args)>::write_tuple(t);
    std::cout << ")\n";
}

//================================================================================================================================

template<class Function, class... Args>
void async_wrapper(Function&& f, Args&&... args, std::future<void>& future,
                   std::future<void>&& is_valid, std::promise<void>&& is_moved) {
    is_valid.wait(); // Wait until the return value of std::async is written to "future"
    auto our_future = std::move(future); // Move "future" to a local variable
    is_moved.set_value(); // Only now we can leave void_async in the main thread

    // This is also used by std::async so that member function pointers work transparently
    auto functor = std::bind(f, std::forward<Args>(args)...);
    functor();
}

template<class Function, class... Args> // This is what you call instead of std::async
void void_async(Function&& f, Args&&... args) {
    std::future<void> future; // This is for std::async return value
    // This is for our synchronization of moving "future" between threads
    std::promise<void> valid;
    std::promise<void> is_moved;
    auto valid_future = valid.get_future();
    auto moved_future = is_moved.get_future();

    // Here we pass "future" as a reference, so that async_wrapper
    // can later work with std::async's return value
    future = std::async(
        async_wrapper<Function, Args...>,
        std::forward<Function>(f), std::forward<Args>(args)...,
        std::ref(future), std::move(valid_future), std::move(is_moved)
    );
    valid.set_value(); // Unblock async_wrapper waiting for "future" to become valid
    moved_future.wait(); // Wait for "future" to actually be moved
}


int long_running_task(int target, const std::atomic_bool& cancelled)
{
    // simulate a long running task for target*100ms, 
    // the task should check for cancelled often enough!
    while(target-- && !cancelled)
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    // return results to the future or raise an error 
    // in case of cancellation
    return cancelled ? 1 : 0;
}

//================================================================================================================================



template<std::size_t N, class TupleT, class NewT>
constexpr auto replace_tuple_element( const TupleT& t, const NewT& n )
{
    constexpr auto tail_size = std::tuple_size<TupleT>::value - N - 1;

    return [&]<std::size_t... I_head, std::size_t... I_tail>
        ( std::index_sequence<I_head...>, std::index_sequence<I_tail...> )
        {
            return std::tuple{
                std::get<I_head>( t )...,
                n,
                std::get<I_tail + N + 1>( t )...
            };
        }(  
           std::make_index_sequence<N>{}, 
           std::make_index_sequence<tail_size>{} 
          );
}

