///////////////////////////////////////////////////////////////////////////////
// foreach.hpp header file
//
// Copyright 2004 Eric Niebler.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_FOREACH
#include <cstddef>
#include <utility>  // for std::pair

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

// Some compilers allow temporaries to be bound to non-const references.
// These compilers make it impossible to for BOOST_FOREACH to detect
// temporaries and avoid reevaluation of the collection expression.
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)                                                       \
 || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))                                      \
 || (BOOST_WORKAROUND(BOOST_INTEL_CXX_VERSION, <= 700) && defined(_MSC_VER))
# define BOOST_FOREACH_NO_RVALUE_DETECTION
#endif

// Some compilers do not correctly implement the L-value/R-value conversion
// rules of the ternary conditional operator.
#if defined(BOOST_FOREACH_NO_RVALUE_DETECTION)                                                  \
 || BOOST_WORKAROUND(BOOST_MSVC, BOOST_TESTED_AT(1400))                                         \
 || BOOST_WORKAROUND(BOOST_INTEL_WIN, BOOST_TESTED_AT(800))                                     \
 || BOOST_WORKAROUND(__GNUC__, < 3)                                                             \
 || (BOOST_WORKAROUND(__GNUC__, == 3) && (__GNUC_MINOR__ <= 2))
# define BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
#endif

#include <boost/mpl/bool.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/range/end.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/result_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/utility/addressof.hpp>

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
# include <new>
# include <boost/aligned_storage.hpp>
# include <boost/utility/enable_if.hpp>
# include <boost/type_traits/is_array.hpp>
#endif

namespace boost
{

// forward declarations for iterator_range
template<typename T>
class iterator_range;

// forward declarations for sub_range
template<typename T>
class sub_range;

namespace foreach
{

///////////////////////////////////////////////////////////////////////////////
// in_range
//
template<typename T>
inline std::pair<T, T> in_range(T begin, T end)
{
    return std::make_pair(begin, end);
}

} // namespace foreach

namespace foreach_detail_
{

///////////////////////////////////////////////////////////////////////////////
// adl_begin/adl_end
//
template<typename T>
inline BOOST_DEDUCED_TYPENAME range_result_iterator<T>::type adl_begin(T &t)
{
    #if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))                                  \
     || BOOST_WORKAROUND(__GNUC__, < 3)
    return boost::begin(t);
    #else
    using boost::begin;
    typedef BOOST_DEDUCED_TYPENAME range_result_iterator<T>::type type;
    return type(begin(t));
    #endif
}

template<typename T>
inline BOOST_DEDUCED_TYPENAME range_result_iterator<T>::type adl_end(T &t)
{
    #if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))                                  \
     || BOOST_WORKAROUND(__GNUC__, < 3)
    return boost::end(t);
    #else
    using boost::end;
    typedef BOOST_DEDUCED_TYPENAME range_result_iterator<T>::type type;
    return type(end(t));
    #endif
}

///////////////////////////////////////////////////////////////////////////////
// auto_any_t/auto_any
//
struct auto_any_base
{
    // auto_any_base must evaluate to false in boolean context so that
    // they can be declared in if() statements.
    operator bool() const
    {
        return false;
    }
};

template<typename T>
struct auto_any : auto_any_base
{
    auto_any(T const &t)
        : item(t)
    {
    }

    // temporaries of type auto_any will be bound to const auto_any_base
    // references, but we still want to be able to mutate the stored
    // data, so declare it as mutable.
    mutable T item;
};

typedef auto_any_base const &auto_any_t;

template<typename T, typename C>
inline BOOST_DEDUCED_TYPENAME boost::mpl::if_<C, T const, T>::type &auto_any_cast(auto_any_t a)
{
    return static_cast<auto_any<T> const &>(a).item;
}

typedef boost::mpl::true_ const_;

///////////////////////////////////////////////////////////////////////////////
// type2type
//
template<typename T, typename C = boost::mpl::false_>
struct type2type
    : boost::mpl::if_<C, T const, T>
{
};

template<typename T, typename C = boost::mpl::false_>
struct foreach_iterator
{
    // If there is no function template ordering, then it may
    // be impossible to strip cv-modifiers from T, so use
    // range_result_iterator. Otherwise, use range_const_iterator
    // and range_iterator.
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    typedef BOOST_DEDUCED_TYPENAME boost::mpl::eval_if<
        C
      , range_const_iterator<T>
      , range_iterator<T>
    >::type type;
#else
    typedef BOOST_DEDUCED_TYPENAME boost::mpl::eval_if<
        C
      , range_result_iterator<T const>
      , range_result_iterator<T>
    >::type type;
#endif
};

template<typename T, typename C = boost::mpl::false_>
struct foreach_reference
    : iterator_reference<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
{
};

///////////////////////////////////////////////////////////////////////////////
// encode_type
//
template<typename T>
inline type2type<T> *encode_type(T &)
{
    return 0;
}

#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
template<typename T>
inline type2type<T, const_> *encode_type(T const &)
{
    return 0;
}
#endif

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION

///////////////////////////////////////////////////////////////////////////////
// rvalue_probe
//
struct rvalue_probe
{
    template<typename T>
    rvalue_probe(T const &t, bool &b)
        : ptemp(const_cast<T *>(&t))
        , rvalue(b)
    {
    }

    template<typename U>
    operator U()
    {
        rvalue = true;
        return *static_cast<U *>(ptemp);
    }

    template<typename V>
    operator V &() const
    {
        return *static_cast<V *>(ptemp);
    }

    void *ptemp;
    bool &rvalue;
};

///////////////////////////////////////////////////////////////////////////////
// simple_variant
//  holds either a T or a T*
template<typename T>
struct simple_variant
{
    simple_variant(T *t)
        : rvalue(false)
    {
        *static_cast<T **>(data.address()) = t;
    }

    simple_variant(T const &t)
        : rvalue(true)
    {
        ::new(data.address()) T(t);
    }

    simple_variant(simple_variant const &that)
        : rvalue(that.rvalue)
    {
        if(rvalue)
            ::new(data.address()) T(*that.get());
        else
            *static_cast<T **>(data.address()) = that.get();
    }

    ~simple_variant()
    {
        if(rvalue)
            get()->~T();
    }

    T *get() const
    {
        if(rvalue)
            return static_cast<T *>(data.address());
        else
            return *static_cast<T **>(data.address());
    }

private:
    enum { size = sizeof(T) > sizeof(T*) ? sizeof(T) : sizeof(T*) };
    simple_variant &operator =(simple_variant const &); 
    bool const                      rvalue;
    mutable aligned_storage<size>   data;
};

#elif !defined(BOOST_FOREACH_NO_RVALUE_DETECTION)

///////////////////////////////////////////////////////////////////////////////
// is_rvalue
//
template<typename T>
inline mpl::false_ *is_rvalue(T &, int)
{
	return 0;
}

template<typename T>
inline mpl::true_ *is_rvalue(T const &, ...)
{
	return 0;
}

#endif // BOOST_FOREACH_NO_CONST_RVALUE_DETECTION

///////////////////////////////////////////////////////////////////////////////
// set_false
//
inline bool set_false(bool &b)
{
    return b = false;
}

///////////////////////////////////////////////////////////////////////////////
// to_ptr
//
template<typename T>
inline T *to_ptr(T const &/*t*/) { return 0; }

// Borland needs a little extra help with arrays
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
template<typename T,std::size_t N>
inline T (*to_ptr(T (&t)[N]))[N] { return 0; }
#endif

///////////////////////////////////////////////////////////////////////////////
// cheap_copy
//   Overload this for user-defined collection types if they are inexpensive to copy.
//   This tells BOOST_FOREACH it can avoid the r-value/l-value detection stuff.
inline mpl::false_ *cheap_copy(...) { return 0; }

template<typename T>
inline mpl::true_ *cheap_copy(std::pair<T, T> *) { return 0; }

template<typename T>
inline mpl::true_ *cheap_copy(iterator_range<T> *) { return 0; }

template<typename T>
inline mpl::true_ *cheap_copy(sub_range<T> *) { return 0; }

template<typename T>
inline mpl::true_ *cheap_copy(T **) { return 0; }

template<typename T,std::size_t N>
inline mpl::false_ *cheap_copy(T (*)[N]) { return 0; }

///////////////////////////////////////////////////////////////////////////////
// derefof
//
template<typename T>
inline T &derefof(T *t)
{
    // This is a work-around for a compiler bug in Borland. If T* is a pointer to array type U(*)[N],
    // then dereferencing it results in a U* instead of U(&)[N]. The cast forces the issue.
    return reinterpret_cast<T &>(
        *const_cast<char *>(
            reinterpret_cast<char const volatile *>(t)
        )
    );
}

///////////////////////////////////////////////////////////////////////////////
// contain
//
template<typename T>
inline auto_any<T> contain(T const &t, void *, boost::mpl::true_ *)
{
    return t;
}

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
template<typename T>
inline auto_any<T *> contain(T &t, bool *, boost::mpl::false_ *)
{
    return boost::addressof(t);
}

template<typename T>
inline BOOST_DEDUCED_TYPENAME disable_if<
    is_array<T>
  , auto_any<simple_variant<T const> >
>::type
contain(T const &t, bool *rvalue, boost::mpl::false_ *)
{
    return *rvalue ? simple_variant<T const>(t) : simple_variant<T const>(&t);
}
#else
template<typename T>
inline auto_any<T *> contain(T &t, boost::mpl::false_ *, boost::mpl::false_ *) // l-value
{
    return boost::addressof(t);
}

template<typename T>
inline auto_any<T> contain(T const &t, boost::mpl::true_ *, boost::mpl::false_ *) // r-value
{
    return t;
}
#endif

/////////////////////////////////////////////////////////////////////////////
// begin
//
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
begin(auto_any_t col, type2type<T, C> *, void *, boost::mpl::true_ *)
{
    return foreach_detail_::adl_begin(auto_any_cast<T, C>(col));
}

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
begin(auto_any_t col, type2type<T, C> *, bool *, boost::mpl::false_ *)
{
    typedef BOOST_DEDUCED_TYPENAME type2type<T, C>::type type;
    return foreach_detail_::adl_begin(derefof(auto_any_cast<type *, boost::mpl::false_>(col)));
}

template<typename T>
inline BOOST_DEDUCED_TYPENAME disable_if<
    is_array<T>
  , auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, boost::mpl::true_>::type>
>::type
begin(auto_any_t col, type2type<T, const_> *, bool *, boost::mpl::false_ *)
{
    return foreach_detail_::adl_begin(*auto_any_cast<simple_variant<T const>, boost::mpl::false_>(col).get());
}
#else
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
begin(auto_any_t col, type2type<T, C> *, boost::mpl::false_ *, boost::mpl::false_ *) // l-value
{
    typedef BOOST_DEDUCED_TYPENAME type2type<T, C>::type type;
    typedef BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type iterator;
    return iterator(foreach_detail_::adl_begin(derefof(auto_any_cast<type *, boost::mpl::false_>(col))));
}

template<typename T>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, boost::mpl::true_>::type>
begin(auto_any_t col, type2type<T, const_> *, boost::mpl::true_ *, boost::mpl::false_ *) // r-value
{
    return foreach_detail_::adl_begin(auto_any_cast<T, boost::mpl::true_>(col));
}
#endif

///////////////////////////////////////////////////////////////////////////////
// end
//
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
end(auto_any_t col, type2type<T, C> *, void *, boost::mpl::true_ *)
{
    return foreach_detail_::adl_end(auto_any_cast<T, C>(col));
}

#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
template<typename T, typename C>
inline auto_any<int>
end(auto_any_t col, type2type<T *, C> *, void *, boost::mpl::true_ *)
{
    return 0; // not used
}
#endif

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
end(auto_any_t col, type2type<T, C> *, bool *, boost::mpl::false_ *)
{
    typedef BOOST_DEDUCED_TYPENAME type2type<T, C>::type type;
    return foreach_detail_::adl_end(derefof(auto_any_cast<type *, boost::mpl::false_>(col)));
}

template<typename T>
inline BOOST_DEDUCED_TYPENAME disable_if<
    is_array<T>
  , auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, boost::mpl::true_>::type>
>::type
end(auto_any_t col, type2type<T, const_> *, bool *, boost::mpl::false_ *)
{
    return foreach_detail_::adl_end(*auto_any_cast<simple_variant<T const>, boost::mpl::false_>(col).get());
}
#else
template<typename T, typename C>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type>
end(auto_any_t col, type2type<T, C> *, boost::mpl::false_ *, boost::mpl::false_ *) // l-value
{
    typedef BOOST_DEDUCED_TYPENAME type2type<T, C>::type type;
    typedef BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type iterator;
    return iterator(foreach_detail_::adl_end(derefof(auto_any_cast<type *, boost::mpl::false_>(col))));
}

template<typename T>
inline auto_any<BOOST_DEDUCED_TYPENAME foreach_iterator<T, boost::mpl::true_>::type>
end(auto_any_t col, type2type<T, const_> *, boost::mpl::true_ *, boost::mpl::false_ *) // r-value
{
    return foreach_detail_::adl_end(auto_any_cast<T, boost::mpl::true_>(col));
}
#endif

///////////////////////////////////////////////////////////////////////////////
// done
//
template<typename T, typename C>
inline bool done(auto_any_t cur, auto_any_t end, type2type<T, C> *)
{
    typedef BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type iter_t;
    return auto_any_cast<iter_t, boost::mpl::false_>(cur) == auto_any_cast<iter_t, boost::mpl::false_>(end);
}

#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
template<typename T, typename C>
inline bool done(auto_any_t cur, auto_any_t, type2type<T *, C> *)
{
    return ! *auto_any_cast<T *, boost::mpl::false_>(cur);
}
#endif

///////////////////////////////////////////////////////////////////////////////
// next
//
template<typename T, typename C>
inline void next(auto_any_t cur, type2type<T, C> *)
{
    typedef BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type iter_t;
    ++auto_any_cast<iter_t, boost::mpl::false_>(cur);
}

///////////////////////////////////////////////////////////////////////////////
// deref
//
template<typename T, typename C>
inline BOOST_DEDUCED_TYPENAME foreach_reference<T, C>::type
deref(auto_any_t cur, type2type<T, C> *)
{
    typedef BOOST_DEDUCED_TYPENAME foreach_iterator<T, C>::type iter_t;
    return *auto_any_cast<iter_t, boost::mpl::false_>(cur);
}

} // namespace foreach_detail_
} // namespace boost

#ifndef BOOST_FOREACH_NO_CONST_RVALUE_DETECTION
///////////////////////////////////////////////////////////////////////////////
// R-values and const R-values supported here
///////////////////////////////////////////////////////////////////////////////

// A sneaky way to get the type of the collection without evaluating the expression
# define BOOST_FOREACH_TYPEOF(COL)                                                              \
    (true ? 0 : boost::foreach_detail_::encode_type(COL))

// Evaluate the collection expression, and detect if it is an l-value or and r-value
# define BOOST_FOREACH_EVAL(COL)                                                                \
    (true ? boost::foreach_detail_::rvalue_probe((COL), _foreach_rvalue) : (COL))

// The R-value/L-value-ness of the collection expression is determined dynamically
# define BOOST_FOREACH_RVALUE(COL)                                                              \
    (&_foreach_rvalue)

# define BOOST_FOREACH_CHEAP_COPY(COL)                                                          \
    (true ? 0 : boost::foreach_detail_::cheap_copy(boost::foreach_detail_::to_ptr(COL)))

# define BOOST_FOREACH_NOOP(COL)                                                                \
    ((void)0)

#elif !defined(BOOST_FOREACH_NO_RVALUE_DETECTION)
///////////////////////////////////////////////////////////////////////////////
// R-values supported here, const R-values NOT supported here
///////////////////////////////////////////////////////////////////////////////

// A sneaky way to get the type of the collection without evaluating the expression
# define BOOST_FOREACH_TYPEOF(COL)                                                              \
    (true ? 0 : boost::foreach_detail_::encode_type(COL))

// Evaluate the collection expression
# define BOOST_FOREACH_EVAL(COL)                                                                \
    (COL)

// Determine whether the collection expression is an l-value or an r-value.
// NOTE: this gets the answer for const R-values wrong.
# define BOOST_FOREACH_RVALUE(COL)                                                              \
    (true ? 0 : boost::foreach_detail_::is_rvalue((COL), 0))

# define BOOST_FOREACH_CHEAP_COPY(COL)                                                          \
    (true ? 0 : boost::foreach_detail_::cheap_copy(boost::foreach_detail_::to_ptr(COL)))

# define BOOST_FOREACH_NOOP(COL)                                                                \
    ((void)0)

#else
///////////////////////////////////////////////////////////////////////////////
// R-values NOT supported here
///////////////////////////////////////////////////////////////////////////////

// A sneaky way to get the type of the collection without evaluating the expression
# define BOOST_FOREACH_TYPEOF(COL)                                                              \
    (true ? 0 : boost::foreach_detail_::encode_type(COL))

// Evaluate the collection expression
# define BOOST_FOREACH_EVAL(COL)                                                                \
    (COL)

// Can't use R-values with BOOST_FOREACH
# define BOOST_FOREACH_RVALUE(COL)                                                              \
    (static_cast<boost::mpl::false_ *>(0))

# define BOOST_FOREACH_CHEAP_COPY(COL)                                                          \
    (true ? 0 : boost::foreach_detail_::cheap_copy(boost::foreach_detail_::to_ptr(COL)))

// Attempt to make uses of BOOST_FOREACH with non-lvalues fail to compile
// BUGBUG but cheap-to-copy containers *would* be handled correctly. Hrm.
# define BOOST_FOREACH_NOOP(COL)                                                                \
    ((void)&(COL))

#endif


#define BOOST_FOREACH_CONTAIN(COL)                                                              \
    boost::foreach_detail_::contain(                                                            \
        BOOST_FOREACH_EVAL(COL)                                                                 \
      , BOOST_FOREACH_RVALUE(COL)                                                               \
      , BOOST_FOREACH_CHEAP_COPY(COL))

#define BOOST_FOREACH_BEGIN(COL)                                                                \
    boost::foreach_detail_::begin(                                                              \
        _foreach_col                                                                            \
      , BOOST_FOREACH_TYPEOF(COL)                                                               \
      , BOOST_FOREACH_RVALUE(COL)                                                               \
      , BOOST_FOREACH_CHEAP_COPY(COL))

#define BOOST_FOREACH_END(COL)                                                                  \
    boost::foreach_detail_::end(                                                                \
        _foreach_col                                                                            \
      , BOOST_FOREACH_TYPEOF(COL)                                                               \
      , BOOST_FOREACH_RVALUE(COL)                                                               \
      , BOOST_FOREACH_CHEAP_COPY(COL))

#define BOOST_FOREACH_DONE(COL)                                                                 \
    boost::foreach_detail_::done(                                                               \
        _foreach_cur                                                                            \
      , _foreach_end                                                                            \
      , BOOST_FOREACH_TYPEOF(COL))

#define BOOST_FOREACH_NEXT(COL)                                                                 \
    boost::foreach_detail_::next(                                                               \
        _foreach_cur                                                                            \
      , BOOST_FOREACH_TYPEOF(COL))

#define BOOST_FOREACH_DEREF(COL)                                                                \
    boost::foreach_detail_::deref(                                                              \
        _foreach_cur                                                                            \
      , BOOST_FOREACH_TYPEOF(COL))

///////////////////////////////////////////////////////////////////////////////
// BOOST_FOREACH
//
//   For iterating over collections. Collections can be
//   arrays, null-terminated strings, or STL containers.
//   The loop variable can be a value or reference. For
//   example:
//
//   std::list<int> int_list(/*stuff*/);
//   BOOST_FOREACH(int &i, int_list)
//   {
//       /* 
//        * loop body goes here.
//        * i is a reference to the int in int_list.
//        */
//   }
//
//   Alternately, you can declare the loop variable first,
//   so you can access it after the loop finishes. Obviously,
//   if you do it this way, then the loop variable cannot be
//   a reference.
//
//   int i;
//   BOOST_FOREACH(i, int_list)
//       { ... }
//
#define BOOST_FOREACH(VAR, COL)                                                                 \
    if (bool _foreach_rvalue = false) {} else                                                   \
    if (boost::foreach_detail_::auto_any_t _foreach_col = BOOST_FOREACH_CONTAIN(COL)) {} else   \
    if (boost::foreach_detail_::auto_any_t _foreach_cur = BOOST_FOREACH_BEGIN(COL)) {} else     \
    if (boost::foreach_detail_::auto_any_t _foreach_end = BOOST_FOREACH_END(COL)) {} else       \
    for (bool _foreach_continue = true;                                                         \
              _foreach_continue && !BOOST_FOREACH_DONE(COL);                                    \
              _foreach_continue ? BOOST_FOREACH_NEXT(COL) : BOOST_FOREACH_NOOP(COL))            \
        if  (boost::foreach_detail_::set_false(_foreach_continue)) {} else                      \
        for (VAR = BOOST_FOREACH_DEREF(COL); !_foreach_continue; _foreach_continue = true)

#endif
