#ifndef _MODELCOREUTILS_HPP
#define _MODELCOREUTILS_HPP 1

namespace Feel {
namespace FeelModels {

namespace detail {

template<uint16_type, typename T> struct ChangeBasisOrder;

template<
    uint16_type NewOrder,
    template<uint16_type, template<uint16_type> class, typename, template<class, int, class> class, uint16_type > class BasisType,
    uint16_type Order,
    template<uint16_type> class PolySetType,
    typename ContinuityType,
    template<class, int, class> class Pts,
    uint16_type Tag
        >
struct ChangeBasisOrder<NewOrder, BasisType<Order, PolySetType, ContinuityType, Pts, Tag>>
{
    typedef BasisType<NewOrder, PolySetType, ContinuityType, Pts, Tag> type;
};

template<template<uint16_type> class, typename T> struct ChangeBasisPolySet;

template<
    template<uint16_type> class NewPolySetType,
    template<uint16_type, template<uint16_type> class, typename, template<class, int, class> class, uint16_type > class BasisType,
    uint16_type Order,
    template<uint16_type> class PolySetType,
    typename ContinuityType,
    template<class, int, class> class Pts,
    uint16_type Tag
    >
struct ChangeBasisPolySet<NewPolySetType, BasisType<Order, PolySetType, ContinuityType, Pts, Tag>>
{
    typedef BasisType<Order, NewPolySetType, ContinuityType, Pts, Tag> type;
};



template <typename T>
std::set<T>
set_difference( std::set<T> const& a,  std::set<T> const& b )
{
    std::set<T> res;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(res, res.begin()));
    return res;
}

template<typename C, typename T> struct ChangeBasisContinuity;

template<
    typename NewContinuityType,
    template<uint16_type, template<uint16_type> class, typename, template<class, int, class> class, uint16_type > class BasisType,
    uint16_type Order,
    template<uint16_type> class PolySetType,
    typename ContinuityType,
    template<class, int, class> class Pts,
    uint16_type Tag
        >
struct ChangeBasisContinuity<NewContinuityType, BasisType<Order, PolySetType, ContinuityType, Pts, Tag>>
{
    typedef BasisType<Order, PolySetType, NewContinuityType, Pts, Tag> type;
};


} // namespace detail
} // namespace FeelModels
} // namespace Feel

#endif // _MODELCOREUTILS_HPP
