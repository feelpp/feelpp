//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file geometryconceptwrappers.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _GEOMETRY_CONCEPT_WRAPPERS_HPP
#define _GEOMETRY_CONCEPT_WRAPPERS_HPP 1

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace Feel {
namespace detail {
namespace geometry {

struct RowMajor {};
struct ColumnMajor {};

template< typename MatrixType, uint16_type Dim >
struct pointWrapper: public std::reference_wrapper< MatrixType >
{
    typedef std::reference_wrapper< MatrixType > super_type;
    typedef std::remove_cv_t<std::remove_reference_t<MatrixType>> matrix_type;
    typedef typename matrix_type::value_type value_type;
    static inline const uint16_type nDim = Dim;

    // Inherit constructors
    using super_type::reference_wrapper;
    // Inherit operators
    using super_type::operator=;
    using super_type::operator MatrixType&;
    using super_type::get;
};

template< uint16_type Dim, typename MatrixType >
pointWrapper<MatrixType, Dim> pointWrap( MatrixType & m )
{
    return pointWrapper<MatrixType, Dim>( m );
}
template< uint16_type Dim, typename MatrixType >
pointWrapper<const MatrixType, Dim> pointCWrap( MatrixType const& m )
{
    return pointWrapper<const MatrixType, Dim>( m );
}

template< typename MatrixType, uint16_type Dim, typename Ordering = ColumnMajor >
struct segmentWrapper: public std::reference_wrapper< MatrixType >
{
    typedef std::reference_wrapper< MatrixType > super_type;
    typedef std::remove_cv_t<std::remove_reference_t<MatrixType>> matrix_type;
    typedef typename matrix_type::value_type value_type;
    static inline const uint16_type nDim = Dim;

    // Inherit constructors
    using super_type::reference_wrapper;
    // Inherit operators
    using super_type::operator=;
    using super_type::operator MatrixType&;
    using super_type::get;
};

template< uint16_type Dim, typename MatrixType, typename Ordering = ColumnMajor >
segmentWrapper<MatrixType, Dim, Ordering> segmentWrap( MatrixType & m )
{
    return segmentWrapper<MatrixType, Dim, Ordering>( m );
}
template< uint16_type Dim, typename MatrixType, typename Ordering = ColumnMajor >
segmentWrapper<const MatrixType, Dim, Ordering> segmentCWrap( MatrixType const& m )
{
    return segmentWrapper<const MatrixType, Dim, Ordering>( m );
}

template< typename MatrixType, uint16_type Dim, typename Ordering = ColumnMajor >
struct ringWrapper: public std::reference_wrapper< MatrixType >
{
    typedef std::reference_wrapper< MatrixType > super_type;
    typedef std::remove_cv_t<std::remove_reference_t<MatrixType>> matrix_type;
    typedef typename matrix_type::value_type value_type;
    static inline const uint16_type nDim = Dim;

    // Inherit constructors
    using super_type::reference_wrapper;
    // Inherit operators
    using super_type::operator=;
    using super_type::operator MatrixType&;
    using super_type::get;
};

template< uint16_type Dim, typename MatrixType, typename Ordering = ColumnMajor >
ringWrapper<MatrixType, Dim, Ordering> ringWrap( MatrixType & m )
{
    return ringWrapper<MatrixType, Dim, Ordering>( m );
}
template< uint16_type Dim, typename MatrixType, typename Ordering = ColumnMajor >
ringWrapper<const MatrixType, Dim, Ordering> ringCWrap( MatrixType const& m )
{
    return ringWrapper<const MatrixType, Dim, Ordering>( m );
}

} // namespace geometry
} // namespace detail
} // namespace Feel


namespace boost {
namespace geometry {

namespace traits {

// pointWrapper traits
template< typename MatrixType, Feel::uint16_type Dim >
struct tag< Feel::detail::geometry::pointWrapper<MatrixType, Dim> >
{
    typedef point_tag type;
};
template< typename MatrixType, Feel::uint16_type Dim >
struct coordinate_type< Feel::detail::geometry::pointWrapper<MatrixType, Dim> >
{
    typedef typename Feel::detail::geometry::pointWrapper<MatrixType, Dim>::value_type type;
};
template< typename MatrixType, Feel::uint16_type Dim >
struct coordinate_system< Feel::detail::geometry::pointWrapper<MatrixType, Dim> >
{
    typedef boost::geometry::cs::cartesian type;
};
template< typename MatrixType, Feel::uint16_type Dim >
struct dimension< Feel::detail::geometry::pointWrapper<MatrixType, Dim> > :
    boost::mpl::int_<Dim>
{};
template< typename MatrixType, Feel::uint16_type Dim, std::size_t AccessDim >
struct access< Feel::detail::geometry::pointWrapper<MatrixType, Dim>, AccessDim >
{
    typedef typename coordinate_type<Feel::detail::geometry::pointWrapper<MatrixType, Dim>>::type CoordinateType;

    static inline CoordinateType get( Feel::detail::geometry::pointWrapper<MatrixType, Dim> const& p )
    {
        return p.get().operator()(AccessDim);
    }

    static inline void set(
        Feel::detail::geometry::pointWrapper<MatrixType, Dim>& p,
        CoordinateType const& value )
    {
        p.get().operator()(AccessDim) = value;
    }
};

// segmentWrapper traits
template< typename MatrixType, Feel::uint16_type Dim, typename Ordering >
struct tag< Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Ordering> >
{
    typedef segment_tag type;
};
template< typename MatrixType, Feel::uint16_type Dim, typename Ordering >
struct point_type< Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Ordering> >
{
    typedef typename Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Ordering>::value_type value_type;
    typedef boost::geometry::model::point<value_type, Dim, boost::geometry::cs::cartesian> type;
};
// RowMajor access case
template< typename MatrixType, Feel::uint16_type Dim, std::size_t AccessDim >
struct indexed_access<Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::RowMajor>, 0, AccessDim>
{
    typedef Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::RowMajor> segment_type;
    typedef typename geometry::coordinate_type<segment_type>::type coordinate_type;

    static inline coordinate_type get(segment_type const& s)
    {
        return s.get().operator()(0, AccessDim);
    }

    static inline void set(segment_type& s, coordinate_type const& value)
    {
        s.get().operator()(0, AccessDim) = value;
    }
};
template <typename MatrixType, Feel::uint16_type Dim, std::size_t AccessDim>
struct indexed_access<Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::RowMajor>, 1, AccessDim>
{
    typedef Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::RowMajor> segment_type;
    typedef typename geometry::coordinate_type<segment_type>::type coordinate_type;

    static inline coordinate_type get(segment_type const& s)
    {
        return s.get().operator()(1, AccessDim);
    }

    static inline void set(segment_type& s, coordinate_type const& value)
    {
        s.get().operator()(1, AccessDim) = value;
    }
};
// ColumnMajor access case
template< typename MatrixType, Feel::uint16_type Dim, std::size_t AccessDim >
struct indexed_access<Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::ColumnMajor>, 0, AccessDim>
{
    typedef Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::ColumnMajor> segment_type;
    typedef typename geometry::coordinate_type<segment_type>::type coordinate_type;

    static inline coordinate_type get(segment_type const& s)
    {
        return s.get().operator()(AccessDim, 0);
    }

    static inline void set(segment_type& s, coordinate_type const& value)
    {
        s.get().operator()(AccessDim, 0) = value;
    }
};
template <typename MatrixType, Feel::uint16_type Dim, std::size_t AccessDim>
struct indexed_access<Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::ColumnMajor>, 1, AccessDim>
{
    typedef Feel::detail::geometry::segmentWrapper<MatrixType, Dim, Feel::detail::geometry::ColumnMajor> segment_type;
    typedef typename geometry::coordinate_type<segment_type>::type coordinate_type;

    static inline coordinate_type get(segment_type const& s)
    {
        return s.get().operator()(AccessDim, 1);
    }

    static inline void set(segment_type& s, coordinate_type const& value)
    {
        s.get().operator()(AccessDim, 1) = value;
    }
};

// ringWrapper traits
template< typename MatrixType, Feel::uint16_type Dim, typename Ordering >
struct tag< Feel::detail::geometry::ringWrapper<MatrixType, Dim, Ordering> >
{
    typedef ring_tag type;
};

} // namespace traits

} // namespace geometry
} // namespace boost

#endif // _GEOMETRY_CONCEPT_WRAPPERS_HPP
