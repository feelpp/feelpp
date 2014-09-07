/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-02-24

  Copyright (C) 2006 EPFL
  Copyright (C) 2008 Université de Grenoble 1 (Joseph Fourier)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file region.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-02-24
 */
#ifndef FEELPP_REGION_HPP
#define FEELPP_REGION_HPP 1

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <feel/feelpoly/lagrange.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>

/**/
namespace Feel
{
namespace lambda = boost::lambda;
/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType, typename Expr>
typename SpaceType::element_type
region( boost::shared_ptr<SpaceType> const& space,
        Expr const& expr )
{
    BOOST_STATIC_ASSERT( SpaceType::fe_type::nOrder == 0 );
    typedef typename SpaceType::element_type element_type;
    element_type v( space, "field" );
    VLOG(1) << "[region] saving region with pid:\n";
    int pid = space->mesh()->comm().rank();
    VLOG(1) << "[region] saving region with pid: " << pid << "\n";
    auto it = space->mesh()->beginElementWithProcessId( pid );
    auto en = space->mesh()->endElementWithProcessId( pid );
    VLOG(1) << "[region] nb elements in region: " << std::distance( it, en ) << "\n";
    for ( ; it != en; ++it )
    {
        size_type dof_id = boost::get<0>( space->dof()->localToGlobal( it->id(),0, 0 ) );

        if ( dof_id >= v.firstLocalIndex() &&
                dof_id < v.lastLocalIndex() )
            v ( dof_id ) = expr( *it );
    }

    return v;
}
struct Region
{
    virtual ~Region() {}
};

/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType, typename Expr>
typename SpaceType::element_type
regionv( boost::shared_ptr<SpaceType> const& space,
         Expr const& expr )
{
    BOOST_STATIC_ASSERT( SpaceType::fe_type::nOrder == 0 );
    typedef typename SpaceType::element_type element_type;
    element_type v( space, "field" );

    int pid = space->mesh()->comm().rank();
    typename SpaceType::mesh_type::element_const_iterator it = space->mesh()->beginElementWithProcessId( pid );
    typename SpaceType::mesh_type::element_const_iterator en = space->mesh()->endElementWithProcessId( pid );

    for ( ; it != en; ++it )
    {
        size_type dof_id = boost::get<0>( space->dof()->localToGlobal( it->id(),0, 0 ) );

        if ( dof_id >= v.firstLocalIndex() &&
                dof_id < v.lastLocalIndex() )
            v ( dof_id ) = expr( *it ).value();
    }

    return v;
}

/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType>
typename SpaceType::element_type
regionProcess( boost::shared_ptr<SpaceType> const& space )
{
    return region( space, lambda::bind( &SpaceType::mesh_type::element_type::processId, lambda::_1 ) );

}

/**
 * class for RegionProcess
 */
struct RegionProcess : public Region
{
    template<typename SpaceType>
    typename SpaceType::element_type
    apply( boost::shared_ptr<SpaceType> const& space )
    {
        return regionProcess( space );
    }
};

/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType>
typename SpaceType::element_type
regionMarker( boost::shared_ptr<SpaceType> const& space )
{
    return region( space, lambda::bind( &SpaceType::mesh_type::element_type::marker, lambda::_1 ) );

}
/**
 * class for RegionMarker
 */
struct RegionMarkre : public Region
{
    template<typename SpaceType>
    typename SpaceType::element_type
    apply( boost::shared_ptr<SpaceType> const& space )
    {
        return regionMarker( space );
    }
};

/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType>
typename SpaceType::element_type
regionMarker2( boost::shared_ptr<SpaceType> const& space )
{
    return region( space, lambda::bind( &SpaceType::mesh_type::element_type::marker2, lambda::_1 ) );

}

/**
 * class for RegionMarker2
 */
struct RegionMarker2 : public Region
{
    template<typename SpaceType>
    typename SpaceType::element_type
    apply( boost::shared_ptr<SpaceType> const& space )
    {
        return regionMarker2( space );
    }
};

/**
 *
 * \ingroup SpaceTime
 */
template<typename SpaceType>
typename SpaceType::element_type
regionMarker3( boost::shared_ptr<SpaceType> const& space )
{
    return region( space, lambda::bind( &SpaceType::mesh_type::element_type::marker3, lambda::_1 ) );

}

/**
 * class for RegionMarker3
 */
struct RegionMarker3 : public Region
{
    template<typename SpaceType>
    typename SpaceType::element_type
    apply( boost::shared_ptr<SpaceType> const& space )
    {
        return regionMarker3( space );
    }
};

} // Feel

#endif /* FEELPP_REGION_HPP */
