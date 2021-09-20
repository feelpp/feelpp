/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-09

  Copyright (C) 2005,2006 EPFL

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
   \file meshbase.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-09
 */
#define FEELPP_INSTANTIATE_NOEXTERN_MESHBASE 1
#include <feel/feelmesh/meshbase.hpp>

namespace Feel
{

template <typename IndexT>
MeshBase<IndexT>::MeshBase( std::string const& name, uint16_type topoDim, uint16_type realDim, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    super2( "Mesh", name ),
    M_topodim( topoDim ),
    M_realdim( realDim ),
    M_components( MESH_ALL_COMPONENTS ),
    M_is_updated( false ),
    M_is_parametric( false ),
    M_n_vertices( 0 ),
    M_n_parts( 1 )
{
    DVLOG(2) << "[MeshBase] constructor...\n";
    DVLOG(2) << "[MeshBase] worldcomm:" << worldComm->globalRank() ;
    CHECK( worldComm ) << "invalid mesh worldcomm";
}


template <typename IndexT>
MeshBase<IndexT>::~MeshBase()
{
    std::vector<std::shared_ptr<MeshBase<IndexT>>> mwns;
    meshesWithNodesSharedImpl( mwns );
    for ( auto m : mwns )
        m->removeMeshWithNodesShared( this );
}

template <typename IndexT>
void
MeshBase<IndexT>::removeMeshWithNodesShared( MeshBase<IndexT> * m )
{
    auto newEnd = std::remove_if(M_meshesWithNodesShared.begin(), M_meshesWithNodesShared.end(),
                                 [&m](auto const& currentPtr) { return currentPtr.lock().get() == m; });
    M_meshesWithNodesShared.erase( newEnd, M_meshesWithNodesShared.end() );
}


template <typename IndexT>
void
MeshBase<IndexT>::clear()
{
    M_is_updated = false;

    M_n_vertices = 0;

    // Reset the number of partitions
    M_n_parts = 1;

    M_components = MESH_ALL_COMPONENTS;

    std::vector<std::shared_ptr<MeshBase<IndexT>>> mwns;
    meshesWithNodesSharedImpl( mwns );
    for ( auto m : mwns )
        m->removeMeshWithNodesShared( this );
}
template <typename IndexT>
void
MeshBase<IndexT>::updateForUse( size_type components )
{
    this->setComponents( components );
    this->updateForUse();
}

template <typename IndexT>
bool
MeshBase<IndexT>::isPartitioned() const
{
    if ( mpi::environment::initialized() )
        return M_n_parts == this->worldComm().localSize();

    else
        return M_n_parts == 1;
}


template <typename IndexT>
flag_type
MeshBase<IndexT>::markerId ( boost::any const& __marker ) const
{
    flag_type theflag = -1;
    if ( boost::any_cast<flag_type>( &__marker ) )
    {
        theflag = boost::any_cast<flag_type>( __marker);
    }
    else if ( boost::any_cast<int>( &__marker ) )
    {
        theflag = boost::any_cast<int>( __marker);
    }
    else if ( boost::any_cast<size_type>( &__marker ) )
    {
        theflag = boost::any_cast<size_type>( __marker);
    }
    else if ( boost::any_cast<uint16_type>( &__marker ) )
    {
        theflag = boost::any_cast<uint16_type>( __marker);
    }
    else if ( boost::any_cast<std::string>( &__marker ) )
    {
        theflag = this->markerName( boost::any_cast<std::string>( __marker) );
    }
    else if ( boost::any_cast<const char*>( &__marker ) )
    {
        theflag = this->markerName( std::string(boost::any_cast<const char*>( __marker) ) );
    }
    else
        CHECK( theflag != -1 ) << "invalid flag type\n";
    return theflag;
}

template <typename IndexT>
std::set<flag_type>
MeshBase<IndexT>::markersId( boost::any const& markerAny ) const
{
    std::set<flag_type> theflags;
    if ( boost::any_cast<std::vector<flag_type> >( &markerAny ) )
    {
        for ( flag_type const& marker : boost::any_cast< std::vector<flag_type> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::list<flag_type> >( &markerAny ) )
    {
        for ( flag_type const& marker : boost::any_cast< std::list<flag_type> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::set<flag_type> >( &markerAny ) )
    {
        for ( flag_type const& marker : boost::any_cast< std::set<flag_type> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::vector<int> >( &markerAny ) )
    {
        for ( int const& marker : boost::any_cast< std::vector<int> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::list<int> >( &markerAny ) )
    {
        for ( int const& marker : boost::any_cast< std::list<int> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::set<int> >( &markerAny ) )
    {
        for ( int const& marker : boost::any_cast< std::set<int> >( markerAny) )
            theflags.insert( marker );
    }
    else if ( boost::any_cast<std::vector<std::string> >( &markerAny ) )
    {
        for ( std::string const& marker : boost::any_cast< std::vector<std::string> >( markerAny) )
            theflags.insert( this->markerId( marker ) );
    }
    else if ( boost::any_cast<std::list<std::string> >( &markerAny ) )
    {
        for ( std::string const& marker : boost::any_cast< std::list<std::string> >( markerAny) )
            theflags.insert( this->markerId( marker ) );
    }
    else if ( boost::any_cast<std::set<std::string> >( &markerAny ) )
    {
        for ( std::string const& marker : boost::any_cast< std::set<std::string> >( markerAny) )
            theflags.insert( this->markerId( marker ) );
    }
    else
    {
        theflags.insert( this->markerId( markerAny ) );
    }
    return theflags;
}

template class MeshBase<uint32_type>;
} // Feel
