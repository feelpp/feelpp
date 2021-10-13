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
MeshBase<IndexT>::setSubMeshData( smd_ptrtype smd )
{
    M_smd = smd;
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
MeshBase<IndexT>::attachMeshSupport( MeshSupportBase* ms ) const
{
    auto itFoundMeshSupport = std::find_if( M_meshSupportsAttached.begin(), M_meshSupportsAttached.end(),
                                            [&ms]( MeshSupportBase* ms2) { return ms == ms2; });
    if ( itFoundMeshSupport != M_meshSupportsAttached.end() )
        return;
    M_meshSupportsAttached.push_back( ms );
}

template <typename IndexT>
void
MeshBase<IndexT>::detachMeshSupport( MeshSupportBase* ms ) const
{
    auto itFoundMeshSupport = std::find_if( M_meshSupportsAttached.begin(), M_meshSupportsAttached.end(),
                                            [&ms]( MeshSupportBase* ms2) { return ms == ms2; });
    if ( itFoundMeshSupport == M_meshSupportsAttached.end() )
        return;
    M_meshSupportsAttached.erase( itFoundMeshSupport );
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


template <typename IndexT>
typename MeshBase<IndexT>::size_type 
MeshBase<IndexT>::subMeshToMesh( size_type id ) const
{
    CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
    return M_smd->leftFind( id ).value();
}

template <typename IndexT>
typename MeshBase<IndexT>::size_type 
MeshBase<IndexT>::meshToSubMesh( size_type id ) const
{
    CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
    // if the submesh element id has not been found, return invalid value
    return M_smd->rightFind( id ).value_or( invalid_v<size_type> );
}
template <typename IndexT>
std::pair<std::vector<typename MeshBase<IndexT>::size_type>,bool> 
MeshBase<IndexT>::meshToSubMesh( std::vector<size_type> const& p, bool add_invalid_indices ) const
{
    CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
    std::vector<size_type> sid;//(std::distance(p.first,p.second) );
    bool has_invalid_values = false;
    std::for_each( p.begin(), p.end(), [&]( auto const& id ){
            if ( auto pid = M_smd->rightFind( id ); pid )
                sid.push_back( *pid );
            else
            {
                if ( add_invalid_indices )
                    sid.push_back( invalid_v<size_type> );
                has_invalid_values = true;
            }
            // the submesh element id has not been found, return invalid value
            //return invalid_v<size_type>;
        });
    return std::make_pair(sid,has_invalid_values);
}

template <typename IndexT>
typename MeshBase<IndexT>::size_type 
MeshBase<IndexT>::subMeshToMesh( std::shared_ptr<MeshBase> const& m, size_type id ) const
{
    if ( this == m.get() )
        return id;
    if ( isRelatedTo( m ) )
    {
        if ( this->isSubMeshFrom( m ) )
        {
            CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
            return M_smd->leftFind( id ).value();
        }
        else if ( this->isSiblingOf( m ) )
        {
            size_type id_in_parent =  M_smd->leftFind( id ).value();
            size_type id_in_sibling =  m->meshToSubMesh( id_in_parent );
            return id_in_sibling;
        }
    }
    return invalid_v<size_type>;
}

template <typename IndexT>
typename MeshBase<IndexT>::size_type 
MeshBase<IndexT>::meshToSubMesh( std::shared_ptr<MeshBase> const& m, size_type id ) const
{
        if ( this == m.get() )
            return id;
        if ( isRelatedTo( m ) )
        {
            if ( this->isSubMeshFrom( m ) )
            {
                CHECK( M_smd ) << "mesh doesn't have any submesh data\n";
                if ( auto p = M_smd->rightFind( id ); p )
                    return *p;
            }
            else if ( this->isSiblingOf( m ) )
            {
                size_type id_in_parent =  m->subMeshToMesh( id );
                size_type id_in_sibling =  this->meshToSubMesh( id_in_parent );
                return id_in_sibling;
            }
            // the submesh element id has not been found, return invalid value
            // will return invalid_v<size_type>
        }
        return invalid_v<size_type>;
}
template <typename IndexT>
void
MeshBase<IndexT>::meshesWithNodesSharedImpl( std::vector<std::shared_ptr<MeshBase<IndexT>>> & ret ) const
    {
        for ( std::weak_ptr<MeshBase<IndexT>> m : M_meshesWithNodesShared )
        {
            if ( m.expired() )
                continue;
            auto otherPtr = m.lock();
            auto foundOther = std::find_if(ret.begin(), ret.end(),
                                           [&otherPtr](auto const& ptr) { return otherPtr == ptr; });
            if ( foundOther == ret.end() )
            {
                ret.push_back( otherPtr );
                otherPtr->meshesWithNodesSharedImpl( ret );
            }
        }
    }

template class MeshBase<uint32_type>;
} // Feel
