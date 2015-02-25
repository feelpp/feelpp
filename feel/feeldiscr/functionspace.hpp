/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2004-11-22

   Copyright (C) 2004 EPFL
   Copyright (C) 2006-2012 Université Joseph Fourier (Grenoble I)

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
   \file FunctionSpace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-11-22
*/
#ifndef __FunctionSpace_H
#define __FunctionSpace_H 1

#include <type_traits>

#include <boost/static_assert.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/mpl/range_c.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/version.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
// #include <boost/numeric/ublas/vector_serialize.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/optional.hpp>
#include <boost/preprocessor/control/if.hpp>

#include <boost/smart_ptr/enable_shared_from_this.hpp>


//#include<boost/filesystem.hpp>


#include <stdexcept>
#include <sstream>
#include <limits>

#include <feel/feelcore/parameter.hpp>
#include <feel/feelpoly/operations.hpp>

#include <feel/feelalg/boundingbox.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelalg/vectorublas.hpp>

#include <feel/feelmesh/regiontree.hpp>
#include <feel/feelpoly/geomap.hpp>


#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/periodic.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feeldiscr/doftable.hpp>
#include <feel/feeldiscr/dofcomposite.hpp>
#include <feel/feeldiscr/parameter.hpp>
#include <feel/feeldiscr/bases.hpp>
#include <feel/feeldiscr/functionspacebase.hpp>
#include <feel/feeldiscr/mortar.hpp>

#include <feel/feeldiscr/region.hpp>
#include <feel/feelvf/exprbase.hpp>
#include <feel/feelvf/detail/gmc.hpp>

namespace Feel
{
namespace fusion = boost::fusion;
namespace parameter = boost::parameter;

namespace detail
{

template<typename T,int M, int N>
struct ID
{
    friend class boost::serialization::access;
    typedef T value_type;
    typedef Eigen::Matrix<value_type,M,N> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename  array_type::index_range range;

    struct result
    {

        typedef typename array_type::template array_view<1>::type type;
    };

    ID()
        :
        M_id()
    {}

    ID( ID const& id )
        :
        M_id( id.M_id )
    {}
    ID( array_type const& id )
        :
        M_id( id )
    {}


    template<typename Elem, typename ContextType>
    ID( Elem const& elem, ContextType const & context )
        :
        M_id( elem.idExtents( context ) )
    {
        for(int k = 0;k < M_id.shape()[0]; k++ )
            M_id[k].setZero();
        elem.id_( context, M_id );
    }

    ID& operator=( ID const& id )
    {
        if ( this != &id )
        {
            //VLOG(1) "[ID] extent = " <<

            M_id = id.M_id;
        }

        return *this;

    }
    m_type const& operator[]( uint16_type q  ) const
    {
        return M_id[q];
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
    {
        return M_id[q]( c1,c2 );
    }

    template <typename ExtentList>
    void resize( const ExtentList& extents )
    {
        M_id.resize( extents );
    }

    array_type M_id;

    template<class Archive>
    void save( Archive & ar, const unsigned int /*version*/ ) const
    {
        size_type e1 = M_id.shape()[0];
        DVLOG(2) << "saving in archive e1= " << e1 << "\n";
        ar  & e1;
        DVLOG(2) << "saving in archive array of size = " << M_id.num_elements() << "\n";
        ar  & boost::serialization::make_array( M_id.data(), M_id.num_elements() );
        DVLOG(2) << "saving in archive done\n";
    }
    template<class Archive>
    void load( Archive & ar, const unsigned int /*version*/ )
    {
        size_type e1, e2, e3;
        ar  & e1;
        DVLOG(2) << "loading from archive e1= " << e1 << "\n";
        M_id.resize( boost::extents[e1] );
        DVLOG(2) << "loading from archive array of size = " << M_id.num_elements() << "\n";
        ar  & boost::serialization::make_array( M_id.data(), M_id.num_elements() );
        DVLOG(2) << "loading from archive done\n";
        DVLOG(2) << "creating view interpolation context done\n";
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};
template<typename T,int M,int N>
struct DD
{
    typedef T value_type;
    friend class boost::serialization::access;
    typedef Eigen::Matrix<value_type,M,N> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename array_type::index_range range;
    struct result
    {

        typedef array_type type;
    };

    DD()
        :
        M_grad()
    {}

    template<typename Elem, typename ContextType>
    DD( Elem const& elem, ContextType const & context )
        :
        M_grad( elem.gradExtents( context ) )
    {
        for(int k = 0;k < M_grad.shape()[0]; k++ )
            M_grad[k].setZero();
        elem.grad_( context, M_grad );
    }

    m_type const& operator[]( uint16_type q  ) const
    {
        return M_grad[q];
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
    {
        return M_grad[q]( c1,c2 );
    }
    array_type M_grad;

    template<class Archive>
    void save( Archive & ar, const unsigned int /*version*/ ) const
    {
        size_type e1 = M_grad.shape()[0];
        DVLOG(2) << "saving in archive e1= " << e1 << "\n";
        ar  & e1;
        DVLOG(2) << "saving in archive array of size = " << M_grad.num_elements() << "\n";
        ar  & boost::serialization::make_array( M_grad.data(), M_grad.num_elements() );
        DVLOG(2) << "saving in archive done\n";
    }
    template<class Archive>
    void load( Archive & ar, const unsigned int /*version*/ )
    {
        size_type e1, e2, e3;
        ar  & e1;
        DVLOG(2) << "loading from archive e1= " << e1 << "\n";
        M_grad.resize( boost::extents[e1] );
        DVLOG(2) << "loading from archive array of size = " << M_grad.num_elements() << "\n";
        ar  & boost::serialization::make_array( M_grad.data(), M_grad.num_elements() );
        DVLOG(2) << "loading from archive done\n";
        DVLOG(2) << "creating view interpolation context done\n";
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

};
template<typename T,int N, int M, int P>
struct D
{
    typedef T value_type;
    typedef Eigen::Matrix<value_type,M,P> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename array_type::index_range range;

    struct result
    {

        typedef typename array_type::template array_view<1>::type type;
    };

    D()
        :
        M_grad()
    {}

    template<typename Elem, typename ContextType>
    D( Elem const& elem, ContextType const & context )
        :
        M_grad( elem.dExtents( context ) )
    {
        for(int k = 0;k < M_grad.shape()[0]; k++ )
            M_grad[k].setZero();
        elem.d_( N, context, M_grad );
    }

    m_type const& operator[]( uint16_type q  ) const
    {
        return M_grad[q];
    }
    value_type operator()( uint16_type c1, uint16_type /*c2*/, uint16_type q  ) const
    {
        return M_grad[q]( c1,0 );
    }
    array_type M_grad;
};

template<typename T>
struct Div
{
    typedef T value_type;
    typedef Eigen::Matrix<value_type,1,1> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename array_type::index_range range;
    struct result
    {

        typedef typename array_type::template array_view<1>::type type;
    };

    Div()
        :
        M_div()
    {}

    template<typename Elem, typename ContextType>
    Div( Elem const& elem, ContextType const & context )
        :
        M_div( elem.divExtents( context ) )
    {
        for(int k = 0;k < M_div.shape()[0]; k++ )
            M_div[k].setZero();
        elem.div_( context, M_div );
    }

    m_type const& operator[]( uint16_type q  ) const
    {
        return M_div[q];
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
    {
        return M_div[q]( c1,c2 );
    }
    array_type M_div;
};
template<typename T, int N, int D>
struct Curl
{
    typedef T value_type;
    typedef Eigen::Matrix<value_type,D,1> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename array_type::index_range range;
    struct result
    {

        typedef typename array_type::template array_view<1>::type type;
    };

    Curl()
        :
        M_curl()
    {}

    template<typename Elem, typename ContextType>
    Curl( Elem const& elem, ContextType const & context )
        :
        M_curl( elem.curlExtents( context ) )
    {
        for(int k = 0;k < M_curl.shape()[0]; k++ )
            M_curl[k].setZero();
        init( elem, context, boost::is_same<mpl::int_<N>, mpl::int_<-1> >() );
    }
    template<typename Elem, typename ContextType>
    void
    init( Elem const& elem, ContextType const& context, mpl::bool_<true> )
    {
        elem.curl_( context, M_curl );
    }
    template<typename Elem, typename ContextType>
    void
    init( Elem const& elem, ContextType const& context, mpl::bool_<false> )
    {
        elem.curl_( context, M_curl, N );
    }

    m_type const& operator[]( uint16_type q  ) const
    {
        return M_curl[q];
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
    {
        return this->operator()( c1, c2, q, mpl::int_<N>() );
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<-1>  ) const
    {
        return M_curl[q]( c1,c2 );
    }
    value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<0>  ) const
    {
        return M_curl[q]( 0,0 );
    }
    value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<1>  ) const
    {
        return M_curl[q]( 0,0 );
    }
    value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<2>  ) const
    {
        return M_curl[q]( 0,0 );
    }
    array_type M_curl;
};

template<typename T,int M, int N>
struct H
{
    friend class boost::serialization::access;
    typedef T value_type;
    typedef Eigen::Matrix<value_type,M,N> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename  array_type::index_range range;

    struct result
    {

        typedef typename array_type::template array_view<1>::type type;
    };

    H()
        :
        M_hess()
    {}

    template<typename Elem, typename ContextType>
    H( Elem const& elem, ContextType const & context )
        :
        M_hess( elem.hessExtents( context ) )
    {
        for(int k = 0;k < M_hess.shape()[0]; k++ )
            M_hess[k].setZero();
        elem.hess_( context, M_hess );
    }
    m_type const& operator[]( uint16_type q  ) const
    {
        return M_hess[q];
    }
    value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
    {
        return M_hess[q]( c1,c2 );
    }
    array_type M_hess;
};

template<typename MeshPtrType, typename PeriodicityType = Periodicity<NoPeriodicity> >
struct InitializeSpace
{
    InitializeSpace( MeshPtrType const& mesh,
                     PeriodicityType const& periodicity,
                     std::vector<Dof> const& dofindices,
                     std::vector<WorldComm> const & worldsComm,
                     std::vector<bool> extendedDofTable )
        :
        M_cursor( 0 ),
        M_worldsComm( worldsComm ),
        M_mesh( mesh ),
        M_dofindices( dofindices ),
        M_periodicity( periodicity ),
        M_extendedDofTable( extendedDofTable )
    {}
    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        operator()( x, is_shared_ptr<MeshPtrType>() );
    }
    template <typename T>
    void operator()( boost::shared_ptr<T> & x, mpl::bool_<true> ) const
    {
        auto p = *fusion::find<typename T::periodicity_0_type>(M_periodicity);
        x = boost::shared_ptr<T>( new T( M_mesh, M_dofindices, p,
                                         std::vector<WorldComm>( 1,M_worldsComm[M_cursor] ),
                                         std::vector<bool>( 1,M_extendedDofTable[M_cursor] ) ) );
        FEELPP_ASSERT( x ).error( "invalid function space" );

        ++M_cursor;// warning M_cursor < nb color
    }
    template <typename T>
    void operator()( boost::shared_ptr<T> & x, mpl::bool_<false> ) const
    {
        auto p = *fusion::find<typename T::periodicity_0_type>(M_periodicity);

        // look for T::mesh_ptrtype in MeshPtrType
        auto m = *fusion::find<typename T::mesh_ptrtype>(M_mesh);
        x = boost::shared_ptr<T>( new T( m, M_dofindices, p,
                                         std::vector<WorldComm>( 1,M_worldsComm[M_cursor] ),
                                         std::vector<bool>( 1,M_extendedDofTable[M_cursor] ) ) );
        FEELPP_ASSERT( x ).error( "invalid function space" );

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    std::vector<WorldComm> M_worldsComm;
    MeshPtrType M_mesh;
    std::vector<Dof> const& M_dofindices;
    PeriodicityType M_periodicity;
    std::vector<bool> M_extendedDofTable;
};
template<typename DofType>
struct updateDataMapProcess
{
    updateDataMapProcess( std::vector<WorldComm> const & worldsComm,
                          WorldComm const& worldCommFusion,
                          uint16_type lastCursor )
        :
        M_cursor( 0 ),
        M_start_index( 0 ),
        M_lastCursor( lastCursor ),
        M_worldsComm( worldsComm ),
        M_dm( new DofType( worldCommFusion ) ),
        M_dmOnOff( new DofType( worldCommFusion ) )
    {}

    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {

        if ( M_worldsComm[M_cursor].isActive() )
        {
            size_type nLocWithGhost=x->nLocalDofWithGhost();
            size_type nLocWithoutGhost=x->nLocalDofWithoutGhost();
            M_dm->setFirstDof( M_dm->worldComm().globalRank(), x->dof()->firstDof() );
            M_dm->setLastDof( M_dm->worldComm().globalRank(), x->dof()->lastDof() );
            M_dm->setFirstDofGlobalCluster( M_dm->worldComm().globalRank(), M_start_index + x->dof()->firstDofGlobalCluster() );
            M_dm->setLastDofGlobalCluster( M_dm->worldComm().globalRank(), M_start_index + x->dof()->lastDofGlobalCluster() );
            M_dm->setNLocalDofWithoutGhost( M_dm->worldComm().globalRank(), x->dof()->nLocalDofWithoutGhost() );
            M_dm->setNLocalDofWithGhost( M_dm->worldComm().globalRank(), x->dof()->nLocalDofWithGhost() );

            M_dm->resizeMapGlobalProcessToGlobalCluster( nLocWithGhost );
            M_dm->resizeMapGlobalClusterToGlobalProcess( nLocWithoutGhost );

            for ( size_type i=0; i<nLocWithGhost; ++i )
                M_dm->setMapGlobalProcessToGlobalCluster( i, M_start_index + x->dof()->mapGlobalProcessToGlobalCluster( i ) );

            for ( size_type i=0; i<nLocWithoutGhost; ++i )
                M_dm->setMapGlobalClusterToGlobalProcess( i, x->dof()->mapGlobalClusterToGlobalProcess( i ) );
        }


        if ( M_cursor==0 )
        {
            M_dmOnOff->setFirstDof( M_dmOnOff->worldComm().globalRank(), x->dof()->firstDof() );
            M_dmOnOff->setFirstDofGlobalCluster( M_dmOnOff->worldComm().globalRank(),
                                                  M_start_index + x->dof()->firstDofGlobalCluster() );
            M_dmOnOff->setNLocalDofWithoutGhost( M_dmOnOff->worldComm().globalRank(),
                                                  0 );
            M_dmOnOff->setNLocalDofWithGhost( M_dmOnOff->worldComm().globalRank(),
                                               0 );
        }

        if ( M_cursor==M_lastCursor )
        {
            M_dmOnOff->setLastDof( M_dmOnOff->worldComm().globalRank(),
                                    M_start_index + x->dof()->lastDof() );
            M_dmOnOff->setLastDofGlobalCluster( M_dmOnOff->worldComm().globalRank(),
                                                 M_start_index + x->dof()->lastDofGlobalCluster() );
        }

        // update nLoc
        size_type nLocWithoutGhostOnOff= M_dmOnOff->nLocalDofWithoutGhost() + x->dof()->nLocalDofWithoutGhost();
        size_type nLocWithGhostOnOff= M_dmOnOff->nLocalDofWithGhost() + x->dof()->nLocalDofWithGhost();

        M_dmOnOff->setNLocalDofWithoutGhost( M_dmOnOff->worldComm().globalRank(),
                                              nLocWithoutGhostOnOff );
        M_dmOnOff->setNLocalDofWithGhost( M_dmOnOff->worldComm().globalRank(),
                                           nLocWithGhostOnOff );

        // update map
        M_dmOnOff->resizeMapGlobalProcessToGlobalCluster( nLocWithGhostOnOff );
        M_dmOnOff->resizeMapGlobalClusterToGlobalProcess( nLocWithoutGhostOnOff );

        size_type startGlobClusterDof = M_dmOnOff->nLocalDofWithoutGhost() - x->dof()->nLocalDofWithoutGhost();
        size_type startGlobProcessDof = M_dmOnOff->nLocalDofWithGhost() - x->dof()->nLocalDofWithGhost();

        for ( size_type i=0; i<x->dof()->nLocalDofWithGhost(); ++i )
        {
            M_dmOnOff->setMapGlobalProcessToGlobalCluster( startGlobProcessDof + i, M_start_index + x->dof()->mapGlobalProcessToGlobalCluster( i ) );
        }

        for ( size_type i=0; i<x->dof()->nLocalDofWithoutGhost(); ++i )
        {
            M_dmOnOff->setMapGlobalClusterToGlobalProcess( startGlobClusterDof + i, x->dof()->mapGlobalClusterToGlobalProcess( i ) );
        }


        M_start_index+=x->nDof();

        ++M_cursor;// warning M_cursor < nb color
    }

    boost::shared_ptr<DofType> dataMap() const
    {
        return M_dm;
    }
    boost::shared_ptr<DofType> dataMapOnOff() const
    {
        return M_dmOnOff;
    }

    mutable uint16_type M_cursor;
    mutable size_type M_start_index;
    uint16_type M_lastCursor;
    std::vector<WorldComm> M_worldsComm;
    mutable boost::shared_ptr<DofType> M_dm;
    mutable boost::shared_ptr<DofType> M_dmOnOff;
}; // updateDataMapProcess



template<typename DofType>
struct updateDataMapProcessStandard
{
    updateDataMapProcessStandard( std::vector<WorldComm> const & worldsComm,
                                  WorldComm const& worldCommFusion,
                                  uint16_type lastCursor,
                                  std::vector<size_type> const& startDofGlobalCluster,
                                  std::vector<size_type> nLocWithoutGhost, std::vector<size_type> nLocWithGhost)
        :
        M_cursor( 0 ),
        M_start_index( worldCommFusion.globalSize(), 0 ),//  M_startDofGlobalCluster[worldCommFusion.globalRank()] ),
        M_startIndexWithGhost( 0 ),
        M_lastCursor( lastCursor ),
        M_worldsComm( worldsComm ),
        M_dm( new DofType( worldCommFusion ) ),
        M_startDofGlobalCluster(startDofGlobalCluster),
        M_nLocWithoutGhost(nLocWithoutGhost),
        M_nLocWithGhost(nLocWithGhost)
    {

        const int myrank = M_dm->worldComm().globalRank();
        const int worldsize = M_dm->worldComm().globalSize();
        for (int proc = 0 ; proc < worldsize ; ++proc)
        {
            M_dm->setNLocalDofWithoutGhost( proc, nLocWithoutGhost[proc] );
            M_dm->setNLocalDofWithGhost( proc, nLocWithGhost[proc] );
            M_start_index[proc] = M_startDofGlobalCluster[proc];// nLocWithoutGhost;
        }

        M_dm->resizeMapGlobalProcessToGlobalCluster( nLocWithGhost[myrank] );
        M_dm->resizeMapGlobalClusterToGlobalProcess( nLocWithoutGhost[myrank] );
    }

    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        const int myrank = M_dm->worldComm().globalRank();
        const int worldsize = M_dm->worldComm().globalSize();
        const size_type nLocWithGhost=x->dof()->nLocalDofWithGhost(myrank);
        const size_type nLocWithoutGhost=x->dof()->nLocalDofWithoutGhost(myrank);

        M_dm->addNeighborSubdomains( x->dof()->neighborSubdomains() );

        for (int proc = 0 ; proc < worldsize ; ++proc)
        {
            if ( M_cursor==0 )
            {
                M_dm->setFirstDof( proc,  x->dof()->firstDof(proc) );
                M_dm->setFirstDofGlobalCluster( proc, M_startDofGlobalCluster[proc] );
                M_dm->setLastDof( proc, x->dof()->lastDof(proc) );
                if (x->dof()->nLocalDofWithoutGhost(proc) == 0)
                    M_dm->setLastDofGlobalCluster( proc, M_startDofGlobalCluster[proc] );
                else
                    M_dm->setLastDofGlobalCluster( proc, M_startDofGlobalCluster[proc] + x->dof()->nLocalDofWithoutGhost(proc) - 1 );
            }
            else
            {
                M_dm->setLastDof( proc, M_dm->lastDof(proc) + x->dof()->nLocalDofWithGhost(proc)  );
                M_dm->setLastDofGlobalCluster( proc, M_dm->lastDofGlobalCluster(proc) + x->dof()->nLocalDofWithoutGhost(proc) );
            }
        }

        const size_type firstDofGC = M_dm->firstDofGlobalCluster();
        const size_type firstDofGCSubSpace = x->dof()->firstDofGlobalCluster();

        for (size_type gdof = x->dof()->firstDof(myrank) ; gdof < x->dof()->nLocalDofWithGhost(myrank) ; ++gdof )
        {
            const size_type localDof = M_startIndexWithGhost+gdof;//nLocalDofStart+gdof;
            const size_type gdofGC = x->dof()->mapGlobalProcessToGlobalCluster(gdof);

            if ( x->dof()->dofGlobalClusterIsOnProc( gdofGC ) )
            {
                const size_type globalDof = M_start_index[myrank] + ( gdofGC-firstDofGCSubSpace );
                M_dm->setMapGlobalProcessToGlobalCluster( localDof, globalDof );
                M_dm->setMapGlobalClusterToGlobalProcess( globalDof-firstDofGC ,localDof );
            }
            else
            {
                const int realproc = x->dof()->procOnGlobalCluster(gdofGC);
                const size_type globalDof = M_start_index[realproc] + (gdofGC- x->dof()->firstDofGlobalCluster(realproc));
                M_dm->setMapGlobalProcessToGlobalCluster( localDof, globalDof );
            }
        }

        for ( auto const& activeDofShared : x->dof()->activeDofSharedOnCluster() )
        {
            M_dm->setActiveDofSharedOnCluster( M_startIndexWithGhost + activeDofShared.first, activeDofShared.second );
        }



        for (int proc = 0 ; proc < worldsize ; ++proc)
            M_start_index[proc] += x->dof()->nLocalDofWithoutGhost(proc);
        M_startIndexWithGhost+=nLocWithGhost;

        ++M_cursor;
    }

    boost::shared_ptr<DofType> dataMap() const
    {
        return M_dm;
    }
    boost::shared_ptr<DofType> dataMapOnOff() const
    {
        return M_dm;
    }

    mutable uint16_type M_cursor;
    mutable std::vector<size_type> M_start_index;
    mutable size_type M_startIndexWithGhost;
    uint16_type M_lastCursor;
    std::vector<WorldComm> M_worldsComm;
    mutable boost::shared_ptr<DofType> M_dm;
    std::vector<size_type> M_startDofGlobalCluster;
    std::vector<size_type> M_nLocWithoutGhost, M_nLocWithGhost;
}; // updateDataMapProcessStandard





struct NbDof
{
    typedef size_type result_type;
    NbDof( size_type start = 0, size_type size = invalid_size_type_value )
        :
        M_cursor( start ),
        M_finish( size )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NbDof( T,S )>
#else
    struct result<NbDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};
    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( !x )
            return ret;

        if ( M_cursor < M_finish )
            ret += x->nDof();

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
};

#if 0
struct NLocalDof
{
    NLocalDof( size_type start = 0, size_type size = invalid_size_type_value )
        :
        M_cursor( start ),
        M_finish( size )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDof( T,S )>
#else
    struct result<NLocalDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};
    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
            ret += x->nLocalDof();

        ++M_cursor;
        return ret;
    }
    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
};
#else // MPI
template< typename IsWithGhostType>
struct NLocalDof
{

    NLocalDof( std::vector<WorldComm> const & worldsComm = std::vector<WorldComm>( 1,Environment::worldComm() ),
               bool useOffSubSpace = false,
               size_type start = 0, size_type size = invalid_size_type_value )
        :
        M_cursor( start ),
        M_finish( size ),
        M_worldsComm( worldsComm ),
        M_useOffSubSpace( useOffSubSpace )
    {}
    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDof( T,S )>
#else
    struct result<NLocalDof( S,T )>
#endif
:
    boost::remove_reference<S>
    {};

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<true> /**/ ) const
    {
        return x->nLocalDofWithGhost();
    }

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<false> /**/ ) const
    {
        return x->nLocalDofWithoutGhost();
    }

    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
        {
            if ( M_useOffSubSpace )
            {
                ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }

            else
            {
                if ( M_worldsComm[M_cursor].isActive() )
                    ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }
        }

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type M_cursor;
    size_type M_finish;
    std::vector<WorldComm> const& M_worldsComm;
    bool M_useOffSubSpace;
};
#endif // end MPI


template< typename IsWithGhostType>
struct NLocalDofOnProc
{

    NLocalDofOnProc( const int proc,
                     std::vector<WorldComm> const & worldsComm = std::vector<WorldComm>( 1,Environment::worldComm() ),
                     bool useOffSubSpace = false,
                     size_type start = 0, size_type size = invalid_size_type_value )
        :
        M_proc(proc),
        M_cursor( start ),
        M_finish( size ),
        M_worldsComm( worldsComm ),
        M_useOffSubSpace( useOffSubSpace )
    {}

    template<typename Sig>
    struct result;

    template<typename T, typename S>
#if BOOST_VERSION < 104200
    struct result<NLocalDofOnProc( T,S )>
#else
    struct result<NLocalDofOnProc( S,T )>
#endif
:
    boost::remove_reference<S>
    {};

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<true> /**/ ) const
    {
        return x->nLocalDofWithGhostOnProc(M_proc);
    }

    template <typename T>
    size_type
    nLocalDof( T const& x, mpl::bool_<false> /**/ ) const
    {
        return x->nLocalDofWithoutGhostOnProc(M_proc);
    }

    template <typename T>
    size_type
    operator()( T const& x, size_type s ) const
    {
        size_type ret = s;

        if ( M_cursor < M_finish )
        {
            if ( M_useOffSubSpace )
            {
                ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }

            else
            {
                if ( M_worldsComm[M_cursor].isActive() )
                    ret += nLocalDof( x, mpl::bool_<IsWithGhostType::value>() );
            }
        }

        ++M_cursor;
        return ret;
    }

    template <typename T>
    size_type
    operator()( size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    int M_proc;
    mutable size_type M_cursor;
    size_type M_finish;
    std::vector<WorldComm> const& M_worldsComm;
    bool M_useOffSubSpace;
}; // NLocalDofOnProc


template<int i,typename SpaceCompositeType>
struct InitializeContainersOff
{
    InitializeContainersOff( boost::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        if ( M_cursor==i && !x )
            x = boost::shared_ptr<T>( new T( M_space->template functionSpace<i>()->dof() ) );

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    boost::shared_ptr<SpaceCompositeType> M_space;
};


template<int i,typename SpaceCompositeType>
struct SendContainersOn
{
    SendContainersOn( boost::shared_ptr<SpaceCompositeType> const& _space,
                      std::vector<double> const& _dataToSend )
        :
        M_cursor( 0 ),
        M_space( _space ),
        M_dataToSend( _dataToSend )
    {}
    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        if ( M_cursor!=i )
        {
            int locRank=M_space->worldComm().localRank();
            int globRank=M_space->worldComm().localColorToGlobalRank( M_cursor,locRank );
            int tag = 0;
            //std::cout << "\n I am proc " << M_space->worldComm().globalRank()
            //          << " I send to proc " << globRank << std::endl;
            M_space->worldComm().globalComm().send( globRank,tag,M_dataToSend );
        }

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    boost::shared_ptr<SpaceCompositeType> M_space;
    std::vector<double> M_dataToSend;
};


template<int i,typename SpaceCompositeType>
struct RecvContainersOff
{
    RecvContainersOff( boost::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        if ( M_cursor==i )
        {
            std::vector<double> dataToRecv( M_space->template functionSpace<i>()->nLocalDof() );
            int locRank=M_space->worldComm().localRank();
            int globRank=M_space->worldComm().localColorToGlobalRank( i,locRank );
            int tag = 0;//locRank;
            //std::cout << "\n I am proc " << M_space->worldComm().globalRank()
            //          << " I recv to proc " << globRank << std::endl;
            M_space->worldComm().globalComm().recv( globRank,tag,dataToRecv );
            std::copy( dataToRecv.begin(), dataToRecv.end(), x->begin() );
        }

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    boost::shared_ptr<SpaceCompositeType> M_space;
};



template< typename map_type >
struct searchIndicesBySpace
{
    searchIndicesBySpace()
    {}

    searchIndicesBySpace( map_type& /*u*/ )
    {}

    template<typename T>
    searchIndicesBySpace( T const& fspace, map_type& u )
    {
        u = getIndicesFromSpace( fspace,u );
    }

    template<typename Sig>
    struct result;

    template<typename T, typename M>
#if BOOST_VERSION < 104200
    struct result<searchIndicesBySpace( T,M )>
#else
    struct result<searchIndicesBySpace( M,T )>
#endif
:
    boost::remove_reference<M>
    {};

    template < typename T >
    map_type getIndicesFromSpace( T const& fspace, map_type t ) const
    {
        if ( fspace->mesh()->numElements() == 0 )
            return t;

        size_type nProc = fspace->dof()->nProcessors();

        //search for the biggest index already in t; this will give the shift for the dofs
        std::vector< size_type > max_per_space;

        for ( size_type j=0; j<t.size(); j++ )
        {
            size_type _end = t[j].size();

            if ( _end )
                max_per_space.push_back( t[j][_end-1] );
        }

        //from all max indices found, determine the biggest
        size_type max_index = 0;

        if ( t.size() )
            max_index = *max_element( max_per_space.begin(), max_per_space.end() ) + 1;

        //std::cout << "maximum index " << max_index << "\n";

        //loop in all processors
        for ( size_type i=0; i<nProc; i++ )
        {
            /*
              std::cout << "Processor " << i << " has dofs"
              << " from " << fspace->dof()->firstDof(i)
              << " to " << fspace->dof()->lastDof(i) << "\n";
            */

            size_type _first = fspace->dof()->firstDof( i );
            size_type _last  = fspace->dof()->lastDof( i );

            //the dofs numbering for the current space start at max_index+1
            for ( size_type j=_first; j<=_last; j++ )
                t[i].push_back( max_index + j );
        }

        return t;
    }
    template <typename T>
    map_type
#if BOOST_VERSION < 104200
    operator()( T const& fspace, map_type t ) const
#else
    operator()( map_type t, T const& fspace ) const
#endif
    {
        return getIndicesFromSpace( fspace,t );
    }
};

// get start for each proc ->( proc0 : 0 ), (proc1 : sumdofproc0 ), (proc2 : sumdofproc0+sumdofproc1 ) ....
struct computeStartOfFieldSplit
{
    typedef boost::tuple< uint16_type , size_type > result_type;

    template<typename T>
    result_type operator()( result_type const &  previousRes, T const& t )
    {
        auto cptSpaces = previousRes.get<0>();
        auto start = previousRes.get<1>();

        for (int proc=0;proc<t->dof()->worldComm().globalSize();++proc)
            {
                if (proc < t->dof()->worldComm().globalRank())
                    start+=t->dof()->nLocalDofWithoutGhost(proc);
            }
        return boost::make_tuple( ++cptSpaces, start );
    }

};

// compute split
struct computeNDofForEachSpace
{

    computeNDofForEachSpace(size_type startSplit)
    :
        M_indexSplit( new IndexSplit() ),
        M_startSplit(startSplit)
    {}

    boost::shared_ptr<IndexSplit> const& indexSplit() const { return M_indexSplit; }

    typedef boost::tuple< uint16_type, size_type, IndexSplit > result_type;

#if 0
    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        const size_type nDofWithoutGhost = t->dof()->nLocalDofWithoutGhost();
        const size_type nDofWithGhost = t->dof()->nLocalDofWithGhost();
        const size_type firstDof = t->dof()->firstDofGlobalCluster();
        uint16_type cptSpaces = previousRes.get<0>();
        const size_type start = previousRes.get<1>();
        auto is = previousRes.get<2>();

        is[cptSpaces].resize(nDofWithoutGhost);

        for ( size_type i=0; i<nDofWithGhost; ++i )
        {
            if ( t->dof()->dofGlobalProcessIsGhost(i) ) continue;
            const size_type globalDof = t->dof()->mapGlobalProcessToGlobalCluster(i);
            M_is[cptSpaces][globalDof - firstDof ] = M_startSplit + start + (globalDof - firstDof);
        }

        return boost::make_tuple( ++cptSpaces, ( start+nDofWithoutGhost ), is );
    }
#else
    template<typename T>
    void operator()( T const& t ) const
    {
        M_indexSplit->addSplit( M_startSplit, t->map().indexSplit() );
    }
#endif

    mutable boost::shared_ptr<IndexSplit> M_indexSplit;
    size_type M_startSplit;
};

struct rebuildDofPointsTool
{

    template <typename T>
    void operator()( boost::shared_ptr<T> & x ) const
    {
        x->dof()->rebuildDofPoints( *x->mesh() );
    }
};

struct BasisName
{
    typedef std::string result_type;

    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        std::ostringstream os;

        if ( previousRes.size() )
            os << previousRes << "_" << t->basis()->familyName();

        else
            os << t->basis()->familyName();

        return os.str();
    }
};

struct BasisOrder
{
    typedef std::vector<int> result_type;

    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        std::vector<int> res( previousRes );
        res.push_back( t->nSubFunctionSpace() );
        return res;
    }
};


} // detail

enum ComponentType
{
    NO_COMPONENT = -1,
    X = 0,
    Y,
    Z,
    NX,
    NY,
    NZ,
    TX,
    TY,
    TZ
};
enum FunctionSpaceType
{
    SCALAR = 0,
    VECTORIAL = 1
};

template<uint16_type PN,
         uint16_type GN = 1>
struct Order
{
    static const uint16_type PolynomialOrder = PN;
    static const uint16_type GeometricOrder = GN;

    static const bool is_isoparametric = ( PN == GN );
    static const bool is_subparametric = ( PN > GN );
    static const bool is_surparametric = ( PN < GN );
};

typedef parameter::parameters<
    parameter::required<tag::mesh_type, boost::is_base_and_derived<MeshBase,_> >
    , parameter::optional<parameter::deduced<tag::bases_list>, boost::is_base_and_derived<Feel::detail::bases_base,_> >
    , parameter::optional<parameter::deduced<tag::value_type>, boost::is_floating_point<_> >
    , parameter::optional<parameter::deduced<tag::mortar_type>, boost::is_base_and_derived<Feel::detail::mortar_base,_> >
    , parameter::optional<parameter::deduced<tag::periodicity_type>, boost::is_base_and_derived<Feel::detail::periodicity_base,_> >
    > functionspace_signature;


/**
 * \class FunctionSpace
 * \ingroup SpaceTime
 * @author Function Space Class
 *
 * \c FunctionSpace is a representation of a functional space parametrized by
 * the type of the mesh (\c MeshType)
 *
 * @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
//template<typename MeshType, typename Basis_t, typename T_type = double, typename PeriodicityType = NoPeriodicity>
template<
typename A0,
         typename A1 = parameter::void_,
         typename A2 = parameter::void_,
         typename A3 = parameter::void_,
         typename A4 = parameter::void_>
class FunctionSpace
    :
public FunctionSpaceBase,
public boost::enable_shared_from_this<FunctionSpace<A0,A1,A2,A3,A4> >
{
public:
    template<typename FA0,typename FA1,typename FA2,typename FA3,typename FA4>
    friend class FunctionSpace;

    typedef typename functionspace_signature::bind<A0,A1,A2,A3,A4>::type args;

    typedef typename parameter::binding<args, tag::mesh_type>::type meshes_list;
    typedef typename parameter::binding<args, tag::value_type, double>::type value_type;
    typedef typename parameter::binding<args, tag::mortar_type, mortars<NoMortar> >::type mortar_list;
    typedef typename parameter::binding<args, tag::periodicity_type, Periodicity<NoPeriodicity> >::type periodicity_type;
    typedef typename parameter::binding<args, tag::bases_list, Feel::bases<Lagrange<1,Scalar> > >::type bases_list;

    BOOST_MPL_ASSERT_NOT( ( boost::is_same<mpl::at<bases_list,mpl::int_<0> >, mpl::void_> ) );

public:

    template<typename ThePeriodicityType, int pos>
    struct GetPeriodicity
    {
#if 0
        typedef typename boost::remove_reference<periodicity_type>::type periodicity_list_noref;
        typedef typename fusion::result_of::at_c<periodicity_list_noref, pos>::type _type;
        typedef typename boost::remove_reference<_type>::type type;
#else
        typedef typename mpl::if_<mpl::equal_to<fusion::result_of::size<ThePeriodicityType>,mpl::int_<1> >,
                                  mpl::identity<fusion::result_of::at_c<ThePeriodicityType,0> >,
                                  mpl::identity<fusion::result_of::at_c<ThePeriodicityType,pos> > >::type::type::type _type;
        typedef typename boost::remove_reference<_type>::type type;

#endif
    };
    template<typename TheMortarType, int pos>
    struct GetMortar
    {
        typedef typename mpl::if_<mpl::equal_to<fusion::result_of::size<TheMortarType>,mpl::int_<1> >,
                                  mpl::identity<fusion::result_of::at_c<TheMortarType,0> >,
                                  mpl::identity<fusion::result_of::at_c<TheMortarType,pos> > >::type::type::type _type;
        typedef typename boost::remove_reference<_type>::type type;
    };

    template<typename BasisType>
    struct ChangeMesh
    {
        typedef typename boost::remove_reference<meshes_list>::type meshes_list_noref;
        typedef typename boost::remove_reference<bases_list>::type bases_list_noref;
        typedef typename fusion::result_of::distance<typename fusion::result_of::begin<meshes_list_noref>::type,
                                                     typename fusion::result_of::find<bases_list_noref,BasisType>::type>::type pos;
        typedef typename fusion::result_of::at_c<meshes_list_noref, pos::value >::type _mesh_type;

        typedef FunctionSpace<typename boost::remove_reference<_mesh_type>::type,
                              Feel::detail::bases<BasisType>,value_type,
                              Periodicity<typename GetPeriodicity<periodicity_type,pos::value>::type >,
                              mortar_list> _type;
        typedef boost::shared_ptr<_type> type;


    };
    template<typename BasisType>
    struct ChangeBasis
    {
        //typedef typename mpl::if_<mpl::and_<boost::is_base_of<MeshBase, meshes_list >, boost::is_base_of<Feel::detail::periodic_base, periodicity_type > >,
        typedef typename boost::remove_reference<bases_list>::type bases_list_noref;
        typedef typename fusion::result_of::distance<typename fusion::result_of::begin<bases_list_noref>::type,
                                                     typename fusion::result_of::find<bases_list_noref,BasisType>::type>::type pos;

        typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                  mpl::identity<mpl::identity<boost::shared_ptr<FunctionSpace<meshes_list,Feel::detail::bases<BasisType>,value_type, Periodicity<typename GetPeriodicity<periodicity_type,pos::value>::type>, mortars<typename GetMortar<mortar_list,pos::value>::type > > > > >,
                                  mpl::identity<ChangeMesh<BasisType> > >::type::type::type type;

//mpl::identity<typename mpl::transform<meshes_list, ChangeMesh<mpl::_1,BasisType>, mpl::back_inserter<fusion::vector<> > >::type > >::type::type type;
    };
    typedef typename mpl::transform<bases_list, ChangeBasis<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type functionspace_vector_type;

    template<typename BasisType>
    struct ChangeBasisToComponentBasis
    {
        typedef typename BasisType::component_basis_type component_basis_type;
        //typedef typename mpl::if_<mpl::and_<boost::is_base_of<MeshBase, meshes_list >, boost::is_base_of<Feel::detail::periodic_base, periodicity_type > >,
        typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                  mpl::identity<mpl::identity<boost::shared_ptr<FunctionSpace<meshes_list,Feel::detail::bases<component_basis_type>,value_type, periodicity_type, mortar_list> > > >,
                                  mpl::identity<ChangeMesh<component_basis_type> > >::type::type::type type;
    };

    typedef typename mpl::transform<bases_list,
            ChangeBasisToComponentBasis<mpl::_1>,
            mpl::back_inserter<fusion::vector<> > >::type component_functionspace_vector_type;

    template<typename BasisType>
    struct GetComponentBasis
    {
        typedef typename BasisType::component_basis_type type;
    };

    typedef Feel::detail::bases2<typename mpl::transform<bases_list,
                                                         GetComponentBasis<mpl::_1>,
                                                         mpl::back_inserter<fusion::vector<> > >::type> component_basis_vector_type;



public:

    /** @name Constants
     */
    //@{
    static const bool is_composite = ( mpl::size<bases_list>::type::value > 1 );

    template<typename MeshListType,int N>
    struct GetMesh
    {
        typedef typename mpl::if_<mpl::or_<boost::is_base_of<MeshBase, MeshListType >,
                                           is_shared_ptr<MeshListType> >,
                                  mpl::identity<mpl::identity<MeshListType> >,
                                  mpl::identity<mpl::at_c<MeshListType,N> > >::type::type::type type;
    };
    // mesh
    typedef meshes_list MeshesListType;
    typedef typename GetMesh<meshes_list,0>::type mesh_0_type;
    typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                              mpl::identity<meshes_list>,
                              mpl::identity<mesh_0_type> >::type::type mesh_type;

    template<typename MeshType>
    struct ChangeToMeshPtr
    {
        typedef boost::shared_ptr<MeshType> type;
    };

    typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                              mpl::identity<boost::shared_ptr<mesh_type> >,
                              mpl::identity<typename mpl::transform<meshes_list, ChangeToMeshPtr<mpl::_1>, mpl::back_inserter<meshes<> > >::type  > >::type::type mesh_ptrtype;
    typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
            mpl::identity<typename mesh_type::element_type>,
            mpl::identity<typename mesh_0_type::element_type> >::type::type convex_type;


    template<typename BasisType>
    struct GetNComponents
    {
        //typedef mpl::int_<BasisType::template apply<mesh_type::nDim,value_type,typename mesh_type::element_type>::type::nComponents> type;
        typedef mpl::int_<BasisType::template apply<mesh_type::nDim,
                mesh_type::nRealDim,
                value_type,
                typename mesh_type::element_type>::type::nComponents> type;
    };
    struct nodim { static const int nDim = -1; static const int nRealDim = -1; };
    static constexpr uint16_type nDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                          mpl::identity<meshes_list >,
                                          mpl::identity<nodim> >::type::type::nDim;
    static constexpr uint16_type nRealDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                              mpl::identity<meshes_list>,
                                              mpl::identity<nodim> >::type::type::nRealDim;



    //typedef typename mpl::at_c<bases_list,0>::type::template apply<mesh_type::nDim,value_type,typename mesh_type::element_type>::type basis_0_type;
    typedef typename mpl::at_c<bases_list,0>::type::template apply<GetMesh<meshes_list,0>::type::nDim,
                                                                   GetMesh<meshes_list,0>::type::nRealDim,
                                                                   value_type,
                                                                   typename GetMesh<meshes_list,0>::type::element_type>::type basis_0_type;

    static constexpr uint16_type rank = ( is_composite? invalid_uint16_type_value : basis_0_type::rank );
    static constexpr bool is_scalar = ( is_composite? false : basis_0_type::is_scalar );
    static constexpr bool is_vectorial = ( is_composite? false : basis_0_type::is_vectorial );
    static constexpr bool is_tensor2 = ( is_composite? false : basis_0_type::is_tensor2 );
    static constexpr bool is_continuous = ( is_composite? false : basis_0_type::isContinuous );
    static constexpr bool is_modal = ( is_composite? false : basis_0_type::is_modal );
    static constexpr bool is_nodal = !is_modal;
    static constexpr uint16_type nComponents1 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents1 );
    static constexpr uint16_type nComponents2 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents2 );
    static constexpr bool is_product = ( is_composite? invalid_uint16_type_value : basis_0_type::is_product );
    typedef typename  basis_0_type::continuity_type continuity_type;

    static const uint16_type nComponents = mpl::transform<bases_list,
                             GetNComponents<mpl::_1>,
                             mpl::inserter<mpl::int_<0>,mpl::plus<mpl::_,mpl::_> > >::type::value;
    static const uint16_type N_COMPONENTS = nComponents;
    static const uint16_type nSpaces = mpl::size<bases_list>::type::value;

    typedef typename GetPeriodicity<periodicity_type,0>::type periodicity_0_type;
    static const bool is_periodic = periodicity_0_type::is_periodic;

    typedef typename GetMortar<mortar_list,0>::type mortar_0_type;
    static const bool is_mortar = mortar_0_type::is_mortar;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef typename ublas::type_traits<value_type>::real_type real_type;
    typedef typename node<value_type>::type node_type;


    typedef FunctionSpace<A0,A1,A2,A3,A4> functionspace_type;
    typedef functionspace_type space_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef boost::shared_ptr<functionspace_type> pointer_type;

    typedef FunctionSpace<meshes_list, component_basis_vector_type, value_type, periodicity_type, mortar_list> component_functionspace_type;
    typedef boost::shared_ptr<component_functionspace_type> component_functionspace_ptrtype;


    // basis
    typedef bases_list BasisType;
    template<int N>
    struct Basis
    {
        typedef typename mpl::at_c<bases_list,N>::type type;
    };
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<bases_list>,
            mpl::identity<basis_0_type> >::type::type basis_type;


    typedef boost::shared_ptr<basis_type> basis_ptrtype;
    typedef basis_type reference_element_type;
    typedef boost::shared_ptr<reference_element_type> reference_element_ptrtype;
    typedef reference_element_type fe_type;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<fe_type>,
                              mpl::identity<typename basis_0_type::SSpace::type> >::type::type mortar_fe_type;
    typedef reference_element_ptrtype fe_ptrtype;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<boost::none_t>,
            mpl::identity<typename basis_0_type::PreCompute> >::type pc_type;
    typedef boost::shared_ptr<pc_type> pc_ptrtype;

    // component basis
#if 0
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<component_basis_vector_type>,
            typename mpl::at_c<component_basis_vector_type,0>::type::template apply<nDim,value_type,convex_type> >::type::type component_basis_type;
#else
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<component_basis_vector_type>,
            typename mpl::at_c<component_basis_vector_type,0>::type::template apply<nDim,
            nRealDim,
            value_type,
            convex_type> >::type::type component_basis_type;
#endif
    typedef boost::shared_ptr<component_basis_type> component_basis_ptrtype;

    // trace space
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename mesh_type::trace_mesh_type>,
            mpl::identity<mpl::void_> >::type::type trace_mesh_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename mesh_type::trace_mesh_ptrtype>,
            mpl::identity<mpl::void_> >::type::type trace_mesh_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<FunctionSpace<trace_mesh_type, bases_list> >,
            mpl::identity<mpl::void_> >::type::type trace_functionspace_type;
    typedef typename boost::shared_ptr<trace_functionspace_type> trace_functionspace_ptrtype;

    // wirebasket
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename mesh_type::trace_trace_mesh_type>,
            mpl::identity<mpl::void_> >::type::type trace_trace_mesh_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename mesh_type::trace_trace_mesh_ptrtype>,
            mpl::identity<mpl::void_> >::type::type trace_trace_mesh_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<FunctionSpace<trace_trace_mesh_type, bases_list> >,
            mpl::identity<mpl::void_> >::type::type trace_trace_functionspace_type;
    typedef typename boost::shared_ptr<trace_trace_functionspace_type> trace_trace_functionspace_ptrtype;

#if 0
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename trace_functionspace_type::element_type>,
            mpl::identity<mpl::void_> >::type::type trace_element_type;
#endif

    // geomap
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm_type>, mpl::identity<mpl::void_> >::type::type gm_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm1_type>, mpl::identity<mpl::void_> >::type::type gm1_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::element_type>, mpl::identity<mpl::void_> >::type::type geoelement_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::face_type>, mpl::identity<mpl::void_> >::type::type geoface_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef boost::shared_ptr<gm1_type> gm1_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context<vm::POINT, geoelement_type> >,
            mpl::identity<mpl::void_> >::type::type pts_gmc_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB, geoelement_type> >,
                              mpl::identity<mpl::void_> >::type::type gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_ptrtype>, mpl::identity<mpl::void_> >::type::type geopc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_type>, mpl::identity<mpl::void_> >::type::type geopc_type;

    // basis context
    //typedef typename basis_type::template Context<vm::POINT, basis_type, gm_type, geoelement_type, pts_gmc_type::context> basis_context_type;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<mpl::void_>,
                              mpl::identity<typename basis_0_type::template Context<vm::POINT|vm::GRAD|vm::JACOBIAN|vm::HESSIAN|vm::KB, basis_0_type, gm_type, geoelement_type> > >::type::type basis_context_type;
                              //mpl::identity<typename basis_0_type::template Context<vm::POINT, basis_0_type, gm_type, geoelement_type, pts_gmc_type::context> > >::type::type basis_context_type;
    typedef boost::shared_ptr<basis_context_type> basis_context_ptrtype;

    // dof
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<DofComposite>,
                              mpl::identity<DofTable<mesh_type, basis_type, periodicity_0_type, mortar_0_type> > >::type::type dof_type;

    typedef boost::shared_ptr<dof_type> dof_ptrtype;
    typedef boost::shared_ptr<DataMap> datamap_ptrtype;
    typedef boost::shared_ptr<IndexSplit> indexsplit_ptrtype;

    // return types
    //typedef typename bases_list::polyset_type return_value_type;

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<N_COMPONENTS>, mpl::int_<1> >,
            mpl::identity<value_type>,
            mpl::identity<node_type> >::type::type return_type;

    typedef boost::function<return_type ( node_type const& )> function_type;
    //@}

    template<int i>
    struct sub_functionspace
    {
        typedef typename mpl::at_c<functionspace_vector_type,i>::type ptrtype;
        typedef typename ptrtype::element_type type;
    };

    /**
       @name Subclasses
    */
    //@{

    /**
     * \class Context
     */
    class Context
        :
        //public std::vector<basis_context_ptrtype>
        //the index of point is associated to a basis_context_ptrtype
        public std::map<int,basis_context_ptrtype>
    {
    public:
        static const bool is_rb_context = false;
        typedef std::map<int,basis_context_ptrtype> super;
        typedef typename super::value_type bc_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;
        typedef typename super::iterator iterator;
        Context( functionspace_ptrtype Xh ) : M_Xh( Xh ) {}
        virtual ~Context() {}

        std::pair<iterator, bool>
        add( node_type const& t )
        {
            return add( t, mpl::bool_<is_composite>() );
        }
        std::pair<iterator, bool>
        add( node_type const& t, mpl::bool_<true> )
        {
            return std::make_pair( this->end(), false );
        }
        std::pair<iterator, bool>
        add( node_type const& t, mpl::bool_<false> )
        {
            std::pair<iterator, bool> ret = std::make_pair(this->end(),false);
            //LOG(INFO)<<"add point\n";

            //rank of the current processor
            int proc_number = Environment::worldComm().globalRank();

            //total number of processors
            int nprocs = Environment::worldComm().globalSize();

            // add point t to list of points
            M_t.push_back( t );

            // localise t in space, find geometrical element in which t
            // belongs
            matrix_node_type m( mesh_type::nDim, 1 );
            for(int i = 0; i < mesh_type::nDim; ++i )
                m(i,0) = t(i);
            auto loc =  M_Xh->mesh()->tool_localization();
            loc->setExtrapolation( false );
            auto analysis = loc->run_analysis( m, invalid_size_type_value );
            auto found_points = analysis.template get<0>();
            bool found = found_points[0];

            std::vector<int> found_pt( nprocs, 0 );
            std::vector<int> global_found_pt( nprocs, 0 );

            if( found ) //we are on the proc that have the searched point
            {

                found_pt[proc_number]=1;

                auto it = loc->result_analysis_begin();
                auto en = loc->result_analysis_end();
                DCHECK( boost::next(it) == en ) << "Logic problem in finding one point in the mesh\n";
                auto eid = it->first;
                auto xref = boost::get<1>( *(it->second.begin()) );
                DVLOG(2) << "found point " << t << " in element " << eid << " on proc "<<proc_number<<"\n";
                DVLOG(2) << "  - reference coordinates " << xref << "\n";

                typename basis_type::points_type p(mesh_type::nDim,1);

                ublas::column( p, 0 ) = xref;
                // compute for each basis function in reference element its
                // value at \hat{t} in reference element
                auto basispc = M_Xh->basis()->preCompute( M_Xh->basis(), p );
                DVLOG(2) << "build precompute data structure for basis functions\n";
                auto gmpc = M_Xh->mesh()->gm()->preCompute( M_Xh->mesh()->gm(), p );
                DVLOG(2) << "build precompute data structure for geometric mapping\n";

                // build geometric mapping
                //auto gmc2 = M_Xh->mesh()->gm()->template context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB>( M_Xh->mesh()->element( eid ),
                //
                //                                                                                            gmpc );
                auto gmc = M_Xh->mesh()->gm()->template context<vm::POINT|vm::GRAD|vm::JACOBIAN|vm::HESSIAN|vm::KB>( M_Xh->mesh()->element( eid ), gmpc );
                //auto gmc = M_Xh->mesh()->gm()->template context<vm::POINT>( M_Xh->mesh()->element( eid ), gmpc );
                DVLOG(2) << "build geometric mapping context\n";
                // compute finite element context
                auto ctx = basis_context_ptrtype( new basis_context_type( M_Xh->basis(), gmc, basispc ) );
                DVLOG(2) << "build basis function context\n";

                //this->push_back( ctx );

                int number = M_t.size()-1;
                ret = this->insert( std::pair<int,basis_context_ptrtype>( number , ctx ) );
                //DVLOG(2) << "Context size: " << this->size() << "\n";

                if ( nprocs > 1 )
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );
                else
                    global_found_pt[ 0 ] = found_pt[ 0 ];

            }//if( found )
            else
            {
                if ( nprocs > 1 )
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<int>() );

            }//not found case

            //verify that the point is on a proc
            bool found_on_a_proc = false;
            int i;
            for (i = 0 ; i < global_found_pt.size(); ++i )
            {
                if ( global_found_pt[i] != 0 )
                {
                    DVLOG(2) << "processor " << i << " has the point " << t << "\n";
                    found_on_a_proc = true;
                    M_t_proc.push_back(i);
                    break;
                }
            }
            CHECK( found_on_a_proc ) << "the point " << t << " was not found ! \n";
            return ret;

        }//add ( non composite case )


        void addCtx( basis_context_ptrtype ctx , int proc_having_the_point)
        {
            int position = M_t.size();
            this->insert( std::pair<int,basis_context_ptrtype>( position , ctx ) );
            node_type n;//only to increase M_t ( this function may be called during online step of crb )
            M_t.push_back( n );
            M_t_proc.push_back( proc_having_the_point );
        }
        void removeCtx()
        {
            this->clear();
            M_t.clear();
            M_t_proc.clear();
        }

        int nPoints() const
        {
            return M_t.size();
        }

        int processorHavingPoint( int point_number ) const
        {
            return M_t_proc[point_number];
        }

        functionspace_ptrtype functionSpace() const
        {
            return M_Xh;
        }

        virtual functionspace_type* ptrFunctionSpace() const
        {
            LOG( INFO ) << "fem ptrFunctionSpace()";
            return M_Xh.get();
        }

        node_type const& node(int i) const
        {
            int size = M_t.size();
            CHECK( i < size ) <<" i  = "<<i<<" and the context has "<< size<<" points \n";
            return M_t[i];
        }

    private:

        std::vector<node_type> M_t;
        std::vector<int> M_t_proc;//point number i is on proc M_t_proc[i]
        functionspace_ptrtype M_Xh;


    };
    /**
     * \return function space context
     */
    Context context() { return Context( this->shared_from_this() ); }

    /*virtual*/ basis_context_ptrtype contextBasis( std::pair<int, basis_context_ptrtype> const& p, Context const& c ) { return p.second; }

    /**
     * \class Element
     */
    template<typename T = double,  typename Cont = VectorUblas<T> >
    class Element
        :
        public Cont,boost::addable<Element<T,Cont> >, boost::subtractable<Element<T,Cont> >, FunctionSpaceBase::ElementBase
    {
    public:
        typedef T value_type;

        template<typename BasisType,typename keyType>
        struct ChangeElement
        {
            typedef T value_type;
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef typename ChangeBasis<BasisType>::type::element_type fs_type;
            typedef typename fs_type::template Element<value_type, typename VectorUblas<T>::range::type > element_type;
            typedef std::pair< keyType, boost::shared_ptr<element_type> > the_type;

            typedef typename mpl::if_<mpl::bool_<FunctionSpace<A0,A1,A2,A3,A4>::is_composite>,
                                      mpl::identity<the_type>,
                                      mpl::identity<boost::none_t> >::type::type type;
        };

        typedef mpl::range_c<int,0, FunctionSpace<A0,A1,A2,A3,A4>::nSpaces> rangeElementStorageType;
        typedef typename mpl::transform<bases_list, rangeElementStorageType, ChangeElement<mpl::_1,mpl::_2>, mpl::back_inserter<fusion::vector<> > >::type element_vector_type;

        //typedef typename fusion::result_of::accumulate<bases_list, fusion::vector<>, ChangeElement<> >
        typedef typename VectorUblas<T>::range::type ct_type;

        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> eigen_type;

        /**
         * usefull in // with composite case
         * store the off views
         */
        template<typename BasisType>
        struct AddOffContainer
        {
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef boost::shared_ptr<Cont> type;
        };
        typedef typename mpl::transform<bases_list, AddOffContainer<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type container_vector_type;



        typedef FunctionSpace<A0,A1,A2,A3,A4> functionspace_type;
        typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

        static const uint16_type nDim = mesh_type::nDim;
        static const uint16_type nRealDim = mesh_type::nRealDim;
        static const bool is_composite = functionspace_type::is_composite;
        static const bool is_scalar = functionspace_type::is_scalar;
        static const bool is_vectorial = functionspace_type::is_vectorial;
        static const bool is_tensor2 = functionspace_type::is_tensor2;
        static const bool is_continuous = functionspace_type::is_continuous;
        static const uint16_type nComponents1 = functionspace_type::nComponents1;
        static const uint16_type nComponents2 = functionspace_type::nComponents2;
        static const uint16_type nComponents = functionspace_type::nComponents;
        static const uint16_type nSpaces = functionspace_type::nSpaces;
        static const bool is_mortar = functionspace_type::is_mortar;

        /** @name Typedefs
         */
        //@{


        typedef typename ublas::type_traits<value_type>::real_type real_type;
        typedef T element_type;

        typedef Cont super;
        typedef Cont container_type;
        typedef container_type vector_temporary_type;

        typedef typename mpl::if_<mpl::bool_<is_composite>,
                mpl::identity<boost::none_t>,
                mpl::identity<typename basis_0_type::polyset_type> >::type::type polyset_type;

        typedef typename mpl::if_<mpl::bool_<is_composite>,
                mpl::identity<boost::none_t>,
                mpl::identity<typename basis_0_type::PreCompute> >::type::type pc_type;
        typedef boost::shared_ptr<pc_type> pc_ptrtype;
        //typedef typename basis_type::polyset_type return_value_type;
        typedef typename functionspace_type::return_type return_type;

        typedef typename matrix_node<value_type>::type matrix_node_type;

        typedef typename mpl::if_<mpl::bool_<is_composite>,
                mpl::identity<boost::none_t>,
                mpl::identity<typename basis_0_type::polynomial_type> >::type::type polynomial_view_type;

        /**
         * interpolate type if available
         */
        typedef typename mpl::if_<mpl::bool_<is_modal>,
                                  mpl::identity<boost::none_t>,
                                  mpl::identity<typename basis_0_type::local_interpolant_type> >::type::type local_interpolant_type;

        typedef Element<T,Cont> this_type;
        template<int i>
        struct sub_element
        {
            //typedef typename mpl::at_c<element_vector_type,i>::type type;
            typedef typename mpl::at_c<element_vector_type,i>::type::second_type ptrtype;
            typedef typename ptrtype::element_type type;
        };
        typedef typename functionspace_type::component_functionspace_type component_functionspace_type;
        typedef typename functionspace_type::component_functionspace_ptrtype component_functionspace_ptrtype;
        typedef typename component_functionspace_type::template Element<T,typename VectorUblas<value_type>::slice::type> component_type;

        /**
         * geometry typedef
         */
        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::gm_type gm_type;
        typedef boost::shared_ptr<gm_type> gm_ptrtype;
        typedef typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;

        //@}

        /** @name Constructors, destructor
         */
        //@{

        Element();

        Element( Element const& __e );

        friend class FunctionSpace<A0,A1,A2,A3,A4>;

        Element( functionspace_ptrtype const& __functionspace,
                 std::string const& __name = "unknown",
                 size_type __start = 0,
                 ComponentType __ct = NO_COMPONENT );

        Element( functionspace_ptrtype const& __functionspace,
                 container_type const& __c,
                 std::string const& __name = "unknown",
                 size_type __start = 0,
                 ComponentType __ct = NO_COMPONENT );

        ~Element();

        void initFromSpace( functionspace_ptrtype const& __functionspace,
                            container_type const& __c );

        //@}

        /** @name Operator overloads
         */
        //@{

        Element& operator=( Element const& __e );
#if 0
        template<typename ExprT>
        Element& operator=( vf::Expr<ExprT> const& __expr )
        {
            *this = project( this->functionspace(), elements( this->functionspace()->mesh() ), __expr );
            return *this;
        }
#endif
        template<typename VectorExpr>
        Element& operator=( VectorExpr const& __v );


        /**
         * \return the container read-only
         */
        super const& container() const
        {
            return *this;
        }

        /**
         * \return the container read-write
         */
        super & container()
        {
            return *this;
        }

        /**
         * get the component of the element
         *
         *
         * @return the THECOMP-th component of the element
         */
        template<ComponentType THECOMP>
        component_type
        comp()
        {
            //return comp( typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
            return this->template comp<THECOMP>( typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
        }
        template<ComponentType THECOMP>
        component_type
        comp( mpl::bool_<true> )
        {
            BOOST_STATIC_ASSERT( THECOMP >= X && THECOMP < ( ComponentType )N_COMPONENTS );
            auto s = ublas::slice( THECOMP, N_COMPONENTS, M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( THECOMP );
            return component_type( compSpace(),
                                   typename component_type::container_type( this->vec().data().expression(), s ),
                                   __name,
                                   start()+THECOMP,
                                   THECOMP );
        }
        template<ComponentType THECOMP>
        component_type
        comp( mpl::bool_<false> )
        {
            BOOST_STATIC_ASSERT( THECOMP >= X && THECOMP < ( ComponentType )N_COMPONENTS );
            auto s = ublas::slice( THECOMP, N_COMPONENTS, M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( THECOMP );
            return component_type( compSpace(),
                                   typename component_type::container_type( ( VectorUblas<value_type>& )*this, s ),
                                   __name,
                                   start()+THECOMP,
                                   THECOMP );
        }

        /**
         * get the component of the element
         * const version
         *
         * @param i component id
         * @return the i-th component of the element
         */
        component_type
        comp( ComponentType i ) const
        {
            //return comp( i, mpl::bool_<boost::is_same<>is_composite>() );
            return comp( i, typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
        }
        component_type
        comp( ComponentType i, mpl::bool_<true> ) const
        {
            FEELPP_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( i );

            component_type c( compSpace(),
                              typename component_type::container_type( this->vec().data().expression(), s ),
                              __name,
                              start()+i,
                              i );
            return c;
        }
        component_type
        comp( ComponentType i, mpl::bool_<false> ) const
        {
            FEELPP_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( i );

            component_type c( compSpace(),
                              typename component_type::container_type( ( VectorUblas<value_type>& )*this, s ),
                              //typename component_type::container_type( this->data().expression(), r ),
                              __name,
                              start()+i,
                              i );
            return c;
        }

        /**
         * get the component of the element
         *
         * @param i component id
         * @return the i-th component of the element
         */
        component_type
        comp( ComponentType i )
        {
            //return comp( i, mpl::bool_<is_composite>() );
            return comp( i, typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
        }
        component_type
        comp( ComponentType i, mpl::bool_<true> )
        {
            FEELPP_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, M_functionspace->nDofPerComponent() );

            std::string __name = this->name() + "_" + componentToString( i );
            component_type c( compSpace(),
                              typename component_type::container_type( this->vec().data().expression(), s ),
                              __name,
                              start()+i,
                              i );
            return c;
        }
        component_type
        comp( ComponentType i, mpl::bool_<false> )
        {
            FEELPP_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, M_functionspace->nDofPerComponent() );

            std::string __name = this->name() + "_" + componentToString( i );
            component_type c( compSpace(),
                              typename component_type::container_type( ( VectorUblas<value_type>& )*this, s ),
                              __name,
                              start()+i,
                              i );
            return c;
        }

        value_type localToGlobal( size_type e, size_type l, int c ) const
        {
            size_type index=boost::get<0>( M_functionspace->dof()->localToGlobal( e, l, c ) );
            return super::operator()( index );
        }
#if 0
        value_type& localToGlobal( size_type e, size_type l, int c )
        {
            size_type index=boost::get<0>( M_functionspace->dof()->localToGlobal( e, l, c ) );
            return super::operator()( index );
        }
#endif
        value_type  operator()( size_t i ) const
        {
            return super::operator()( i );
        }
        value_type& operator()( size_t i )
        {
            return super::operator()( i );
        }
        Element& operator+=( Element const& _e )
        {
            for ( int i=0; i < _e.nLocalDof(); ++i )
                this->operator()( i ) += _e( i );

            return *this;
        }
        Element& operator-=( Element const& _e )
        {
            for ( size_type i=0; i < _e.nLocalDof(); ++i )
                this->operator()( i ) -= _e( i );

            return *this;
        }
        /**
         * update element when mesh has been changed
         */
        void operator()( MESH_CHANGES mesh_changes )
        {
            DVLOG(2) << "Update element after a change in the mesh\n";
            switch(mesh_changes)
            {
              case MESH_CHANGES_POINTS_COORDINATES:
                M_dof->rebuildDofPoints( *M_mesh );
              break;
              default:
              break;
            }
        }

        template<typename AE>
        container_type& assign( const ublas::vector_expression<AE> &ae )
        {
            return super::assign( ae );
        }
        void assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index = M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            this->operator[]( index ) = __v;
        }
        void plus_assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index = M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            this->operator[]( index ) += __v;
        }

        void assign( geoelement_type const& e, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.id() );
            for( auto const& ldof : M_functionspace->dof()->localDof( e.id() ) )
            {
                size_type index = ldof.second.index();
                this->operator[]( index ) = s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
            }
        }
        void assign( geoface_type const& e, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.element(0).id() );
            for( auto const& ldof : M_functionspace->dof()->faceLocalDof( e.id() ) )
            {
                size_type index = ldof.index();
                this->operator[]( index ) = s(ldof.localDof())*Ihloc( ldof.localDofInFace() );
            }
        }
        void plus_assign( geoelement_type const& e, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.id() );
            for( auto const& ldof : M_functionspace->dof()->localDof( e.id() ) )
            {
                size_type index = ldof.second.index();
                this->operator[]( index ) += s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
            }
        }
        void plus_assign( geoface_type const& e, local_interpolant_type const& Ihloc )
         {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.element(0).id() );
            for( auto const& ldof : M_functionspace->dof()->faceLocalDof( e.id() ) )
            {
                size_type index = ldof.index();
                this->operator[]( index ) += s(ldof.localDof())*Ihloc( ldof.localDofInFace() );
            }
        }

        //@}

        /** @name Accessors
         */
        //@{

        typedef boost::multi_array<value_type,3> array_type;
        typedef Eigen::Matrix<value_type,nComponents1,1> _id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nRealDim> _grad_type;
        typedef Eigen::Matrix<value_type,nRealDim,nRealDim> _hess_type;
        typedef Eigen::Matrix<value_type,1,1> _div_type;
        typedef Eigen::Matrix<value_type,nRealDim,1> _curl_type;
        typedef boost::multi_array<_id_type,1> id_array_type;
        typedef boost::multi_array<_grad_type,1> grad_array_type;
        typedef boost::multi_array<_hess_type,1> hess_array_type;
        typedef boost::multi_array<_div_type,1> div_array_type;
        typedef boost::multi_array<_curl_type,1> curl_array_type;
        typedef boost::multi_array<_div_type,1> comp_curl_array_type;

        /**
         * \return the map
         */
        DataMap const& map() const
        {
            return M_functionspace->map();
        }

        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh()
        {
            return M_functionspace->mesh();
        }

        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh() const
        {
            return M_functionspace->mesh();
        }

        /**
         * \return the element-wise square root
         */
        Element<T,Cont> sqrt() const
        {
            Element<T,Cont> _e( M_functionspace );

            for ( int i=0; i < _e.nLocalDof(); ++i ) _e( i )  = math::sqrt( this->operator()( i ) );

            return _e;
        }

        /**
         * \return the element-wise power to n
         */
        Element<T,Cont> pow( int n ) const
        {
            Element<T,Cont> _e( M_functionspace );

            for ( int i=0; i < _e.nLocalDof(); ++i ) _e( i )  = math::pow( this->operator()( i ),n );

            return _e;
        }

        // Only works for scalar fields
        template < typename p0_space_type >
        typename p0_space_type::element_type extremeValue( boost::shared_ptr<p0_space_type> const& P0h, std::string extreme )
        {
            // check if the mesh coming from P0h and the class elements is the same
            FEELPP_ASSERT( P0h->mesh() == this->mesh() ).error( "mesh is not the same" );
            FEELPP_ASSERT( is_scalar ).error( "only works for scalar fields" );

            typename p0_space_type::element_type p0Element( P0h );

            for ( auto elt_it = P0h->mesh()->beginElementWithProcessId(); elt_it != P0h->mesh()->endElementWithProcessId(); ++elt_it )
            {
                size_type eid = elt_it->id();

                size_type dofp0 = boost::get<0>( P0h->dof()->localToGlobal( eid, 0, 0 ) );
                std::vector<value_type> values ( functionspace_type::fe_type::nLocalDof );

                size_type dofpn = 0;

                for ( uint16_type local_id=0; local_id < functionspace_type::fe_type::nLocalDof; ++local_id )
                {
                    dofpn = boost::get<0>( this->functionSpace()->dof()->localToGlobal( eid, local_id, 0 ) );
                    values[local_id] = this->operator()( dofpn );
                }

                if ( extreme == "max" )
                    p0Element.assign( eid, 0, 0, *( std::max_element( values.begin(), values.end() ) ) );

                else
                    p0Element.assign( eid, 0, 0, *( std::min_element( values.begin(), values.end() ) ) );
            }

            return p0Element;
        }

        value_type max() const
        {
            return this->max( true );
        }
        value_type max( bool parallel ) const
        {
            return super::max( parallel );
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type max( boost::shared_ptr<p0_space_type> const& P0h )
        {
            return this->extremeValue( P0h, "max" );
        }

        value_type min() const
        {
            return this->min( true );
        }
        value_type min( bool parallel ) const
        {
            return super::min( parallel );
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type min( boost::shared_ptr<p0_space_type> const& P0h )
        {
            return this->extremeValue( P0h, "min" );
        }

        //! Interpolation at a set of points
        //@{
        /**
         * data structure that stores the interpolated values of the
         * element at a set of points
         */
        typedef Feel::detail::ID<value_type,nComponents1,nComponents2> id_type;


        /**
         * \return the extents of the interpolation of the function at
         * a set of points
         */
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        idExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nComponents1;
            //shape[2] = nComponents2;
            return shape;
        }

        template<typename Context_t>
        id_type
        id( Context_t const & context ) const
        {
            return id_type( *this, context );
        }

        /**
         * interpolate the function at a set of points(defined in the
         * reference element) in an element.
         *
         * \warning this operation assumes that the interpolation is
         * done on the mesh associated with the discretization space
         *
         * \todo handle the case where the mesh/element is not the
         * same as the one associated with the discretization space.
         *
         */
        template<typename Context_t>
        void
        id_( Context_t const & context, id_array_type& v ) const;

        template<typename Context_t>
        void
        id( Context_t const & context, id_array_type& v ) const
        {
            id_( context, v );
        }


        /*
         * evaluate the function at all points added to functionspace_type::Context
         */
        Eigen::Matrix<value_type, Eigen::Dynamic, 1>
        evaluate( functionspace_type::Context const & context , bool do_communications=true) const
        {
            const int npoints = context.nPoints();

            //number of component
            const int ncdof  = is_product?nComponents:1;

            //rank of the current processor
            int proc_number = this->worldComm().globalRank();

            //total number of processors
            int nprocs = this->worldComm().globalSize();

            auto it = context.begin();
            auto en = context.end();

            eigen_type __globalr( npoints*ncdof );
            __globalr.setZero();
            eigen_type __localr( npoints*ncdof );
            __localr.setZero();

            boost::array<typename array_type::index, 1> shape;
            shape[0] = 1;

            id_array_type v( shape );
            if( context.size() > 0 )
            {
                for( int i = 0 ; it != en; ++it, ++i )
                {
                    v[0].setZero(1);
                    auto basis = it->second;
                    id( *basis, v );
                    int global_index = it->first;
                    for(int comp=0; comp<ncdof; comp++)
                    {
                        __localr( global_index*ncdof+comp ) = v[0]( comp, 0 );
                    }
                }
            }

            if( do_communications )
                mpi::all_reduce( this->worldComm() , __localr, __globalr, std::plus< eigen_type >() );
            else
                __globalr = __localr;

            return __globalr;
        }

        value_type
        evaluate( functionspace_type::Context const & context, int i , bool do_communications=true) const
        {
            //rank of the current processor
            int proc_number = Environment::worldComm().globalRank();

            //total number of processors
            int nprocs = Environment::worldComm().globalSize();

            int npoints = context.nPoints();
            CHECK( i >= 0 && i < npoints ) << "the index " << i << " of the point where you want to evaluate the element is out of range\n";
            boost::array<typename array_type::index, 1> shape;
            shape[0] = 1;
            id_array_type v( shape );

            value_type result=0;
            int proc_having_the_point = context.processorHavingPoint( i );
            if( proc_number == proc_having_the_point )
            {
                auto basis = context.at( i );
                id( *basis , v );
                result = v[0](0,0);
            }

            if( do_communications )
                boost::mpi::broadcast( Environment::worldComm() , result , proc_having_the_point );

            return result;
        }

        void
        idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;


        /*
         * Get the reals points matrix in a context
         * 1 : Element
         *
         * \todo store a geometric mapping context to evaluate the real points
         * from a set of point in the referene element, should probably done in
         * the real element (geond)
         */
        template<typename Context_t>
        matrix_node_type
        ptsInContext( Context_t const & context, mpl::int_<1> ) const
        {
            //new context for evaluate the points
            typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
            typedef boost::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

            gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),  context.pc() ) );

            return __c_interp->xReal();
        }

        /*
         * Get the real point matrix in a context
         * 2 : Face
         * \todo see above
         */
        template<typename Context_t>
        matrix_node_type
        ptsInContext( Context_t const & context,  mpl::int_<2> ) const
        {
            //new context for the interpolation
            typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
            typedef boost::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

            typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::permutation_type permutation_type;
            typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::precompute_ptrtype precompute_ptrtype;

            //not good because ?
            //gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),context.pcFaces(), context.faceId()) );
            //good with this
            std::vector<std::map<permutation_type, precompute_ptrtype> > __geo_pcfaces = context.pcFaces();
            gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId() ) );

            return __c_interp->xReal();
        }

        /**
         * \return a geometric mapping context associated with the element
         * \todo should be in the element of the mesh
         */
        gmc_ptrtype geomapPtr( geoelement_type const& elt ) const
        {
            gm_ptrtype __gm = functionSpace()->gm();
            geopc_ptrtype __geopc( new geopc_type( __gm, elt.G() ) );
            return gmc_ptrtype( new gmc_type( __gm, elt, __geopc ) );
        }
        /**
         * \return a precomputation of the basis functions
         * \warning this seems quite buggy (where is it used?)
         */
        pc_ptrtype pcPtr( geoelement_type const& elt ) const
        {
            return pc_ptrtype( new pc_type( functionSpace()->fe(), elt.G() ) );
        }

        id_type operator()( Eigen::Matrix<value_type,nDim,1> const& __x, bool extrapolate = false ) const
            {
                node_type n( nDim );
                for(int i = 0; i < nDim; ++i ) n[i]=__x[i];
                return operator()( n, extrapolate );
            }

        /**
         * interpolate the function at node (real coordinate) x
         *
         * @return the interpolated value of the function at the real point x
         */
        id_type operator()( node_type const& __x, bool extrapolate = false ) const;

        //@}

        //! gradient interpolation tool
        //@{

        typedef Feel::detail::DD<value_type, nComponents1, nRealDim> grad_type;
        typedef Feel::detail::D<value_type,0,nComponents1, 1> dx_type;
        typedef Feel::detail::D<value_type,1,nComponents1, 1> dy_type;
        typedef Feel::detail::D<value_type,2,nComponents1, 1> dz_type;

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        gradExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nComponents1;
            //shape[2] = nRealDim;

            return shape;
        }
        template<typename ContextType>
        grad_type
        grad( ContextType const & context ) const
        {
            return grad_type( *this, context );
        }
        template<typename ElementType>
        void
        grad( ElementType const& elt, grad_array_type& v )
        {
            gmc_ptrtype gmc( geomapPtr( elt ) );
            v.resize( gradExtents(  *gmc ) );
            grad_( *gmc, *pcPtr( elt ), v );
        }

        /**
         * interpolate the function derivative in the first component
         * at a set of points(defined in the reference element) in an
         * element.
         *
         * \warning this operation assumes that the interpolation is
         * done on the mesh associated with the discretization space
         *
         * \todo handle the case where the mesh/element is not the
         * same as the one associated with the discretization space.
         *
         * \todo need to introduce the geometric transformation to
         * incorporate the pseudo-inverse of the jacobian for the derivative
         */
        template<typename ContextType>
        void grad_( ContextType const & context, grad_array_type& v ) const;

        template<typename ContextType>
        void grad( ContextType const & context, grad_array_type& v ) const
        {
            grad_( context, v );
        }

        void
        gradInterpolate( matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        /**
         * interpolate the gradient of the function at node (real coordinate) x
         *
         * @return the interpolated value of the gradient function at the real point x
         */
        grad_type grad( node_type const& __x ) const;

        template<typename ContextType>
        dx_type
        dx( ContextType const & context ) const
        {
            return dx_type( *this, context );
        }

        template<typename ContextType>
        dy_type
        dy( ContextType const & context ) const
        {
            return dy_type( *this, context );
        }
        template<typename ContextType>
        dz_type
        dz( ContextType const & context ) const
        {
            return dz_type( *this, context );
        }

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dxExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nComponents1;
            //shape[2] = nComponents2;
            return shape;
        }
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dyExtents( ContextType const & context ) const
        {
            return dxExtents( context );
        }
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dzExtents( ContextType const & context ) const
        {
            return dxExtents( context );
        }
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dExtents( ContextType const & context ) const
        {
            return dxExtents( context );
        }

        template<typename ContextType>
        void d_( int N, ContextType const & context, id_array_type& v ) const;

        template<typename ContextType>
        void dx( ContextType const & context, id_array_type& v ) const
        {
            d_( 0, context, v );
        }
        template<typename ContextType>
        void dy( ContextType const & context, id_array_type& v ) const
        {
            d_( 1, context, v );
        }
        template<typename ContextType>
        void dz( ContextType const & context, id_array_type& v ) const
        {
            d_( 2, context, v );
        }

        void
        dxInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        void
        dyInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        void
        dzInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;


        //@}

        typedef Feel::detail::Div<value_type> div_type;

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        divExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = 1;
            //shape[2] = 1;
            return shape;
        }
        template<typename ContextType>
        div_type
        div( ContextType const & context ) const
        {
            //BOOST_STATIC_ASSERT( (rank == 1) || (rank==2) );
            return div_type( *this, context );
        }


        template<typename ContextType>
        void
        div( ContextType const & context, div_array_type& v ) const
        {
            //BOOST_STATIC_ASSERT( (rank == 1) || (rank==2) );
            div_( context, v );
        }
        template<typename ContextType>
        void div_( ContextType const & context, div_array_type& v ) const;

        void
        divInterpolate( matrix_node_type __ptsReal, div_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        typedef Feel::detail::Curl<value_type,-1,nRealDim> curl_type;
        typedef Feel::detail::Curl<value_type,0,1> curlx_type;
        typedef Feel::detail::Curl<value_type,1,1> curly_type;
        typedef Feel::detail::Curl<value_type,2,1> curlz_type;

        template<typename ContextType>
        curl_type
        curl( ContextType const & context ) const
        {
            return curl_type( *this, context );
        }

        template<typename ContextType>
        curlx_type
        curlx( ContextType const & context ) const
        {
            return curlx_type( *this, context );
        }

        template<typename ContextType>
        curly_type
        curly( ContextType const & context ) const
        {
            return curly_type( *this, context );
        }

        template<typename ContextType>
        curlz_type
        curlz( ContextType const & context ) const
        {
            return curlz_type( *this, context );
        }

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        curlExtents( ContextType const & context ) const
        {
            //BOOST_MPL_ASSERT_MSG( ( rank == 1 ), INVALID_RANK_FOR_CURL, (mpl::int_<rank>) );

            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nComponents1;
            //shape[2] = 1;
            return shape;
        }
        template<typename ContextType>
        void
        curl( ContextType const & context, curl_array_type& v ) const
        {
            //BOOST_MPL_ASSERT_MSG( ( rank == 1 ), INVALID_RANK_FOR_CURL, (mpl::int_<rank>) );

            curl_( context, v );
        }

        void
        curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        template<typename ContextType>
        void curl_( ContextType const & context, curl_array_type& v ) const;

        template<typename ContextType>
        void curl_( ContextType const & context, comp_curl_array_type& v, int comp ) const;


        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        curlxExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nComponents1;
            //shape[2] = nComponents2;
            return shape;
        }
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        curlyExtents( ContextType const & context ) const
        {
            return curlxExtents( context );
        }
        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        curlzExtents( ContextType const & context ) const
        {
            return curlxExtents( context );
        }

        template<typename ContextType>
        void
        curlx( ContextType const & context, comp_curl_array_type& v ) const
        {
            //BOOST_STATIC_ASSERT( rank == 1 );
            curl_( context, v, 0 );
        }
        template<typename ContextType>
        void
        curly( ContextType const & context, comp_curl_array_type& v ) const
        {
            //BOOST_STATIC_ASSERT( rank == 1 );
            curl_( context, v, 1 );
        }
        template<typename ContextType>
        void
        curlz( ContextType const & context, comp_curl_array_type& v ) const
        {
            //BOOST_STATIC_ASSERT( rank == 1 );
            curl_( context, v, 2 );
        }

        void
        curlxInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        void
        curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        void
        curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        hessExtents( ContextType const & context ) const
        {
            BOOST_STATIC_ASSERT( rank == 0 );

            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            //shape[1] = nRealDim;
            //shape[2] = nRealDim;
            return shape;
        }
        /**
         * interpolate the function second order derivatives at a set
         * of points(defined in the reference element) in an element.
         *
         * \warning this operation assumes that the interpolation is
         * done on the mesh associated with the discretization space
         *
         * \todo handle the case where the mesh/element is not the
         * same as the one associated with the discretization space.
         *
         * \todo need to introduce the geometric transformation to
         * incorporate the pseudo-inverse of the jacobian for the derivative
         */
        template<typename ContextType>
        void
        hess_( ContextType const & context, hess_array_type& v ) const;

        template<typename ContextType>
        void
        hess_( ContextType const & context, hess_array_type& v, mpl::int_<0> ) const;

        void
        hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        typedef Feel::detail::H<value_type,nRealDim,nRealDim> hess_type;

        template<typename ContextType>
        hess_type
        hess( ContextType const & context ) const
        {
            return hess_type( *this, context );
        }

        template<typename ContextType>
        void
        hess( ContextType const & context, hess_array_type& v ) const
        {
            hess_( context, v );
        }
        template<typename ElementType>
        void
        hess( ElementType const& elt, hess_array_type& v )
        {
            gmc_ptrtype gmc( geomapPtr( elt ) );
            v.resize( hessExtents(  *gmc ) );
            hess_( *gmc, *pcPtr( elt ), v );
        }

        /**
           \return the finite element space
        */
        functionspace_ptrtype const& functionSpace() const
        {
            return M_functionspace;
        }

        /**
         * get the \c FunctionSpace vector
         */
        //functionspace_vector_type const&
        //functionSpaces() const { return M_functionspace->functionSpaces(); }

        functionspace_vector_type const&
        functionSpaces() const
        {
            return M_functionspace->functionSpaces();
        }

        /**
         * get the \p i -th \c FunctionSpace out the list
         */
        template<int i>
        typename mpl::at_c<functionspace_vector_type,i>::type
        functionSpace() const
        {
            return M_functionspace->template functionSpace<i>();
        }

        template<typename BackendType>
        typename BackendType::vector_ptrtype
        extractValuesWithMarker( std::string const& m,
                                 boost::shared_ptr<BackendType> backend )
            {
                auto r =  functionSpace()->dof()->markerToDof( m );
                size_type s = std::distance( r.first, r.second );
                vector_ptrtype res = backend->newVector(s, s);
                size_type i = 0;
                for( auto it = r.first, en = r.second; it != en; ++it, ++i )
                    res->set( i, this->operator[]( it->second ) );
                return res;
            }

        template<typename BackendType>
        typename BackendType::vector_ptrtype
        extractValuesWithoutMarker( std::string const& m,
                                    boost::shared_ptr<BackendType> backend )
            {
                auto r1 =  functionSpace()->dof()->markerToDofLessThan( m );
                auto r2 =  functionSpace()->dof()->markerToDofGreaterThan( m );
                size_type s = std::distance( r1.first, r1.second ) + std::distance( r2.first, r2.second );
                vector_ptrtype res = backend->newVector(s, s);
                size_type i = 0;
                for( auto it = r1.first, en = r1.second; it != en; ++it, ++i )
                    res->set( i, this->operator[]( it->second ) );
                for( auto it = r2.first, en = r2.second; it != en; ++it, ++i )
                    res->set( i, this->operator[]( it->second ) );
                return res;
            }


        template<int i>
        typename sub_element<i>::type
        elementImpl( std::string const& name ="u", bool updateOffViews=true );

        template<int i,typename ExprT>
        typename sub_element<i>::type &
        element( ExprT e, std::string const& name = "u",
                 bool updateOffViews = true,
                 typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 );

        template<int i>
        typename sub_element<i>::type
        elementImpl( std::string const& name ="u", bool updateOffViews=true ) const;

        template<int i>
        typename sub_element<i>::type &
        element( std::string const& name ="u", bool updateOffViews=true )
        {
            CHECK ( fusion::at_c<i>( M_elements ).second ) << " has not element \n";
            return *(fusion::at_c<i>( M_elements ).second);
        }
        template<int i>
        typename sub_element<i>::type const&
        element( std::string const& name ="u", bool updateOffViews=true ) const
        {
            CHECK ( fusion::at_c<i>( M_elements ).second ) << " has not element \n";
            return *(fusion::at_c<i>( M_elements ).second);
        }
        template<int i>
        typename sub_element<i>::ptrtype &
        elementPtr()
        {
            CHECK ( fusion::at_c<i>( M_elements ).second ) << " has not element \n";
            return fusion::at_c<i>( M_elements ).second;
        }
        template<int i>
        typename sub_element<i>::ptrtype const&
        elementPtr() const
        {
            CHECK ( fusion::at_c<i>( M_elements ).second ) << " has not element \n";
            return fusion::at_c<i>( M_elements ).second;
        }

        /**
         *
         * @return the finite element space associated with the n-th component
         */
        component_functionspace_ptrtype const& compSpace() const
        {
            return M_functionspace->compSpace();
        }

        /**
         * subWorlds : vector of WorldComm ( >1 if composite)
         */
        std::vector<WorldComm> const& worldsComm() const
        {
            return M_functionspace->worldsComm();
        }

        /**
         * world communicator
         */
        WorldComm const& worldComm() const
        {
            return M_functionspace->worldComm();
        }

        /**
         * \return the number of dof
         */
        size_type nDof() const
        {
            return M_functionspace->nDof();
        }

        /**
         * \return the number of local dof (dof per processor)
         */
        size_type nLocalDof() const
        {
            return M_functionspace->nLocalDof();
        }

        /**
           \return the number of dof per component
        */
        size_type nDofPerComponent() const
        {
            return M_functionspace->nDofPerComponent();
        }

        /**
           \return the name of the element
        */
        std::string const& name() const
        {
            return M_name;
        }

        size_type start() const
        {
            return M_start;
        }

        bool isAComponent() const
        {
            return M_ct >= X && M_ct <= Z;
        }

        ComponentType component() const
        {
            if (  M_ct < X || M_ct > Z )
            {
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << M_ct;
                throw std::logic_error( __str.str() );
            }

            return M_ct;
        }

        std::string componentToString( ComponentType ct ) const
        {
            if (  ct < X || ct > Z )
            {
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << ct;
                throw std::logic_error( __str.str() );
            }

            switch ( ct )
            {
            case X:
                return "X";

            case Y:
                return "Y";

            case Z:
                return "Z";

            default:
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << ct;
                throw std::logic_error( __str.str() );
            }
        }

        //@}

        /** @name Mutators
         */
        //@{

        /**
         * set the name of the field
         * @param __name name of the field
         */
        void setName( std::string const & __name )
        {
            M_name = __name;
        }

        void setFunctionSpace( functionspace_ptrtype space )
        {
            M_functionspace = space;
            super::init( M_functionspace->dof() );
            //super::init( M_functionspace->nDof(),  M_functionspace->nLocalDof() );
        }

        BOOST_PARAMETER_MEMBER_FUNCTION( ( void ),
                                         save,
                                         tag,
                                         ( required
                                           ( path,* ) )
                                         ( optional
                                           ( type,( std::string ),std::string( "binary" ) )
                                           ( suffix,( std::string ),std::string( "" ) )
                                           ( sep,( std::string ),std::string( "" ) )
                                         ) )
        {
            Feel::detail::ignore_unused_variable_warning( args );

            if ( !fs::exists( fs::path( path ) ) )
            {
                fs::create_directories( fs::path( path ) );
            }

            std::ostringstream os1;
            os1 << M_name << sep << suffix << "-" << this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
            fs::path p = fs::path( path ) / os1.str();
            LOG(INFO) << "saving "  << p << "\n";
            fs::ofstream ofs( p );

            if ( type == "binary" )
            {
                boost::archive::binary_oarchive oa( ofs );
                oa << *this;
            }

            else if ( type == "text" )
            {
                boost::archive::text_oarchive oa( ofs );
                oa << *this;
            }

            else if ( type == "xml" )
            {
                //boost::archive::xml_oarchive oa(ofs);
                //oa << *this;
            }
        }
        BOOST_PARAMETER_MEMBER_FUNCTION(
            ( bool ),
            load,
            tag,
            ( required
              ( path,* ) )
            ( optional
              ( type,( std::string ),std::string( "binary" ) )
              ( suffix,( std::string ),std::string( "" ) )
              ( sep,( std::string ),std::string( "" ) )
            )
        )
        {
            Feel::detail::ignore_unused_variable_warning( args );
            std::ostringstream os1;
            os1 << M_name << sep << suffix << "-" <<  this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
            fs::path p = fs::path( path ) / os1.str();
            fs::path partial_path = fs::path(path);

            fs::path full_path_dir_sol(fs::current_path());
            full_path_dir_sol = full_path_dir_sol/partial_path;
            //std::cout << " In load the first full path is " << p << std::endl;
            if ( !fs::exists( p ) )
            {
                std::ostringstream os2;
                os2 << M_name << sep << suffix<< "-" <<  this->worldComm().globalSize() << "." << this->worldComm().globalRank();
                p = fs::path( path ) / os2.str();

                if ( !fs::exists( p ) )
                {
                    LOG(WARNING)  << "[load] :" <<  full_path_dir_sol << "  FILE : " << os1.str() << " OR " << os2.str() << " DO NOT EXIST" << std::endl ;
                    //std::cerr << "ATTENTION :  p does not exist
                    return 0;
                }
            }
            LOG(INFO) << p << " exists, is is a regular file : " << fs::is_regular_file( p ) << "\n";
            if ( !fs::is_regular_file( p ) )
            {

                LOG(WARNING) << "[load] : " << full_path_dir_sol << p << " is not a  regular_file !" << std::endl;
                return 0;
            }

            fs::ifstream ifs( p );

            if ( type == "binary" )
            {
                boost::archive::binary_iarchive ia( ifs );
                ia >> *this;
            }

            else if ( type == "text" )
            {
                boost::archive::text_iarchive ia( ifs );
                ia >> *this;
            }

            else if ( type == "xml" )
            {
                //boost::archive::xml_iarchive ia(ifs);
                //ia >> *this;
                return false;
            }
            return true;
        }

        void printMatlab( std::string fname, bool gmsh = false ) const
            {
                container_type m( *this );
                if ( gmsh )
                {
                    auto relation = this->functionSpace()->dof()->pointIdToDofRelation();
                    for( size_type i = 0; i < this->localSize(); ++i )
                    {
                        m[relation.first[i]-1] = this->operator[]( i );
                    }
                }
                m.printMatlab( fname );
            }


        BOOST_PARAMETER_MEMBER_FUNCTION( (void),
                                         on,
                                         tag,
                                         ( required
                                           ( range, *  )
                                           ( expr,   * )
                                             ) // 4. one required parameter, and

                                         ( optional
                                           ( prefix,   ( std::string ), "" )
                                           ( geomap,         *, GeomapStrategyType::GEOMAP_OPT )
                                           ( accumulate,     *( boost::is_integral<mpl::_> ), false )
                                           ( verbose,   ( bool ), boption(_prefix=prefix,_name="on.verbose") )))
            {
                return onImpl( range, expr, prefix, geomap, accumulate, verbose );
            }


        //@}
    private:

        /*
         *
         */
        void initSubElementView( mpl::true_ );

        void initSubElementView( mpl::false_ ) {}


        friend class boost::serialization::access;

        template<class Archive>
        void serialize( Archive & ar, const unsigned int version )
        {
            //ar & BOOST_SERIALIZATION_NVP( boost::serialization::base_object<super>(*this) );
            ar & boost::serialization::make_nvp( "name", M_name );
            DVLOG(2) << "got name " << M_name << "\n";

            if ( Archive::is_saving::value )
            {
                //std::cout << "saving in version " << version << "\n";
                size_type s = this->functionSpace()->nLocalDofWithGhost();
                ar & boost::serialization::make_nvp( "size", s );

                std::vector<int> no = M_functionspace->basisOrder();
                //for( int i = 0; i < no.size(); ++i ) std::cout << no[i] << std::endl;
                std::string family = M_functionspace->basisName();
                //std::cout << "family name = " << family << std::endl;

                ar & boost::serialization::make_nvp( "order",  no );
                //std::cout << "saving order done" << std::endl;
                ar & boost::serialization::make_nvp( "family", family );
                //std::cout << "saving family done" << std::endl;

                typename container_type::const_iterator it = this->begin();
                typename container_type::const_iterator en = this->end();

                for ( size_type i = 0; it != en; ++it, ++i )
                {
                    T value = this->operator[]( i );
                    std::ostringstream v_str;
                    v_str << "value_" << i;

                    //                     DVLOG(2) << "save value " << value << " at " << v_str.str() << "\n";

                    ar & boost::serialization::make_nvp( v_str.str().c_str(), value );
                }
            }

            if ( Archive::is_loading::value )
            {
                //std::cout << "loading in version " << version << "\n";

                size_type s( 0 );
                ar & boost::serialization::make_nvp( "size", s );

                // verify number of degree of freedom
                DVLOG(2) << "loading ublas::vector of size " << s << "\n";

                if ( s != this->functionSpace()->nLocalDofWithGhost() )
                    throw std::logic_error( ( boost::format( "load function: invalid number of degrees of freedom, read %1% but has %2%" ) % s % this->functionSpace()->nLocalDofWithGhost() ).str() );


                std::vector<int> order;
                std::string family;
                ar & boost::serialization::make_nvp( "order", order );
                //for( int i = 0; i < order.size(); ++i ) std::cout << order[i] << std::endl;

                ar & boost::serialization::make_nvp( "family", family );
                //std::cout << "family name = " << family << std::endl;
#if 0
                auto orders = M_functionspace->basisOrder();

                if ( order !=  orders )
                    throw std::logic_error( ( boost::format( "load function: invalid polynomial order, read %1% but has %2%" ) % order % orders ).str() );

                std::string bname = M_functionspace->basisName();

                if ( family !=  bname )
                    throw std::logic_error( ( boost::format( "load function: invalid polynomial family, read %1% but has %2%" ) % family % bname ).str() );

#endif

                for ( size_type i = 0; i < s ; ++i )
                {
                    value_type value(  0 );
                    std::ostringstream v_str;
                    v_str << "value_" << i;
                    //                     DVLOG(2) << "load value at " << v_str.str() << "\n";
                    ar & boost::serialization::make_nvp( v_str.str().c_str(), value );
                    //                     DVLOG(2) << "got value " << value << " at index " << i << "\n";
                    this->operator[]( i ) = value;
                }
            }


        }
    private:

        template<typename RangeType, typename ExprType>
        void onImpl( RangeType const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate = true, bool verbose = false )
        {
            onImplBase( r, e, prefix, geomap_strategy, accumulate, verbose, boost::is_std_list<RangeType>()  );
        }

        template<typename RangeType, typename ExprType>
        void onImplBase( RangeType const& rList, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::true_ )
        {
            const int iDim = boost::tuples::template element<0, typename RangeType::value_type>::type::value;
            for ( auto const& r : rList )
                onImpl( std::make_pair( r.template get<1>(), r.template get<2>()), e, prefix, geomap_strategy, accumulate, verbose, mpl::int_<iDim>() );
        }
        template<typename RangeType, typename ExprType>
        void onImplBase( RangeType const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::false_ )
        {
            const int iDim = boost::tuples::template element<0, RangeType>::type::value;
            onImpl( std::make_pair( r.template get<1>(), r.template get<2>()), e, prefix, geomap_strategy, accumulate, verbose, mpl::int_<iDim>() );
        }

        template<typename IteratorType, typename ExprType>
        void onImpl( std::pair<IteratorType,IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_ELEMENTS>  );
        template<typename IteratorType, typename ExprType>
        void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_FACES>  );

    private:

        /**
           Finite Element space
        */
        functionspace_ptrtype M_functionspace;

        std::string M_name;

        size_type M_start;

        //! type of the component
        ComponentType M_ct;

        // only init in // with composite case : ex p = U.element<1>()
        mutable boost::optional<container_vector_type> M_containersOffProcess;

        // views on subelement
        element_vector_type M_elements;


    }; // Element

    //@}
    /** @name Typedefs
     */
    //@{
    typedef Element<value_type> element_type;
    typedef Element<value_type> real_element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef std::map< size_type, std::vector< size_type > > proc_dist_map_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * \c FunctionSpace constructor
     *
     * \param mesh a mesh data structure
     */
    FunctionSpace( mesh_ptrtype const& mesh,
                   size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
                   periodicity_type  periodicity = periodicity_type(),
                   std::vector<WorldComm> const& _worldsComm = Environment::worldsComm(nSpaces),
                   std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false) )
        :
        M_worldsComm( _worldsComm ),
        M_worldComm( new WorldComm( _worldsComm[0] ) ),
        M_extendedDofTableComposite( extendedDofTable ),
        M_extendedDofTable( extendedDofTable[0] )
    {
        this->init( mesh, mesh_components, periodicity );
    }

    FunctionSpace( mesh_ptrtype const& mesh,
                   std::vector<Dof > const& dofindices,
                   periodicity_type periodicity = periodicity_type(),
                   std::vector<WorldComm> const& _worldsComm = Environment::worldsComm(nSpaces),
                   std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false) )
        :
        M_worldsComm( _worldsComm ),
        M_worldComm( new WorldComm( _worldsComm[0] ) ),
        M_extendedDofTableComposite( extendedDofTable ),
        M_extendedDofTable( extendedDofTable[0] )
    {
        this->init( mesh, 0, dofindices, periodicity );
    }

    FunctionSpace() {}
    // template<typename... FSpaceList>
    // FunctionSpace( FSpaceList... space_list )
    //     :
    //     M_functionspaces( fusion::make_vector( space_list... ) )
    // {
    //     this->initList( space_list... );
    // }

    /**
     * helper static function to create a boost::shared_ptr<> out of
     * the \c FunctionSpace
     */
#if 0 // ambiguous call with new below
    static pointer_type New( mesh_ptrtype const& __m, size_type mesh_components = MESH_RENUMBER | MESH_CHECK )
    {
        return pointer_type( new functionspace_type( __m, mesh_components ) );
    }
#endif // 0

    static pointer_type New( mesh_ptrtype const& __m, std::vector<Dof > const& dofindices )
    {
        return pointer_type( new functionspace_type( __m, dofindices ) );
    }
    BOOST_PARAMETER_MEMBER_FUNCTION( ( pointer_type ),
                                     static New,
                                     tag,
                                     ( required
                                       ( mesh,* )
                                     )
                                     ( optional
                                       ( worldscomm, *, Environment::worldsComm(nSpaces) )
                                       ( components, ( size_type ), MESH_RENUMBER | MESH_CHECK )
                                       ( periodicity,*,periodicity_type() )
                                       ( extended_doftable,*,std::vector<bool>(nSpaces,false) )
                                     )
                                   )
    {
        return NewImpl( mesh, worldscomm, components, periodicity, extended_doftable );
    }

    static pointer_type NewImpl( mesh_ptrtype const& __m,
                                 std::vector<WorldComm> const& worldsComm = Environment::worldsComm(nSpaces),
                                 size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
                                 periodicity_type periodicity = periodicity_type(),
                                 std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false) )
    {

        return pointer_type( new functionspace_type( __m, mesh_components, periodicity, worldsComm, extendedDofTable ) );
    }

    template<typename ...FSpaceList>
    static pointer_type
    NewFromList( FSpaceList... fspacelist )
        {
            auto X = pointer_type( new functionspace_type );
            X->M_functionspaces = fusion::make_vector(fspacelist...);
            X->initList( fspacelist... );
            return X;
        }


    /**
     * initialize the function space
     */
    void init( mesh_ptrtype const& mesh,
               size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
               periodicity_type periodicity = periodicity_type() )
    {
        Feel::Context ctx( mesh_components );
        DVLOG(2) << "component     MESH_RENUMBER: " <<  ctx.test( MESH_RENUMBER ) << "\n";
        DVLOG(2) << "component MESH_UPDATE_EDGES: " <<  ctx.test( MESH_UPDATE_EDGES ) << "\n";
        DVLOG(2) << "component MESH_UPDATE_FACES: " <<  ctx.test( MESH_UPDATE_FACES ) << "\n";
        DVLOG(2) << "component    MESH_PARTITION: " <<  ctx.test( MESH_PARTITION ) << "\n";

        this->init( mesh, mesh_components, periodicity, std::vector<Dof >(), mpl::bool_<is_composite>() );
        //mesh->addObserver( *this );
    }

    void init( mesh_ptrtype const& mesh,
               size_type mesh_components,
               std::vector<Dof > const& dofindices,
               periodicity_type periodicity = periodicity_type() )
    {

        Feel::Context ctx( mesh_components );
        DVLOG(2) << "component     MESH_RENUMBER: " <<  ctx.test( MESH_RENUMBER ) << "\n";
        DVLOG(2) << "component MESH_UPDATE_EDGES: " <<  ctx.test( MESH_UPDATE_EDGES ) << "\n";
        DVLOG(2) << "component MESH_UPDATE_FACES: " <<  ctx.test( MESH_UPDATE_FACES ) << "\n";
        DVLOG(2) << "component    MESH_PARTITION: " <<  ctx.test( MESH_PARTITION ) << "\n";

        this->init( mesh, mesh_components, periodicity, dofindices, mpl::bool_<is_composite>() );
        //mesh->addObserver( *this );
#if !defined(__INTEL_COMPILER ) 
        if(boption( "connect"))
          mesh->addObserver( *this );
#endif
    }

    //! destructor: do nothing thanks to shared_ptr<>
    ~FunctionSpace()
        {
            VLOG(1) << "FunctionSpace Destructor...";
            M_dof.reset();
            CHECK( M_dof.use_count() == 0 ) << "Invalid Dof Table shared_ptr";
            M_dofOnOff.reset();
            CHECK( M_dofOnOff.use_count() == 0 ) << "Invalid Dof OnOff shared_ptr";
            M_ref_fe.reset();
            CHECK( M_ref_fe.use_count() == 0 ) << "Invalid reffe shared_ptr";
        }


    void setWorldsComm( std::vector<WorldComm> const& _worldsComm )
    {
        M_worldsComm=_worldsComm;
    }
    void setWorldComm( WorldComm const& _worldComm )
    {
        M_worldComm.reset( new WorldComm( _worldComm ) );
    }

    std::vector<WorldComm> const& worldsComm() const
    {
        return M_worldsComm;
    }
    WorldComm const& worldComm() const
    {
        return *M_worldComm;
    }

    bool hasEntriesForAllSpaces()
    {
        return (this->template mesh<0>()->worldComm().localSize() == this->template mesh<0>()->worldComm().globalSize() );
    }
    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * update element when mesh has been changed
     */
    void operator()( MESH_CHANGES mesh_changes )
    {
        DVLOG(2) << "Update function space after a change in the mesh\n";
        switch(mesh_changes)
        {
          case MESH_CHANGES_POINTS_COORDINATES:
            rebuildDofPoints( );
          break;
          default:
          break;
        }
    }


    //@}

    /** @name Accessors
     */
    //@{


    /*
     * Get the real point matrix in a context
     * 2 : Face
     * \todo see above
     */
    template<typename Context_t>
    typename matrix_node<value_type>::type
    ptsInContext( Context_t const & context,  mpl::int_<2> ) const
    {
        //new context for the interpolation
        typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
        typedef boost::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

        typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::permutation_type permutation_type;
        typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::precompute_ptrtype precompute_ptrtype;

        //not good because ?
        //gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),context.pcFaces(), context.faceId()) );
        //good with this
        std::vector<std::map<permutation_type, precompute_ptrtype> > __geo_pcfaces = context.pcFaces();
        gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId() ) );

        return __c_interp->xReal();
    }


    /**
     * return 1 if scalar space or the number of components if vectorial space
     */
    uint16_type qDim() const
    {
        return N_COMPONENTS;
    }

    /**
     * \return the number of degrees of freedom for each space overall the entire domain
     */
    size_type nDof() const
    {
        return this->nDof( mpl::bool_<is_composite>() );
    }

    /**
     * \return the number of degrees of freedom for each space on the current subdomain
     */
    size_type nLocalDof() const
    {
        return this->nLocalDof( mpl::bool_<is_composite>() );
    }

    size_type nLocalDofWithGhost() const
    {
        return this->nLocalDofWithGhost( mpl::bool_<is_composite>() );
    }

    size_type nLocalDofWithoutGhost() const
    {
        return this->nLocalDofWithoutGhost( mpl::bool_<is_composite>() );
    }

    size_type nLocalDofWithGhostOnProc( const int proc ) const
    {
        return this->nLocalDofWithGhostOnProc( proc, mpl::bool_<is_composite>() );
    }

    size_type nLocalDofWithoutGhostOnProc(const int proc) const
    {
        return this->nLocalDofWithoutGhostOnProc( proc, mpl::bool_<is_composite>() );
    }

    /**
     * \return the distribution of the dofs among the processors
     */
    proc_dist_map_type getProcDistMap() const
    {
        return procDistMap;
    }

    /**
     * \return the starting value of the global dof numbering
     */
    size_type nDofStart( size_type i = /*invalid_size_type_value*/0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NbDof( 0, i ) );
        return start;
    }

    size_type nLocalDofStart( size_type i = 0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(),true,0,i ) );
        return start;
    }

    size_type nLocalDofWithGhostStart( size_type i = 0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(),true,0,i ) );
        return start;
    }

    size_type nLocalDofWithoutGhostStart( size_type i = 0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<false> >( this->worldsComm(),true,0,i ) );
        return start;
    }

    size_type nLocalDofWithGhostOnProcStart( const int proc, size_type i = 0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<true> >( proc, this->worldsComm(),true,0,i ) );
        return start;
    }

    size_type nLocalDofWithoutGhostOnProcStart( const int proc, size_type i = 0 ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<false> >( proc, this->worldsComm(),true,0,i ) );
        return start;
    }

    uint16_type nSubFunctionSpace() const
    {
        return nSpaces;
    }

    /**
     * \return the number of degrees of freedom per dim
     */
    size_type nDofPerComponent() const
    {
        return this->nDof()/qDim();
    }


    /**
       \return the mesh associated with the approximation
    */
    mesh_ptrtype& mesh()
    {
        return M_mesh;
    }

    /**
       \return the mesh associated with the approximation
    */
    mesh_ptrtype const& mesh() const
    {
        return M_mesh;
    }

    template<int i>
    typename GetMesh<mesh_ptrtype,i>::type
    mesh() const
    {
        return mesh<i>( is_shared_ptr<mesh_ptrtype>() );
    }
    template<int i>
    typename GetMesh<mesh_ptrtype,i>::type
    mesh( mpl::bool_<true> ) const
    {
        return M_mesh;
    }
    template<int i>
    typename GetMesh<mesh_ptrtype,i>::type
    mesh( mpl::bool_<false> ) const
    {
        return fusion::at_c<i>(M_mesh);
    }

    /**
       \return the reference finite element
    */
    basis_ptrtype const& basis() const
    {
        DCHECK( M_ref_fe ) << "Invalid reference element\n";
        return M_ref_fe;
    }

    /**
     * \return the basis name
     *
     * \note in the case of product of space, this is the concatenation of the
     * basis name of each function space
     */
    std::string basisName() const
    {
        return basisName( mpl::bool_<( nSpaces>1 )>() );
    }
    std::string basisName( mpl::bool_<true> ) const
    {
        return  fusion::accumulate( this->functionSpaces(), std::string(), Feel::detail::BasisName() );
    }
    std::string basisName( mpl::bool_<false> ) const
    {
        return this->basis()->familyName();
    }

    /**
     * \return the basis order
     *
     * \note in the case of product of space, this is the concatenation of the
     * orders in a vector
     */
    std::vector<int> basisOrder() const
    {
        return basisOrder( mpl::bool_<( nSpaces>1 )>() );
    }
    std::vector<int> basisOrder( mpl::bool_<true> ) const
    {
        return  fusion::accumulate( this->functionSpaces(), std::vector<int>(), Feel::detail::BasisOrder() );
    }
    std::vector<int> basisOrder( mpl::bool_<false> ) const
    {
        std::vector<int> o( 1 );
        o[0]=basis_type::nOrder;
        return o;
    }

    /**
       \return the reference finite element
    */
    reference_element_ptrtype const& fe() const
    {
        DCHECK( M_ref_fe ) << "Invalid reference element\n";
        return M_ref_fe;
    }

    /**
       \return the geometric mapping
    */
    gm_ptrtype const& gm() const
    {
        return M_mesh->gm();
    }

    /**
       \return the geometric mapping
    */
    gm1_ptrtype const& gm1() const
    {
        return M_mesh->gm1();
    }

    /**
       \return the degrees of freedom
    */
    DataMap const& map() const
    {
        return *M_dof;
    }
    /**
     \return the degrees of freedom
     */
    datamap_ptrtype mapPtr() const
        {
            return M_dof;
        }

    /**
       \return the degrees of freedom
    */
    DataMap const& mapOn() const
    {
        return *M_dof;
    }

    /**
       \return the degrees of freedom
    */
    DataMap const& mapOnOff() const
    {
        return *M_dofOnOff;
    }

    /**
       \return the degrees of freedom
    */
    dof_ptrtype dof()
    {
        return M_dof;
    }

    /**
       \return the degrees of freedom (ON processor)
    */
    dof_ptrtype dofOn()
    {
        return M_dof;
    }

    /**
       \return the degrees of freedom (ON and OFF processor)
    */
    dof_ptrtype dofOnOff()
    {
        return M_dofOnOff;
    }

    /**
       \return the degrees of freedom
    */
    dof_ptrtype const& dof() const
    {
        return M_dof;
    }

    /**
       \return the degrees of freedom (ON processor)
    */
    dof_ptrtype const& dofOn() const
    {
        return M_dof;
    }

    /**
       \return the degrees of freedom (ON and OFF processor)
    */
    dof_ptrtype const& dofOnOff() const
    {
        return M_dofOnOff;
    }

    /**
     * \return true if need to build extended DofTable
     */
    std::vector<bool> extendedDofTableComposite() const { return M_extendedDofTableComposite; }
    bool extendedDofTable() const { return M_extendedDofTable; }


    //! \return true if mortar, false otherwise
    constexpr bool isMortar() const { return is_mortar; }

    /**
     * get the \c FunctionSpace vector
     */
    //functionspace_vector_type const&
    //functionSpaces() const { return M_functionspaces; }

    functionspace_vector_type const&
    functionSpaces() const
    {
        return M_functionspaces;
    }


    boost::shared_ptr<IndexSplit>
    buildDofIndexSplit()
    {
        auto startSplit = boost::fusion::fold( functionSpaces(), boost::make_tuple(0,0), Feel::detail::computeStartOfFieldSplit() ).template get<1>();
        auto computeSplit = Feel::detail::computeNDofForEachSpace(startSplit);
        boost::fusion::for_each( functionSpaces(), computeSplit );
        return computeSplit.indexSplit();

    }

    boost::shared_ptr<IndexSplit> const&
    dofIndexSplit() const
    {
        return this->dof()->indexSplit();
    }


    /**
     * \return an element of the function space
     */
    element_type
    element( std::string const& name = "u" )
    {
        element_type u( this->shared_from_this(), name );
        u.zero();
        return u;
    }

    /**
     * \param e expression to initialize the element
     * \param u name of the element
     * \return an element initialized with expression \p e
     */
    template<typename ExprT>
    element_type
    element( ExprT e, std::string const& name = "u", typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 )
    {
        element_type u( this->shared_from_this(), name );
        bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
        EntityProcessType entityProcess = (addExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        u.on( _range=elements(M_mesh,entityProcess), _expr=e );
        return u;
    }

    /**
     * \return an element of the function space
     */
    element_ptrtype
    elementPtr( std::string const& name = "u" )
    {
        element_ptrtype u( new element_type( this->shared_from_this(), name ) );
        u->zero();
        return u;
    }

    /**
     * \param e expression to initialize the element
     * \param u name of the element
     * \return a pointer to an element initialized with expression \p e
     */
    template<typename ExprT>
    element_ptrtype
    elementPtr( ExprT e, std::string const& name = "u", typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 )
    {
        element_ptrtype u( new element_type( this->shared_from_this(), name ) );
        bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
        u->on( _range=elements(M_mesh,addExtendedElt), _expr=e );
        return u;
    }

    typedef std::vector<element_ptrtype> elements_ptrtype;

    /**
     * \param u name of the elements
     * \return a pointer to an element initialized with expression \p e
     */
    elements_ptrtype
    elementsPtr( int r, std::string const& name = "u" )
    {
        elements_ptrtype u( r );
        for( int i = 0; i <  r; ++i )
        {
            u[i] = element_ptrtype(new element_type( this->shared_from_this(), name ));
            u[i]->zero();
        }
        return u;
    }
    /**
     * \param e expression to initialize the element
     * \param u name of the element
     * \return a pointer to an element initialized with expression \p e
     */
    template<typename ExprT>
    element_ptrtype
    elementPtr( int r, ExprT e, std::string const& name = "u", typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 )
    {
        elements_ptrtype u( r );
        for( int i = 0; i <  r; ++i )
        {
            u[i] = element_ptrtype(new element_type( this->shared_from_this(), name ));
            bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
            u[i]->on( _range=elements(M_mesh,addExtendedElt), _expr=e );
        }
        return u;
    }

    /**
     * get the \p i -th \c FunctionSpace out the list
     */
    template<int i>
    typename mpl::at_c<functionspace_vector_type,i>::type
    functionSpace()
    {
        return fusion::at_c<i>( M_functionspaces );
    }

    /**
     *
     * @return the finite element space associated with the n-th component
     */
    component_functionspace_ptrtype const& compSpace() const
    {
        buildComponentSpace();
        return M_comp_space;
    }

    /**
     * \return trace space
     */
    trace_functionspace_ptrtype
    trace()  const
    {
        //return trace( mpl::greater<mpl::int_<nDim>,mpl::int_<1> >::type() )
        return trace_functionspace_type::New( mesh()->trace( boundaryfaces( mesh() ) ) );
    }
    template<typename RangeT>
    trace_functionspace_ptrtype
    trace( RangeT range  )  const
    {
        return trace_functionspace_type::New( mesh()->trace( range ) );
    }


    /**
     * \return trace trace space
     */
    trace_trace_functionspace_ptrtype
    wireBasket()  const
    {
        //return trace( mpl::greater<mpl::int_<nDim>,mpl::int_<1> >::type() )
        return trace_trace_functionspace_type::New( _mesh=mesh()->wireBasket( markededges( mesh(),"WireBasket" ) ), _worldscomm=this->worldsComm() );
    }
    template<typename RangeT>
    trace_trace_functionspace_ptrtype
    wireBasket( RangeT range  )  const
    {
        return trace_trace_functionspace_type::New( mesh()->wireBasket( range ) );
    }


    /**
       \return true if Space has a region tree to localize points
    */
    bool hasRegionTree() const
    {
        return M_rt != boost::none;
    }

    /**
       \return the Region Tree
    */
    region_tree_ptrtype const& regionTree() const;

    /**
       \return update the Region Tree
    */
    void updateRegionTree() const;

    /**
       find a point(real coordinate) \c pt in the mesh associated with
       the space and provide: (i) the reference coordinates \c ptr and
       (ii) the index of the element containing it.

       \param pt point to be looked for (real coordinates)
       \param ptr reference coordinates of the point
       \param cv index of the element containing the point \c pt

    */
    bool findPoint( node_type const& pt, size_type &cv, node_type& ptr ) const;

    /**
     * build component space in vectorial case
     */
    void buildComponentSpace() const;

    /**
     * rebuild dof points after a mesh mover for example
     */
    void rebuildDofPoints() { rebuildDofPoints(mpl::bool_<is_composite>()); }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void printInfo() const
    {
        LOG(INFO) << " number of components : " << qDim() << "\n";
        LOG(INFO) << "         n Global Dof : " << nDof() << "\n";
        LOG(INFO) << "         n Local  Dof : " << nLocalDof() << "\n";
    }

    //@}


    FunctionSpace( FunctionSpace const& __fe )
        :
        M_worldsComm( __fe.M_worldsComm ),
        M_worldComm( __fe.M_worldComm ),
        M_mesh( __fe.M_mesh ),
        M_ref_fe( __fe.M_ref_fe ),
        M_comp_space( __fe.M_comp_space ),
        M_dof( __fe.M_dof ),
        M_dofOnOff( __fe.M_dofOnOff ),
        M_extendedDofTableComposite( __fe.M_extendedDofTableComposite ),
        M_extendedDofTable( __fe.M_extendedDofTable ),
        M_rt( __fe.M_rt )
    {
        DVLOG(2) << "copying FunctionSpace\n";
    }

protected:

private:

    void init( mesh_ptrtype const& mesh,
               size_type mesh_components,
               periodicity_type const& periodicity,
               std::vector<Dof > const& dofindices,
               mpl::bool_<false> );
    void init( mesh_ptrtype const& mesh,
               size_type mesh_components,
               periodicity_type const& periodicity,
               std::vector<Dof > const& dofindices,
               mpl::bool_<true> );
    template<typename FSpaceHead, typename... FSpaceTail>
    void initList( FSpaceHead& fspacehead, FSpaceTail... fspacetail );

    void initList();

    template<typename FSpaceHead>
    void initHead( FSpaceHead& fspacehead );

    size_type nDof( mpl::bool_<false> ) const;
    size_type nDof( mpl::bool_<true> ) const;

    size_type nLocalDof( mpl::bool_<false> ) const;
    size_type nLocalDof( mpl::bool_<true> ) const;

    size_type nLocalDofWithGhost( mpl::bool_<false> ) const;
    size_type nLocalDofWithGhost( mpl::bool_<true> ) const;
    size_type nLocalDofWithoutGhost( mpl::bool_<false> ) const;
    size_type nLocalDofWithoutGhost( mpl::bool_<true> ) const;

    size_type nLocalDofWithGhostOnProc( const int proc, mpl::bool_<false> ) const;
    size_type nLocalDofWithGhostOnProc( const int proc, mpl::bool_<true> ) const;
    size_type nLocalDofWithoutGhostOnProc( const int proc, mpl::bool_<false> ) const;
    size_type nLocalDofWithoutGhostOnProc( const int proc, mpl::bool_<true> ) const;

    void rebuildDofPoints( mpl::bool_<false> );
    void rebuildDofPoints( mpl::bool_<true> );

    friend class ComponentSpace;
    class ComponentSpace
    {
    public:
        typedef FunctionSpace<A0,A1,A2,A3,A4> functionspace_type;
        typedef functionspace_type* functionspace_ptrtype;
        typedef functionspace_type const* functionspace_cptrtype;

        typedef typename FunctionSpace<A0,A1,A2,A3,A4>::component_functionspace_type component_functionspace_type;
        typedef typename FunctionSpace<A0,A1,A2,A3,A4>::component_functionspace_ptrtype component_functionspace_ptrtype;
        typedef component_functionspace_type const* component_functionspace_cptrtype;

        ComponentSpace( FunctionSpace<A0,A1,A2,A3,A4> * __functionspace,
                        mesh_ptrtype __m )
            :
            M_functionspace( __functionspace ),
            M_mesh( __m )
        {}
        component_functionspace_ptrtype operator()()
        {
            return operator()( mpl::bool_<functionspace_type::is_scalar>() );
        }
        component_functionspace_ptrtype operator()( mpl::bool_<true> )
        {
            FEELPP_ASSERT( 0 ).error( "invalid call for component space extraction" );
            //return M_functionspace;
            return component_functionspace_ptrtype();
        }
        component_functionspace_ptrtype operator()( mpl::bool_<false> )
        {
            return component_functionspace_type::New( M_mesh );
        }

    private:

        FunctionSpace<A0,A1,A2,A3,A4> * M_functionspace;
        mesh_ptrtype M_mesh;
    };


protected:

    //friend class FunctionSpace<mesh_type, typename bases_list::component_basis_type, value_type>;
    //friend class FunctionSpace<mesh_type, bases_list, value_type>;

    std::vector<WorldComm> M_worldsComm;
    boost::shared_ptr<WorldComm> M_worldComm;

    // finite element mesh
    mutable mesh_ptrtype M_mesh;
    mutable periodicity_type M_periodicity;
    //! finite element reference type
    reference_element_ptrtype M_ref_fe;

    //! component fe space
    mutable component_functionspace_ptrtype M_comp_space;

    //! Degrees of freedom
    dof_ptrtype M_dof;

    //! Degrees of freedom (only init wiht mpi)
    dof_ptrtype M_dofOnOff;

    //! build the extended dof table in //
    std::vector<bool> M_extendedDofTableComposite;
    bool M_extendedDofTable;

    /** region tree associated with the mesh */
    mutable boost::optional<region_tree_ptrtype> M_rt;

    //mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, real_type, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, real_type> > > M_prof_find_points;
private:

    //! disable default constructor
    //FunctionSpace();

    functionspace_vector_type M_functionspaces;

    proc_dist_map_type procDistMap;
#if 0
    std::list<Element> M_constraints;
#endif
}; // FunctionSpace

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename T,  typename Cont>
const bool FunctionSpace<A0,A1,A2,A3, A4>::Element<T,Cont>::is_scalar;

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename T,  typename Cont>
const bool FunctionSpace<A0,A1,A2,A3,A4>::Element<T,Cont>::is_vectorial;

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename T,  typename Cont>
const bool FunctionSpace<A0,A1,A2,A3,A4>::Element<T,Cont>::is_tensor2;

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename T,  typename Cont>
const uint16_type FunctionSpace<A0,A1,A2,A3,A4>::Element<T,Cont>::nComponents;

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::init( mesh_ptrtype const& __m,
        size_type mesh_components,
        periodicity_type const& periodicity,
        std::vector<Dof> const& dofindices,
        mpl::bool_<false> )
{
    DVLOG(2) << "calling init(<space>) begin\n";
    DVLOG(2) << "calling init(<space>) is_periodic: " << is_periodic << "\n";

    M_mesh = __m;
    M_periodicity = periodicity;
    VLOG(1) << "FunctionSpace init begin mesh use_count : " << M_mesh.use_count();


    if ( basis_type::nDofPerEdge || nDim >= 3 )
        mesh_components |= MESH_UPDATE_EDGES;

    /*
     * update faces info in mesh only if dofs exists on faces or the
     * expansion is continuous between elements. This case handles strong
     * Dirichlet imposition
     */
    if ( basis_type::nDofPerFace || is_continuous  || nDim >= 3 )
        mesh_components |= MESH_UPDATE_FACES;

    M_mesh->components().set( mesh_components );

    M_mesh->updateForUse();

    if ( is_periodic )
    {
        M_mesh->removeFacesFromBoundary( { periodicity.tag1(), periodicity.tag2() } );
    }

    M_ref_fe = basis_ptrtype( new basis_type );

    M_dof = dof_ptrtype( new dof_type( M_ref_fe, fusion::at_c<0>(periodicity), this->worldsComm()[0] ) );

    M_dof->setBuildDofTableMPIExtended( this->extendedDofTable() );

    DVLOG(2) << "[functionspace] Dof indices is empty ? " << dofindices.empty() << "\n";
    M_dof->setDofIndices( dofindices );
    DVLOG(2) << "[functionspace] is_periodic = " << is_periodic << "\n";

    M_dof->build( M_mesh );

    M_dofOnOff = M_dof;

    

    DVLOG(2) << "nb dim : " << qDim() << "\n";
    DVLOG(2) << "nb dof : " << nDof() << "\n";
    DVLOG(2) << "nb dof per component: " << nDofPerComponent() << "\n";


    //detail::searchIndicesBySpace<proc_dist_map_type>( this, procDistMap);

    DVLOG(2) << "calling init(<space>) end\n";
    VLOG(1) << "FunctionSpace init begin mesh use_count : " << M_mesh.use_count();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::init( mesh_ptrtype const& __m,
                                         size_type mesh_components,
                                         periodicity_type const& periodicity,
                                         std::vector<Dof> const& dofindices,
                                         mpl::bool_<true> )
{
    DVLOG(2) << "calling init(<composite>) begin\n";
    M_mesh = __m;

    // todo : check worldsComm size and M_functionspaces are the same!
    fusion::for_each( M_functionspaces,
                      Feel::detail::InitializeSpace<mesh_ptrtype,periodicity_type>( __m, periodicity,
                                                                                    dofindices,
                                                                                    this->worldsComm(),
                                                                                    this->extendedDofTableComposite() ) );

    this->initList();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::initList()
{
    if ( !M_worldComm )
        M_worldComm = boost::shared_ptr<WorldComm>( new WorldComm( M_worldsComm[0] ));

    if ( this->worldComm().globalSize()>1 )
    {
        if ( this->hasEntriesForAllSpaces() )
            {
                // construction with same partionment for all subspaces
                // and each processors has entries for all subspaces
                DVLOG(2) << "init(<composite>) type hasEntriesForAllSpaces\n";
                // build usefull data for Feel::detail::updateDataMapProcessStandard
                const int worldsize = this->worldComm().globalSize();
                std::vector<size_type> startDofGlobalCluster(worldsize);
                std::vector<size_type> nLocalDofWithoutGhostWorld(worldsize),nLocalDofWithGhostWorld(worldsize);
                for (int proc = 0; proc<worldsize ; ++proc )
                {
                    nLocalDofWithoutGhostWorld[proc] = this->nLocalDofWithoutGhostOnProc(proc);
                    nLocalDofWithGhostWorld[proc] = this->nLocalDofWithGhostOnProc(proc);
                }

                startDofGlobalCluster[0]=0;
                for (int p=1;p<worldsize;++p)
                    {
                        startDofGlobalCluster[p] = startDofGlobalCluster[p-1] + nLocalDofWithoutGhostWorld[p-1];
                    }

                // build datamap
                auto dofInitTool=Feel::detail::updateDataMapProcessStandard<dof_type>( this->worldsComm(), this->worldComm(),
                                                                                 this->nSubFunctionSpace()-1,
                                                                                 startDofGlobalCluster,
                                                                                 nLocalDofWithoutGhostWorld,
                                                                                 nLocalDofWithGhostWorld
                                                                                 );
                fusion::for_each( M_functionspaces, dofInitTool );
                // finish update datamap
                M_dof = dofInitTool.dataMap();
                M_dof->setNDof( this->nDof() );
                M_dofOnOff = M_dof;
            }
        else
            {
                // construction with same partionment for all subspaces
                // and one processor has entries for only one subspace
                DVLOG(2) << "init(<composite>) type Not hasEntriesForAllSpaces\n";

                // build the WorldComm associated to mix space
                WorldComm mixSpaceWorldComm = this->worldsComm()[0];

                if ( this->worldsComm().size()>1 )
                    for ( int i=1; i<( int )this->worldsComm().size(); ++i )
                        {
                            mixSpaceWorldComm = mixSpaceWorldComm + this->worldsComm()[i];
                        }

                this->setWorldComm( mixSpaceWorldComm );
                //mixSpaceWorldComm.showMe();

                // update DofTable for the mixedSpace (we have 2 dofTables : On and OnOff)
                auto dofInitTool=Feel::detail::updateDataMapProcess<dof_type>( this->worldsComm(), mixSpaceWorldComm, this->nSubFunctionSpace()-1 );
                fusion::for_each( M_functionspaces, dofInitTool );
                // finish update datamap
                M_dof = dofInitTool.dataMap();
                M_dof->setNDof( this->nDof() );
                M_dof->updateDataInWorld();
                M_dofOnOff = dofInitTool.dataMapOnOff();
                M_dofOnOff->setNDof( this->nDof() );
                M_dofOnOff->updateDataInWorld();
            }
    }

    else // sequential
    {
        // update DofTable for the mixedSpace (here On is not build properly but OnOff yes and On=OnOff, see Feel::detail::updateDataMapProcess)
        auto dofInitTool=Feel::detail::updateDataMapProcess<dof_type>( this->worldsComm(), this->worldComm(), this->nSubFunctionSpace()-1 );
        fusion::for_each( M_functionspaces, dofInitTool );
        M_dof = dofInitTool.dataMapOnOff();
        M_dof->setNDof( this->nDof() );
        M_dofOnOff = dofInitTool.dataMapOnOff();
        M_dofOnOff->setNDof( this->nDof() );
    }
    M_dof->setIndexSplit( this->buildDofIndexSplit() );
    //M_dof->indexSplit().showMe();

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename FSpaceHead>
void
FunctionSpace<A0, A1, A2, A3, A4>::initHead( FSpaceHead& head )
{
    DVLOG(2) << "calling initHead(<composite>) begin\n";
    M_mesh = head->mesh();
    M_worldsComm.push_back( head->worldComm() );
    M_extendedDofTableComposite.push_back( head->M_extendedDofTable );
    M_extendedDofTable = head->M_extendedDofTable;
}
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename FSpaceHead, typename... FSpaceTail>
void
FunctionSpace<A0, A1, A2, A3, A4>::initList( FSpaceHead& head, FSpaceTail... tail )
{
    initHead( head );
    initList( tail... );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nDof( mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NbDof() );
    DVLOG(2) << "calling nDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nDof( mpl::bool_<false> ) const
{
    return M_dof->nDof();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDof( mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm() ) );
    DVLOG(2) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDof( mpl::bool_<false> ) const
{
    //return M_dof->nLocalDof();
    return M_dof->nLocalDofWithGhost();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithGhost( mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm() ) );
    DVLOG(2) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithGhost( mpl::bool_<false> ) const
{
    return M_dof->nLocalDofWithGhost();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithGhostOnProc( const int proc, mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<true> >( proc, this->worldsComm() ) );
    DVLOG(2) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithGhostOnProc( const int proc, mpl::bool_<false> ) const
{
    return M_dof->nLocalDofWithGhost(proc);
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithoutGhost( mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<false> >( this->worldsComm() ) );
    DVLOG(2) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithoutGhost( mpl::bool_<false> ) const
{
    return M_dof->nLocalDofWithoutGhost();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithoutGhostOnProc( const int proc, mpl::bool_<true> ) const
{
    DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<false> >( proc, this->worldsComm() ) );
    DVLOG(2) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDofWithoutGhostOnProc( const int proc, mpl::bool_<false> ) const
{
    return M_dof->nLocalDofWithoutGhost(proc);
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::buildComponentSpace() const
{
    if ( is_vectorial && !M_comp_space )
    {
        // Warning: this works regarding the communicator . for the component space
        // it will use in mixed spaces only numberofSudomains/numberofspace processors
        //
        M_comp_space = component_functionspace_ptrtype( new component_functionspace_type( M_mesh,
                                                                                          MESH_COMPONENTS_DEFAULTS,
                                                                                          M_periodicity,
                                                                                          std::vector<WorldComm>( 1,this->worldsComm()[0] ) ) );
        VLOG(2) << " - component space :: nb dim : " << M_comp_space->qDim() << "\n";
        VLOG(2) << " - component space :: nb dof : " << M_comp_space->nDof() << "\n";
        VLOG(2) << " - component space :: nb dof per component: " << M_comp_space->nDofPerComponent() << "\n";
    }
}
template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::rebuildDofPoints( mpl::bool_<false> )
{
    M_dof->rebuildDofPoints( *M_mesh );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::rebuildDofPoints( mpl::bool_<true> )
{
    fusion::for_each( M_functionspaces, Feel::detail::rebuildDofPointsTool() );
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::updateRegionTree() const
{
    scalar_type EPS=1E-13;

    region_tree_ptrtype __rt( new region_tree_type );

    __rt->clear();
    BoundingBox<> __bb( M_mesh->gm()->isLinear() );

    typedef typename mesh_type::element_iterator mesh_element_iterator;
    mesh_element_iterator it = M_mesh->beginElementWithProcessId( M_mesh->comm().rank() );
    mesh_element_iterator en = M_mesh->endElementWithProcessId( M_mesh->comm().rank() );

    for ( size_type __i = 0; it != en; ++__i, ++it )
    {
        __bb.make( it->G() );

        for ( unsigned k=0; k < __bb.min.size(); ++k )
        {
            __bb.min[k]-=EPS;
            __bb.max[k]+=EPS;
        }

        __rt->addBox( __bb.min, __bb.max, it->id() );
    }

    //__rt->dump();
    M_rt = __rt;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
region_tree_ptrtype const&
FunctionSpace<A0, A1, A2, A3, A4>::regionTree() const
{
    if ( !M_rt )
    {
        scalar_type EPS=1E-13;

        region_tree_ptrtype __rt( new region_tree_type );

        __rt->clear();
        BoundingBox<> __bb( M_mesh->gm()->isLinear() );

        typedef typename mesh_type::element_iterator mesh_element_iterator;
        mesh_element_iterator it = M_mesh->beginElementWithProcessId( M_mesh->comm().rank() );
        mesh_element_iterator en = M_mesh->endElementWithProcessId( M_mesh->comm().rank() );

        for ( size_type __i = 0; it != en; ++__i, ++it )
        {
            __bb.make( it->G() );

            for ( unsigned k=0; k < __bb.min.size(); ++k )
            {
                __bb.min[k]-=EPS;
                __bb.max[k]+=EPS;
            }

            __rt->addBox( __bb.min, __bb.max, it->id() );
        }

        //__rt->dump();
        M_rt = __rt;
    }

    return M_rt.get();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
bool
FunctionSpace<A0, A1, A2, A3, A4>::findPoint( node_type const& pt,size_type &cv , node_type &ptr ) const
{
    if ( !hasRegionTree() )
        regionTree();

    //M_prof_find_points.resume();

    region_tree_type* __rt = M_rt.get().get();

    region_tree_type::pbox_set_type boxlst;
    __rt->findBoxesAtPoint( pt, boxlst );

    typedef typename gm_type::Inverse inv_trans_type;
    typename gm_type::reference_convex_type refelem;

    std::pair<size_type,value_type> closest = std::make_pair( invalid_size_type_value, -1e30 );

    region_tree_type::pbox_set_type::const_iterator it = boxlst.begin();
    region_tree_type::pbox_set_type::const_iterator ite = boxlst.end();

    for ( ; it != ite; ++it )
    {
        inv_trans_type __git( M_mesh->gm(), M_mesh->element( ( *it )->id ), this->worldComm().subWorldCommSeq() );

        size_type cv_stored = ( *it )->id;


        DVLOG(2) << "[FunctionSpace::findPoint] id : " << cv_stored << "\n";

        __git.setXReal( pt );
        ptr = __git.xRef();


        bool isin;
        value_type dmin;
        boost::tie( isin, dmin ) = refelem.isIn( ptr );
        DVLOG(2) << "[FunctionSpace::findPoint] isin: " << isin << " dmin: " << dmin << "\n";

        closest =  ( dmin > closest.second )?std::make_pair( cv_stored, dmin ):closest;

        if ( isin )
        {
            DVLOG(2) << "[FunctionSpace::findPoint] id of the convex where " << pt << " belongs : " << cv_stored << "\n";
            DVLOG(2) << "[FunctionSpace::findPoint] ref coordinate: " << ptr << "\n";
            cv = ( *it )->id;
            //M_prof_find_points.pause();
            return true;
        }
    }

    cv=closest.first;
    //M_prof_find_points.pause();
    return false;
}

template<typename T,int M,int N>
std::ostream&
operator<<( std::ostream& os, Feel::detail::ID<T,M,N> const& id )
{
    const size_type* shape =  id.M_id.shape();

    for ( size_type i = 0; i < shape[0]; ++i )
    {
        os << id[i] << std::endl;
    }
    os << std::endl;

    return os;
}

/**
   iostream operator to print some information about the function space \p Xh
   \code
   auto Xh = Pch<2>(mesh);
   // Xh is a pointer
   std::cout << "Xh:" << *Xh << "\n";
   \endcode
 */
template<typename A0, typename A1, typename A2, typename A3, typename A4>
std::ostream&
operator<<( std::ostream& os, FunctionSpace<A0, A1, A2, A3, A4> const& Xh )
{
    os << "Number of Dof: Global=" << Xh.nDof() << " , Local=" << Xh.nLocalDof();
    return os;
}


namespace detail
{

template<typename FuncSpaceType>
struct is_function_space_ptr : mpl::false_ {};

template<typename FuncSpaceType>
struct is_function_space_ptr<boost::shared_ptr<FuncSpaceType> > : mpl::true_ {};
} // detail





} // Feel

#include <feel/feeldiscr/detail/element_impl.hpp>



#endif /* __FunctionSpace_H */
