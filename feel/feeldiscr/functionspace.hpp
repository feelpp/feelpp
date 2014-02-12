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

#include <feel/feelvf/detail/gmc.hpp>

namespace Feel
{
namespace fusion = boost::fusion;
namespace parameter = boost::parameter;

namespace detail
{


template<typename T>
struct vector_plus
{
    std::vector<T> operator()( std::vector<T> const& v1, std::vector<T> const& v2 ) const
    {
        FEELPP_ASSERT( v1.size() == v2.size() )( v1.size() )( v2.size() ).error( "invalid vector size for vector_plus<>" );
        std::vector<T> res( v1.size() );

        for ( size_type i = 0; i < v1.size(); ++i )
            res[i]=v1[i]+v2[i];

        return res;
    }
};
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
        size_type e2 = M_id.shape()[1];
        DVLOG(2) << "saving in archive e2= " << e2 << "\n";
        ar  & e2;
        size_type e3 = M_id.shape()[2];
        DVLOG(2) << "saving in archive e3= " << e3 << "\n";
        ar  & e3;
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
        ar  & e2;
        DVLOG(2) << "loading from archive e2= " << e2 << "\n";
        ar  & e3;
        DVLOG(2) << "loading from archive e3= " << e3 << "\n";
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
        elem.grad_( context, M_grad );
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
        size_type e2 = M_grad.shape()[1];
        DVLOG(2) << "saving in archive e2= " << e2 << "\n";
        ar  & e2;
        size_type e3 = M_grad.shape()[2];
        DVLOG(2) << "saving in archive e3= " << e3 << "\n";
        ar  & e3;
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
        ar  & e2;
        DVLOG(2) << "loading from archive e2= " << e2 << "\n";
        ar  & e3;
        DVLOG(2) << "loading from archive e3= " << e3 << "\n";
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
        elem.d_( N, context, M_grad );
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
        elem.div_( context, M_div );
#if 0
        uint16_type nComponents1 = elem.nComponents1;
        std::fill( M_div.data(), M_div.data()+M_div.num_elements(), value_type( 0 ) );

        M_grad( elem.div_( context, pc, M_div ) ),


                 const uint16_type nq = context.xRefs().size2();

        for ( int c1 = 0; c1 < nComponents1; ++c1 )
            for ( uint16_type q = 0; q < nq ; ++q )
            {
                M_div[q]( 0,0 ) += M_grad[q]( c1,c1 );
            }

#endif
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
        elem.hess_( context, M_hess );
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
        M_startSplit(startSplit)
    {}

    typedef boost::tuple< uint16_type, size_type, std::vector<std::vector<size_type> > > result_type;

    template<typename T>
    result_type operator()( result_type const & previousRes, T const& t )
    {
        const size_type nDofWithoutGhost = t->dof()->nLocalDofWithoutGhost();
        const size_type nDofWithGhost = t->dof()->nLocalDofWithGhost();
        const size_type firstDof = t->dof()->firstDofGlobalCluster();
        uint16_type cptSpaces = previousRes.get<0>();
        const size_type start = previousRes.get<1>();
        auto is = previousRes.get<2>();

        //std::cout << "compo " << cptSpaces << " start " << start <<" split nDofWithout " << nDof << " with ghost "<< nDofWithGhost<< " M_startSplit " << M_startSplit <<std::endl;

        //is.push_back( std::vector<size_type>( nDofWithoutGhost ) );
        is[cptSpaces].resize(nDofWithoutGhost);

        for ( size_type i=0; i<nDofWithGhost; ++i )
        {
            if ( t->dof()->dofGlobalProcessIsGhost(i) ) continue;
            const size_type globalDof = t->dof()->mapGlobalProcessToGlobalCluster(i);
            is[cptSpaces][globalDof - firstDof ] = M_startSplit + start + (globalDof - firstDof);
        }

        return boost::make_tuple( ++cptSpaces, ( start+nDofWithoutGhost ), is );
    }

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

template<typename ElementType>
struct CreateElementVector
{
public:
    template<typename Sig>
    struct result;

    template<typename Lhs, typename Rhs>
    struct result<CreateElementVector( Lhs,Rhs )>
    {
	    typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
	    typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type::element_type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<ElementType>::type>::type ElementType_noref_type;
	    typedef typename boost::fusion::result_of::size<lhs_noref_type>::type index;
        typedef typename ElementType_noref_type::template sub_element<index::value>::type elt_type;
        BOOST_MPL_ASSERT( ( boost::is_same<typename elt_type::functionspace_type,rhs_noref_type> ) );
        typedef typename boost::fusion::result_of::make_vector<elt_type>::type v_elt_type;

	    typedef typename boost::fusion::result_of::push_back<lhs_noref_type, elt_type>::type ptype;
        typedef typename boost::fusion::result_of::as_vector<ptype>::type type;
    };
    CreateElementVector( ElementType const& e ) : M_e( e ), M_names() {}
    CreateElementVector( ElementType const& e, std::vector<std::string> const& names ) : M_e( e ), M_names( names ) {}
    ElementType const& M_e;
    std::vector<std::string> M_names;

    template<typename Lhs, typename Rhs>
    typename result<CreateElementVector( Lhs,Rhs )>::type
    operator()( Lhs const&  lhs, Rhs const& rhs ) const
    {
        typedef typename boost::remove_const<typename boost::remove_reference<Lhs>::type>::type lhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<Rhs>::type>::type rhs_noref_type;
        typedef typename boost::remove_const<typename boost::remove_reference<ElementType>::type>::type ElementType_noref_type;
	    typedef typename boost::fusion::result_of::size<lhs_noref_type>::type index;
        typename ElementType_noref_type::template sub_element<index::value>::type elt = M_e.template element<index::value>();
        static const int s = mpl::size<typename ElementType::functionspace_type::bases_list>::type::value;
        BOOST_STATIC_ASSERT( (boost::is_same<decltype(elt), typename ElementType::template sub_element<index::value>::type>::value ) );
        if ( !M_names.empty() && M_names.size() > index::value )
        {

            FEELPP_ASSERT( M_names.size() == s  )
                ( M_names.size() )( s ).error( "incompatible number of function names and functions");
            elt.setName( M_names[index::value] );
        }
        else if  ( ( M_names.size() == 1 )  && s > 1 )
        {
            elt.setName( (boost::format( "%1%-%2%" ) % M_names[0] % index::value ).str() );
        }
        else
        {
            elt.setName( (boost::format( "%1%-%2%" ) % M_e.name() % index::value ).str() );
        }
        return boost::fusion::as_vector( boost::fusion::push_back( lhs, elt ) );
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
    typedef typename functionspace_signature::bind<A0,A1,A2,A3,A4>::type args;

    typedef typename parameter::binding<args, tag::mesh_type>::type meshes_list;
    typedef typename parameter::binding<args, tag::value_type, double>::type value_type;
    typedef typename parameter::binding<args, tag::mortar_type, mortars<NoMortar> >::type mortar_list;
    typedef typename parameter::binding<args, tag::periodicity_type, Periodicity<NoPeriodicity> >::type periodicity_type;
    typedef typename parameter::binding<args, tag::bases_list, Feel::detail::bases<Lagrange<1,Scalar> > >::type bases_list;

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
    static const uint16_type nDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                             mpl::identity<meshes_list >,
                                             mpl::identity<nodim> >::type::type::nDim;
    static const uint16_type nRealDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                                 mpl::identity<meshes_list>,
                                                 mpl::identity<nodim> >::type::type::nRealDim;



    //typedef typename mpl::at_c<bases_list,0>::type::template apply<mesh_type::nDim,value_type,typename mesh_type::element_type>::type basis_0_type;
    typedef typename mpl::at_c<bases_list,0>::type::template apply<GetMesh<meshes_list,0>::type::nDim,
                                                                   GetMesh<meshes_list,0>::type::nRealDim,
                                                                   value_type,
                                                                   typename GetMesh<meshes_list,0>::type::element_type>::type basis_0_type;

    static const uint16_type rank = ( is_composite? invalid_uint16_type_value : basis_0_type::rank );
    static const bool is_scalar = ( is_composite? false : basis_0_type::is_scalar );
    static const bool is_vectorial = ( is_composite? false : basis_0_type::is_vectorial );
    static const bool is_tensor2 = ( is_composite? false : basis_0_type::is_tensor2 );
    static const bool is_continuous = ( is_composite? false : basis_0_type::isContinuous );
    static const bool is_modal = ( is_composite? false : basis_0_type::is_modal );
    static const uint16_type nComponents1 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents1 );
    static const uint16_type nComponents2 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents2 );
    static const bool is_product = ( is_composite? invalid_uint16_type_value : basis_0_type::is_product );
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
        typedef typename mpl::at_c<functionspace_vector_type,i>::type type;
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
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );
                else
                    global_found_pt[ 0 ] = found_pt[ 0 ];

            }//if( found )
            else
            {
                if ( nprocs > 1 )
                    mpi::all_reduce( M_Xh->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );

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

        template<typename BasisType>
        struct ChangeElement
        {
            typedef T value_type;
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef typename ChangeBasis<BasisType>::type::element_type fs_type;
            typedef typename fs_type::template Element<value_type, typename VectorUblas<T>::range::type > element_type;
            typedef element_type type;
        };
        typedef typename mpl::transform<bases_list, ChangeElement<mpl::_1>, mpl::back_inserter<mpl::vector<> > >::type element_vector_type;

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

        typedef Element<T,Cont> this_type;
        template<int i>
        struct sub_element
        {
            typedef typename mpl::at_c<element_vector_type,i>::type type;
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
            size_type index=start()+boost::get<0>( M_functionspace->dof()->localToGlobal( e, l, c ) );
            return super::operator()( index );
        }
#if 0
        value_type& localToGlobal( size_type e, size_type l, int c )
        {
            size_type index=start()+boost::get<0>( M_functionspace->dof()->localToGlobal( e, l, c ) );
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
        }

        template<typename AE>
        container_type& assign( const ublas::vector_expression<AE> &ae )
        {
            return super::assign( ae );
        }
        void assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index=start()+ M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            this->operator[]( index ) = __v;
        }
        void plus_assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index=start()+ M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            this->operator[]( index ) += __v;
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
            return super::max();
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type max( boost::shared_ptr<p0_space_type> const& P0h )
        {
            return this->extremeValue( P0h, "max" );
        }

        value_type min() const
        {
            return super::min();
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
            int npoints = context.nPoints();

            //rank of the current processor
            int proc_number = Environment::worldComm().globalRank();

            //total number of processors
            int nprocs = Environment::worldComm().globalSize();

            auto it = context.begin();
            auto en = context.end();

            eigen_type __globalr( npoints );
            __globalr.setZero();
            eigen_type __localr( npoints );
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
                    __localr( global_index ) = v[0]( 0, 0 );
                }
            }

            if( do_communications )
                mpi::all_reduce( Environment::worldComm() , __localr, __globalr, std::plus< eigen_type >() );
            else
                __globalr = __localr;

            return __globalr;
        }

        double
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

            double result=0;
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
        typename mpl::at_c<element_vector_type,i>::type
        element( std::string const& name ="u", bool updateOffViews=true )
        {
            if ( this->worldComm().globalSize()>1 ) this->worldComm().globalComm().barrier();

            size_type nbdof_start =  fusion::accumulate( this->functionSpaces(),
                                     size_type( 0 ),
                                     Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(), false, 0, i ) );

            typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );
            DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
            DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
            DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
            DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

            if ( this->functionSpace()->template functionSpace<i>()->worldComm().isActive() )
            {
                ct_type ct( *this, ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ),
                            M_functionspace->template functionSpace<i>()->dof() );

                // update M_containersOffProcess<i> : send
                if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
                {
                    std::vector<double> dataToSend( space->nLocalDof() );
                    std::copy( ct.begin(), ct.end(), dataToSend.begin() );

                    if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

                    fusion::for_each( *M_containersOffProcess, Feel::detail::SendContainersOn<i,functionspace_type>( this->functionSpace(), dataToSend ) );
                }

                DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
                DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
                return typename mpl::at_c<element_vector_type,i>::type( space, ct, name );
            }

            else
            {
                // initialize if not the case
                if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

                fusion::for_each( *M_containersOffProcess, Feel::detail::InitializeContainersOff<i,functionspace_type>( this->functionSpace() ) );

                // update M_containersOffProcess<i> : recv
                if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
                {
                    fusion::for_each( *M_containersOffProcess, Feel::detail::RecvContainersOff<i,functionspace_type>( this->functionSpace() ) );
                }

                // build a subrange view identical
                ct_type ct( *fusion::at_c<i>( *M_containersOffProcess ),
                            ublas::range( 0, space->nLocalDof() ),
                            M_functionspace->template functionSpace<i>()->dof() );

                DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
                DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";

                return typename mpl::at_c<element_vector_type,i>::type( space, ct, name );
            }
        }

        template<int i>
        typename mpl::at_c<element_vector_type,i>::type
        element( std::string const& name ="u", bool updateOffViews=true ) const
        {
            size_type nbdof_start =  fusion::accumulate( M_functionspace->functionSpaces(),
                                     size_type( 0 ),
                                     Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm(), false, 0, i ) );
            typename mpl::at_c<functionspace_vector_type,i>::type space( M_functionspace->template functionSpace<i>() );

            DVLOG(2) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
            DVLOG(2) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
            DVLOG(2) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
            DVLOG(2) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

            if ( this->functionSpace()->worldsComm()[i].isActive() )
            {
                ct_type ct( const_cast<VectorUblas<value_type>&>( dynamic_cast<VectorUblas<value_type> const&>( *this ) ),
                            ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ),
                            M_functionspace->template functionSpace<i>()->dof() );

                // update M_containersOffProcess<i> : send
                if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
                {
                    std::vector<double> dataToSend( space->nLocalDof() );
                    std::copy( ct.begin(), ct.end(), dataToSend.begin() );

                    if ( !M_containersOffProcess ) M_containersOffProcess = boost::in_place();

                    fusion::for_each( *M_containersOffProcess, Feel::detail::SendContainersOn<i,functionspace_type>( this->functionSpace(), dataToSend ) );
                }

                DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
                DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
                return typename mpl::at_c<element_vector_type,i>::type( space, ct, name );
            }

            else
            {
                // update M_containersOffProcess<i> : recv
                if ( this->worldComm().globalSize()>1 && updateOffViews && !this->functionSpace()->hasEntriesForAllSpaces() )
                {
                    fusion::for_each( *M_containersOffProcess, Feel::detail::RecvContainersOff<i,functionspace_type>( this->functionSpace() ) );
                }

                // build a subrange view identical
                ct_type ct( *fusion::at_c<i>( *M_containersOffProcess ),
                            ublas::range( 0, space->nLocalDof() ),
                            M_functionspace->template functionSpace<i>()->dof() );

                DVLOG(2) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
                DVLOG(2) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
                return typename mpl::at_c<element_vector_type,i>::type( space, ct, name );
            }

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
                                           ( verbose,   ( bool ), option(_prefix=prefix,_name="on.verbose").template as<bool>() )))
            {
                return onImpl( range, expr, prefix, geomap, accumulate, verbose );
            }


        //@}
    private:

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

    std::vector<std::vector<size_type> > dofIndexSplit()
    {
        if ( nSpaces > 1 )
        {
            uint16_type cptSpaces=0;
            size_type start=0;
            std::vector<std::vector<size_type> > is(nSpaces);
            auto initial = boost::make_tuple( cptSpaces,start,is );

            // get start for each proc ->( proc0 : 0 ), (proc1 : sumdofproc0 ), (proc2 : sumdofproc0+sumdofproc1 ) ....
            auto startSplit = boost::fusion::fold( functionSpaces(), boost::make_tuple(0,0), Feel::detail::computeStartOfFieldSplit() ).template get<1>();
            // compute split
            auto result = boost::fusion::fold( functionSpaces(), initial,  Feel::detail::computeNDofForEachSpace(startSplit) );
            is = result.template get<2>();


#if 0
            for ( int proc = 0; proc<this->worldComm().globalSize(); ++proc )
            {
                this->worldComm().globalComm().barrier();
                if ( proc==this->worldComm().globalRank() )
                {
                    std::cout << "proc " << proc << "\n";
                    std::cout << "split size=" << result.template get<2>().size() << " nspace=" << nSpaces << "\n";
                    std::cout << "split:\n";

                    //std::cout << "\n\n";


                    for ( int s= 0; s < nSpaces; ++s )
                    {
                        std::cout << "space: " << is[s].size() << "\n";

                        for ( int i = 0; i < is[s].size(); ++i )
                        {
                            std::cout << is[s][i] << " ";
                        }

                        std::cout << "\n\n";
                    }
                }
            }

#endif
            return is;
        }

        std::vector<std::vector<size_type> > is;
        is.push_back( std::vector<size_type>( nLocalDof() ) );
        int index = 0;
        //for( int& i : is[0] ) { i = index++; }
        BOOST_FOREACH( auto& i, is[0] )
        {
            i = index++;
        }
        return is;

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
        return trace_trace_functionspace_type::New( mesh()->wireBasket( markededges( mesh(),"WireBasket" ) ) );
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
        return M_rt;
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
    mesh_ptrtype M_mesh;

    //! finite element reference type
    reference_element_ptrtype M_ref_fe;

    //! component fe space
    component_functionspace_ptrtype M_comp_space;

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
    FunctionSpace();

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

    if ( is_vectorial )
    {
        // Warning: this works regarding the communicator . for the component space
        // it will use in mixed spaces only numberofSudomains/numberofspace processors
        //
        M_comp_space = component_functionspace_ptrtype( new component_functionspace_type( M_mesh,
                                                                                           MESH_COMPONENTS_DEFAULTS,
                                                                                           periodicity,
                                                                                           std::vector<WorldComm>( 1,this->worldsComm()[0] ) ) );
    }

    DVLOG(2) << "nb dim : " << qDim() << "\n";
    DVLOG(2) << "nb dof : " << nDof() << "\n";
    DVLOG(2) << "nb dof per component: " << nDofPerComponent() << "\n";

    if ( is_vectorial )
    {
        DVLOG(2) << "component space :: nb dim : " << M_comp_space->qDim() << "\n";
        DVLOG(2) << "component space :: nb dof : " << M_comp_space->nDof() << "\n";
        DVLOG(2) << "component space :: nb dof per component: " << M_comp_space->nDofPerComponent() << "\n";
    }

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

//
// Element implementation
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element()
    :
    super(),
    M_start( 0 ),
    M_ct( NO_COMPONENT ),
    M_containersOffProcess( boost::none )
{}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( Element const& __e )
    :
    super( __e ),
    M_functionspace( __e.M_functionspace ),
    M_name( __e.M_name ),
    M_start( __e.M_start ),
    M_ct( __e.M_ct ),
    M_containersOffProcess( __e.M_containersOffProcess )
{
    DVLOG(2) << "Element<copy>::range::start = " << this->start() << "\n";
    DVLOG(2) << "Element<copy>::range::size = " << this->size() << "\n";

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
        std::string const& __name,
        size_type __start,
        ComponentType __ct )
    :
    super( __functionspace->dof() ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_start( __start ),
    M_ct( __ct ),
    M_containersOffProcess( boost::none )
{
    DVLOG(2) << "Element::start = " << this->start() << "\n";
    DVLOG(2) << "Element::size = " << this->size() << "\n";
    DVLOG(2) << "Element::ndof = " << this->nDof() << "\n";
    DVLOG(2) << "Element::nlocaldof = " << this->nLocalDof() << "\n";
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
        container_type const& __c,
        std::string const& __name,
        size_type __start,
        ComponentType __ct )
    :
    super( __c ),
    M_functionspace( __functionspace ),
    M_name( __name ),
    M_start( __start ),
    M_ct( __ct ),
    M_containersOffProcess( boost::none )
{
    DVLOG(2) << "Element<range>::range::start = " << __c.start() << "\n";
    DVLOG(2) << "Element<range>::range::size = " << __c.size() << "\n";
    DVLOG(2) << "Element<range>::start = " << this->start() << "\n";
    DVLOG(2) << "Element<range>::size = " << this->size() << "\n";
    DVLOG(2) << "Element<range>::ndof = " << this->nDof() << "\n";
    DVLOG(2) << "Element<range>::nlocaldof = " << this->nLocalDof() << "\n";
    M_start = __c.start();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::~Element()
{
    VLOG(1) << "Element destructor...";
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::initFromSpace( functionspace_ptrtype const& __functionspace,
        container_type const& __c )
{
    M_functionspace = __functionspace;
    ( container_type )*this = __c;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element<Y,Cont> const& __e )
{
    if (  this != &__e )
    {
        M_functionspace = __e.M_functionspace;

        if ( __e.M_name != "unknown" )
            M_name = __e.M_name;

        M_start = __e.M_start;
        M_ct = __e.M_ct;
        M_containersOffProcess = __e.M_containersOffProcess;
        this->resize( M_functionspace->nLocalDof() );
        super::operator=( __e );
        this->outdateGlobalValues();
    }

    return *this;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename VectorExpr>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( VectorExpr const& __v )
{
    if (  __v.size() != this->size() )
    {
        std::ostringstream __err;
        __err << "Invalid vector size this->size()=" << this->size()
              << " and v.size()=" << __v.size();
        throw std::logic_error( __err.str() );
    }

    this->outdateGlobalValues();
    super::operator=( __v );
    return *this;
}


//
// Interpolation tools
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::id_type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator()( node_type const& __x, bool extrapolate ) const
{
    this->updateGlobalValues();

    node_type __x_ref;
    size_type __cv_id;
    int rank = functionSpace()->mesh()->comm().rank();
    int nprocs = functionSpace()->mesh()->comm().size();
    std::vector<int> found_pt( nprocs, 0 );
    std::vector<int> global_found_pt( nprocs, 0 );

    if ( functionSpace()->findPoint( __x, __cv_id, __x_ref ) || extrapolate )
    {
#if !defined( NDEBUG )
        DVLOG(2) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";
#endif
        gm_ptrtype __gm = functionSpace()->gm();
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        typename matrix_node<value_type>::type pts( __x_ref.size(), 1 );
        ublas::column( pts, 0 ) = __x_ref;
        geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );


        typedef typename gm_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        gmc_ptrtype __c( new gmc_type( __gm,
                                       functionSpace()->mesh()->element( __cv_id ),
                                       __geopc ) );
        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );

        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::fe_type fe_type;
        typedef typename fe_type::template Context<vm::POINT|vm::GRAD|vm::KB|vm::JACOBIAN, fe_type, gm_type, geoelement_type> fectx_type;
        typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
        fectx_ptrtype fectx( new fectx_type( this->functionSpace()->fe(),
                                             __c,
                                             pc ) );
        found_pt[ rank ] = 1;

#if defined(FEELPP_HAS_MPI)

        if ( nprocs > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );
        }

#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* FEELPP_HAS_MPI */

        id_type __id( this->id( *fectx ) );

        //DVLOG(2) << "[interpolation]  id = " << __id << "\n";
#if defined(FEELPP_HAS_MPI)
        DVLOG(2) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), __id, functionSpace()->mesh()->comm().rank() );
        }

        //DVLOG(2) << "[interpolation] after broadcast id = " << __id << "\n";
#endif /* FEELPP_HAS_MPI */
        return __id;
    }

    else
    {
#if defined(FEELPP_HAS_MPI)

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );
        }

#endif /* FEELPP_HAS_MPI */
        bool found = false;
        size_type i = 0;

        for ( ; i < global_found_pt.size(); ++i )
            if ( global_found_pt[i] != 0 )
            {
                DVLOG(2) << "processor " << i << " has the point " << __x << "\n";
                found = true;
                break;
            }

        id_type __id;

        if ( found )
        {
            DVLOG(2) << "receiving interpolation context from processor " << i << "\n";
#if defined(FEELPP_HAS_MPI)

            if ( functionSpace()->mesh()->comm().size() > 1 )
                mpi::broadcast( functionSpace()->mesh()->comm(), __id, i );

#endif /* FEELPP_HAS_MPI */
        }

        else
        {
            LOG(WARNING) << "no processor seems to have the point " << __x << "\n";
        }

        return __id;
    }

} // operator()
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename Context_t>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    const uint16_type nq = context.xRefs().size2();

    //double vsum=0;

    //array_type v( boost::extents[nComponents1][nComponents2][context.xRefs().size2()] );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( typename array_type::index c1 = 0; c1 < ncdof; ++c1 )
        {
            typename array_type::index ldof = basis_type::nDof*c1+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            //std::cout << "ldof = " << ldof << "\n";
            //std::cout << "gdof = " << gdof << "\n";

            FEELPP_ASSERT( gdof < this->size() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .warn( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            //std::cout << "v_ =" << v_ << "\n";
            //for( typename array_type::index c2 = 0; c2 < nComponents2; ++c2 )
            for ( uint16_type q = 0; q < nq; ++q )
            {
                for ( typename array_type::index i = 0; i < nComponents1; ++i )
                    //for( typename array_type::index j = 0; j < nComponents2; ++j )
                {
                    v[q]( i,0 ) += v_*context.id( ldof, i, 0, q );
                    //vsum +=v_*context.id( ldof, i, 0, q );
                    //v[q](i,0) += v_*context.gmc()->J(*)*context.pc()->phi( ldof, i, 0, q );
                }
            }
        }
    }

    //LOG( INFO ) << "vsum : "<<vsum;
    //return v;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );

    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );

        //update precompute of basis functions
        __pc->update( pts );

        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        id_type __id( this->id( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
        {
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __id( i,0,k );
            }
        }
    }

}

//
// Grad
//
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::grad_type
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::grad( node_type const& __x ) const
{
    this->updateGlobalValues();

    node_type __x_ref;
    size_type __cv_id;
    std::vector<int> found_pt( functionSpace()->mesh()->comm().size(), 0 );
    std::vector<int> global_found_pt( functionSpace()->mesh()->comm().size(), 0 );

    if ( functionSpace()->findPoint( __x, __cv_id, __x_ref ) )
    {
#if !defined( NDEBUG )
        DVLOG(2) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";
#endif
        gm_ptrtype __gm = functionSpace()->gm();
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        typename matrix_node<value_type>::type pts( __x_ref.size(), 1 );
        ublas::column( pts, 0 ) = __x_ref;
        geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );


        typedef typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::KB, geoelement_type> gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        gmc_ptrtype __c( new gmc_type( __gm,
                                       functionSpace()->mesh()->element( __cv_id ),
                                       __geopc ) );
        pc_ptrtype pc( new pc_type( this->functionSpace()->fe(), pts ) );
        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::fe_type fe_type;
        typedef typename fe_type::template Context<vm::GRAD, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
        typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
        fectx_ptrtype fectx( new fectx_type( this->functionSpace()->fe(),
                                             __c,
                                             pc ) );

        found_pt[ functionSpace()->mesh()->comm().rank() ] = 1;

#if defined(FEELPP_HAS_MPI)

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );
        }

#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* FEELPP_HAS_MPI */

        grad_type g_( this->grad( *fectx ) );
        //DVLOG(2) << "[interpolation]  id = " << v << "\n";
#if defined(FEELPP_HAS_MPI)
        DVLOG(2) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), g_, functionSpace()->mesh()->comm().rank() );
        }

        //DVLOG(2) << "[interpolation] after broadcast g_ = " << g_ << "\n";
#endif /* FEELPP_HAS_MPI */
        return g_;
    }

    else
    {
#if defined(FEELPP_HAS_MPI)

        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, Feel::detail::vector_plus<int>() );
        }

#endif /* FEELPP_HAS_MPI */
        bool found = false;
        size_type i = 0;

        for ( ; i < global_found_pt.size(); ++i )
            if ( global_found_pt[i] != 0 )
            {
                DVLOG(2) << "processor " << i << " has the point " << __x << "\n";
                found = true;
                break;
            }

        grad_type g_;

        if ( found )
        {
            DVLOG(2) << "receiving interpolation context from processor " << i << "\n";
#if defined(FEELPP_HAS_MPI)

            if ( functionSpace()->mesh()->comm().size() > 1 )
                mpi::broadcast( functionSpace()->mesh()->comm(), g_, i );

#endif /* FEELPP_HAS_MPI */

            //DVLOG(2) << "[interpolation] after broadcast id = " << v << "\n";
        }

        else
        {
            Warning() << "no processor seems to have the point " << __x << "\n";
        }

        return g_;
    }

} // grad
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::grad_( ContextType const & context, grad_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    //std::cout << "coeff=" << coeff << "\n";
    //array_type v( boost::extents[nComponents1][nRealDim][context.xRefs().size2()] );
    //std::fill( v.data(), v.data()+v.num_elements(), value_type( 0 ) );
    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( int k = 0; k < nComponents1; ++k )
                    for ( int j = 0; j < nRealDim; ++j )
                    {
                        v[q]( k,j ) += v_*context.grad( ldof, k, j, q );
                    }
            }
        }
    }

#if 0
    std::cout << "xref = " << context.xRefs() << "\n";
    std::cout << "xreal = " << context.xReal() << "\n";
    std::cout << "pts = " << context.G() << "\n";

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        for ( typename array_type::index i = 0; i < nDim; ++i )
        {
            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                std::cout << "Gref[ " << l << "," << i << "," << q << "] = " << pc.grad( l, 0, i, q ) << "\n";

        }
    }

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        for ( int j = 0; j < nRealDim; ++j )
            for ( typename array_type::index i = 0; i < nDim; ++i )
            {
                for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                    std::cout << "G[ " << l << "," << j << "," << i << "," << q << "] = " << context.B( q )( j, i )*pc.grad( l, 0, i, q ) << "\n";

            }
    }

    for ( int k = 0; k < nComponents1; ++k )
        for ( int j = 0; j < nRealDim; ++j )
        {
            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                std::cout << "v[ " << k << "," << j << "," << q << "]= " << v[q]( k,j )  << "\n";
                std::cout << "B(" << q << ")=" << context.B( q ) << "\n";
                std::cout << "J(" << q << ")=" << context.J( q ) << "\n";
                std::cout << "K(" << q << ")=" << context.K( q ) << "\n";
            }
        }

#endif
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::gradInterpolate(  matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );

    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );


        //update precompute of basis functions
        __pc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        grad_type __grad( this->grad( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
        {
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
                for ( uint j = 0; j < nRealDim; ++j )
                {
                    v[boost::get<0>( *itL )]( i,j ) = __grad( i,j,k );
                }
        }

    }

}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::div_( ContextType const & context, div_array_type& v ) const
{
#if 1

    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    const size_type Q = context.xRefs().size2();

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );
            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < Q; ++q )
            {
#if 0
                std::cout << "v(" << gdof << ")=" << v_ << "\n";
                std::cout << "context.div(" << ldof << "," << q << ")="
                          << context.div( ldof, 0, 0, q ) << "\n" ;
#endif
                v[q]( 0,0 ) += v_*context.div( ldof, 0, 0, q );
            }
        }
    }

#else

    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );
            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );

            for ( int k = 0; k < nComponents1; ++k )
            {
                for ( typename array_type::index i = 0; i < nDim; ++i )
                {
                    for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                    {
                        v[q]( 0,0 ) += v_*context.gmContext()->B( q )( k, i )*context.pc()->grad( ldof, k, i, q );
                    }
                }
            }
        }
    }

#endif
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::divInterpolate( matrix_node_type __ptsReal, div_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::DIV|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        div_type __div( this->div( *__ctx ) );

        //update the output data
        //for( typename array_type::index i = 0; i < nComponents1; ++i )
        //   {
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            v[boost::get<0>( *itL )]( 0,0 ) =  __div( 0,0,k );

        //   }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curl_( ContextType const & context, curl_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;

    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );
            const uint16_type nq = context.xRefs().size2();

            for ( uint16_type q = 0; q < nq ; ++q )
            {
                if ( nDim == 3 )
                {
                    for ( typename array_type::index i = 0; i < nDim; ++i )
                    {
                        v[q]( i,0 ) += v_*context.curl( ldof, i, 0, q );
                    }
                }

                else if ( nDim == 2 )
                {
                    v[q]( 0,0 ) += v_*context.curl( ldof, 0, 0, q );
                }

            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curl_( ContextType const & context, comp_curl_array_type& v, int comp ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_size_type_value )
        return;


    for ( int l = 0; l < basis_type::nDof; ++l )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            int ldof = c1*basis_type::nDof+l;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( elt_id, l, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( l )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            //value_type v_ = (*this)( gdof );
            value_type v_ = this->globalValue( gdof );
            const uint16_type nq = context.xRefs().size2();

            for ( uint16_type q = 0; q < nq ; ++q )
            {
                if ( nDim == 3 )
                {
                    v[q]( 0,0 ) += v_*context.curl( ldof, comp, 0, q );
                }

                else if ( nDim == 2 )
                {
                    v[q]( 0,0 ) += v_*context.curl( ldof, 2, 0, q );
                }

            }
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );


    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curl_type __curl( this->curl( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curl( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curl( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curl( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curl( 0,0,k );
        }
    }

}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlxInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );


    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlx_type __curlx( this->curlx( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlx( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlx( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlx( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlx( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curly_type __curly( this->curly( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curly( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curly( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curly( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curly( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::CURL|vm::JACOBIAN|vm::KB|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        curlz_type __curlz( this->curlz( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        if ( nDim == 3 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
            {
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlz( 0,0,k );
                v[boost::get<0>( *itL )]( 1,0 ) =  __curlz( 1,0,k );
                v[boost::get<0>( *itL )]( 2,0 ) =  __curlz( 2,0,k );
            }
        }

        else if ( nDim == 2 )
        {
            itL=it->second.begin();

            for ( uint k=0; k<nbPtsElt; ++k,++itL )
                v[boost::get<0>( *itL )]( 0,0 ) =  __curlz( 0,0,k );
        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::d_( int N, ContextType const & context, id_array_type& v ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    for ( int i = 0; i < basis_type::nDof; ++i )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            size_type ldof = basis_type::nDof*c1 + i;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( i )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[q]( i,0 ) += v_*context.d( ldof, i, N, q );
                }
            }

        }
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dxInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        dx_type __dx( this->dx( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dx( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dyInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        dy_type __dy( this->dy( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dy( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dzInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::GRAD|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );
        //evaluate element for these points
        dz_type __dz( this->dz( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( typename array_type::index i = 0; i < nComponents1; ++i )
            {
                v[boost::get<0>( *itL )]( i,0 ) =  __dz( i,0,k );
            }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hess_( ContextType const & context, hess_array_type& v ) const
{
    hess_( context, v, mpl::int_<rank>() );
}
template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename ContextType>
//typename FunctionSpace<A0, A1, A2, A3, A4>::template Element<Y,Cont>::array_type
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hess_( ContextType const & context, hess_array_type& v, mpl::int_<0> ) const
{
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    //FEELPP_ASSERT( comp < nRealDim )( comp )( nRealDim ).error( "[FunctionSpace::Element] grad: invalid component" );
    for ( int i = 0; i < basis_type::nDof; ++i )
    {
        const int ncdof = is_product?nComponents1:1;

        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            size_type ldof = basis_type::nDof*c1 + i;
            size_type gdof = boost::get<0>( M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
            FEELPP_ASSERT( gdof >= this->firstLocalIndex() &&
                           gdof < this->lastLocalIndex() )
            ( context.eId() )
            ( i )( c1 )( ldof )( gdof )
            ( this->size() )( this->localSize() )
            ( this->firstLocalIndex() )( this->lastLocalIndex() )
            .error( "FunctionSpace::Element invalid access index" );

            value_type v_ = this->globalValue( gdof );

            for ( size_type q = 0; q < context.xRefs().size2(); ++q )
            {
                for ( int i = 0; i < nRealDim; ++i )
                    for ( int j = 0; j < nRealDim; ++j )
                    {
                        v[q]( i,j ) += v_*context.hess( ldof, i, j, q );
                    } // i,j
            } // q

        }
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    if ( conformalEval )
        __loc->run_analysis( __ptsReal,invalid_size_type_value, setPointsConf, mpl::int_<1>() );
    else
        __loc->run_analysis( __ptsReal,invalid_size_type_value );
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if ( it==it_end ) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>( *( it->second.begin() ) ).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ) );
    gmc_ptrtype __c( new gmc_type( __gm,
                                   this->functionSpace()->mesh()->element( it->first ),
                                   __geopc ) );
    typedef typename mesh_type::element_type geoelement_type;
    typedef typename functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::JACOBIAN|vm::KB|vm::HESSIAN|vm::FIRST_DERIVATIVE|vm::POINT, fe_type, gm_type, geoelement_type,gmc_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;
    fectx_ptrtype __ctx( new fectx_type( this->functionSpace()->fe(),
                                         __c,
                                         __pc ) );

    for ( ; it!=it_end; ++it )
    {
        nbPtsElt = it->second.size();

        //iterate in the list pt for a element
        itL=it->second.begin();
        itL_end=it->second.end();

        //compute a point matrix with the list of point
        pts= matrix_node_type( nbCoord, nbPtsElt );

        for ( size_type i=0; i<nbPtsElt; ++i,++itL )
            ublas::column( pts, i ) = boost::get<1>( *itL );

        //update geomap context
        __geopc->update( pts );
        __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
        //update precompute of basis functions
        __pc->update( pts );
        __ctx->update( __c, __pc );

        //evaluate element for these points
        hess_type __hess( this->hess( *__ctx ) );

        //update the output data
        itL=it->second.begin();

        for ( uint k=0; k<nbPtsElt; ++k,++itL )
            for ( int i = 0; i < nRealDim; ++i )
                for ( int j = 0; j < nRealDim; ++j )
                    v[boost::get<0>( *itL )]( i,j ) =  __hess( i,j,k );
    }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
                                                            ExprType const& ex,
                                                            std::string const& prefix,
                                                            GeomapStrategyType geomap_strategy,
                                                            bool accumulate,
                                                            bool verbose,
                                                            mpl::int_<MESH_ELEMENTS> )
{
    const size_type context = ExprType::context|vm::POINT;
    typedef ExprType expression_type;
    typedef Element<Y,Cont> element_type;
    // mesh element
    typedef typename functionspace_type::mesh_type::element_type geoelement_type;

    // geometric mapping context
    typedef typename functionspace_type::mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;
    typedef typename gm_type::template Context<context, geoelement_type> gm_context_type;
    typedef typename geoelement_type::gm1_type gm1_type;
    typedef typename gm1_type::template Context<context, geoelement_type> gm1_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;
    typedef boost::shared_ptr<gm1_context_type> gm1_context_ptrtype;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm_context_ptrtype> > map_gmc_type;
    typedef fusion::map<fusion::pair<vf::detail::gmc<0>, gm1_context_ptrtype> > map_gmc1_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc1_type> t_expr1_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto it = r.first;
    auto en = r.second;

    //
    // Precompute some data in the reference element for
    // geometric mapping and reference finite element
    //
    auto* __fe = this->functionSpace()->fe().get();
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( this->functionSpace()->gm(),
                                                                                         __fe->points() ) );
    typename gm1_type::precompute_ptrtype __geopc1( new typename gm1_type::precompute_type( this->mesh()->gm1(),
                                                                                            __fe->points() ) );



    const uint16_type ndofv = functionspace_type::fe_type::nDof;

    // return if no elements
    if ( it == en )
        return;

    gm_context_ptrtype __c( new gm_context_type( this->functionSpace()->gm(),*it,__geopc ) );
    gm1_context_ptrtype __c1( new gm1_context_type( this->mesh()->gm1(),*it,__geopc1 ) );

    typedef typename t_expr_type::shape shape;
    static const bool is_rank_ok = ( shape::M == nComponents1 &&
                                     shape::N == nComponents2 );

    BOOST_MPL_ASSERT_MSG( mpl::bool_<is_rank_ok>::value,
                          INVALID_TENSOR_RANK,
                          ( mpl::int_<shape::M>, mpl::int_<nComponents>, shape ) );

    map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
    t_expr_type tensor_expr( basis_type::isomorphism( ex ), mapgmc );

    map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );

    t_expr1_type tensor_expr1( basis_type::isomorphism( ex ), mapgmc1 );

    std::vector<bool> points_done( this->functionSpace()->dof()->nLocalDof()/this->nComponents );
    std::fill( points_done.begin(), points_done.end(),false );

    const int ncdof  = fe_type::is_product?nComponents:1;
    for ( ; it!=en ; ++it )
    {
        geoelement_type const& curElt = boost::unwrap_ref(*it);
        switch ( geomap_strategy )
        {
        case GeomapStrategyType::GEOMAP_HO:
        {
            __c->update( *it );
            map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
            tensor_expr.update( mapgmc );

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {

                for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                {
                    if ( accumulate )
                        this->plus_assign( curElt.id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );

                    else
                        this->assign( curElt.id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_O1:
        {
            __c1->update( *it );
            map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
            tensor_expr1.update( mapgmc1 );

            for ( uint16_type __j = 0; __j < ndofv; ++__j )
            {
                for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                    //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                {
                    if ( accumulate )
                        this->plus_assign( curElt.id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );

                    else
                        this->assign( curElt.id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );
                }
            }
        }
        break;

        case GeomapStrategyType::GEOMAP_OPT:
        {
            if ( curElt.isOnBoundary() )
            {
                // HO if on boundary
                __c->update( *it );
                map_gmc_type mapgmc( fusion::make_pair<vf::detail::gmc<0> >( __c ) );
                tensor_expr.update( mapgmc );

                for ( uint16_type __j = 0; __j < ndofv; ++__j )
                {
                    for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                        //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        if ( accumulate )
                            this->plus_assign( curElt.id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );

                        else
                            this->assign( curElt.id(), __j, c1, tensor_expr.evalq( c1, 0, __j ) );
                    }
                }
            }

            else
            {
                __c1->update( *it );
                map_gmc1_type mapgmc1( fusion::make_pair<vf::detail::gmc<0> >( __c1 ) );
                tensor_expr1.update( mapgmc1 );

                for ( uint16_type __j = 0; __j < ndofv; ++__j )
                {
                    for ( uint16_type c1 = 0; c1 < ncdof; ++c1 )
                        //for ( uint16_type c2 = 0; c2 < shape::N;++c2 )
                    {
                        if ( accumulate )
                            this->plus_assign( curElt.id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );

                        else
                            this->assign( curElt.id(), __j, c1, tensor_expr1.evalq( c1, 0, __j ) );
                    }
                }
            }
        }
        break;
        }

        //if P0 continuous finish loop here
        if ( isP0Continuous<fe_type>::result )
        {
            break;
        }
    }
} // onImpl (MESH_ELEMENTS)

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
template<typename IteratorType,  typename ExprType>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::onImpl( std::pair<IteratorType,IteratorType> const& r,
                                                            ExprType const& ex,
                                                            std::string const& prefix,
                                                            GeomapStrategyType geomap_strategy,
                                                            bool accumulate,
                                                            bool verbose,
                                                            mpl::int_<MESH_FACES> )
{
    typedef ExprType expression_type;
    // mesh element
    typedef typename element_type::functionspace_type::mesh_type::element_type geoelement_type;
    typedef typename geoelement_type::face_type face_type;

    // geometric mapping context
    typedef typename element_type::functionspace_type::mesh_type::gm_type gm_type;
    typedef boost::shared_ptr<gm_type> gm_ptrtype;
    typedef typename gm_type::template Context<context, geoelement_type> gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, gmc_ptrtype> > map_gmc_type;


    // dof
    typedef typename element_type::functionspace_type::dof_type dof_type;

    // basis
    typedef typename element_type::functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context< context, fe_type, gm_type, geoelement_type> fecontext_type;
    typedef boost::shared_ptr<fecontext_type> fecontext_ptrtype;
    //typedef fusion::map<fusion::pair<Feel::vf::detail::gmc<0>, fecontext_ptrtype> > map_gmc_type;

    // expression
    //typedef typename expression_type::template tensor<map_gmc_type,fecontext_type> t_expr_type;
    typedef typename expression_type::template tensor<map_gmc_type> t_expr_type;
    typedef typename t_expr_type::shape shape;

    //
    // start
    //
    boost::timer __timer;

    std::vector<int> dofs;
    std::vector<value_type> values;
    auto __face_it = r.first;
    auto __face_en = r.second;
    if ( __face_it != __face_en )
    {
        // get the first face properly connected
        for( ; __face_it != __face_en; ++__face_it )
            if ( __face_it->isConnectedTo0() )
                break;


        dof_type const* __dof = this->functionSpace()->dof().get();

        fe_type const* __fe = this->functionSpace()->fe().get();

        gm_ptrtype __gm( new gm_type );


        //
        // Precompute some data in the reference element for
        // geometric mapping and reference finite element
        //
        typedef typename geoelement_type::permutation_type permutation_type;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;
        DVLOG(2)  << "[elementon] numTopologicalFaces = " << geoelement_type::numTopologicalFaces << "\n";
        std::vector<std::map<permutation_type, geopc_ptrtype> > __geopc( geoelement_type::numTopologicalFaces );

        for ( uint16_type __f = 0; __f < geoelement_type::numTopologicalFaces; ++__f )
        {
            for ( permutation_type __p( permutation_type::IDENTITY );
                  __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
            {
                __geopc[__f][__p] = geopc_ptrtype(  new geopc_type( __gm, __fe->points( __f ) ) );
                //DVLOG(2) << "[geopc] FACE_ID = " << __f << " ref pts=" << __fe->dual().points( __f ) << "\n";
                FEELPP_ASSERT( __geopc[__f][__p]->nPoints() ).error( "invalid number of points" );
            }
        }

        uint16_type __face_id = __face_it->pos_first();
        gmc_ptrtype __c( new gmc_type( __gm, __face_it->element( 0 ), __geopc, __face_id ) );

        map_gmc_type mapgmc( fusion::make_pair<Feel::vf::detail::gmc<0> >( __c ) );

        DVLOG(2)  << "face_type::numVertices = " << face_type::numVertices << ", fe_type::nDofPerVertex = " << fe_type::nDofPerVertex << "\n"
                  << "face_type::numEdges = " << face_type::numEdges << ", fe_type::nDofPerEdge = " << fe_type::nDofPerEdge << "\n"
                  << "face_type::numFaces = " << face_type::numFaces << ", fe_type::nDofPerFace = " << fe_type::nDofPerFace << "\n";

        size_type nbFaceDof = invalid_size_type_value;

        if ( !fe_type::is_modal )
            nbFaceDof = ( face_type::numVertices * fe_type::nDofPerVertex +
                          face_type::numEdges * fe_type::nDofPerEdge +
                          face_type::numFaces * fe_type::nDofPerFace );

        else
            nbFaceDof = face_type::numVertices * fe_type::nDofPerVertex;

        DVLOG(2)  << "nbFaceDof = " << nbFaceDof << "\n";
        //const size_type nbFaceDof = __fe->boundaryFE()->points().size2();

        for ( ;
              __face_it != __face_en;
              ++__face_it )
        {
            if ( !__face_it->isConnectedTo0() )
            {
                LOG( WARNING ) << "face not connected" << *__face_it;

                continue;
            }
            // do not process the face if it is a ghost face: belonging to two
            // processes and being in a process id greater than the one
            // corresponding face
            if ( __face_it->isGhostFace() )
            {
                LOG(WARNING) << "face id : " << __face_it->id() << " is a ghost face";
                continue;
            }

            DVLOG(2) << "FACE_ID = " << __face_it->id()
                     << " element id= " << __face_it->ad_first()
                     << " pos in elt= " << __face_it->pos_first()
                     << " marker: " << __face_it->marker() << "\n";
            DVLOG(2) << "FACE_ID = " << __face_it->id() << " face pts=" << __face_it->G() << "\n";

            uint16_type __face_id = __face_it->pos_first();
            __c->update( __face_it->element( 0 ), __face_id );

            DVLOG(2) << "FACE_ID = " << __face_it->id() << "  ref pts=" << __c->xRefs() << "\n";
            DVLOG(2) << "FACE_ID = " << __face_it->id() << " real pts=" << __c->xReal() << "\n";

            map_gmc_type mapgmc( fusion::make_pair<Feel::vf::detail::gmc<0> >( __c ) );

            t_expr_type expr( ex, mapgmc );
            expr.update( mapgmc );

            std::pair<size_type,size_type> range_dof( std::make_pair( this->start(),
                                                                      this->functionSpace()->nDof() ) );
            DVLOG(2)  << "[elementon] dof start = " << range_dof.first << "\n";
            DVLOG(2)  << "[elementon] dof range = " << range_dof.second << "\n";

            for ( uint16_type c1 = 0; c1 < shape::M; ++c1 )
                for ( uint16_type c2 = 0; c2 < shape::N; ++c2 )
                {
                    for ( uint16_type l = 0; l < nbFaceDof; ++l )
                    {
                        DVLOG(2) << "[elementonexpr] local dof=" << l
                                 << " |comp1=" << c1 << " comp 2= " << c2 << " | pt = " <<  __c->xReal( l ) << "\n";
                        typename expression_type::value_type __value = expr.evalq( c1, c2, l );
                        DVLOG(2) << "[elementonexpr] value=" << __value << "\n";

                        // global Dof
                        size_type thedof =  this->start() +
                            boost::get<0>( __dof->faceLocalToGlobal( __face_it->id(), l, c1 ) );

                        //size_type thedof_nproc = __dof->dofNProc( thedof );
                        if ( std::find( dofs.begin(),
                                        dofs.end(),
                                        thedof ) != dofs.end() )
                            continue;

                        this->operator()( thedof ) = __value;
                    } // loop on space components

                } // loop on face dof
        }

    } // __face_it != __face_en
}
template<typename T,int M,int N>
std::ostream&
operator<<( std::ostream& os, Feel::detail::ID<T,M,N> const& id )
{
    const size_type* shape =  id.M_id.shape();

    for ( size_type i = 0; i < shape[0]; ++i )
    {
        for ( size_type j = 0; j < shape[1]; ++j )
        {
            for ( size_type k = 0; k < shape[2]; ++k )
            {
                os << id( i, j, k ) << ' ';
            }

            os << std::endl;
        }

        os << std::endl;
    }

    return os;
}



/**
   TODO: write the documentation
 */
template<typename EltType>
typename fusion::result_of::accumulate<typename EltType::functionspace_type::functionspace_vector_type,
                                       fusion::vector<>,
                                       Feel::detail::CreateElementVector<EltType> >::type
subelements( EltType const& e, std::vector<std::string> const& n )
{
    return fusion::accumulate( e.functionSpaces(), fusion::vector<>(), Feel::detail::CreateElementVector<EltType>( e, n ) );
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




#endif /* __FunctionSpace_H */
