/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2004-11-22

   Copyright (C) 2004 EPFL
   Copyright (C) 2006-2012 Universite Joseph Fourier (Grenoble I)
   Copyright (C) 2011-2021 Feel++ Consortium

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
#ifndef FEELPP_DISCR_FUNCTIONSPACE_H
#define FEELPP_DISCR_FUNCTIONSPACE_H 1

#include <type_traits>
#include <variant>
#include <boost/static_assert.hpp>

#include <boost/version.hpp>
#if BOOST_VERSION >= 106700 && BOOST_VERSION < 107100
#include <contrib/boost/fusion/include/boost/fusion/container/vector/vector.hpp>
#else
#include <boost/fusion/container/vector.hpp>
#endif
//#include <boost/fusion/container/generation/make_vector.hpp>

//#include <boost/mpl/at.hpp>
//#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/fusion/support/pair.hpp>
//#include <boost/fusion/support/is_sequence.hpp>
//#include <boost/fusion/sequence.hpp>
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


#include <feel/feeldiscr/enums.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/periodic.hpp>
#include <feel/feelpoly/expansiontypes.hpp>
#include <feel/feeldiscr/doftable.hpp>
#include <feel/feeldiscr/dofcomposite.hpp>
#include <feel/feeldiscr/parameter.hpp>
#include <feel/feeldiscr/bases.hpp>
#include <feel/feeldiscr/functionspacebase.hpp>
#include <feel/feeldiscr/mortar.hpp>
#include <feel/feeldiscr/traits.hpp>

#include <feel/feeldiscr/region.hpp>
#include <feel/feelvf/exprbase.hpp>
#include <feel/feelvf/ginac.hpp>
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
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<M,N>> m_type;
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
        {
            M_id[k].setZero();
        }
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
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<M,N>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
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
        {
            M_grad[k].setZero();
        }
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

template<typename T,int M>
struct SymmetricDD
{
    typedef T value_type;
    friend class boost::serialization::access;
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<M,M>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
    typedef boost::multi_array<m_type,1> array_type;
    typedef typename array_type::index_range range;
    struct result
    {

        typedef array_type type;
    };

    SymmetricDD()
        :
        M_grad()
    {}

    template<typename Elem, typename ContextType>
    SymmetricDD( Elem const& elem, ContextType const & context )
        :
        M_grad( elem.gradExtents( context ) )
    {
        for(int k = 0;k < M_grad.shape()[0]; k++ )
        {
            M_grad[k].setZero();
        }
        elem.symmetricGradient( context, M_grad );
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
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<M,P>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
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
        {
            //M_grad[k].resize(M,P);
            M_grad[k].setZero();
        }
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

template<typename T, int D = 1>
struct Div
{
    typedef T value_type;
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<D,1>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
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
        {
            //M_div[k].resize(D,1);
            M_div[k].setZero();
        }
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
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<(D==3)?D:1,1>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
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
        {
            //M_curl[k].resize(D,1);
            M_curl[k].setZero();
        }
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
    typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<M,N>> m_type;
    //typedef Eigen::Tensor<value_type,2> m_type;
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
        {
            //M_hess[k].resize(M,N);
            M_hess[k].setZero();
        }
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

template<typename SpaceType>
struct InitializeSpace
{
    typedef typename SpaceType::functionspace_vector_type functionspace_vector_type;
    typedef typename SpaceType::mesh_ptrtype MeshPtrType;
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;
    typedef typename SpaceType::periodicity_type PeriodicityType;
    using globaldof_type = Dof<typename SpaceType::mesh_type::size_type>;
    InitializeSpace( functionspace_vector_type & functionspaces,
                     MeshPtrType const& mesh,
                     mesh_support_vector_type const& meshSupport,
                     PeriodicityType const& periodicity,
                     std::vector<globaldof_type> const& dofindices,
                     worldscomm_ptr_t const & worldsComm,
                     std::vector<bool> extendedDofTable )
        :
        M_functionspaces( functionspaces ),
        M_cursor( 0 ),
        M_worldsComm( worldsComm ),
        M_mesh( mesh ),
        M_meshSupport( meshSupport ),
        M_dofindices( dofindices ),
        M_periodicity( periodicity ),
        M_extendedDofTable( extendedDofTable )
    {}
    template <typename T>
    void operator()( T const& t ) const
        {
            if constexpr ( is_shared_ptr<MeshPtrType>() )
            {
                typedef typename fusion::result_of::at_c<functionspace_vector_type,T::value>::type _subspace_ptrtype;
                typedef typename boost::remove_reference<_subspace_ptrtype>::type subspace_ptrtype;
                typedef typename subspace_ptrtype::element_type subspace_type;

                auto & subSpace = boost::fusion::at_c<T::value>( M_functionspaces );
                auto subMeshSupport = typename subspace_type::mesh_support_vector_type( boost::fusion::at_c<T::value>( M_meshSupport ) );
                auto p = *fusion::find<typename subspace_type::periodicity_0_type>(M_periodicity);
                subSpace = subspace_ptrtype( new subspace_type( M_mesh, subMeshSupport, M_dofindices, p,
                                                                makeWorldsComm( 1,M_worldsComm[M_cursor] ),
                                                                std::vector<bool>( 1,M_extendedDofTable[M_cursor] ) ) );
                FEELPP_ASSERT( subSpace ).error( "invalid function space" );

                ++M_cursor;// warning M_cursor < nb color
            }
            else
            {
                typedef typename fusion::result_of::at_c<functionspace_vector_type,T::value>::type _subspace_ptrtype;
                typedef typename boost::remove_reference<_subspace_ptrtype>::type subspace_ptrtype;
                typedef typename subspace_ptrtype::element_type subspace_type;

                auto & subSpace = boost::fusion::at_c<T::value>( M_functionspaces );
                auto p = *fusion::find<typename subspace_type::periodicity_0_type>(M_periodicity);
                // look for T::mesh_ptrtype in MeshPtrType
                //auto m = *fusion::find<typename subspace_type::mesh_ptrtype>(M_mesh);
                auto m = boost::fusion::at_c<T::value>( M_mesh );
                auto subMeshSupport = typename subspace_type::mesh_support_vector_type( boost::fusion::at_c<T::value>( M_meshSupport ) );
                subSpace = subspace_ptrtype( new subspace_type( m, subMeshSupport, M_dofindices, p,
                                                                makeWorldsComm( 1,M_worldsComm[M_cursor] ),
                                                                std::vector<bool>( 1,M_extendedDofTable[M_cursor] ) ) );
                FEELPP_ASSERT( subSpace ).error( "invalid function space" );

                ++M_cursor;// warning M_cursor < nb color
            }
        }
    functionspace_vector_type & M_functionspaces;
    mutable uint16_type M_cursor;
    worldscomm_ptr_t M_worldsComm;
    MeshPtrType M_mesh;
    mesh_support_vector_type const& M_meshSupport;
    std::vector<globaldof_type> const& M_dofindices;
    PeriodicityType M_periodicity;
    std::vector<bool> M_extendedDofTable;
};

template<typename DofType>
struct updateDataMapProcess
{
    updateDataMapProcess( worldscomm_ptr_t const & worldsComm,
                          worldcomm_ptr_t const& worldCommFusion,
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
    void operator()( std::shared_ptr<T> & x ) const
    {

        if ( M_worldsComm[M_cursor]->isActive() )
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

            for ( size_type i=0; i<nLocWithGhost; ++i )
                M_dm->setMapGlobalProcessToGlobalCluster( i, M_start_index + x->dof()->mapGlobalProcessToGlobalCluster( i ) );

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

        size_type startGlobClusterDof = M_dmOnOff->nLocalDofWithoutGhost() - x->dof()->nLocalDofWithoutGhost();
        size_type startGlobProcessDof = M_dmOnOff->nLocalDofWithGhost() - x->dof()->nLocalDofWithGhost();

        for ( size_type i=0; i<x->dof()->nLocalDofWithGhost(); ++i )
        {
            M_dmOnOff->setMapGlobalProcessToGlobalCluster( startGlobProcessDof + i, M_start_index + x->dof()->mapGlobalProcessToGlobalCluster( i ) );
        }


        M_start_index+=x->nDof();

        ++M_cursor;// warning M_cursor < nb color
    }

    std::shared_ptr<DofType> dataMap() const
    {
        return M_dm;
    }
    std::shared_ptr<DofType> dataMapOnOff() const
    {
        return M_dmOnOff;
    }

    mutable uint16_type M_cursor;
    mutable size_type M_start_index;
    uint16_type M_lastCursor;
    worldscomm_ptr_t M_worldsComm;
    mutable std::shared_ptr<DofType> M_dm;
    mutable std::shared_ptr<DofType> M_dmOnOff;
}; // updateDataMapProcess

template<typename DofType>
struct updateDataMapProcessStandard
{
    typedef std::shared_ptr<DofType> result_type;

    updateDataMapProcessStandard( worldcomm_ptr_t const& worldComm,
                                  uint16_type nSpaces )
        :
        M_worldComm( worldComm ),
        M_cursor( 0 ),
        M_lastCursor( nSpaces-1 )
    {}

    template <typename T>
    result_type operator()( result_type const& r, std::shared_ptr<T> & x ) const
    {
        M_subdm.push_back( x->mapPtr() );
        if ( M_cursor == M_lastCursor )
        {
            result_type dm = std::make_shared<DofType>( M_subdm,M_worldComm );
            return dm;
        }
        ++M_cursor;
        return r;
    }

    worldcomm_ptr_t M_worldComm;
    mutable uint16_type M_cursor;
    uint16_type M_lastCursor;
    mutable std::vector<datamap_ptrtype<>> M_subdm;
};





struct NbDof
{
    typedef size_type result_type;
    NbDof( size_type start = 0, size_type size = invalid_v<size_type> )
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
    NLocalDof( size_type start = 0, size_type size = invalid_v<size_type> )
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

    NLocalDof( worldscomm_ptr_t const & worldsComm = Environment::worldsComm(1),
               bool useOffSubSpace = false,
               size_type start = 0, size_type size = invalid_v<size_type> )
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
                if ( M_worldsComm[M_cursor]->isActive() )
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
    worldscomm_ptr_t M_worldsComm;
    bool M_useOffSubSpace;
};
#endif // end MPI


template< typename IsWithGhostType>
struct NLocalDofOnProc
{

    NLocalDofOnProc( const int proc,
                     worldscomm_ptr_t const & worldsComm = Environment::worldsComm(1),
                     bool useOffSubSpace = false,
                     size_type start = 0, size_type size = invalid_v<size_type> )
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
                if ( M_worldsComm[M_cursor]->isActive() )
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
    worldscomm_ptr_t M_worldsComm;
    bool M_useOffSubSpace;
}; // NLocalDofOnProc


template<int i,typename SpaceCompositeType>
struct InitializeContainersOff
{
    explicit InitializeContainersOff( std::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
    {
        if ( M_cursor==i && !x )
            x = std::shared_ptr<T>( new T( M_space->template functionSpace<i>()->dof() ) );

        ++M_cursor;// warning M_cursor < nb color
    }
    mutable uint16_type M_cursor;
    std::shared_ptr<SpaceCompositeType> M_space;
};


template<int i,typename SpaceCompositeType>
struct SendContainersOn
{
    SendContainersOn( std::shared_ptr<SpaceCompositeType> const& _space,
                      std::vector<double> const& _dataToSend )
        :
        M_cursor( 0 ),
        M_space( _space ),
        M_dataToSend( _dataToSend )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
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
    std::shared_ptr<SpaceCompositeType> M_space;
    std::vector<double> M_dataToSend;
};


template<int i,typename SpaceCompositeType>
struct RecvContainersOff
{
    explicit RecvContainersOff( std::shared_ptr<SpaceCompositeType> const& _space )
        :
        M_cursor( 0 ),
        M_space( _space )
    {}
    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
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
    std::shared_ptr<SpaceCompositeType> M_space;
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

struct hasSubSpaceWithComponentsSplit
{
    typedef bool result_type;
    template<typename T>
    result_type operator()( result_type const &  previousRes, T const& t )
    {
        //return ( T::element_type::dof_type::is_product && T::element_type::dof_type::nComponents > 1 ) || previousRes;
        return t->map().hasIndexSplitWithComponents() || previousRes;
    }
};

// compute split
template<bool UseComponentsSplit>
struct computeNDofForEachSpace
{
    computeNDofForEachSpace(size_type startSplit)
        :
        M_indexSplit( new IndexSplit() ),
        M_startSplit(startSplit)
    {}

    std::shared_ptr<IndexSplit> const& indexSplit() const { return M_indexSplit; }

    typedef boost::tuple< uint16_type, size_type, IndexSplit > result_type;

    template<typename T>
    void operator()( T const& t ) const
    {
        this->operator()( t, mpl::bool_<UseComponentsSplit>() );
    }
    template<typename T>
    void operator()( T const& t, mpl::false_ ) const
    {
        M_indexSplit->addSplit( M_startSplit, t->map().indexSplit() );
    }
    template<typename T>
    void operator()( T const& t, mpl::true_ ) const
    {
        M_indexSplit->addSplit( M_startSplit, t->map().indexSplitWithComponents() );
    }

    mutable std::shared_ptr<IndexSplit> M_indexSplit;
    size_type M_startSplit;
};

struct rebuildDofPointsTool
{

    template <typename T>
    void operator()( std::shared_ptr<T> & x ) const
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

template<typename SpaceType>
struct createWorldsComm
{

    typedef typename SpaceType::mesh_ptrtype mesh_ptrtype;
    typedef typename SpaceType::meshes_list meshes_list;
    static const bool useMeshesList = !boost::is_base_of<MeshBase<>, meshes_list >::value;

    struct UpdateWorldsComm
    {
        UpdateWorldsComm( createWorldsComm<SpaceType> & cwc )
            :
            M_cwc( cwc )
            {}
        template<typename T>
        void operator()( T const& t) const
            {
                M_cwc.M_worldsComm.push_back( t->worldComm().shared_from_this() );
            }
        createWorldsComm<SpaceType> & M_cwc;
    };

    createWorldsComm( mesh_ptrtype const& mesh )
        {
            this->init<useMeshesList>( mesh );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< !_UseMeshesList >::type* = nullptr )
        {
            M_worldsComm.resize( SpaceType::nSpaces, mesh->worldComm().shared_from_this() );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< _UseMeshesList >::type* = nullptr )
        {
            boost::fusion::for_each( mesh, UpdateWorldsComm( *this ) );
        }
    worldscomm_ptr_t worldsComm()       { return M_worldsComm; }
    worldscomm_ptr_t worldsComm() const { return M_worldsComm; }

    worldscomm_ptr_t M_worldsComm;
};

template<typename SpaceType>
std::vector<bool>
createInfoExtendedDofTable( bool b )
{
    return std::vector<bool>( SpaceType::nSpaces,b );
}
template<typename SpaceType>
std::vector<bool>
createInfoExtendedDofTable( std::vector<bool> const& b )
{
    CHECK( b.size() == SpaceType::nSpaces ) << "invalid extended doftable info vector size : " << b.size() << " should be : " << SpaceType::nSpaces;
    return b;
}

template<typename SpaceType>
struct createMeshSupport
{
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;
    typedef typename fusion::result_of::at_c<mesh_support_vector_type,0>::type _mesh_support_ptrtype;
    typedef typename boost::remove_reference<_mesh_support_ptrtype>::type mesh_support_ptrtype;
    typedef typename mesh_support_ptrtype::element_type mesh_support_type;
    typedef typename SpaceType::mesh_ptrtype mesh_ptrtype;
    typedef typename mesh_support_type::range_elements_type range_elements_type;

    typedef typename SpaceType::meshes_list meshes_list;
    static const bool useMeshesList = !boost::is_base_of<MeshBase<>, meshes_list >::value;

    struct HasAllMeshSupportDefined
    {
        typedef bool result_type;
        template<typename T>
        result_type operator()( result_type const& r, T const& t) const
            {
                if ( !t )
                    return false;
                else
                    return r;
            }
    };

    struct UpdateMeshSupport
    {
        UpdateMeshSupport( createMeshSupport<SpaceType> & cms )
            :
            M_cms( cms )
            {}
        template<typename T>
        void operator()( T const& t) const
            {
                this->updateImpl<T,useMeshesList>( t );
            }
        template<typename T,bool _UseMeshesList >
        void updateImpl( T const& t, std::enable_if_t< !_UseMeshesList >* = nullptr ) const
            {
                auto & meshSupport = boost::fusion::at_c<T::value>( M_cms.M_meshSupportVector );
                if ( meshSupport )
                    return;
                CHECK( M_cms.M_meshSupport0 ) << "no mesh support defined";
                meshSupport = M_cms.M_meshSupport0;
            }
        template<typename T,bool _UseMeshesList >
        void updateImpl( T const& t, std::enable_if_t< _UseMeshesList >* = nullptr ) const
            {
                auto & meshSupport = boost::fusion::at_c<T::value>( M_cms.M_meshSupportVector );
                if ( meshSupport )
                    return;
                typedef typename fusion::result_of::at_c<mesh_support_vector_type,T::value>::type _submesh_support_ptrtype;
                typedef typename boost::remove_reference<_submesh_support_ptrtype>::type submesh_support_ptrtype;
                typedef typename submesh_support_ptrtype::element_type submesh_support_type;
                auto const& mesh = boost::fusion::at_c<T::value>( M_cms.M_mesh );
                meshSupport.reset( new submesh_support_type(mesh) );
            }

        createMeshSupport<SpaceType> & M_cms;
    };

    createMeshSupport( mesh_ptrtype const& mesh, mesh_support_vector_type const& meshSupport )
        :
        M_mesh( mesh ),
        M_meshSupportVector( meshSupport )
        {
            this->init<useMeshesList>(mesh);
        }
    template<typename RangeType, typename std::enable_if_t<is_range_v<RangeType>,int> = 0 >
    createMeshSupport( mesh_ptrtype const& mesh, RangeType && rangeMeshElt )
        :
        M_mesh( mesh )
        {
            this->init2<useMeshesList>(mesh,std::forward<RangeType>(rangeMeshElt));
        }
    createMeshSupport( mesh_ptrtype const& mesh, mesh_support_ptrtype const& meshSupport )
        :
        M_mesh( mesh )
        {
            M_meshSupport0 = meshSupport;
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }

    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< !_UseMeshesList >::type* = nullptr )
        {
            HasAllMeshSupportDefined hasMSFunctor;
            bool hasMS = boost::fusion::fold( M_meshSupportVector, true, hasMSFunctor );
            if ( !hasMS )
                M_meshSupport0.reset( new mesh_support_type(mesh) );

            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }
    template<bool _UseMeshesList >
    void init( mesh_ptrtype const& mesh, typename std::enable_if< _UseMeshesList >::type* = nullptr )
        {
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }
    template<bool _UseMeshesList, typename RangeType >
    void init2( mesh_ptrtype const& mesh, RangeType && rangeMeshElt )
        {
            if constexpr ( _UseMeshesList )
                CHECK( false ) << fmt::format( "MeshSupport not allowed in Mesh List" );
            else
            {
                if ( std::forward<RangeType>( rangeMeshElt ).container() )
                    M_meshSupport0.reset( new mesh_support_type(mesh,std::forward<RangeType>(rangeMeshElt) ) );
                else
                    M_meshSupport0.reset( new mesh_support_type(mesh) );

                mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
                boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
            }
        }

    mesh_ptrtype const& M_mesh;
    mesh_support_vector_type M_meshSupportVector;
    mesh_support_ptrtype M_meshSupport0;
};

template<typename SpaceType>
struct FunctionSpaceMeshSupport
{
    typedef typename SpaceType::mesh_support_vector_type mesh_support_vector_type;

    struct UpdateMeshSupport
    {
        UpdateMeshSupport( FunctionSpaceMeshSupport<SpaceType> & fsms )
            :
            M_fsms( fsms )
            {}

        template<typename T>
        void operator()( T const& t) const
            {
                this->updateImpl<T,SpaceType::is_composite>( t );
            }

        template<typename T,bool _IsComposite >
        void updateImpl( T const& t, typename std::enable_if< !_IsComposite >::type* = nullptr ) const
            {
                auto doftable = M_fsms.M_space.dof();
                if ( !doftable )
                    return;
                if ( doftable->hasMeshSupport() )
                {
                    auto & meshSupport = boost::fusion::at_c<T::value>( M_fsms.M_meshSupportVector );
                    meshSupport = doftable->meshSupport();
                }
            }
        template<typename T,bool _IsComposite >
        void updateImpl( T const& t, typename std::enable_if< _IsComposite >::type* = nullptr ) const
            {
                auto subspace = M_fsms.M_space.template functionSpace<T::value>();
                if ( !subspace )
                    return;
                auto doftable = subspace->dof();
                if ( !doftable )
                    return;
                if ( doftable->hasMeshSupport() )
                {
                    auto & meshSupport = boost::fusion::at_c<T::value>( M_fsms.M_meshSupportVector );
                    meshSupport = doftable->meshSupport();
                }
            }
        FunctionSpaceMeshSupport<SpaceType> & M_fsms;
    };

    FunctionSpaceMeshSupport( SpaceType const& space )
        :
        M_space( space )
        {
            mpl::range_c<int,0,SpaceType::nSpaces> keySpaces;
            boost::fusion::for_each( keySpaces, UpdateMeshSupport( *this ) );
        }

    SpaceType const& M_space;
    mesh_support_vector_type M_meshSupportVector;
};

} // detail


template<uint16_type PN,
         uint16_type GN = 1>
struct Order
{
    static inline const uint16_type PolynomialOrder = PN;
    static inline const uint16_type GeometricOrder = GN;

    static const bool is_isoparametric = ( PN == GN );
    static const bool is_subparametric = ( PN > GN );
    static const bool is_surparametric = ( PN < GN );
};

typedef parameter::parameters<
    parameter::required<tag::mesh_type, boost::is_base_and_derived<MeshBase<>,_> >
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
    typename MeshTypes,
    typename BasisTypes = Feel::bases<Lagrange<1,Scalar> >,
    typename T = double,
    typename PeriodicityType = Periodicity<NoPeriodicity>,
    typename MortarType = mortars<NoMortar>>
class FunctionSpace
    :
    public FunctionSpaceBase,
    public std::enable_shared_from_this<FunctionSpace<MeshTypes,BasisTypes,T,PeriodicityType,MortarType> >
{
public:

    using meshes_list = MeshTypes;
    using bases_list = BasisTypes;
    using value_type = T;
    using periodicity_type = PeriodicityType;
    using mortar_list = MortarType;
    using mortar_type = mortar_list;

    static_assert(!mp11::mp_same<mp11::mp_at_c<bases_list, 0>, void>::value, "The first type in bases_list should not be void");

public:

    using super = FunctionSpaceBase;

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
        typedef typename fusion::result_of::distance<typename fusion::result_of::begin<bases_list_noref>::type,
                                                     typename fusion::result_of::find<bases_list_noref,BasisType>::type>::type pos;
        typedef typename fusion::result_of::at_c<meshes_list_noref, pos::value >::type _mesh_type;

        typedef FunctionSpace<typename boost::remove_reference<_mesh_type>::type,
                              Feel::detail::bases<BasisType>,value_type,
                              Periodicity<typename GetPeriodicity<periodicity_type,pos::value>::type >,
                              mortar_list> _type;
        typedef std::shared_ptr<_type> type;
    };
    template<typename BasisType>
    struct ChangeBasis
    {
        //typedef typename mpl::if_<mpl::and_<boost::is_base_of<MeshBase<>, meshes_list >, boost::is_base_of<Feel::detail::periodic_base, periodicity_type > >,
        typedef typename boost::remove_reference<bases_list>::type bases_list_noref;
        typedef typename fusion::result_of::distance<typename fusion::result_of::begin<bases_list_noref>::type,
                                                     typename fusion::result_of::find<bases_list_noref,BasisType>::type>::type pos;

        typedef typename mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                                  mpl::identity<mpl::identity<std::shared_ptr<FunctionSpace<meshes_list,Feel::detail::bases<BasisType>,value_type, Periodicity<typename GetPeriodicity<periodicity_type,pos::value>::type>, mortars<typename GetMortar<mortar_list,pos::value>::type > > > > >,
                                  mpl::identity<ChangeMesh<BasisType> > >::type::type::type type;

//mpl::identity<typename mpl::transform<meshes_list, ChangeMesh<mpl::_1,BasisType>, mpl::back_inserter<fusion::vector<> > >::type > >::type::type type;
    };
    typedef typename mpl::transform<bases_list, ChangeBasis<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type functionspace_vector_type;

    template<typename BasisType>
    struct ChangeMeshToComponentBasis
    {
        typedef typename boost::remove_reference<meshes_list>::type meshes_list_noref;
        typedef typename boost::remove_reference<bases_list>::type bases_list_noref;
        typedef typename fusion::result_of::distance<typename fusion::result_of::begin<bases_list_noref>::type,
                                                     typename fusion::result_of::find<bases_list_noref,BasisType>::type>::type pos;
        typedef typename fusion::result_of::at_c<meshes_list_noref, pos::value >::type _mesh_type;

        typedef typename BasisType::component_basis_type component_basis_type;

        typedef FunctionSpace<typename boost::remove_reference<_mesh_type>::type,
                              Feel::detail::bases<component_basis_type>,value_type,
                              Periodicity<typename GetPeriodicity<periodicity_type,pos::value>::type >,
                              mortar_list> _type;
        typedef std::shared_ptr<_type> type;
    };

    template<typename BasisType>
    struct ChangeBasisToComponentBasis
    {
        typedef typename BasisType::component_basis_type component_basis_type;
        //typedef typename mpl::if_<mpl::and_<boost::is_base_of<MeshBase<>, meshes_list >, boost::is_base_of<Feel::detail::periodic_base, periodicity_type > >,
        typedef typename mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                                  mpl::identity<mpl::identity<std::shared_ptr<FunctionSpace<meshes_list,Feel::detail::bases<component_basis_type>,value_type, periodicity_type, mortar_list> > > >,
                                  mpl::identity<ChangeMeshToComponentBasis<BasisType> > >::type::type::type type;
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
        typedef typename mpl::if_<mpl::or_<boost::is_base_of<MeshBase<>, MeshListType >,
                                           is_shared_ptr<MeshListType> >,
                                  mpl::identity<mpl::identity<MeshListType> >,
                                  mpl::identity<mpl::at_c<MeshListType,N> > >::type::type::type type;
    };
    template<typename MeshListType,int N>
    struct GetMeshSupport
    {
        typedef typename GetMesh<MeshListType,N>::type mesh_ptrtype;
        typedef MeshSupport<typename mesh_ptrtype::element_type> type;
        typedef typename type::range_elements_type range_type;
        typedef std::shared_ptr<type> ptrtype;
    };
    // mesh
    typedef meshes_list MeshesListType;
    typedef typename GetMesh<meshes_list,0>::type mesh_0_type;
    using index_type = typename mesh_0_type::index_type;
    using size_type = typename mesh_0_type::size_type;
    typedef typename mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                              mpl::identity<meshes_list>,
                              mpl::identity<mesh_0_type> >::type::type mesh_type;

    template<typename MeshType>
    struct ChangeToMeshPtr
    {
        typedef std::shared_ptr<MeshType> type;
    };

    typedef typename mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                              mpl::identity<std::shared_ptr<mesh_type> >,
                              mpl::identity<typename mpl::transform<meshes_list, ChangeToMeshPtr<mpl::_1>, mpl::back_inserter<meshes<> > >::type  > >::type::type mesh_ptrtype;
    typedef typename mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
            mpl::identity<typename mesh_type::element_type>,
            mpl::identity<typename mesh_0_type::element_type> >::type::type convex_type;

    template<typename SpaceType>
    struct ChangeMeshSupport
    {
        typedef typename SpaceType::element_type::mesh_type mesh_type;
        typedef MeshSupport<mesh_type> mesh_support_type;
        typedef std::shared_ptr<mesh_support_type> mesh_support_ptrtype;
        typedef mesh_support_ptrtype type;
    };
    typedef typename mpl::transform<functionspace_vector_type, ChangeMeshSupport<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type mesh_support_vector_type;


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
    static constexpr uint16_type nDim = mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                                          mpl::identity<meshes_list >,
                                          mpl::identity<nodim> >::type::type::nDim;
    static constexpr uint16_type nRealDim = mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
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
    static constexpr bool is_tensor2 = ( is_composite? false : basis_0_type::is_tensor2  && !is_symm_v<basis_0_type> );
    static constexpr bool is_tensor2symm = ( is_composite? false : basis_0_type::is_tensor2 && is_symm_v<basis_0_type> );
    static constexpr bool is_continuous = ( is_composite? false : basis_0_type::isContinuous );
    static constexpr bool is_modal = ( is_composite? false : basis_0_type::is_modal );
    static constexpr bool is_nodal = !is_modal;


    static constexpr uint16_type nComponents1 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents1 );
    static constexpr uint16_type nComponents2 = ( is_composite? invalid_uint16_type_value : basis_0_type::nComponents2 );
    static constexpr bool is_product = ( is_composite? false : basis_0_type::is_product );
    typedef typename  basis_0_type::continuity_type continuity_type;

    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<boost::none_t>,
                              mpl::identity<typename basis_0_type::polyset_type> >::type::type polyset_type;

    static constexpr uint16_type nComponents = mpl::transform<bases_list,
                             GetNComponents<mpl::_1>,
                             mpl::inserter<mpl::int_<0>,mpl::plus<mpl::_,mpl::_> > >::type::value;
    static constexpr uint16_type N_COMPONENTS = nComponents;
    static constexpr uint16_type nSpaces = mpl::size<bases_list>::type::value;
    static constexpr uint16_type nRealComponents = is_tensor2symm?basis_0_type::nComponents1*(basis_0_type::nComponents1+1)/2:nComponents;
    typedef typename GetPeriodicity<periodicity_type,0>::type periodicity_0_type;
    static constexpr bool is_periodic = periodicity_0_type::is_periodic;

    typedef typename GetMortar<mortar_list,0>::type mortar_0_type;
    static const bool is_mortar = mortar_0_type::is_mortar;
    static const bool is_hdiv_conforming = Feel::is_hdiv_conforming<basis_0_type>::value;
    static const bool is_hcurl_conforming = Feel::is_hcurl_conforming<basis_0_type>::value;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef typename ublas::type_traits<value_type>::real_type real_type;
    typedef typename node<value_type>::type node_type;


    typedef FunctionSpace<MeshTypes,BasisTypes,T,PeriodicityType,MortarType> functionspace_type;
    typedef functionspace_type space_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef std::shared_ptr<functionspace_type> pointer_type;

    typedef FunctionSpace<meshes_list, component_basis_vector_type, value_type, periodicity_type, mortar_list> component_functionspace_type;
    typedef std::shared_ptr<component_functionspace_type> component_functionspace_ptrtype;


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


    typedef std::shared_ptr<basis_type> basis_ptrtype;
    typedef basis_type reference_element_type;
    typedef std::shared_ptr<reference_element_type> reference_element_ptrtype;
    typedef reference_element_type fe_type;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<fe_type>,
                              mpl::identity<typename basis_0_type::SSpace::type> >::type::type mortar_fe_type;
    typedef reference_element_ptrtype fe_ptrtype;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
            mpl::identity<boost::none_t>,
            mpl::identity<typename basis_0_type::PreCompute> >::type pc_type;
    typedef std::shared_ptr<pc_type> pc_ptrtype;

    /**
     * interpolate type if available
     */
    using local_interpolant_type = local_interpolant_t<basis_0_type>;

    using local_interpolants_type = local_interpolants_t<basis_0_type>;

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
    typedef std::shared_ptr<component_basis_type> component_basis_ptrtype;

    // trace space
    using is_trace_applicable = mp11::mp_bool<(nDim > 1)>;

    using trace_mesh_type = mp11::mp_if<is_trace_applicable, typename mesh_type::template trace_mesh_type<>, void>;
    using trace_mesh_ptrtype = mp11::mp_if<is_trace_applicable, typename mesh_type::template trace_mesh_ptrtype<>, void>;
    using trace_functionspace_type = mp11::mp_if<is_trace_applicable, FunctionSpace<trace_mesh_type, bases_list>, void>;
    using trace_functionspace_ptrtype = std::shared_ptr<trace_functionspace_type>;


    // wirebasket
    using is_trace_trace_applicable = mp11::mp_bool<(nDim > 2)>;

    using trace_trace_mesh_type = mp11::mp_if<is_trace_trace_applicable, typename mesh_type::template trace_trace_mesh_type<>, void>;
    using trace_trace_mesh_ptrtype = mp11::mp_if<is_trace_trace_applicable, typename mesh_type::template trace_trace_mesh_ptrtype<>, void>;
    using trace_trace_functionspace_type = mp11::mp_if<is_trace_trace_applicable, FunctionSpace<trace_trace_mesh_type, bases_list>, void>;
    using trace_trace_functionspace_ptrtype = std::shared_ptr<trace_trace_functionspace_type>;

#if 0
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<1> >,
            mpl::identity<typename trace_functionspace_type::element_type>,
            mpl::identity<mpl::void_> >::type::type trace_element_type;
#endif

    // geomap
    //typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm_type>, mpl::identity<mpl::void_> >::type::type gm_type;
    //typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm1_type>, mpl::identity<mpl::void_> >::type::type gm1_type;

    using gm_type = typename mesh_type::gm_type;
    using gm1_type = typename mesh_type::gm1_type;
    //typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::element_type>, mpl::identity<mpl::void_> >::type::type geoelement_type;
    using geoelement_type = typename mesh_type::element_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::face_type>, mpl::identity<mpl::void_> >::type::type geoface_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::edge_type>, mpl::identity<mpl::void_> >::type::type geoedge_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::point_type>, mpl::identity<mpl::void_> >::type::type geopoint_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;
    //typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context<vm::POINT, geoelement_type> >,
    //mpl::identity<mpl::void_> >::type::type pts_gmc_type;

    using pts_gmc_type = typename gm_type::template Context</*vm::POINT,*/ geoelement_type>;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context</*vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB,*/ geoelement_type> >,
                              mpl::identity<mpl::void_> >::type::type gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_ptrtype>, mpl::identity<mpl::void_> >::type::type geopc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_type>, mpl::identity<mpl::void_> >::type::type geopc_type;

    // basis context
    //typedef typename basis_type::template Context<vm::POINT, basis_type, gm_type, geoelement_type, pts_gmc_type::context> basis_context_type;
#if 0
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<mpl::void_>,
                              mpl::identity<typename basis_0_type::template Context<vm::POINT|vm::GRAD|vm::JACOBIAN|vm::Hessian|vm::KB, basis_0_type, gm_type, geoelement_type> > >::type::type basis_context_type;
                              //mpl::identity<typename basis_0_type::template Context<vm::POINT, basis_0_type, gm_type, geoelement_type, pts_gmc_type::context> > >::type::type basis_context_type;
#else
    static const size_type basis_context_value = ( basis_0_type::is_product && nComponents > 1 )?
        vm::POINT|vm::GRAD|vm::JACOBIAN|vm::KB :
        vm::POINT|vm::GRAD|vm::JACOBIAN|vm::HESSIAN|vm::KB;
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<mpl::void_>,
                              mpl::identity<typename basis_0_type::template Context<basis_context_value, basis_0_type, gm_type, geoelement_type> > >::type::type basis_context_type;
#endif
    typedef std::shared_ptr<basis_context_type> basis_context_ptrtype;

    // dof
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<DofComposite>,
                              mpl::identity<DofTable<mesh_type, basis_type, periodicity_0_type, mortar_0_type> > >::type::type dof_type;

    typedef std::shared_ptr<dof_type> dof_ptrtype;
    typedef std::shared_ptr<DataMap<>> datamap_ptrtype;
    typedef std::shared_ptr<IndexSplit> indexsplit_ptrtype;

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
    template<int i> using sub_functionspace_type = typename sub_functionspace<i>::type;
    template<int i> using sub_functionspace_ptrtype = typename sub_functionspace<i>::ptrtype;
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
        public std::map<int,std::pair<basis_context_ptrtype,std::vector<index_type>>>
    {
    public:
        static const bool is_rb_context = false;
        //typedef std::map<int,basis_context_ptrtype> super;
        using super = std::map<int,std::pair<basis_context_ptrtype,std::vector<index_type>>>;
        typedef typename super::value_type bc_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;
        typedef typename super::iterator iterator;
        Context( functionspace_ptrtype const& Xh ) : M_Xh( Xh ) {}
        Context() = default;
        Context( Context const& ) = default;
        Context& operator=( Context const& ) = default;
        virtual ~Context() {}

        std::pair<iterator, bool>
        add( node_type const& t )
        {
            return update( t, -1, mpl::bool_<is_composite>() );
        }
        std::pair<iterator, bool>
        replace( int ptIdInCtx, node_type const& t )
        {
            return update( t, ptIdInCtx, mpl::bool_<is_composite>() );
        }
    private :
        std::pair<iterator, bool>
        update( node_type const& t, int ptIdInCtx, mpl::bool_<true> )
        {
            return std::make_pair( this->end(), false );
        }
        std::pair<iterator, bool>
        update( node_type const& t, int ptIdInCtx, mpl::bool_<false> )
        {
            std::pair<iterator, bool> ret = std::make_pair(this->end(),false);
            //LOG(INFO)<<"add point\n";

            //rank of the current processor
            rank_type procId = this->functionSpace()->worldComm().globalRank();
            //total number of processors
            rank_type nprocs = this->functionSpace()->worldComm().globalSize();

            // add point t to list of points
            if ( ptIdInCtx < 0 )
                M_t.push_back( t );
            else
            {
                CHECK( ptIdInCtx <= (M_t.size()-1) ) << "ptIdInCtx invalid " << ptIdInCtx << " and " << M_t.size();
                auto const& ptRegister = M_t[ptIdInCtx];
                bool ptsAreIdentical = true;
                for (uint16_type d=0;d<mesh_type::nRealDim;++d)
                    ptsAreIdentical = ptsAreIdentical && (std::abs( ptRegister[d]-t[d] )<1e-14);
                // if pt are identical, do nothing and keep the context
                auto itFindCtx = this->find( ptIdInCtx );
                if ( ptsAreIdentical )
                    return std::make_pair( itFindCtx, true );
                else
                    M_t[ptIdInCtx] = t;
                // erase previous context stored
                if ( itFindCtx != this->end() )
                    this->erase( ptIdInCtx );
            }

            // localise t in space, find geometrical element in which t
            // belongs
            matrix_node_type m( mesh_type::nRealDim, 1 );
            for(int i = 0; i < mesh_type::nRealDim; ++i )
                m(i,0) = t(i);
            auto loc =  M_Xh->mesh()->tool_localization();
            loc->setExtrapolation( false );
            auto analysis = loc->run_analysis( m, invalid_v<index_type> );
            auto found_points = analysis.template get<0>();
            bool found = found_points[0];

            std::vector<uint8_type> found_pt( nprocs, 0 );
            if ( found )
                found_pt[procId] = 1;

            std::vector<uint8_type> global_found_pt( nprocs, 0 );
            if ( nprocs > 1 )
                mpi::all_reduce( M_Xh->mesh()->comm(), found_pt.data(), found_pt.size(), global_found_pt.data(), std::plus<uint8_type>() );
            else
                global_found_pt[ procId ] = found_pt[ procId ];
            // only one proc has the point
            bool findOneProcess = false;
            for ( rank_type p=0;p<nprocs;++p )
            {
                if ( findOneProcess )
                    global_found_pt[p] = 0;
                if ( global_found_pt[p] == 1 )
                    findOneProcess=true;
            }

            if ( global_found_pt[ procId ] == 1 ) //we are on the proc that have the searched point
            {
                auto it = loc->result_analysis_begin();
                auto en = loc->result_analysis_end();
                DCHECK( boost::next(it) == en ) << "Logic problem in finding one point in the mesh\n";
                auto eid = it->first;
                auto xref = boost::get<1>( *(it->second.begin()) );
                DVLOG(2) << "found point " << t << " in element " << eid << " on proc "<< procId <<"\n";
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
                auto gmc = M_Xh->mesh()->gm()->template context<basis_context_value/*basis_context_type::gmc_type::context*/>( M_Xh->mesh()->element( eid ), gmpc );
                DVLOG(2) << "build geometric mapping context\n";
                // compute finite element context
                auto ctx = basis_context_ptrtype( new basis_context_type( M_Xh->basis(), gmc, basispc ) );
                DVLOG(2) << "build basis function context\n";

                int number = (ptIdInCtx < 0)? M_t.size()-1 : ptIdInCtx;
                std::vector<index_type> ptIds;ptIds.push_back( number );
                ret = this->insert( std::make_pair( number , std::make_pair( ctx, std::move(ptIds) ) ) );
                //DVLOG(2) << "Context size: " << this->size() << "\n";

            }//if( found )

            //verify that the point is on a proc
            bool found_on_a_proc = false;
            for (int i = 0 ; i < global_found_pt.size(); ++i )
            {
                if ( global_found_pt[i] == 1 )
                {
                    DVLOG(2) << "processor " << i << " has the point " << t << "\n";
                    found_on_a_proc = true;
                    if ( ptIdInCtx < 0 )
                        M_t_proc.push_back(i);
                    else
                        M_t_proc[ptIdInCtx] = i;
                    break;
                }
            }
            CHECK( found_on_a_proc ) << "the point " << t << " was not found ! \n";
            return ret;

        }//add ( non composite case )

    public :
        void addCtx( std::pair<basis_context_ptrtype,std::vector<index_type>> const& /*basis_context_ptrtype*/ ctx , int proc_having_the_point)
        {
            int position = M_t.size();
            this->insert( std::make_pair( position, ctx ) );
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

        constexpr bool ctxHaveBeenMpiBroadcasted() const { return false; }

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
    template<typename TT = double,  typename Cont = VectorUblas<TT> >
    class Element
        :
        public Cont, boost::addable<Element<TT,Cont> >, boost::subtractable<Element<TT,Cont> >, FunctionSpaceBase::ElementBase, basis_0_type::polyset_type
    {
    public:
        typedef TT value_type;

        using functionspace_type = FunctionSpace<MeshTypes,BasisTypes,T,PeriodicityType,MortarType>;
        friend class FunctionSpace<MeshTypes,BasisTypes,T,PeriodicityType,MortarType>;

        template<typename BasisType,typename keyType>
        struct ChangeElement
        {
            typedef TT value_type;
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef typename ChangeBasis<BasisType>::type::element_type fs_type;
            //typedef typename fs_type::template Element<value_type, typename Cont::range::type > element_type;
            typedef typename fs_type::template Element<value_type, Cont > element_type;
            typedef std::pair< keyType, std::shared_ptr<element_type> > the_type;

            typedef typename mpl::if_<mpl::bool_<functionspace_type::is_composite>,
                                      mpl::identity<the_type>,
                                      mpl::identity<boost::none_t> >::type::type type;
        };

        typedef mpl::range_c<int,0, functionspace_type::nSpaces> rangeElementStorageType;
        typedef typename mpl::transform<bases_list, rangeElementStorageType, ChangeElement<mpl::_1,mpl::_2>, mpl::back_inserter<fusion::vector<> > >::type element_vector_type;

        //typedef typename fusion::result_of::accumulate<bases_list, fusion::vector<>, ChangeElement<> >
        //typedef typename Cont[>VectorUblas<T><]::range::type ct_type;
        typedef Cont ct_type;

        typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> eigen_type;

        /**
         * useful in // with composite case
         * store the off views
         */
        template<typename BasisType>
        struct AddOffContainer
        {
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef std::shared_ptr<Cont> type;
        };
        typedef typename mpl::transform<bases_list, AddOffContainer<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type container_vector_type;


        using mesh_type = typename functionspace_type::mesh_type;
        using mesh_ptrtype = typename functionspace_type::mesh_ptrtype;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
        using index_type = typename mesh_type::index_type;

        static constexpr uint16_type nDim = mesh_type::nDim;
        static constexpr uint16_type nRealDim = mesh_type::nRealDim;
        static constexpr bool is_composite = functionspace_type::is_composite;
        static constexpr bool is_scalar = functionspace_type::is_scalar;
        static constexpr bool is_vectorial = functionspace_type::is_vectorial;
        static constexpr bool is_tensor2 = functionspace_type::is_tensor2;
        static constexpr bool is_tensor2symm = functionspace_type::is_tensor2symm;
        static constexpr bool is_continuous = functionspace_type::is_continuous;
        static constexpr uint16_type nComponents1 = functionspace_type::nComponents1;
        static constexpr uint16_type nComponents2 = functionspace_type::nComponents2;
        static constexpr uint16_type nComponents = functionspace_type::nComponents;
        static constexpr uint16_type nRealComponents = functionspace_type::nRealComponents;
        static constexpr uint16_type nSpaces = functionspace_type::nSpaces;
        static constexpr bool is_mortar = functionspace_type::is_mortar;
        static constexpr int rank = functionspace_type::rank;
        static constexpr bool is_hcurl_conforming = functionspace_type::is_hcurl_conforming;
        static constexpr bool is_hdiv_conforming = functionspace_type::is_hdiv_conforming;

        /** @name Typedefs
         */
        //@{


        typedef typename ublas::type_traits<value_type>::real_type real_type;
        typedef T element_type;

        typedef Cont super;
        typedef Cont container_type;
        typedef container_type vector_temporary_type;

        using polyset_type = mp11::mp_if_c<is_composite, boost::none_t, typename basis_0_type::polyset_type>;

        using pc_type = mp11::mp_if_c<is_composite, boost::none_t, typename basis_0_type::PreCompute>;
        using pc_ptrtype = std::shared_ptr<pc_type>;

        //typedef typename basis_type::polyset_type return_value_type;
        typedef typename functionspace_type::return_type return_type;

        typedef typename matrix_node<value_type>::type matrix_node_type;

        using polynomial_view_type = mp11::mp_if_c<is_composite, boost::none_t, typename basis_0_type::polynomial_type>;

        using local_interpolant_type = mp11::mp_if_c<is_modal, boost::none_t, local_interpolant_t<basis_0_type>>;

        using local_interpolants_type = mp11::mp_if_c<is_modal, boost::none_t, local_interpolants_t<basis_0_type>>;


        typedef Element<T,Cont> this_type;
        using self_t = this_type;

        template<int i>
        struct sub_element
        {
            //typedef typename mpl::at_c<element_vector_type,i>::type type;
            typedef typename mpl::at_c<element_vector_type,i>::type::second_type ptrtype;
            typedef typename ptrtype::element_type type;
        };
        template<int i> using sub_element_ptrtype = typename mpl::at_c<element_vector_type,i>::type::second_type;
        template<int i> using sub_element_type = typename mpl::at_c<element_vector_type,i>::type::second_type::element_type;

        typedef typename functionspace_type::component_functionspace_type component_functionspace_type;
        typedef typename functionspace_type::component_functionspace_ptrtype component_functionspace_ptrtype;
        //typedef typename component_functionspace_type::template Element<T,typename Cont[>VectorUblas<value_type><]::slice::type> component_type;
        typedef typename component_functionspace_type::template Element<T, Cont> component_type;

        /**
         * geometry typedef
         */
        typedef typename mesh_type::element_type geoelement_type;
        typedef typename functionspace_type::gm_type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;
        typedef typename gm_type::template Context</*vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB,*/ geoelement_type> gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;

        static constexpr size_type DefaultCTX = (is_hdiv_conforming?vm::POINT|vm::JACOBIAN:(is_hcurl_conforming?vm::POINT|vm::KB:vm::POINT));
#if 0
        template <size_type CTX=DefaultCTX>
        using gmc_t = typename gm_type::template Context<CTX, geoelement_type>;
        template <size_type CTX=DefaultCTX>
        using gmc_ptr_t = std::shared_ptr<gmc_t<CTX>>;
#endif

        using fe_type = typename functionspace_type::fe_type;
        template <size_type FECTX=DefaultCTX, size_type GEOCTX=DefaultCTX>
        using fec_t = typename basis_0_type::template Context<FECTX, fe_type, gm_type, geoelement_type,GEOCTX>;
        template <size_type FECTX=DefaultCTX, size_type GEOCTX=DefaultCTX>
        using fec_ptr_t = std::shared_ptr<fec_t<FECTX,GEOCTX>>;

        //@}

        /** @name Constructors, destructor
         */
        //@{

        Element();
        Element( Element&& ) = default;
        Element( Element const& __e );




        Element( functionspace_ptrtype const& __functionspace,
                 std::string const& __name,
                 std::string const& __desc,
                 size_type __start = 0,
                 ComponentType __ct = ComponentType::NO_COMPONENT);

        Element( functionspace_ptrtype const& __functionspace,
                 std::string const& __name = "unknown",
                 size_type __start = 0,
                 ComponentType __ct = ComponentType::NO_COMPONENT );

        Element( functionspace_ptrtype const& __functionspace,
                 container_type const& __c,
                 std::string const& __name,
                 std::string const& __desc,
                 size_type __start = 0,
                 ComponentType __ct = ComponentType::NO_COMPONENT,
                 ComponentType __ct2 = ComponentType::NO_COMPONENT );

        Element( functionspace_ptrtype const& __functionspace,
                 container_type const& __c,
                 std::string const& __name = "unknown",
                 size_type __start = 0,
                 ComponentType __ct = ComponentType::NO_COMPONENT,
                 ComponentType __ct2 = ComponentType::NO_COMPONENT );

        Element( functionspace_ptrtype const& __functionspace,
                 size_type nActiveDof, value_type* arrayActiveDof,
                 size_type nGhostDof, value_type* arrayGhostDof );


        ~Element() override;

        void initFromSpace( functionspace_ptrtype const& __functionspace,
                            container_type const& __c );

        //@}

        /** @name Operator overloads
         */
        //@{
        Element& operator=( Element && __e );
        Element& operator=( Element const& __e );
#if 0
        template<typename ExprT>
        Element& operator=( vf::Expr<ExprT> const& __expr )
        {
            *this = project( this->functionspace(), elements( this->functionspace()->mesh() ), __expr );
            return *this;
        }
#endif

        template<typename ContOtherType>
        Element& operator=( Element<TT,ContOtherType> const& v );

        template<typename VectorExpr>
        Element& operator=( VectorExpr const& v );


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

        value_type globalValue( size_type i ) const
        {
            return this->operator()( i );
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
            //return this->template comp<THECOMP>( typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
            return this->comp( THECOMP );
        }
#if 0
        template<ComponentType THECOMP>
        component_type
        comp( mpl::bool_<true> )
        {
            CHECK( THECOMP >= ComponentType::X && (int)THECOMP < N_COMPONENTS ) << "Invalid component " << (int)THECOMP;
            auto s = ublas::slice( (int)THECOMP, N_COMPONENTS, M_functionspace->nLocalDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( THECOMP );
            return component_type( compSpace(),
                                   typename component_type::container_type( this->vec().data().expression(), s, this->compSpace()->dof() ),
                                   __name,
                                   start()+(int)THECOMP,
                                   THECOMP );
        }
        template<ComponentType THECOMP>
        component_type
        comp( mpl::bool_<false> )
        {
            CHECK( THECOMP >= ComponentType::X && (int)THECOMP < N_COMPONENTS ) << "Invalid component " << (int)THECOMP;
            auto s = ublas::slice( (int)THECOMP, N_COMPONENTS, M_functionspace->nLocalDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( THECOMP );
            return component_type( compSpace(),
                                   typename component_type::container_type( ( VectorUblas<value_type>& )*this, s, this->compSpace()->dof()  ),
                                   __name,
                                   start()+(int)THECOMP,
                                   THECOMP );
        }
#endif
        /**
         * get the component of the element
         * const version
         *
         * @param i component id
         * @return the i-th component of the element
         */
        component_type
        comp( ComponentType i, ComponentType j = ComponentType::NO_COMPONENT ) const
        {
            CHECK( i >= ComponentType::X && (int)i < nComponents1 ) << "Invalid component " << (int)i;
            int startSlice = ((int)i);
            std::string __name = this->name() + "_" + componentToString( i );
            if ( j != ComponentType::NO_COMPONENT )
            {
                CHECK( j >= ComponentType::X && (int)j < nComponents2 ) << "Invalid component " << (int)j;

                if ( is_tensor2symm )
                {
                    startSlice = Feel::detail::symmetricIndex( (int)i, (int)j, nComponents1 );
                }
                else
                {
                    startSlice = ((int)i)*nComponents2+((int)j);
                }
                __name += "_" + componentToString( j );

            }

            auto sActive = ublas::slice( startSlice, nRealComponents, M_functionspace->nLocalDofWithoutGhostPerComponent() );
            auto sGhost = ublas::slice( startSlice, nRealComponents, M_functionspace->nLocalGhostPerComponent() );

            //std::cout << "extract component " << (int)i << " start+i:" << start()+(int)i << " slice size:" << s.size();

            size_type startContainerIndex = start() + startSlice;
            // Warning: drop const-correctness to avoid redefining full "const"-Elements
            // should be fixed...
            component_type c( compSpace(),
                    const_cast<this_type *>(this)->container().slice( sActive, sGhost, this->compSpace()->dof() ),
                    __name,
                    startContainerIndex,
                    i,j );
            return c;
        }

        component_type
        comp( ComponentType i, ComponentType j = ComponentType::NO_COMPONENT )
        {
            CHECK( i >= ComponentType::X && (int)i < nComponents1 ) << "Invalid component " << (int)i;
            int startSlice = ((int)i);
            std::string __name = this->name() + "_" + componentToString( i );
            if ( j != ComponentType::NO_COMPONENT )
            {
                CHECK( j >= ComponentType::X && (int)j < nComponents2 ) << "Invalid component " << (int)j;

                if ( is_tensor2symm )
                {
                    startSlice = Feel::detail::symmetricIndex( (int)i, (int)j, nComponents1 );
                }
                else
                {
                    startSlice = ((int)i)*nComponents2+((int)j);
                }
                __name += "_" + componentToString( j );

            }

            auto sActive = ublas::slice( startSlice, nRealComponents, M_functionspace->nLocalDofWithoutGhostPerComponent() );
            auto sGhost = ublas::slice( startSlice, nRealComponents, M_functionspace->nLocalGhostPerComponent() );

            //std::cout << "extract component " << (int)i << " start+i:" << start()+(int)i << " slice size:" << s.size();

            size_type startContainerIndex = start() + startSlice;
            component_type c( compSpace(),
                    this->container().slice( sActive, sGhost, this->compSpace()->dof() ),
                    __name,
                    startContainerIndex,
                    i,j );
            return c;
        }

        component_type
        operator[]( ComponentType i )
            {
                //return comp( i, mpl::bool_<boost::is_same<>is_composite>() );
                return comp( i, ComponentType::NO_COMPONENT );//, typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
            }
        component_type
        operator[]( ComponentType i ) const
            {
                //return comp( i, mpl::bool_<boost::is_same<>is_composite>() );
                return comp( i, ComponentType::NO_COMPONENT );//, typename mpl::not_<boost::is_same<container_type,VectorUblas<value_type> > >::type() );
            }
        value_type&
        operator[]( index_type i )
            {
                return super::operator[]( i );
            }
        value_type
        operator[]( index_type i )  const
            {
                return super::operator[]( i );
            }

        value_type localToGlobal( index_type e, index_type l, int c ) const
        {
            index_type index=M_functionspace->dof()->localToGlobal( e, l, c ).index();
            return super::operator()( index );
        }
#if 0
        value_type& localToGlobal( index_type e, index_type l, int c )
        {
            index_type index=boost::get<0>( M_functionspace->dof()->localToGlobal( e, l, c ) );
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

        using p0dh_t =  FunctionSpace<mesh_type,bases<Lagrange<0,Scalar,Discontinuous>>>;
        using p0dh_element_t =  typename FunctionSpace<mesh_type,bases<Lagrange<0,Scalar,Discontinuous>>>::template Element<value_type>;

        template<typename elt_t = p0dh_element_t,class = std::enable_if_t<!std::is_same<elt_t,this_type>::value>>
        Element& plusAssign( elt_t const& _e, const value_type sign = 1. )
            {
                if ( this->mesh()  != _e.mesh() )
                    throw std::logic_error("Invalid mesh, they should be the same");
                for ( auto const& rangeElt : elements( this->mesh() ) )
                {
                    auto const& meshElt = boost::unwrap_ref( rangeElt );
                    index_type e = meshElt.id();
                    auto p0_eid = _e.functionSpace()->dof()->localDof( e ).first->second.index();
                    auto _v = _e.operator[](p0_eid);
                    for( auto const& ldof : M_functionspace->dof()->localDof( e ) )
                    {
                        index_type index = ldof.second.index();
                        super::operator[]( index ) += sign*_v;
                    }
                }
                return *this;
            }
        template<typename elt_t = p0dh_element_t,class = std::enable_if_t<!std::is_same<elt_t,this_type>::value>>
        Element& operator+=( elt_t const& _e )
            {
                return plusAssign( _e );
            }
        template<typename elt_t = p0dh_element_t,class = std::enable_if_t<!std::is_same<elt_t,this_type>::value>>
        Element& operator-=( elt_t const& _e )
            {
                return plusAssign( _e, -1. );
            }

        Element& operator-=( Element const& _e )
        {
            for ( index_type i=0; i < _e.nLocalDof(); ++i )
                this->operator()( i ) -= _e( i );

            return *this;
        }
        Element operator-() const { Element r_(*this); r_.scale( -1.); return r_; }
        self_t operator+( value_type f ) const { Element r_(*this); r_.setConstant(f); r_.add( 1., *this); return r_; }
        friend self_t operator+( value_type f_, self_t const& e_ ) { Element r_(e_); r_.setConstant(f_); r_.add( 1., e_); return r_; }
        self_t operator-( value_type f ) const { Element r_(*this); r_.setConstant(-f); r_.add( 1., *this); return r_; }
        friend self_t operator-( value_type f_, self_t const& e_ ) { Element r_(e_); r_.setConstant(f_); r_.add( -1., e_); return r_; }
        self_t operator*( value_type f_ ) const { Element r_(*this); r_.scale(f_); return r_; }
        friend self_t operator*( value_type f_, self_t const& e_ ) { Element r_(e_); r_.scale(f_); return r_; }
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
        void assign( index_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            index_type index = M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            super::operator[]( index ) = __v;
        }
        void plus_assign( index_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            index_type index = M_functionspace->dof()->localToGlobal( ie, il, c ).index();
            super::operator[]( index ) += __v;
        }
        local_interpolant_type element( std::vector<index_type> const& e ) const
            {
                local_interpolant_type l( M_functionspace->basis()->localInterpolant(e.size()) ) ;
                element( e, l );
                return l;
            }
        //!
        //! @return the components of the element associated to the list of elements in e
        //! @note the vector of components is already allocated
        //! @code
        //! Xh->element( K, EigenVector )
        //! @endcode
        //!
        template<typename B = basis_0_type>
        void element( std::vector<index_type> const& e, Eigen::Ref<local_interpolant_t<B>> l, std::enable_if_t<!B::is_modal>* = nullptr ) const
            {
                int s = l.size()/e.size();
                int n = 0;
                std::for_each( e.begin(), e.end(), [&]( auto const& id ){

                        for( auto const& ldof : M_functionspace->dof()->localDof( id ) )
                        {
                            index_type index = ldof.second.index();
                            l( n*s+ldof.first.localDof() ) = super::operator[]( index );
                        }
                        ++n;
                    });
            }
        //!
        //! @return the components of the element associated to the list of elements in e
        //! @note the vector of components is already allocated
        //!
        template<typename B = basis_0_type>
        void element( std::vector<index_type> const& e, local_interpolant_type& l, std::enable_if_t<!B::is_modal>* = nullptr ) const
            {
                int s = l.size()/e.size();
                int n = 0;
                std::for_each( e.begin(), e.end(), [&]( auto const& id ){

                        for( auto const& ldof : M_functionspace->dof()->localDof( id ) )
                        {
                            index_type index = ldof.second.index();
                            l( n*s+ldof.first.localDof() ) = super::operator[]( index );
                        }
                        ++n;
                    });
            }
        std::vector<int> dofs( std::vector<index_type> const& e ) const
            {
                std::vector<int> d;
                int N = M_functionspace->dof()->nRealLocalDof();
                d.reserve( e.size()*N );
                std::for_each( e.begin(), e.end(), [&]( auto const& id ){

                        for( auto const& ldof : M_functionspace->dof()->localDof( id ) )
                        {
                            d.push_back( ldof.second.index() );
                        }
                    });
                return d;

            }
        std::vector<int> dofs( std::vector<index_type> const& e, DataMap<size_type> const& dm, int block  ) const
            {
                std::vector<int> d;
                int N = M_functionspace->dof()->nRealLocalDof();
                d.reserve( e.size()*N );
                std::for_each( e.begin(), e.end(), [&]( auto const& id ){

                        for( auto const& ldof : M_functionspace->dof()->localDof( id ) )
                        {
                            d.push_back( dm.dofIdToContainerId( block, ldof.second.index() ) );
                        }
                    });
                return d;

            }
        template<typename Tloc>
        void assignE( index_type e, Tloc&& loc, bool symm = true )
            {
                int N0 = M_functionspace->dof()->nRealLocalDof();
                int N0c = M_functionspace->dof()->nRealLocalDof( true );
                auto const& s = M_functionspace->dof()->localToGlobalSigns( e );
                    for( auto const& ldof : M_functionspace->dof()->localDof( e ) )
                    {
                        index_type index = ldof.second.index();
                        if ( is_tensor2symm )
                        {
                            int i = M_functionspace->dof()->fe().unsymmToSymm(ldof.first.localDof());
                            super::operator[]( index ) = s(i)*loc(i);
                        }
                        else
                        {
                            super::operator[]( index ) = s(ldof.first.localDof())*loc( ldof.first.localDof() );
                        }
                    }
            }
        void assign( index_type e, local_interpolant_type const& Ihloc )
            {
                auto const& s = M_functionspace->dof()->localToGlobalSigns( e );
                for( auto const& ldof : M_functionspace->dof()->localDof( e ) )
                {
                    index_type index = ldof.second.index();
                    super::operator[]( index ) = s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
                }
            }
        void plus_assign( index_type e, local_interpolant_type const& Ihloc )
            {
                auto const& s = M_functionspace->dof()->localToGlobalSigns( e );
                for( auto const& ldof : M_functionspace->dof()->localDof( e ) )
                {
                    index_type index = ldof.second.index();
                    super::operator[]( index ) += s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
                }
            }

        void assign( geoelement_type const& e, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.id() );
            for( auto const& ldof : M_functionspace->dof()->localDof( e.id() ) )
            {
                index_type index = ldof.second.index();
                super::operator[]( index ) = s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
            }
        }
        template<typename EltType>
        void assign( EltType const& e, local_interpolant_type const& Ihloc,
                     typename std::enable_if<is_3d_real<EltType>::value && is_edge<EltType>::value>::type* = nullptr )
            {
                // we assume here that we are in CG
                // TODO : adapt to DG and loop over all element to which the point belongs to
                // TODO : check if the doftable is computed for this eid
                index_type eid = e.elements().begin()->first;
                uint16_type edgeid_in_element = e.elements().begin()->second;
                this->assign( e, Ihloc, std::make_pair( eid,edgeid_in_element ) );
            }
        template<typename EltType>
        void assign( EltType const& e, local_interpolant_type const& Ihloc, std::pair<index_type, uint16_type> const& eltsInfo,
                     typename std::enable_if<is_3d_real<EltType>::value && is_edge<EltType>::value>::type* = nullptr )
            {
                const auto [ eid,edgeid_in_element ]  = eltsInfo;

                auto const& s = M_functionspace->dof()->localToGlobalSigns( eid );
                for( auto const& ldof : M_functionspace->dof()->edgeLocalDof( eid, edgeid_in_element ) )
                {
                    index_type index= ldof.index();
                    //super::operator[]( index ) = s(edgeid_in_element)*Ihloc( ldof.localDofInFace() );
                    super::operator[]( index ) = Ihloc( ldof.localDofInFace() );
                }
            }

        void assign( geopoint_type const& p, local_interpolant_type const& Ihloc )
            {
                // we assume here that we are in CG
                // TODO : adapt to DG and loop over all element to which the point belongs to
                // TODO : check if the doftable is computed for this eid
                index_type eid = p.elements().begin()->first;
                uint16_type ptid_in_element = p.elements().begin()->second;
                this->assign( p, Ihloc, std::make_pair( eid, ptid_in_element ) );
            }
        void assign( geopoint_type const& p, local_interpolant_type const& Ihloc, std::pair<index_type, uint16_type> const& eltsInfo )
            {
                index_type eid = eltsInfo.first;
                uint16_type ptid_in_element = eltsInfo.second;
                auto const& s = M_functionspace->dof()->localToGlobalSigns( eid );
                for( int c = 0; c < (is_product?nComponents:1); ++c )
                {
                    index_type index = M_functionspace->dof()->localToGlobal( eid, ptid_in_element, c ).index();
                    super::operator[]( index ) = s(ptid_in_element)*Ihloc( c );
                }
            }

        void assign( geoface_type const& e, uint16_type faceConnectionId, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.element( faceConnectionId ).id() );
            for( auto const& ldof : M_functionspace->dof()->faceLocalDof( e.id() ) )
            {
                index_type index = ldof.index();
                super::operator[]( index ) = s(ldof.localDof())*Ihloc( ldof.localDofInFace() );
            }
        }
        void plus_assign( geoelement_type const& e, local_interpolant_type const& Ihloc )
        {
            auto const& s = M_functionspace->dof()->localToGlobalSigns( e.id() );
            for( auto const& ldof : M_functionspace->dof()->localDof( e.id() ) )
            {
                index_type index = ldof.second.index();
                super::operator[]( index ) += s(ldof.first.localDof())*Ihloc( ldof.first.localDof() );
            }
        }
        void plus_assign( geoface_type const& e, uint16_type faceConnectionId, local_interpolant_type const& Ihloc )
         {
             auto const& s = M_functionspace->dof()->localToGlobalSigns( e.element( faceConnectionId ).id() );
            for( auto const& ldof : M_functionspace->dof()->faceLocalDof( e.id() ) )
            {
                index_type index = ldof.index();
                super::operator[]( index ) += s(ldof.localDof())*Ihloc( ldof.localDofInFace() );
            }
        }

        template<typename EltType>
        void plus_assign( EltType const& e, local_interpolant_type const& Ihloc,
                     typename std::enable_if<is_3d_real<EltType>::value && is_edge<EltType>::value>::type* = nullptr )
        {
            // we assume here that we are in CG
            // TODO : adapt to DG and loop over all element to which the point belongs to
            // TODO : check if the doftable is computed for this eid
            index_type eid = e.elements().begin()->first;
            uint16_type edgeid_in_element = e.elements().begin()->second;
            this->plus_assign( e, Ihloc, std::make_pair( eid,edgeid_in_element ) );
        }
        template<typename EltType>
        void plus_assign( EltType const& e, local_interpolant_type const& Ihloc, std::pair<index_type, uint16_type> const& eltsInfo,
                     typename std::enable_if<is_3d_real<EltType>::value && is_edge<EltType>::value>::type* = nullptr )
        {
            index_type eid = eltsInfo.first;
            uint16_type edgeid_in_element = eltsInfo.second;
            auto const& s = M_functionspace->dof()->localToGlobalSigns( eid );
            for( auto const& ldof : M_functionspace->dof()->edgeLocalDof( eid, edgeid_in_element ) )
            {
                index_type index= ldof.index();
                super::operator[]( index ) += s(edgeid_in_element)*Ihloc( ldof.localDofInFace() );
            }
        }

        void plus_assign( geopoint_type const& p, local_interpolant_type const& Ihloc )
            {
                // we assume here that we are in CG
                // TODO : adapt to DG and loop over all element to which the point belongs to
                // TODO : check if the doftable is computed for this eid
                index_type eid = p.elements().begin()->first;
                uint16_type ptid_in_element = p.elements().begin()->second;
                this->plus_assign( p, Ihloc, std::make_pair( eid,ptid_in_element ) );
            }

        void plus_assign( geopoint_type const& p, local_interpolant_type const& Ihloc, std::pair<index_type, uint16_type> const& eltsInfo )
            {
                index_type eid = eltsInfo.first;
                uint16_type ptid_in_element = eltsInfo.second;
                auto const& s = M_functionspace->dof()->localToGlobalSigns( eid );
                for( int c = 0; c < (is_product?nComponents:1); ++c )
                {
                    index_type index = M_functionspace->dof()->localToGlobal( eid, ptid_in_element, c ).index();
                    super::operator[]( index ) += s(ptid_in_element)*Ihloc( c );
                }
            }

        //@}

        /** @name Accessors
         */
        //@{

        typedef boost::multi_array<value_type,3> array_type;
#if 0
        typedef Eigen::Tensor<value_type,nComponents1,nComponents2> _id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nRealDim> _grad_type;
        typedef Eigen::Matrix<value_type,nRealDim,nRealDim> _hess_type;
        typedef Eigen::Matrix<value_type,nComponents2,1> _div_type;
        typedef Eigen::Matrix<value_type,nRealDim,1> _curl_type;
#else
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>> _id_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nRealDim>> _grad_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>> _dn_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nRealDim,nRealDim>> _hess_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents2,1>> _div_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<(nRealDim==3)?nRealDim:1,1>> _curl_type;
        //typedef Eigen::Tensor<value_type,2> _dn_type;
        //typedef Eigen::Tensor<value_type,2> _hess_type;
        //typedef Eigen::Tensor<value_type,2> _div_type;
        //typedef Eigen::Tensor<value_type,2> _curl_type;
#endif
        typedef boost::multi_array<_id_type,1> id_array_type;
        typedef boost::multi_array<_grad_type,1> grad_array_type;
        typedef boost::multi_array<_dn_type,1> dn_array_type;
        typedef boost::multi_array<_hess_type,1> hess_array_type;
        typedef boost::multi_array<_div_type,1> div_array_type;
        typedef boost::multi_array<_curl_type,1> curl_array_type;
        typedef boost::multi_array<_div_type,1> comp_curl_array_type;

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

        //!
        //! \return the dof table
        //!
        dof_ptrtype dof() const
            {
                return M_functionspace->dof();
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
        typename p0_space_type::element_type extremeValue( std::shared_ptr<p0_space_type> const& P0h, std::string const& extreme )
        {
            // check if the mesh coming from P0h and the class elements is the same
            FEELPP_ASSERT( P0h->mesh() == this->mesh() ).error( "mesh is not the same" );
            FEELPP_ASSERT( is_scalar ).error( "only works for scalar fields" );

            typename p0_space_type::element_type p0Element( P0h );

            for ( auto const& rangeElt : P0h->template rangeElements<0>() /*elements( P0h->mesh() )*/ )
            {
                auto const& meshElt = boost::unwrap_ref( rangeElt );
                index_type eid = meshElt.id();

                index_type dofp0 = P0h->dof()->localToGlobal( eid, 0, 0 ).index();
                std::vector<value_type> values ( functionspace_type::fe_type::nLocalDof );

                index_type dofpn = 0;

                for ( uint16_type local_id=0; local_id < functionspace_type::fe_type::nLocalDof; ++local_id )
                {
                    dofpn = this->functionSpace()->dof()->localToGlobal( eid, local_id, 0 ).index();
                    values[local_id] = this->operator()( dofpn );
                }

                if ( extreme == "max" )
                    p0Element.assign( eid, 0, 0, *( std::max_element( values.begin(), values.end() ) ) );

                else
                    p0Element.assign( eid, 0, 0, *( std::min_element( values.begin(), values.end() ) ) );
            }

            return p0Element;
        }

        value_type max() const override
        {
            return this->max( true );
        }
        value_type max( bool parallel ) const
        {
            return super::max( parallel );
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type max( std::shared_ptr<p0_space_type> const& P0h )
        {
            return this->extremeValue( P0h, "max" );
        }

        value_type min() const override
        {
            return this->min( true );
        }
        value_type min( bool parallel ) const
        {
            return super::min( parallel );
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type min( std::shared_ptr<p0_space_type> const& P0h )
        {
            return this->extremeValue( P0h, "min" );
        }

        template <typename ... CTX>
        basis_context_ptrtype
        selectContext( CTX const& ... ctx ) const
            {
                basis_context_ptrtype res;
                auto allCtx = hana::make_tuple( ctx... );
                hana::for_each( allCtx, [this,&res]( auto const& e )
                                {
                                    if constexpr ( std::is_same_v<std::decay_t<decltype(e)>,basis_context_ptrtype> )
                                    {
                                        // Maybe TODO : differentiate context between several spaces/meshes
                                        res = e;
                                    }
                                } );
                return res;
            }

        //! Interpolation at a set of points
        //@{
        /**
         * data structure that stores the interpolated values of the
         * element at a set of points
         */
        using id_type =  Feel::detail::ID<value_type,nComponents1,nComponents2>;
        using laplacian_type = id_type;

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
            _id_type idzero;
            std::fill( v.data(), v.data()+v.num_elements(), idzero.constant(0.));

            if( context.size() > 0 )
            {
                for( int i = 0 ; it != en; ++it, ++i )
                {
                    v[0].setZero();
                    auto basis = std::get<0>( it->second );
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
            const int ncdof  = is_product?nComponents:1;
            _id_type idzero;
            std::fill( v.data(), v.data()+v.num_elements(), idzero.constant(0.));

            value_type result=0;
            int proc_having_the_point = context.processorHavingPoint( i );
            if( proc_number == proc_having_the_point )
            {
                auto basis = context.at( i );
                id( *std::get<0>( basis ) , v );
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
         * from a set of point in the reference element, should probably done in
         * the real element (geond)
         */
        template<typename Context_t>
        matrix_node_type
        ptsInContext( Context_t const & context, mpl::int_<1> ) const
        {
            //new context for evaluate the points
            // typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
            // typedef std::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

            // gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),  context.pc() ) );
            auto __c_interp = context.geometricMapping()->template context<vm::POINT>( context.element_c(),  context.pc() );

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
#if 0
            //new context for the interpolation
            typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type, Context_t::subEntityCoDim> gmc_interp_type;
            typedef std::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

            // typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::permutation_type permutation_type;
            // typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::precompute_ptrtype precompute_ptrtype;

            //not good because ?
            //gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),context.pcFaces(), context.faceId()) );
            //good with this
            //std::vector<std::map<permutation_type, precompute_ptrtype> > __geo_pcfaces = context.pcFaces();
            auto __geo_pcfaces = context.pcFaces();
            gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId() ) );
#endif
            auto __geo_pcfaces = context.pcFaces();
            auto __c_interp = context.geometricMapping()->template context<vm::POINT>( context.element_c(), __geo_pcfaces , context.faceId() );

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
            //return gmc_ptrtype( new gmc_type( __gm, elt, __geopc ) );
            return __gm->template context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB>( elt, __geopc );
        }
        /**
         * \return a precomputation of the basis functions
         * \warning this seems quite buggy (where is it used?)
         */
        pc_ptrtype pcPtr( geoelement_type const& elt ) const
        {
            return pc_ptrtype( new pc_type( functionSpace()->fe(), elt.G() ) );
        }

        id_type operator()( Eigen::Matrix<value_type,nDim,1> const& __x, bool extrapolate = false, bool parallel = true ) const
            {
                node_type n( nDim );
                for(int i = 0; i < nDim; ++i ) n[i]=__x[i];
                return operator()( n, extrapolate, parallel );
            }

        /**
         * interpolate the function at node (real coordinate) x
         *
         * @return the interpolated value of the function at the real point x
         */
        id_type operator()( node_type const& __x, bool extrapolate = false, bool parallel = true ) const;

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
        grad( ElementType const& elt, grad_array_type& v ) const
        {
            if constexpr ( std::is_base_of_v<ElementType,ConvexBase> )
            {
                gmc_ptrtype gmc( geomapPtr( elt ) );
                v.resize( gradExtents(  *gmc ) );
                grad_( *gmc, *pcPtr( elt ), v );
            }
            else
            {
                grad_( elt, v );
            }
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

        //!
        //! compute Symmetric Gradient only in the vectorial case
        //!
        template<typename ContextType,typename EType = this_type>
        void symmetricGradient( ContextType const & context, grad_array_type& v, std::enable_if_t<EType::is_vectorial>* = nullptr ) const;

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


        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dnExtents( ContextType const & context ) const
        {
            boost::array<typename array_type::index, 1> shape;
            shape[0] = context.xRefs().size2();
            return shape;
        }
        template<typename ContextType>
        void
        dn( ContextType const & context, dn_array_type& v ) const
        {
            dn_( context, v );
        }
        template<typename ContextType>
        void dn_( ContextType const & context, dn_array_type& v ) const;

        void
        dnInterpolate( matrix_node_type __ptsReal, dn_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;


        //@}

        typedef Feel::detail::Div<value_type,nComponents2> div_type;

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

        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        laplacianExtents( ContextType const & context ) const
            {
                BOOST_STATIC_ASSERT( rank <= 1 );

                boost::array<typename array_type::index, 1> shape;
                shape[0] = context.xRefs().size2();
                //shape[1] = nRealDim;
                //shape[2] = nRealDim;
                return shape;
            }

        template<typename ContextType>
        void
        laplacian_( ContextType const & context, id_array_type& v ) const;

        template<typename ContextType>
        void
        laplacian_( ContextType const & context, id_array_type& v, mpl::int_<0> ) const;
        template<typename ContextType>
        void
        laplacian_( ContextType const & context, id_array_type& v, mpl::int_<1> ) const;

        void
        laplacianInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        template<typename ContextType>
        laplacian_type
        laplacian( ContextType const & context ) const
        {
            return laplacian_type( *this, context );
        }

        template<typename ContextType>
        void
        laplacian( ContextType const & context, id_array_type& v ) const
        {
            laplacian_( context, v );
        }
        template<typename ElementType>
        void
        laplacian( ElementType const& elt, id_array_type& v )
        {
            gmc_ptrtype gmc( geomapPtr( elt ) );
            v.resize( laplacianExtents(  *gmc ) );
            laplacian_( *gmc, *pcPtr( elt ), v );
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
                                 std::shared_ptr<BackendType> backend )
            {
                auto r =  functionSpace()->dof()->markerToDof( m );
                index_type s = std::distance( r.first, r.second );
                vector_ptrtype res = backend->newVector(s, s);
                index_type i = 0;
                for( auto it = r.first, en = r.second; it != en; ++it, ++i )
                    res->set( i, super::operator[]( it->second ) );
                return res;
            }

        template<typename BackendType>
        typename BackendType::vector_ptrtype
        extractValuesWithoutMarker( std::string const& m,
                                    std::shared_ptr<BackendType> backend )
            {
                auto r1 =  functionSpace()->dof()->markerToDofLessThan( m );
                auto r2 =  functionSpace()->dof()->markerToDofGreaterThan( m );
                index_type s = std::distance( r1.first, r1.second ) + std::distance( r2.first, r2.second );
                vector_ptrtype res = backend->newVector(s, s);
                index_type i = 0;
                for( auto it = r1.first, en = r1.second; it != en; ++it, ++i )
                    res->set( i, super::operator[]( it->second ) );
                for( auto it = r2.first, en = r2.second; it != en; ++it, ++i )
                    res->set( i, super::operator[]( it->second ) );
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
#if 0
        template<int i>
        typename sub_element<i>::type
        elementImpl( std::string const& name ="u", bool updateOffViews=true ) const;
#endif
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
        worldscomm_ptr_t const& worldsComm() const
        {
            return M_functionspace->worldsComm();
        }

        /**
         * world communicator
         */
        WorldComm & worldComm()
        {
            return M_functionspace->worldComm();
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

        /**
         \return the description of the element
         */
        std::string const& description() const
        {
            return M_desc;
        }

        size_type start() const
        {
            return M_start;
        }

        bool isAComponent() const
        {
            return M_ct >= ComponentType::X && M_ct <= ComponentType::Z;
        }

        ComponentType component() const
        {
            if (  M_ct < ComponentType::X || M_ct > ComponentType::Z )
            {
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << (int)M_ct;
                throw std::logic_error( __str.str() );
            }

            return M_ct;
        }

        std::string componentToString( ComponentType ct ) const
        {
            if (  ct < ComponentType::X || ct > ComponentType::Z )
            {
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << (int)ct;
                throw std::logic_error( __str.str() );
            }

            switch ( ct )
            {
            case ComponentType::X:
                return "X";

            case ComponentType::Y:
                return "Y";

            case ComponentType::Z:
                return "Z";

            default:
                std::ostringstream __str;
                __str << "invalid component extraction (should be 0 or 1 or 2): " << (int)ct;
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

        /**
         * set the description of the field (e.g. a formula)
         * @param __desc description of the field
         */
        void seDescription( std::string const & __desc )
            {
                M_desc = __desc;
            }

        void setFunctionSpace( functionspace_ptrtype space )
        {
            M_functionspace = space;
            super::init( M_functionspace->dof() );
            //super::init( M_functionspace->nDof(),  M_functionspace->nLocalDof() );
        }

        template <typename ... Ts>
        void save( Ts && ... v ) const
            {
                auto args = NA::make_arguments( std::forward<Ts>(v)... );
                auto && path = args.get(_path);
                std::string const& name = args.get_else(_name,M_name);
                std::string const& type = args.get_else(_type,"default");
                std::string const& suffix = args.get_else(_suffix,"");
                std::string const& sep = args.get_else(_sep,"");
                saveImpl( Environment::expand( path ), name, type, suffix, sep );
            }

        //!
        //! save function space element in file
        //! @param path path to files
        //! @param type file type binary, ascii, hdf5, xml
        //! @param suffix filename suffix to use
        //! @param sep separator to use in filenames
        //!
        void saveImpl( std::string const& path, std::string const& name, std::string const& type = "binary", std::string const& suffix = "", std::string const & sep = "") const
        {
            std::string typeUsed = type;
            if ( typeUsed == "default" )
            {
#ifdef FEELPP_HAS_HDF5
                typeUsed = "hdf5";
#else
                typeUsed = "binary";
#endif
            }

            if ( typeUsed != "binary" && typeUsed != "text" && typeUsed != "xml" &&  typeUsed != "hdf5" )
            {
                LOG(WARNING)  << "[save] : invalid format " << typeUsed << " (type available : binary,text,xml,hdf5)";
                return;
            }

            // if directory does not exist, create it only by one process
            if ( this->worldComm().isMasterRank() && !fs::exists( fs::path( path ) ) )
            {
                fs::create_directories( fs::path( path ) );
            }
            // wait creating directory
            this->worldComm().barrier();

            if ( typeUsed == "binary" || typeUsed == "text" || typeUsed == "xml" )
            {
                std::ostringstream os1;
                os1 << name << sep << suffix << "-" << this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
                fs::path p = fs::path( path ) / os1.str();
                LOG(INFO) << "saving "  << p << "\n";

                if ( typeUsed == "binary" )
                {
                    std::ofstream ofs( p );
                    boost::archive::binary_oarchive oa( ofs );
                    oa << *this;
                }

                else if ( typeUsed == "text" )
                {
                    std::ofstream ofs( p );
                    boost::archive::text_oarchive oa( ofs );
                    oa << *this;
                }

                else if ( typeUsed == "xml" )
                {
                    //boost::archive::xml_oarchive oa(ofs);
                    //oa << *this;
                }
            }
            else if ( typeUsed == "hdf5" )
            {
                std::ostringstream os2;
                os2 << name << sep << suffix << ".h5";
                fs::path filename = fs::path( path ) / os2.str();
#ifdef FEELPP_HAS_HDF5
                this->saveHDF5( filename.string() );
#else
                CHECK( false ) << "Feel++ is not compiled with hdf5";
#endif
            }
        }

        template <typename ... Ts>
        bool load( Ts && ... v )
            {
                auto args = NA::make_arguments( std::forward<Ts>(v)... );
                auto && path = args.get(_path);
                std::string const& name = args.get_else(_name,M_name);
                std::string const& type = args.get_else(_type,"default");
                std::string const& suffix = args.get_else(_suffix,"");
                std::string const& sep = args.get_else(_sep,"");
                std::string const& space_path = args.get_else(_space_path,"");
                return loadImpl( Environment::expand( path ), name, type, suffix, sep, space_path );
            }

        //!
        //! load function space element from file
        //! @param path path to file
        //! @param type file type binary, ascii, hdf5, xml
        //! @param suffix filename suffix to use
        //! @param sep separator to use in filename
        //! @param space_path path to space file related to input file path
        //!
        bool loadImpl( std::string const& path, std::string const& name, std::string const& type = "binary", std::string const& suffix = "", std::string const& sep = "", std::string const& space_path = "" )
        {
            std::ostringstream oss;
            fs::path p;

            if ( !fs::exists( path ) )
            {
                LOG(WARNING) << "[load] : directory or file " << p << " not exists";
                return false;
            }

            std::string typeUsed = type;
            if ( typeUsed == "default" )
            {
#ifdef FEELPP_HAS_HDF5
                typeUsed = "hdf5";
#else
                typeUsed = "binary";
#endif
            }

            if ( fs::is_directory(path) )
            {
                if ( typeUsed == "hdf5" )
                    oss << name << sep << suffix << ".h5";
                else
                    oss << name << sep << suffix << "-" <<  this->worldComm().globalSize() << "." << this->worldComm().globalRank() << ".fdb";
                p = fs::path(path)/oss.str();

                if ( !fs::exists( p ) )
                {
                    LOG(WARNING) << "[load] : file " << p << " not exists";
                    return 0;
                }
                if ( !fs::is_regular_file( p ) )
                {
                    LOG(WARNING) << "[load] : file " << p << " is not a  regular_file";
                    return 0;
                }
            }
            else if ( fs::is_regular_file( path ) )
            {
                p = path;
            }
            else
            {
                LOG(WARNING) << "[load] : file " << path << " is not a directory or a regular_file";
                return 0;
            }

            std::optional<std::vector<index_type>> spacesRelation;
            if ( !space_path.empty() )
                spacesRelation = this->functionSpace()->relationFromFile( space_path );

            if ( typeUsed == "binary" || typeUsed == "text" || typeUsed == "xml" )
            {
                std::ifstream ifs( p );
                if ( typeUsed == "binary" )
                {
                    boost::archive::binary_iarchive ia( ifs );
                    ia >> *this;
                }

                else if ( typeUsed == "text" )
                {
                    boost::archive::text_iarchive ia( ifs );
                    ia >> *this;
                }
                else if ( typeUsed == "xml" )
                {
                    //boost::archive::xml_iarchive ia(ifs);
                    //ia >> *this;
                    return false;
                }
            }
            else if ( typeUsed == "hdf5" )
            {
#ifdef FEELPP_HAS_HDF5
                this->loadHDF5( p.string(), spacesRelation );
#else
                CHECK( false ) << "Feel++ is not compiled with hdf5";
#endif
            }
            else
            {
                LOG(WARNING)  << "[load] : invalid format " << typeUsed << " (type available : binary,text,xml,hdf5)";
                return false;
            }

            return true;

        }

        void printMatlab( std::string fname, bool gmsh = false ) const override
            {
                container_type m( *this );
                if ( gmsh )
                {
                    auto relation = this->functionSpace()->dof()->pointIdToDofRelation();
                    for( size_type i = 0; i < this->localSize(); ++i )
                    {
                        m[relation.first[i]-1] = super::operator[]( i );
                    }
                }
                m.printMatlab( fname );
            }

        //!
        //! compute the element wise mean value of an element
        //!
        p0dh_element_t ewiseMean( std::shared_ptr<p0dh_t> & P0dh ) const
            {
                p0dh_element_t v = P0dh->element();
                v.zero();
                typename basis_type::points_type p(mesh_type::nDim,1);
                auto basispc = this->functionSpace()->basis()->preCompute( this->functionSpace()->basis(), p );
                auto gmpc = this->mesh()->gm()->preCompute( this->mesh()->gm(), p );
                auto range_elts = elements( this->mesh() );
                auto elt_beg = begin( range_elts );
                auto elt_end = end( range_elts );
                if ( elt_beg != elt_end )
                {
                    auto gmc = this->mesh()->gm()->template context<vm::JACOBIAN>( boost::unwrap_ref(*elt_beg), gmpc );
                    for ( auto const& rangeElt : elements( this->mesh() ) )
                    {
                        auto const& meshElt = boost::unwrap_ref( rangeElt );
                        size_type e = meshElt.id();
                        gmc->template update<vm::JACOBIAN>( meshElt );
                        auto p0_eid = v.functionSpace()->dof()->localDof( e ).first->second.index();
                        for( auto const& ldof : M_functionspace->dof()->localDof( e ) )
                        {
                            size_type index = ldof.second.index();
                            v.operator[](p0_eid) += basispc->firstMoment(ldof.first.localDof())*super::operator[]( index );
                        }
                        v.operator[](p0_eid) *= gmc->J(0)/meshElt.measure();
                    }
                }
                return v;
            }

        template <typename ... Ts>
        void on( Ts && ... v )
            {
                auto args = NA::make_arguments( std::forward<Ts>(v)... );
                auto && expr = args.get(_expr);
                auto && range = args.get_else_invocable(_range, [this]() { return elements(this->functionSpace()->template meshSupport<0>()); } );
                std::string const& prefix = args.get_else(_prefix,"");
                GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT);
                bool accumulate = args.get_else(_accumulate,false);
                bool close = args.get_else(_close,false);
                bool verbose = args.get_else_invocable(_verbose,[&prefix](){ return boption(_prefix=prefix,_name="on.verbose"); } );

                onImpl( range, expr, prefix, Feel::detail::geomapStrategy(range,geomap), accumulate, verbose );
                if ( close )
                {
                    std::string opUsed = ( accumulate )? "+" : "=";
                    sync( *this, opUsed, this->functionSpace()->dofs( range, ComponentType::NO_COMPONENT, true ) );
                }
            }

        template <typename ... Ts>
        void plus( Ts && ... v )
            {
                auto args = NA::make_arguments( std::forward<Ts>(v)... );
                auto && expr = args.get(_expr);
                auto && range = args.get_else_invocable(_range, [this]() { return elements(this->functionSpace()->template meshSupport<0>()); } );
                std::string const& prefix = args.get_else(_prefix,"");
                GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT);
                bool close = args.get_else(_close,false);
                bool verbose = args.get_else_invocable(_verbose,[&prefix](){ return boption(_prefix=prefix,_name="on.verbose"); } );

                this->on(_range=range,_expr=expr,_prefix=prefix,_geomap=geomap,_accumulate=true,_close=close, _verbose=verbose );
            }


        //@}
    private:

        /*
         *
         */
        FEELPP_NO_EXPORT void initSubElementView( mpl::true_ );

        FEELPP_NO_EXPORT void initSubElementView( mpl::false_ ) {}


        friend class boost::serialization::access;

        template<class Archive>
        void serialize( Archive & ar, const unsigned int version )
        {
            //ar & BOOST_SERIALIZATION_NVP( boost::serialization::base_object<super>(*this) );
            ar & boost::serialization::make_nvp( "name", M_name );
            DVLOG(2) << "got name " << M_name << "\n";
            //ar & boost::serialization::make_nvp( "description", M_desc );
            //DVLOG(2) << "got description " << M_desc << "\n";

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

                auto it = this->begin();
                auto en = this->end();

                for ( size_type i = 0; it != en; ++it, ++i )
                {
                    T value = super::operator[]( i );
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
                    super::operator[]( i ) = value;
                }
            }


        }
    public:
        template<typename RangeType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( RangeType const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate = true, bool verbose = false )
        {
            onImplBase( r, e, prefix, geomap_strategy, accumulate, verbose, boost::is_std_list<RangeType>()  );
        }
    private:
        template<typename RangeType, typename ExprType>
        FEELPP_NO_EXPORT void onImplBase( RangeType const& rList, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::true_ )
        {
            for ( auto const& r : rList )
                onImpl( std::make_pair( r.template get<1>(), r.template get<2>()), e, prefix, geomap_strategy, accumulate, verbose, mpl::int_<RangeType::iDim()>() );
        }
        template<typename RangeType, typename ExprType>
        FEELPP_NO_EXPORT void onImplBase( RangeType const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::false_ )
        {
            onImpl( std::make_pair( r.begin(), r.end()), e, prefix, geomap_strategy, accumulate, verbose, mpl::int_<RangeType::iDim()>() );
        }

        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType,IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_ELEMENTS>  );

        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_FACES>  );
        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_FACES>, mpl::true_  );
        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_FACES>, mpl::false_  );

        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_EDGES>  );

        template<typename IteratorType, typename ExprType>
        FEELPP_NO_EXPORT void onImpl( std::pair<IteratorType, IteratorType> const& r, ExprType const& e, std::string const& prefix, GeomapStrategyType geomap_strategy, bool accumulate, bool verbose, mpl::int_<MESH_POINTS>  );

    private:

        /**
           Finite Element space
        */
        functionspace_ptrtype M_functionspace;

        std::string M_name;
        std::string M_desc;

        size_type M_start;

        //! type of the component
        ComponentType M_ct, M_ct2;

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
    typedef std::shared_ptr<element_type> element_ptrtype;
    typedef Element<value_type> real_element_type;

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
                   mesh_support_vector_type const& meshSupport = mesh_support_vector_type(),
                   size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
                   periodicity_type  periodicity = periodicity_type(),
                   worldscomm_ptr_t const& _worldsComm = Environment::worldsComm(nSpaces),
                   std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false),
                   const std::string& name = "" )
        :
        super( name, _worldsComm[0]->clone() ),
        M_worldsComm( _worldsComm ),
        M_extendedDofTableComposite( extendedDofTable ),
        M_extendedDofTable( extendedDofTable[0] )
    {
        this->init( mesh, meshSupport, mesh_components, periodicity );
    }

    // Constructor
    FunctionSpace( mesh_ptrtype const& mesh,
                   mesh_support_vector_type const& meshSupport,
                   std::vector<Dof<typename mesh_type::size_type> > const& dofindices,
                   periodicity_type periodicity = periodicity_type(),
                   worldscomm_ptr_t const& _worldsComm = Environment::worldsComm(nSpaces),
                   std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false),
                   const std::string& name = "" )
        :
        super( name, _worldsComm[0]->clone() ),
        M_worldsComm( _worldsComm ),
        M_extendedDofTableComposite( extendedDofTable ),
        M_extendedDofTable( extendedDofTable[0] )
    {
        this->init( mesh, meshSupport, 0, dofindices, periodicity );
    }

    explicit FunctionSpace( worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr(),
                            const std::string& name = "" )
        :
        super( name, worldcomm ),
        M_worldsComm( makeWorldsComm( nSpaces, worldcomm ) ),
        M_extendedDofTableComposite( std::vector<bool>(nSpaces,false) ),
        M_extendedDofTable( false )
    {}

    /**
     * helper static function to create a std::shared_ptr<> out of
     * the \c FunctionSpace
     */
    template <typename ... Ts,typename  = typename std::enable_if_t< sizeof...(Ts) != 0 && ( NA::is_named_argument_v<Ts> && ...) > >
    static pointer_type New( Ts && ... v )
    {
        auto args = NA::make_arguments( std::forward<Ts>(v)... );
        auto && mesh = args.get(_mesh);
        worldscomm_ptr_t worldscomm = args.get_else_invocable(_worldscomm,[&mesh](){ return Feel::detail::createWorldsComm<functionspace_type>(mesh).worldsComm(); } );
        size_type components = args.get_else(_components, MESH_RENUMBER | MESH_CHECK);
        auto && periodicity = args.get_else(_periodicity,periodicity_type());
        auto && extended_doftable = args.get_else(_extended_doftable,std::vector<bool>(nSpaces,false) );
        auto && range = args.get_else(_range,mesh_support_vector_type());

        auto cms = Feel::detail::createMeshSupport<functionspace_type>( mesh, range );
        std::vector<bool> edt = Feel::detail::createInfoExtendedDofTable<functionspace_type>( extended_doftable );
        return NewImpl( mesh, cms.M_meshSupportVector, worldscomm, components, periodicity, edt );
    }

    static pointer_type New( mesh_ptrtype const& m ) { return New(_mesh=m); }

    static pointer_type NewImpl( mesh_ptrtype const& __m,
                                 mesh_support_vector_type const& meshSupport,
                                 worldscomm_ptr_t const& worldscomm = Environment::worldsComm(nSpaces),
                                 size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
                                 periodicity_type periodicity = periodicity_type(),
                                 std::vector<bool> extendedDofTable = std::vector<bool>(nSpaces,false) )
    {

        return pointer_type( new functionspace_type( __m, meshSupport, mesh_components, periodicity, worldscomm, extendedDofTable ) );
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

#if 1
    /**
     * copy shared_ptr pointers of a FunctionSpace
     * but not totally a view : must be use we caution
     */
    template<typename SimilarSpaceType>
    void shallowCopy( std::shared_ptr<SimilarSpaceType> const& simSpace )
        {
            BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename SimilarSpaceType::mesh_type, mesh_type>::value ),
                                  INVALID_MESH_TYPE_COMPATIBILITY, (typename SimilarSpaceType::mesh_type, mesh_type ) );
            BOOST_MPL_ASSERT_MSG( ( boost::is_same<typename SimilarSpaceType::basis_type, basis_type>::value ),
                                  INVALID_BASIS_TYPE_COMPATIBILITY, (typename SimilarSpaceType::basis_type, basis_type ) );

            M_worldsComm = simSpace->worldsComm();
            this->setWorldComm( simSpace->worldCommPtr() );
            M_mesh = simSpace->mesh();
            M_periodicity = simSpace->periodicity();
            M_ref_fe=simSpace->fe();
            if ( simSpace->hasCompSpace() )
                M_comp_space = simSpace->compSpace();
            M_dof = simSpace->dof();
            M_dofOnOff = simSpace->dofOnOff();
            M_extendedDofTableComposite = simSpace->extendedDofTableComposite();
            M_extendedDofTable = simSpace->extendedDofTable();
            M_functionspaces = simSpace->functionSpaces();
            if ( simSpace->hasRegionTree() )
                M_rt = simSpace->regionTree();
        }
#endif


    /**
     * initialize the function space
     * \param mesh mesh data
     * \param meshSupport
     * \param mesh_components
     */
    void init( mesh_ptrtype const& mesh,
               mesh_support_vector_type const& meshSupport,
               size_type mesh_components,
               std::vector<Dof<typename mesh_type::size_type> > const& dofindices,
               periodicity_type periodicity = periodicity_type() );

    void init( mesh_ptrtype const& mesh,
               mesh_support_vector_type const& meshSupport,
               size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
               periodicity_type periodicity = periodicity_type() )
    {
        this->init( mesh, meshSupport, mesh_components, std::vector<Dof<typename mesh_type::size_type> >(), periodicity );
    }




    //! destructor: do nothing thanks to shared_ptr<>
    ~FunctionSpace() override
        {
            VLOG(1) << "FunctionSpace Destructor...";
            M_dof.reset();
            CHECK( M_dof.use_count() == 0 ) << "Invalid Dof Table shared_ptr";
            M_dofOnOff.reset();
            CHECK( M_dofOnOff.use_count() == 0 ) << "Invalid Dof OnOff shared_ptr";
            M_ref_fe.reset();
            CHECK( M_ref_fe.use_count() == 0 ) << "Invalid reffe shared_ptr";
        }


    void setWorldsComm( worldscomm_ptr_t const& _worldsComm )
    {
        M_worldsComm=_worldsComm;
    }
    worldscomm_ptr_t & worldsComm()
    {
        return M_worldsComm;
    }
    worldscomm_ptr_t const& worldsComm() const
    {
        return M_worldsComm;
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
        //typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
        //typedef std::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

        typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::permutation_type permutation_type;
        typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::precompute_ptrtype precompute_ptrtype;

        //not good because ?
        //gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),context.pcFaces(), context.faceId()) );
        //good with this
        std::vector<std::map<permutation_type, precompute_ptrtype> > __geo_pcfaces = context.pcFaces();
        //gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId() ) );
        auto __c_interp = context.geometricMapping()->template context<Context_t::context|vm::POINT>( context.element_c(), __geo_pcfaces , context.faceId() );

        return __c_interp->xReal();
    }

    periodicity_type const& periodicity() const { return M_periodicity; }

    /**
     * return 1 if scalar space or the number of components if vectorial space
     */
    uint16_type qDim() const
    {
        return is_tensor2symm?nRealComponents:N_COMPONENTS;
    }

    /**
     * \return the number of degrees of freedom for each space overall the entire domain
     */
    size_type nDof() const
    {
        if constexpr ( is_composite )
        {
            DVLOG(2) << "calling nDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NbDof() );
            DVLOG(2) << "calling nDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            return M_dof->nDof();
        }
    }

    /**
     * \return the number of degrees of freedom for each space on the current subdomain
     */
    size_type nLocalDof() const
    {
        if constexpr( is_composite )
        {
            DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm() ) );
            DVLOG(2) << "calling nLocalDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            //return M_dof->nLocalDof();
            return M_dof->nLocalDofWithGhost();
        }
    }

    size_type nLocalDofWithGhost() const
    {
        if constexpr ( is_composite )
        {
            DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<true> >( this->worldsComm() ) );
            DVLOG(2) << "calling nLocalDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            return M_dof->nLocalDofWithGhost();
        }
    }

    size_type nLocalDofWithoutGhost() const
    {
        if constexpr ( is_composite )
        {
            DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDof<mpl::bool_<false> >( this->worldsComm() ) );
            DVLOG(2) << "calling nLocalDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            return M_dof->nLocalDofWithoutGhost();
        }
    }

    size_type nLocalDofWithGhostOnProc( const int proc ) const
    {
        if constexpr ( is_composite )
        {
            DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<true> >( proc, this->worldsComm() ) );
            DVLOG(2) << "calling nLocalDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            return M_dof->nLocalDofWithGhost(proc);
        }
     }

    size_type nLocalDofWithoutGhostOnProc(const int proc) const
    {
        if constexpr ( is_composite )
        {
            DVLOG(2) << "calling nLocalDof(<composite>) begin\n";
            size_type ndof =  fusion::accumulate( M_functionspaces, size_type( 0 ), Feel::detail::NLocalDofOnProc<mpl::bool_<false> >( proc, this->worldsComm() ) );
            DVLOG(2) << "calling nLocalDof(<composite>) end\n";
            return ndof;
        }
        else
        {
            return M_dof->nLocalDofWithoutGhost(proc);
        }
    }

    /**
     * \return the distribution of the dofs among the processors
     */
    FEELPP_DEPRECATED
    proc_dist_map_type getProcDistMap() const
    {
        return procDistMap;
    }

    /**
     * \return the starting value of the global dof numbering
     */
    size_type nDofStart( size_type i = /*invalid_v<size_type>*/0 ) const
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
     * \return the number of degrees of freedom per dim
     */
    size_type nLocalDofPerComponent() const
    {
        return this->nLocalDof()/qDim();
    }

    /**
     * \return the number of degrees of freedom per dim
     */
    size_type nLocalDofWithoutGhostPerComponent() const
    {
        return this->nLocalDofWithoutGhost()/qDim();
    }

    /**
     * \return the number of degrees of freedom per dim
     */
    size_type nLocalGhostPerComponent() const
    {
        return (this->nLocalDofWithGhost()-this->nLocalDofWithoutGhost())/qDim();
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
        if constexpr ( is_shared_ptr<mesh_ptrtype>() )
            return M_mesh;
        else
            return fusion::at_c<i>(M_mesh);
    }

    /**
     * \return the i-th mesh support if it exists, if not create one
     */
    template<int i>
    typename GetMeshSupport<mesh_ptrtype,i>::ptrtype
    meshSupport() const
    {
        auto meshSupportVector = Feel::detail::FunctionSpaceMeshSupport<functionspace_type>( *this ).M_meshSupportVector;
        auto & meshSupport = boost::fusion::at_c<i>( meshSupportVector );
        if( !meshSupport )
            meshSupport = std::make_shared<typename GetMeshSupport<mesh_ptrtype,i>::type>(mesh<i>());
        return meshSupport;
    }

    /**
     * \return the range of elements on which the i-th space is defined
     */
    template<int i>
    typename GetMeshSupport<mesh_ptrtype,i>::range_type
    rangeElements() const
    {
        return meshSupport<i>()->rangeElements();
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
        return { basis_type::nOrder };
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
    DataMap<size_type> const& map() const
    {
        return *M_dof;
    }
    /**
     \return the degrees of freedom
     */
    datamap_ptrtype mapPtr() const override
        {
            return M_dof;
        }

    /**
       \return the degrees of freedom
    */
    DataMap<size_type> const& mapOn() const
    {
        return *M_dof;
    }

    /**
       \return the degrees of freedom
    */
    DataMap<size_type> const& mapOnOff() const
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
    std::vector<bool> const& extendedDofTableComposite() const { return M_extendedDofTableComposite; }
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


    std::shared_ptr<IndexSplit>
    buildDofIndexSplit()
    {
        auto startSplit = boost::fusion::fold( functionSpaces(), boost::make_tuple(0,0), Feel::detail::computeStartOfFieldSplit() ).template get<1>();
        auto computeSplit = Feel::detail::computeNDofForEachSpace<false>(startSplit);
        boost::fusion::for_each( functionSpaces(), computeSplit );
        return computeSplit.indexSplit();
    }
    std::shared_ptr<IndexSplit>
    buildDofIndexSplitWithComponents()
    {
        bool hasCompSplit = boost::fusion::fold( functionSpaces(), false, Feel::detail::hasSubSpaceWithComponentsSplit() );
        if ( hasCompSplit )
        {
            auto startSplit = boost::fusion::fold( functionSpaces(), boost::make_tuple(0,0), Feel::detail::computeStartOfFieldSplit() ).template get<1>();
            auto computeSplit = Feel::detail::computeNDofForEachSpace<true>(startSplit);
            boost::fusion::for_each( functionSpaces(), computeSplit );
            return computeSplit.indexSplit();
        }
        else
            return std::shared_ptr<IndexSplit>();
    }

    std::shared_ptr<IndexSplit> const&
    dofIndexSplit() const
    {
        return this->dof()->indexSplit();
    }


    /**
     * \return the element 0 of the function space
     */
    element_type
    element( std::string const& name = "u", std::string const& desc = "u" )
    {
        element_type u( this->shared_from_this(), name, desc );
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
    element( ExprT e, std::string const& name = "u", std::string const& desc = "u",
             typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 )
    {
        element_type u( this->shared_from_this(), name, desc );
        bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
        EntityProcessType entityProcess = (addExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        u.on( _range=elements(M_mesh,entityProcess), _expr=e );
        return u;
    }

    element_type
    elementFromExpr( std::string const& e, std::string const& name, std::string const& desc = "u" )
    {
        element_type u( this->shared_from_this(), name, desc );
        bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
        EntityProcessType entityProcess = (addExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        u.on( _range=elements(M_mesh,entityProcess), _expr=expr<nComponents1,nComponents2>(e) );
        return u;
    }

    /**
     * \return the element 0 of the function space
     */
    element_ptrtype
    elementPtr( std::string const& name = "u", std::string const& desc = "u" )
    {
        return std::make_shared<element_type>( this->shared_from_this(), name, desc );
    }

    /**
     * \param e expression to initialize the element
     * \param u name of the element
     * \return a pointer to an element initialized with expression \p e
     */
    template<typename ExprT>
    element_ptrtype
    elementPtr( ExprT e, std::string const& name = "u", std::string const& desc = "u",
                typename std::enable_if<std::is_base_of<ExprBase,ExprT>::value >::type* = 0 )
    {
        //return std::make_shared<element_type>( e, name, desc );
        element_ptrtype u = this->elementPtr(name,desc);
        bool addExtendedElt = this->dof()->buildDofTableMPIExtended();
        EntityProcessType entityProcess = (addExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        u->on( _range=elements(M_mesh,entityProcess), _expr=e );
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
            u[i] = std::make_shared<element_type>( this->shared_from_this(), name );
        }
        return u;
    }
    /**
     * map of elements
     * \return a map of elements
     */
    static std::map<std::string,element_type>
    elementsMap()
    {
        return std::map<std::string,element_type>{};
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
            u[i] = std::make_shared<element_type>( e, name );
        }
        return u;
    }

    /**
     * \param vec input vector
     * \param blockIdStart if vec was built from a VectorBlock, need to specify the id of first block
     * \return the element of the function space with value shared by the input vector
     */
    element_type
    element( std::shared_ptr<Vector<value_type> > const& vec, int blockIdStart = 0 )
    {
        return this->element( *vec, blockIdStart );
    }

    /**
     * \param vec input vector
     * \param blockIdStart if vec was built from a VectorBlock, need to specify the id of first block
     * \return the element of the function space which use the storage values of the input vector
     */
    element_type
    element( Vector<value_type> const& vec, int blockIdStart = 0 )
    {
#if FEELPP_HAS_PETSC
        VectorPetsc<value_type> * vecPetsc = const_cast< VectorPetsc<value_type> *>( dynamic_cast< VectorPetsc<value_type> const*>( &vec ) );
        //VectorPetscMPI<value_type> * vecPetsc = const_cast< VectorPetscMPI<value_type> *>( dynamic_cast< VectorPetscMPI<value_type> const*>( &vec ) );
        CHECK( vecPetsc ) << "only petsc vector";

        auto const& dmVec = vec.map();
        CHECK( blockIdStart < dmVec.numberOfDofIdToContainerId() ) << "invalid blockId : " << blockIdStart << " must be less than " << dmVec.numberOfDofIdToContainerId();
        size_type nActiveDof = this->dof()->nLocalDofWithoutGhost();
        value_type* arrayActiveDof = (nActiveDof>0)? std::addressof( (*vecPetsc)( dmVec.dofIdToContainerId(blockIdStart,0) ) ) : nullptr;
        size_type nGhostDof = this->dof()->nLocalGhosts();
        size_type nActiveDofFirstSubSpace = (is_composite)? this->template functionSpace<0>()->dof()->nLocalDofWithoutGhost() : nActiveDof;
        value_type* arrayGhostDof = (nGhostDof>0)? std::addressof( (*vecPetsc)( dmVec.dofIdToContainerId(blockIdStart,nActiveDofFirstSubSpace) ) ) : nullptr;
        element_type u( this->shared_from_this(),
                nActiveDof, arrayActiveDof,
                nGhostDof, arrayGhostDof );
#else
        LOG(WARNING) << "element(Vector<value_type> const& vec, int blockIdStart): This function is disabled when Feel++ is not built with PETSc";
        element_type u;
#endif
        return u;
    }

    /**
     * \param vec input vector
     * \param blockIdStart if vec was built from a VectorBlock, need to specify the id of first block
     * \return the element of the function space which use the storage values of the input vector
     */
    element_ptrtype
    elementPtr( Vector<value_type> const& vec, int blockIdStart = 0 )
    {
#if FEELPP_HAS_PETSC
        VectorPetsc<value_type> * vecPetsc = const_cast< VectorPetsc<value_type> *>( dynamic_cast< VectorPetsc<value_type> const*>( &vec ) );
        //VectorPetscMPI<value_type> * vecPetsc = const_cast< VectorPetscMPI<value_type> *>( dynamic_cast< VectorPetscMPI<value_type> const*>( &vec ) );
        CHECK( vecPetsc ) << "only petsc vector";

        auto const& dmVec = vec.map();
        CHECK( blockIdStart < dmVec.numberOfDofIdToContainerId() ) << "invalid blockId : " << blockIdStart << " must be less than " << dmVec.numberOfDofIdToContainerId();
        size_type nActiveDof = this->dof()->nLocalDofWithoutGhost();
        value_type* arrayActiveDof = (nActiveDof>0)? std::addressof( (*vecPetsc)( dmVec.dofIdToContainerId(blockIdStart,0) ) ) : nullptr;
        size_type nGhostDof = this->dof()->nLocalGhosts();
        size_type nActiveDofFirstSubSpace = (is_composite)? this->template functionSpace<0>()->dof()->nLocalDofWithoutGhost() : nActiveDof;
        value_type* arrayGhostDof = (nGhostDof>0)? std::addressof( (*vecPetsc)( dmVec.dofIdToContainerId(blockIdStart,nActiveDofFirstSubSpace) ) ) : nullptr;
        element_ptrtype u( new element_type(
                    this->shared_from_this(),
                    nActiveDof, arrayActiveDof,
                    nGhostDof, arrayGhostDof )
                );
#else
        LOG(WARNING) << "element(Vector<value_type> const& vec, int blockIdStart): This function is disabled when Feel++ is not built with PETSc";
        element_type u;
#endif
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
     * get the \p i -th \c FunctionSpace out the list
     */
    template<int i>
    typename mpl::at_c<functionspace_vector_type,i>::type const&
    functionSpace() const
    {
        return fusion::at_c<i>( M_functionspaces );
    }


    /**
     * \return true if component Space was built
    */
    bool hasCompSpace() const
    {
        return ( M_comp_space )? true : false;
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
        if constexpr ( nDim > 2 )
            return trace_trace_functionspace_type::New( _mesh=mesh()->wireBasket( markededges( mesh(),"WireBasket" ) ), _worldscomm=this->worldsComm() );
        return trace_trace_functionspace_ptrtype{};
    }
    template<typename RangeT>
    trace_trace_functionspace_ptrtype
    wireBasket( RangeT range  )  const
    {
        if constexpr ( nDim > 2 )
            return trace_trace_functionspace_type::New( mesh()->wireBasket( range ) );
        return trace_trace_functionspace_ptrtype{};
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
    void rebuildDofPoints();

    //! return the index of dof defined on a range
    template <typename RangeType>
    std::set<size_type> dofs( RangeType const& range,
                              ComponentType c1 = ComponentType::NO_COMPONENT,
                              bool onlyMultiProcessDofs = false ) const
        {
            std::set<size_type> res;
            if ( onlyMultiProcessDofs && this->worldComm().localSize() == 1 )
                return res;
            this->dofs( range, c1, onlyMultiProcessDofs, mpl::bool_<is_composite>(), res );
            return res;
        }

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

    /**
     * @brief save on disk some information in json + hdF5 files as doftable (require HDF5 support)
     *
     * @param filepath path of the json files generated (extension can be automatically added if not given)
     */
    template <typename TT=functionspace_type,std::enable_if_t< !TT::is_composite, bool> = true >
    void save( std::string const& filepathstr ) const
        {
            fs::path argfilepath = filepathstr;
            fs::path jsonfilepath = argfilepath.replace_extension("json");
            fs::path h5filepath = argfilepath.replace_extension("h5");

            auto const& mpicomm = this->worldComm().localComm();
            rank_type mpirank = this->worldComm().localRank();

            if ( this->worldComm().isMasterRank() )
            {
                fs::path filedir = (jsonfilepath.is_relative()? fs::absolute(jsonfilepath) : jsonfilepath).parent_path();
                if ( !fs::exists( filedir ) )
                    fs::create_directories( filedir );
            }
            this->worldComm().barrier();

            // save json file
            nl::json jdata;
            this->updateInformationObject( jdata["info"] );
            if ( this->worldComm().isMasterRank() )
            {
                std::ofstream ojson( jsonfilepath );
                ojson << jdata.dump(/*1*/);
            }

            // serialize doftable
            std::vector<uint> doftableSerialization;
            auto meshsupport = this->template meshSupport<0>();
            auto dof = this->dof();
            std::vector<index_type> dd;
            for ( auto const& eltWrap : elements(meshsupport) )
            {
                auto const& elt = unwrap_ref(eltWrap);
                doftableSerialization.push_back( elt.id() );
                for ( int p=0;p<elt.nVertices();++p )
                    doftableSerialization.push_back( elt.point(p).id() );

                auto dofMapping = dof->localDof( elt.id() );
                dd.resize( std::distance( Feel::begin(dofMapping),Feel::end(dofMapping ) ) );
                for( auto const& ldof : dofMapping )
                {
                    index_type thedof = ldof.second.index();
                    uint16_type thelocdof = ldof.first.localDof();
                    index_type thedofgp = dof->mapGlobalProcessToGlobalCluster( thedof );
                    dd[thelocdof] = thedofgp;
                }

                doftableSerialization.insert( doftableSerialization.end(), dd.begin(), dd.end() );
            }


            // get size/position data of hdf5 file
            index_type local_sum = doftableSerialization.size();
            std::vector<index_type> processTolocalSum;
            mpi::all_gather( mpicomm, local_sum,processTolocalSum );
            index_type global_sum = std::accumulate(processTolocalSum.begin(), processTolocalSum.end(), 0);
            index_type offset = std::accumulate(processTolocalSum.begin(), std::next( processTolocalSum.begin(), mpirank) , 0);

            // save hdf5 file
            HDF5 hdf5;
            hdf5.openFile( h5filepath.string(), mpicomm, false );

            bool useTransposedStorage = true;
            const int dimsComp0 = (useTransposedStorage)? 1 : 0;
            const int dimsComp1 = (useTransposedStorage)? 0 : 1;

            std::string tableName = "doftable";
            hsize_t dimsElt[2];
            dimsElt[dimsComp0] = global_sum;
            dimsElt[dimsComp1] = 1;
            hsize_t dimsElt2[2];
            dimsElt2[dimsComp0] = local_sum;
            dimsElt2[dimsComp1] = 1;
            hsize_t offsetElt[2];
            offsetElt[dimsComp0] = offset;
            offsetElt[dimsComp1] = 0;

            hdf5.createTable( tableName, H5T_NATIVE_UINT, dimsElt );
            if ( !doftableSerialization.empty() )
                hdf5.write( tableName, H5T_NATIVE_UINT, dimsElt2, offsetElt, doftableSerialization.data() );
            hdf5.closeTable( tableName );

            hdf5.closeFile();
        }

    /**
     * @brief get relation with current doftable and one stored on disk (mesh points should be same)
     *
     * @param filepath path of the json files on the disk
     */
    template <typename TT=functionspace_type,std::enable_if_t< !TT::is_composite, bool> = true >
    std::vector<index_type> relationFromFile( std::string const& filepathstr ) const
        {
            fs::path argfilepath = filepathstr;
            fs::path jsonfilepath = argfilepath.replace_extension("json");
            fs::path h5filepath = argfilepath.replace_extension("h5");

            // fetch current mapping : (pt ids in elt) -> (dof in process ordered by local dof)
            using doftable_relation_type = std::unordered_map<std::vector /*set*/<index_type>, std::vector<index_type>, Feel::HashTables::HasherContainers<index_type>>;
            doftable_relation_type currentDofTableMapping;
            std::vector<index_type> ptIds;
            std::vector<index_type> dd;
            auto meshsupport = this->template meshSupport<0>();
            auto dof = this->dof();
            for ( auto const& eltWrap : elements(meshsupport) )
            {
                auto const& elt = unwrap_ref(eltWrap);
                ptIds.resize( elt.nVertices() );
                for ( int p=0;p<elt.nVertices();++p )
                    ptIds[p] = elt.point(p).id();

                auto dofMapping = dof->localDof( elt.id() );
                dd.resize( std::distance( Feel::begin(dofMapping),Feel::end(dofMapping ) ) );
                for( auto const& ldof : dofMapping )
                {
                    index_type thedof = ldof.second.index();
                    uint16_type thelocdof = ldof.first.localDof();
                    dd[thelocdof] = thedof;
                }

                currentDofTableMapping.insert( {ptIds, dd} );
            }

            // read hdf5 file
            HDF5 hdf5;
#if 0
            hdf5.openFile( h5filepath.string(), (subComm)? *subComm : this->comm().comm(), true );
#else
            hdf5.openFile( h5filepath.string(), this->worldComm().comm(), true );
#endif
            std::string tableName = "doftable";
            hsize_t dimsGlob[2];
            hsize_t offsetElt[2] = {0,0};
            hdf5.openTable( tableName, dimsGlob );

            std::vector<uint> dataReaded( dimsGlob[0]*dimsGlob[1] );

            hdf5.read( tableName, H5T_NATIVE_UINT, dimsGlob, offsetElt, dataReaded.data() );

            hdf5.closeTable( tableName );
            hdf5.closeFile();

            // update mapping
            uint16_type nVerticesInElt = ptIds.size();
            uint16_type nDofByElt = dd.size();
            std::vector<index_type> mappingWithFile( dof->nLocalDofWithGhost(), invalid_v<index_type> );

            for ( size_type k=0; k<dataReaded.size(); )
            {
                index_type eltId = dataReaded[k++];
                ptIds.resize( nVerticesInElt );
                for ( int p=0;p<nVerticesInElt;++p )
                    ptIds[p] = dataReaded[k++];
                auto const& curentLpDofs = currentDofTableMapping.at( ptIds );

                for ( int ld=0;ld<nDofByElt;++ld )
                    mappingWithFile[ curentLpDofs[ld] ] = dataReaded[k++];
            }

            return mappingWithFile;
        }

    template <typename TT=functionspace_type,std::enable_if_t< TT::is_composite, bool> = true >
    std::vector<index_type> relationFromFile( std::string const& filepathstr ) const
        {
            CHECK( false ) << "composite case not implemented";
            return {};
        }

    //@}


    // Copy constructor.
    FunctionSpace( FunctionSpace const& __fe )
        :
        super( __fe ),
        M_worldsComm( __fe.M_worldsComm ),
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


private:

    template<typename FSpaceHead, typename... FSpaceTail>
    FEELPP_NO_EXPORT void initList( FSpaceHead& fspacehead, FSpaceTail... fspacetail );

    FEELPP_NO_EXPORT void initList();

    template<typename FSpaceHead>
    FEELPP_NO_EXPORT void initHead( FSpaceHead& fspacehead );

    template <typename RangeType>
    void dofs( RangeType const& rangeElt, ComponentType c1, bool onlyMultiProcessDofs, mpl::false_, std::set<size_type> & res,
               std::enable_if_t< RangeType::isOnElements() >* = nullptr ) const
        {
            if ( c1 == ComponentType::NO_COMPONENT )
            {
                for ( auto const& eltWrap : rangeElt )
                {
                    auto const& elt = unwrap_ref( eltWrap );
                    for( auto const& ldof : this->dof()->localDof( elt.id() ) )
                    {
                        size_type index = ldof.second.index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                            res.insert( index );
                    }
                }
            }
            else
            {
                auto XhComp = this->compSpace();
                static const uint16_type nDofComponents = this->dof()->nDofComponents();
                int c1DofShift = ((int)c1);
                for ( auto const& eltWrap : rangeElt )
                {
                    auto const& elt = unwrap_ref( eltWrap );
                    for( auto const& ldof : XhComp->dof()->localDof( elt.id() ) )
                    {
                        size_type index = c1DofShift + nDofComponents*ldof.second.index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                        res.insert( index );
                    }
                }
            }
        }
    template <typename RangeType>
    void dofs( RangeType const& rangeFace, ComponentType c1, bool onlyMultiProcessDofs, mpl::false_, std::set<size_type> & res,
               std::enable_if_t< RangeType::isOnFaces() >* = nullptr ) const
        {
            if ( c1 == ComponentType::NO_COMPONENT )
            {
                for ( auto const& faceWrap : rangeFace )
                {
                    auto const& face = unwrap_ref( faceWrap );
                    auto facedof = this->dof()->faceLocalDof( face.id() );
                    for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                    {
                        size_type index = it->index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                            res.insert( index );
                    }
                }
            }
            else
            {
                auto XhComp = this->compSpace();
                static const uint16_type nDofComponents = this->dof()->nDofComponents();
                int c1DofShift = ((int)c1);
                for ( auto const& faceWrap : rangeFace )
                {
                    auto const& face = unwrap_ref( faceWrap );
                    auto facedof = XhComp->dof()->faceLocalDof( face.id() );
                    for ( auto it= facedof.first, en= facedof.second ; it!=en;++it )
                    {
                        size_type index = c1DofShift + nDofComponents*it->index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                        res.insert( index );
                    }
                }

            }
        }
    template <typename RangeType>
    void dofs( RangeType const& rangeFace, ComponentType c1, bool onlyMultiProcessDofs, mpl::false_, std::set<size_type> & res,
               std::enable_if_t< boost::tuples::template element<0, typename RangeType::super>::type::value == MESH_FACES && 
                                !std::is_same<typename RangeType::super,faces_reference_wrapper_t<mesh_type> >::value >* = nullptr ) const
        {
            CHECK(false) << "TODO";
        }

    template <typename RangeType>
    void dofs( RangeType const& rangeEdge, ComponentType c1, bool onlyMultiProcessDofs, mpl::false_, std::set<size_type> & res,
               std::enable_if_t< RangeType::isOnEdges() >* = nullptr ) const
        {
            size_type eid = invalid_v<size_type>;
            uint16_type edgeid_in_element;
            if ( c1 == ComponentType::NO_COMPONENT )
            {
                for ( auto const& edgeWrap : rangeEdge )
                {
                    auto const& edge = unwrap_ref( edgeWrap );
                    auto itEltInfo = edge.elements().begin();
                    if ( itEltInfo == edge.elements().end() )
                        continue;
                    eid = invalid_v<size_type>;
                    for ( auto const& eltConnectedToEdge : edge.elements() )
                    {
                        size_type eltIdConnected = eltConnectedToEdge.first;
                        if ( this->dof()->isElementDone( eltIdConnected ) )
                        {
                            eid = eltIdConnected;
                            edgeid_in_element = eltConnectedToEdge.second;
                            break;
                        }
                    }
                    if ( eid == invalid_v<size_type> )
                        continue;

                    for( auto const& ldof : this->dof()->edgeLocalDof( eid, edgeid_in_element ) )
                    {
                        size_type index = ldof.index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                        res.insert( index );
                    }
                }
            }
            else
            {
                auto XhComp = this->compSpace();
                static const uint16_type nDofComponents = this->dof()->nDofComponents();
                int c1DofShift = ((int)c1);
                for ( auto const& edgeWrap : rangeEdge )
                {
                    auto const& edge = unwrap_ref( edgeWrap );
                    auto itEltInfo = edge.elements().begin();
                    if ( itEltInfo == edge.elements().end() )
                        continue;
                    eid = invalid_v<size_type>;
                    for ( auto const& eltConnectedToEdge : edge.elements() )
                    {
                        size_type eltIdConnected = eltConnectedToEdge.first;
                        if ( this->dof()->isElementDone( eltIdConnected ) )
                        {
                            eid = eltIdConnected;
                            edgeid_in_element = eltConnectedToEdge.second;
                            break;
                        }
                    }
                    if ( eid == invalid_v<size_type> )
                        continue;

                    for( auto const& ldof : XhComp->dof()->edgeLocalDof( eid, edgeid_in_element ) )
                    {
                        size_type index = c1DofShift + nDofComponents*ldof.index();
                        if ( !onlyMultiProcessDofs ||
                             ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                               this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                        res.insert( index );
                    }
                }
            }
        }
    template <typename RangeType>
    void dofs( RangeType const& rangePoint, ComponentType c1, bool onlyMultiProcessDofs, mpl::false_, std::set<size_type> & res,
               std::enable_if_t< RangeType::isOnPoints()>* = nullptr ) const
        {
            std::vector<uint16_type> compUsed;
            static const uint16_type nDofComponents = this->dof()->nDofComponents();
            if ( c1 == ComponentType::NO_COMPONENT )
            {
                for( uint16_type c = 0; c < nDofComponents; ++c )
                    compUsed.push_back( c );
            }
            else
                compUsed.push_back( (int)c1 );

            size_type eid = invalid_v<size_type>;
            uint16_type ptid_in_element;
            for ( auto const& pointWrap : rangePoint )
            {
                auto const& point = unwrap_ref( pointWrap );
                auto itPointInfo = point.elements().begin();
                if ( itPointInfo == point.elements().end() )
                    continue;
                eid = invalid_v<size_type>;
                for ( auto const& eltConnectedToPoint : point.elements() )
                {
                    size_type eltIdConnected = eltConnectedToPoint.first;
                    if ( this->dof()->isElementDone( eltIdConnected ) )
                    {
                        eid = eltIdConnected;
                        ptid_in_element = eltConnectedToPoint.second;
                        break;
                    }
                }
                if ( eid == invalid_v<size_type> )
                    continue;

                for( uint16_type c : compUsed )
                {
                    size_type index = this->dof()->localToGlobal( eid, ptid_in_element, c ).index();
                    if ( !onlyMultiProcessDofs ||
                         ( this->dof()->dofGlobalProcessIsGhost( index ) ||
                           this->dof()->activeDofSharedOnCluster().find( index ) != this->dof()->activeDofSharedOnCluster().end() ) )
                    res.insert( index );
                }
            }
        }
    template <typename RangeType>
    void dofs( RangeType const& rangeElt, ComponentType c1, bool onlyMultiProcessDofs, mpl::true_, std::set<size_type> & res ) const
        {
            CHECK(false) << "not implemented with composite space";
        }

    friend class ComponentSpace;
    class ComponentSpace
    {
    public:
        using functionspace_type = FunctionSpace<MeshTypes, BasisTypes, T, PeriodicityType, MortarType>;
        using functionspace_ptrtype = functionspace_type*;
        using functionspace_cptrtype = functionspace_type const*;

        using component_functionspace_type = typename functionspace_type::component_functionspace_type;
        using component_functionspace_ptrtype = typename functionspace_type::component_functionspace_ptrtype;
        using component_functionspace_cptrtype = component_functionspace_type const*;


        ComponentSpace( functionspace_ptrtype  __functionspace,
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

        functionspace_ptrtype  M_functionspace;
        mesh_ptrtype M_mesh;
    };
public :
    //! update information for the current object
    void updateInformationObject( nl::json & p ) const override;

protected:

    //friend class FunctionSpace<mesh_type, typename bases_list::component_basis_type, value_type>;
    //friend class FunctionSpace<mesh_type, bases_list, value_type>;

    worldscomm_ptr_t M_worldsComm;

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
void
FunctionSpace<A0, A1, A2, A3, A4>::init( mesh_ptrtype const& __m,
                                         mesh_support_vector_type const& meshSupport,
                                         size_type mesh_components,
                                         std::vector<Dof<typename mesh_type::size_type>> const& dofindices,
                                         periodicity_type periodicity )
{
    Feel::Context ctx( mesh_components );
    DVLOG( 2 ) << "component     MESH_RENUMBER: " << ctx.test( MESH_RENUMBER ) << "\n";
    DVLOG( 2 ) << "component MESH_UPDATE_EDGES: " << ctx.test( MESH_UPDATE_EDGES ) << "\n";
    DVLOG( 2 ) << "component MESH_UPDATE_FACES: " << ctx.test( MESH_UPDATE_FACES ) << "\n";
    DVLOG( 2 ) << "component    MESH_PARTITION: " << ctx.test( MESH_PARTITION ) << "\n";
    if constexpr ( !is_composite )
    {
        DVLOG(2) << "calling init(<space>) begin\n";
        DVLOG(2) << "calling init(<space>) is_periodic: " << is_periodic << "\n";

        M_mesh = __m;
        M_periodicity = periodicity;
        VLOG(1) << "FunctionSpace init begin mesh use_count : " << M_mesh.use_count();

        if ( M_mesh->components().test( MESH_DO_NOT_UPDATE ) )
        {

            if ( basis_type::nDofPerEdge || nDim >= 3 )
                mesh_components |= MESH_UPDATE_EDGES;

            /*
             * update faces info in mesh only if dofs exists on faces or the
             * expansion is continuous between elements. This case handles strong
             * Dirichlet imposition
             */
            if ( basis_type::nDofPerFace || is_continuous || nDim >= 3 )
                mesh_components |= MESH_UPDATE_FACES;

            if ( !M_mesh->isUpdatedForUse() )
            {
                M_mesh->components().set( mesh_components );
                M_mesh->updateForUse();
            }
        }
        if ( is_periodic )
        {
            M_mesh->removeFacesFromBoundary( { periodicity.tag1(), periodicity.tag2() } );
        }

        M_ref_fe = std::make_shared<basis_type>();

        tic();
        tic();
        M_dof = std::make_shared<dof_type>( M_ref_fe, fusion::at_c<0>(periodicity), *this->worldsComm()[0] );
        toc("FunctionSpace dof-1", FLAGS_v>0);
        tic();
        M_dof->setBuildDofTableMPIExtended( this->extendedDofTable() );
        toc("FunctionSpace dof-2", FLAGS_v>0);
        DVLOG(2) << "[functionspace] Dof indices is empty ? " << dofindices.empty() << "\n";
        tic();
        M_dof->setDofIndices( dofindices );
        toc("FunctionSpace dof-3", FLAGS_v>0);
        DVLOG(2) << "[functionspace] is_periodic = " << is_periodic << "\n";
        tic();
        if ( fusion::at_c<0>( meshSupport ) && fusion::at_c<0>( meshSupport )->isPartialSupport() )
            M_dof->setMeshSupport( fusion::at_c<0>( meshSupport ) );
        toc("FunctionSpace dof-4", FLAGS_v>0);
        tic();
        M_dof->build( M_mesh );
        toc("FunctionSpace dof-5", FLAGS_v>0);
        toc("FunctionSpace dof table", FLAGS_v > 0 );
        M_dofOnOff = M_dof;

        this->applyUpdateInformationObject();

        DVLOG(2) << "nb dim : " << qDim() << "\n";
        DVLOG(2) << "nb dof : " << nDof() << "\n";
        DVLOG(2) << "nb dof per component: " << nDofPerComponent() << "\n";


        //detail::searchIndicesBySpace<proc_dist_map_type>( this, procDistMap);

        DVLOG(2) << "calling init(<space>) end\n";
        VLOG(1) << "FunctionSpace init begin mesh use_count : " << M_mesh.use_count();

#if !defined( __INTEL_COMPILER )
        if ( boption( _name="connect" ) )
            M_mesh->addObserver( *this );
#endif
    }
    else if constexpr ( is_composite )
    {
        M_mesh = __m;

        // todo : check worldsComm size and M_functionspaces are the same!
        mpl::range_c<int,0,nSpaces> keySpaces;
        fusion::for_each( keySpaces,
                        Feel::detail::InitializeSpace<functionspace_type>( M_functionspaces,__m, meshSupport, periodicity,
                                                                            dofindices,
                                                                            this->worldsComm(),
                                                                            this->extendedDofTableComposite() ) );

        this->initList();
    }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::initList()
{
    if constexpr ( is_composite )
    {
        if ( !this->hasWorldComm() )
            this->setWorldComm( M_worldsComm[0] );

        if ( true )// this->worldComm().globalSize()>1 )
        {
            if ( this->hasEntriesForAllSpaces() )
                {
                    // construction with same partionment for all subspaces
                    // and each processors has entries for all subspaces
                    DVLOG(2) << "init(<composite>) type hasEntriesForAllSpaces\n";

                    // build datamap
                    auto dofInitTool=Feel::detail::updateDataMapProcessStandard<dof_type>( this->worldCommPtr(),
                                                                                        this->nSubFunctionSpace() );
                    M_dof = fusion::fold( M_functionspaces, M_dof, dofInitTool );
                    // finish update datamap
                    M_dof->setNDof( this->nDof() );
                    M_dofOnOff = M_dof;
                }
            else
                {
                    CHECK( false ) << "deprecated";
                    // construction with same partionment for all subspaces
                    // and one processor has entries for only one subspace
                    DVLOG(2) << "init(<composite>) type Not hasEntriesForAllSpaces\n";

                    // build the WorldComm associated to mix space
                    worldcomm_ptr_t mixSpaceWorldComm = this->worldsComm()[0]->clone();

                    if ( this->worldsComm().size()>1 )
                        for ( int i=1; i<( int )this->worldsComm().size(); ++i )
                            {
                                mixSpaceWorldComm = *mixSpaceWorldComm + *this->worldsComm()[i];
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
    #if 0
        M_dof->setIndexSplit( this->buildDofIndexSplit() );
        M_dof->setIndexSplitWithComponents( this->buildDofIndexSplitWithComponents() );
    #endif
        //M_dof->indexSplit().showMe();
        this->applyUpdateInformationObject();
    }
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
void
FunctionSpace<A0, A1, A2, A3, A4>::buildComponentSpace() const
{
    if ( ( is_vectorial || is_tensor2 || is_tensor2symm ) && !M_comp_space )
    {
        // Warning: this works regarding the communicator . for the component space
        // it will use in mixed spaces only numberofSudomains/numberofspace processors
        //
        auto meshSupport = Feel::detail::FunctionSpaceMeshSupport<functionspace_type>( *this ).M_meshSupportVector;
        M_comp_space = component_functionspace_type::New(_mesh=M_mesh,
                                                         _worldscomm=this->worldsComm(),
                                                         _periodicity=M_periodicity,
                                                         _extended_doftable=this->extendedDofTableComposite(),
                                                         _range=meshSupport );

        VLOG(2) << " - component space :: nb dim : " << M_comp_space->qDim() << "\n";
        VLOG(2) << " - component space :: nb dof : " << M_comp_space->nDof() << "\n";
        VLOG(2) << " - component space :: nb dof per component: " << M_comp_space->nDofPerComponent() << "\n";
    }
}
template <typename A0, typename A1, typename A2, typename A3, typename A4>
void FunctionSpace<A0, A1, A2, A3, A4>::rebuildDofPoints()
{
    if constexpr ( !is_composite )
    {
        M_dof->rebuildDofPoints( *M_mesh );
    }
    else
    {
        fusion::for_each( M_functionspaces, Feel::detail::rebuildDofPointsTool() );
    }
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
#if 0
    mesh_element_iterator it = M_mesh->beginElementWithProcessId( M_mesh->comm().rank() );
    mesh_element_iterator en = M_mesh->endElementWithProcessId( M_mesh->comm().rank() );
#else
    auto it = M_mesh->beginElement();
    auto en = M_mesh->endElement();
#endif

    for ( size_type __i = 0; it != en; ++__i, ++it )
    {
        __bb.make( it->second.G() );

        for ( unsigned k=0; k < __bb.min.size(); ++k )
        {
            __bb.min[k]-=EPS;
            __bb.max[k]+=EPS;
        }

        __rt->addBox( __bb.min, __bb.max, it->second.id() );
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

        auto rangeElements = M_mesh->elementsWithProcessId();
        auto it = std::get<0>( rangeElements );
        auto en = std::get<1>( rangeElements );

        for ( size_type __i = 0; it != en; ++__i, ++it )
        {
            auto const& elt = boost::unwrap_ref( *it );
            __bb.make( elt.G() );

            for ( unsigned k=0; k < __bb.min.size(); ++k )
            {
                __bb.min[k]-=EPS;
                __bb.max[k]+=EPS;
            }

            __rt->addBox( __bb.min, __bb.max, elt.id() );
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

    std::pair<size_type,value_type> closest = std::make_pair( invalid_v<size_type>, -1e30 );

    region_tree_type::pbox_set_type::const_iterator it = boxlst.begin();
    region_tree_type::pbox_set_type::const_iterator ite = boxlst.end();

    for ( ; it != ite; ++it )
    {
        inv_trans_type __git( M_mesh->gm(), M_mesh->element( ( *it )->id ), this->worldComm().subWorldCommSeqPtr() );

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

template <typename SpaceType>
struct UpdateInformationObject
{
    UpdateInformationObject( SpaceType const& space, nl::json& p )
        : M_space( space ), M_p( p ) {}
    template <typename T>
    void operator()( T const& t ) const
    {
        std::string jsname = boost::fusion::at_c<T::value>( M_space.functionSpaces() )->journalSection().to_string();
        M_p.push_back( jsname );
    }
    SpaceType const& M_space;
    nl::json& M_p;
};

template <typename A0, typename A1, typename A2, typename A3, typename A4>
void FunctionSpace<A0, A1, A2, A3, A4>::updateInformationObject( nl::json& p ) const
{
    if constexpr ( !is_composite )
    {
        if ( p.contains( "nSpace" ) )
            return;

        p["nSpace"] = functionspace_type::nSpaces;

        if ( this->mesh() )
            p["mesh"] = this->mesh()->journalSection().to_string();

        if ( M_ref_fe )
        {
            std::string shape;
            if ( is_scalar )
                shape = "scalar";
            else if ( is_vectorial )
                shape = "vectorial";
            else if ( is_tensor2 )
                shape = "tensor2";
            else if ( is_tensor2symm )
                shape = "tensor2symm";
            p.emplace( "basis", nl::json( { { "name", basisName() },
                                            { "order", basisOrder().front() },
                                            { "shape", shape },
                                            { "is_continuous", is_continuous },
                                            { "nComponents", nComponents },
                                            { "nComponents1", nComponents1 },
                                            { "nComponents2", nComponents2 },
                                            { "nLocalDof", fe_type::nLocalDof } } ) );
            if ( is_tensor2symm )
                p["/basis/nRealComponents"_json_pointer] = nRealComponents;
        }

        if ( M_dof ) // sometime not define (example online rbspace)
        {
            nl::json subPt;
            subPt["nDof"] = this->nDof();
            rank_type nProc = this->worldComm().localSize();
            if ( nProc > 1 )
            {
                nl::json::array_t subPt1, subPt2, subPt3;
                for ( rank_type p = 0; p < nProc; ++p )
                {
                    subPt1.push_back( this->dof()->nLocalDofWithGhost( p ) );
                    subPt2.push_back( this->dof()->nLocalDofWithoutGhost( p ) );
                    subPt3.push_back( this->dof()->nLocalGhosts( p ) );
                }
                subPt.emplace( "nLocalDofWithGhost", subPt1 );
                subPt.emplace( "nLocalDofWithoutGhost", subPt2 );
                subPt.emplace( "nLocalGhost", subPt3 );
                subPt.emplace( "extended-doftable", this->dof()->buildDofTableMPIExtended() );
            }
            p.emplace( "doftable", std::move( subPt ) );
        }
    }
    else // composite case
    {
        if ( p.contains( "nSpace" ) )
            return;

        auto subPt = nl::json::array( {} );
        ;
        mpl::range_c<int, 0, functionspace_type::nSpaces> keySpaces;
        boost::fusion::for_each( keySpaces, UpdateInformationObject<functionspace_type>( *this, subPt ) );
        p.emplace( "nSpace", functionspace_type::nSpaces );
        if ( M_dof )
            p.emplace( "nDof", this->nDof() );
        p.emplace( "subfunctionspaces", subPt );
    }
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

} // Feel

#include <feel/feeldiscr/detail/element_impl.hpp>

namespace Feel {
//!
//! @return the support of a function space
//!
template<typename SpaceT, typename = std::enable_if_t<is_functionspace_v<SpaceT>>>
constexpr typename SpaceT::template GetMeshSupport<typename SpaceT::mesh_ptrtype,0>::ptrtype
support( std::shared_ptr<SpaceT> const& X )
{
    return X->template meshSupport<0>();
}

} // Feel

#endif /* FEELPP_DISCR_FUNCTIONSPACE_H */
