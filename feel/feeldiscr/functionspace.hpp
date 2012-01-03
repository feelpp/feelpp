/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2004-11-22

   Copyright (C) 2004 EPFL
   Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2004-11-22
*/
#ifndef __FunctionSpace_H
#define __FunctionSpace_H 1

#include <boost/static_assert.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/support/is_sequence.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/algorithm.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>

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
#include <feel/feelfilters/pointsettomesh.hpp>


namespace Feel
{
    namespace fusion = boost::fusion;
    namespace parameter = boost::parameter;

    namespace detail{

        template<typename T>
        struct vector_plus
        {
            std::vector<T> operator()( std::vector<T> const& v1, std::vector<T> const& v2 ) const
            {
                FEEL_ASSERT( v1.size() == v2.size() )( v1.size() )( v2.size() ).error( "invalid vector size for vector_plus<>");
                std::vector<T> res( v1.size() );
                for( size_type i = 0;i < v1.size(); ++i )
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
                        //Debug() "[ID] extent = " <<

                        M_id = id.M_id;
                    }
                return *this;

            }
            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
            {
                return M_id[q](c1,c2);
            }

            template <typename ExtentList>
            void resize(const ExtentList& extents) { M_id.resize( extents ); }

            array_type M_id;

            template<class Archive>
            void save(Archive & ar, const unsigned int /*version*/) const
            {
                size_type e1 = M_id.shape()[0];
                Debug( 5010 ) << "saving in archive e1= " << e1 << "\n";
                ar  & e1;
                size_type e2 = M_id.shape()[1];
                Debug( 5010 ) << "saving in archive e2= " << e2 << "\n";
                ar  & e2;
                size_type e3 = M_id.shape()[2];
                Debug( 5010 ) << "saving in archive e3= " << e3 << "\n";
                ar  & e3;
                Debug( 5010 ) << "saving in archive array of size = " << M_id.num_elements() << "\n";
                ar  & boost::serialization::make_array( M_id.data(), M_id.num_elements() );
                Debug( 5010 ) << "saving in archive done\n";
            }
            template<class Archive>
            void load(Archive & ar, const unsigned int /*version*/)
            {
                size_type e1, e2, e3;
                ar  & e1;
                Debug( 5010 ) << "loading from archive e1= " << e1 << "\n";
                ar  & e2;
                Debug( 5010 ) << "loading from archive e2= " << e2 << "\n";
                ar  & e3;
                Debug( 5010 ) << "loading from archive e3= " << e3 << "\n";
                M_id.resize( boost::extents[e1] );
                Debug( 5010 ) << "loading from archive array of size = " << M_id.num_elements() << "\n";
                ar  & boost::serialization::make_array( M_id.data(), M_id.num_elements() );
                Debug( 5010 ) << "loading from archive done\n";
                Debug( 5010 ) << "creating view interpolation context done\n";
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
                _M_grad()
            {}

            template<typename Elem, typename ContextType>
            DD( Elem const& elem, ContextType const & context )
                :
                _M_grad( elem.gradExtents( context ) )
            {
                elem.grad_( context, _M_grad );
            }

            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
            {
                return _M_grad[q](c1,c2);
            }
            array_type _M_grad;

            template<class Archive>
            void save(Archive & ar, const unsigned int /*version*/) const
            {
                size_type e1 = _M_grad.shape()[0];
                Debug( 5010 ) << "saving in archive e1= " << e1 << "\n";
                ar  & e1;
                size_type e2 = _M_grad.shape()[1];
                Debug( 5010 ) << "saving in archive e2= " << e2 << "\n";
                ar  & e2;
                size_type e3 = _M_grad.shape()[2];
                Debug( 5010 ) << "saving in archive e3= " << e3 << "\n";
                ar  & e3;
                Debug( 5010 ) << "saving in archive array of size = " << _M_grad.num_elements() << "\n";
                ar  & boost::serialization::make_array( _M_grad.data(), _M_grad.num_elements() );
                Debug( 5010 ) << "saving in archive done\n";
            }
            template<class Archive>
            void load(Archive & ar, const unsigned int /*version*/)
            {
                size_type e1, e2, e3;
                ar  & e1;
                Debug( 5010 ) << "loading from archive e1= " << e1 << "\n";
                ar  & e2;
                Debug( 5010 ) << "loading from archive e2= " << e2 << "\n";
                ar  & e3;
                Debug( 5010 ) << "loading from archive e3= " << e3 << "\n";
                _M_grad.resize( boost::extents[e1] );
                Debug( 5010 ) << "loading from archive array of size = " << _M_grad.num_elements() << "\n";
                ar  & boost::serialization::make_array( _M_grad.data(), _M_grad.num_elements() );
                Debug( 5010 ) << "loading from archive done\n";
                Debug( 5010 ) << "creating view interpolation context done\n";
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
                _M_grad()
            {}

            template<typename Elem, typename ContextType>
            D( Elem const& elem, ContextType const & context )
                :
                _M_grad( elem.dExtents( context ) )
            {
                elem.d_( N, context, _M_grad );
            }

            value_type operator()( uint16_type c1, uint16_type /*c2*/, uint16_type q  ) const
            {
                return _M_grad[q](c1,0);
            }
            array_type _M_grad;
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
                _M_div()
            {}

            template<typename Elem, typename ContextType>
            Div( Elem const& elem, ContextType const & context )
                :
                _M_div( elem.divExtents( context ) )
            {
                elem.div_( context, _M_div );
#if 0
                uint16_type nComponents1 = elem.nComponents1;
                std::fill( _M_div.data(), _M_div.data()+_M_div.num_elements(), value_type( 0 ) );

                _M_grad( elem.div_( context, pc, _M_div ) ),


                    const uint16_type nq = context.xRefs().size2();
                for( int c1 = 0; c1 < nComponents1; ++c1 )
                    for( uint16_type q = 0; q < nq ; ++q )
                        {
                            _M_div[q](0,0) += _M_grad[q](c1,c1);
                        }
#endif
            }

            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
            {
                return _M_div[q](c1,c2);
            }
            array_type _M_div;
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
                _M_curl()
            {}

            template<typename Elem, typename ContextType>
            Curl( Elem const& elem, ContextType const & context )
                :
                _M_curl( elem.curlExtents( context ) )
            {
                init( elem, context, boost::is_same<mpl::int_<N>, mpl::int_<-1> >() );
            }
            template<typename Elem, typename ContextType>
            void
            init( Elem const& elem, ContextType const& context, mpl::bool_<true> )
                {
                    elem.curl_( context, _M_curl );
                }
            template<typename Elem, typename ContextType>
            void
            init( Elem const& elem, ContextType const& context, mpl::bool_<false> )
                {
                    elem.curl_( context, _M_curl, N );
                }

            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
            {
                return this->operator()( c1, c2, q, mpl::int_<N>() );
            }
            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<-1>  ) const
            {
                return _M_curl[q](c1,c2);
            }
            value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<0>  ) const
            {
                return _M_curl[q](0,0);
            }
            value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<1>  ) const
            {
                return _M_curl[q](0,0);
            }
            value_type operator()( uint16_type /*c1*/, uint16_type /*c2*/, uint16_type q, mpl::int_<2>  ) const
            {
                return _M_curl[q](0,0);
            }
            array_type _M_curl;
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
                _M_hess()
            {}

            template<typename Elem, typename ContextType>
            H( Elem const& elem, ContextType const & context )
                :
                _M_hess( elem.hessExtents( context ) )
            {
                elem.hess_( context, _M_hess );
            }

            value_type operator()( uint16_type c1, uint16_type c2, uint16_type q  ) const
            {
                return _M_hess[q](c1,c2);
            }
            array_type _M_hess;
        };

        template<typename MeshType>
        struct InitializeSpace
        {
            InitializeSpace( boost::shared_ptr<MeshType> const& mesh, std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices )
                :
                _M_mesh( mesh ),
                _M_dofindices( dofindices )

            {}
            template <typename T>
            void operator()(boost::shared_ptr<T> & x) const
            {
                x = boost::shared_ptr<T>( new T( _M_mesh, _M_dofindices ) );
                FEEL_ASSERT( x ).error( "invalid function space" );
            }
            boost::shared_ptr<MeshType> _M_mesh;
            std::vector<boost::tuple<size_type, uint16_type, size_type> > const& _M_dofindices;
        };
        struct NbDof
        {
            typedef size_type result_type;
            NbDof( size_type start = 0, size_type size = invalid_size_type_value )
                :
                _M_cursor( start ),
                _M_finish( size )
            {}
            template<typename Sig>
            struct result;

            template<typename T, typename S>
#if BOOST_VERSION < 104200
            struct result<NbDof(T,S)>
#else
            struct result<NbDof(S,T)>
#endif
                :
                boost::remove_reference<S>
            {};
        template <typename T>
        size_type
        operator()(T const& x, size_type s ) const
        {
            size_type ret = s;
            if ( !x )
                return ret;
            if ( _M_cursor < _M_finish )
                ret += x->nDof();
            ++_M_cursor;
            return ret;
        }

        template <typename T>
        size_type
        operator()(size_type s, T const& x ) const
        {
            return this->operator()( x, s );
        }
    private:
        mutable size_type _M_cursor;
        size_type _M_finish;
    };

    struct NLocalDof
    {
        NLocalDof( size_type start = 0, size_type size = invalid_size_type_value )
            :
            _M_cursor( start ),
            _M_finish( size )
        {}
        template<typename Sig>
        struct result;

        template<typename T, typename S>
#if BOOST_VERSION < 104200
        struct result<NLocalDof(T,S)>
#else
        struct result<NLocalDof(S,T)>
#endif
            :
            boost::remove_reference<S>
        {};
    template <typename T>
    size_type
    operator()(T const& x, size_type s ) const
    {
        size_type ret = s;
        if ( _M_cursor < _M_finish )
            ret += x->nLocalDof();
        ++_M_cursor;
        return ret;
    }
    template <typename T>
    size_type
    operator()(size_type s, T const& x ) const
    {
        return this->operator()( x, s );
    }
private:
    mutable size_type _M_cursor;
    size_type _M_finish;
};
template< typename map_type >
struct searchIndicesBySpace
{
    searchIndicesBySpace()
    {}

    searchIndicesBySpace( map_type& /*u*/ )
    {}

    template<typename T>
    searchIndicesBySpace(T const& fspace, map_type& u)
    {
        u = getIndicesFromSpace(fspace,u);
    }

    template<typename Sig>
    struct result;

    template<typename T, typename M>
#if BOOST_VERSION < 104200
    struct result<searchIndicesBySpace(T,M)>
#else
    struct result<searchIndicesBySpace(M,T)>
#endif
        :
        boost::remove_reference<M>
    {};

template < typename T >
map_type getIndicesFromSpace(T const& fspace, map_type t ) const
{
    if ( fspace->mesh()->numElements() == 0 )
        return t;

    size_type nProc = fspace->dof()->nProcessors();

    //search for the biggest index already in t; this will give the shift for the dofs
    std::vector< size_type > max_per_space;

    for ( size_type j=0; j<t.size(); j++)
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
    for ( size_type i=0; i<nProc; i++)
        {
            /*
              std::cout << "Processor " << i << " has dofs"
              << " from " << fspace->dof()->firstDof(i)
              << " to " << fspace->dof()->lastDof(i) << "\n";
            */

            size_type _first = fspace->dof()->firstDof(i);
            size_type _last  = fspace->dof()->lastDof(i);

            //the dofs numbering for the current space start at max_index+1
            for ( size_type j=_first; j<=_last; j++)
                t[i].push_back( max_index + j );
        }

    return t;
}
template <typename T>
map_type
#if BOOST_VERSION < 104200
operator()(T const& fspace, map_type t ) const
#else
    operator()(map_type t, T const& fspace ) const
#endif
{
    return getIndicesFromSpace(fspace,t);
}
};

struct computeNDofForEachSpace
{
    typedef boost::tuple< uint, uint, std::vector<std::vector<int> > > result_type;

    template<typename T>
    result_type operator()(result_type const & previousRes, T const& t)
        {
            auto nDof = t->nDof();

            auto cptSpaces = previousRes.get<0>();
            auto start = previousRes.get<1>();
            auto is = previousRes.get<2>();

            is.push_back( std::vector<int>( nDof ) );
            for (uint i=0;i<nDof;++i) { is[cptSpaces][i] = start+i; }

            return boost::make_tuple( ++cptSpaces, (start+nDof), is );
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

    static const bool is_isoparametric = (PN == GN);
    static const bool is_subparametric = (PN > GN);
    static const bool is_surparametric = (PN < GN);
};

typedef parameter::parameters<
//    parameter::required<tag::mesh_type, mpl::or_<boost::is_base_and_derived<MeshBase,_> >, mpl::or_<fusion::traits::is_sequence<_>, mpl::is_sequence<_> > >
    parameter::required<tag::mesh_type, boost::is_base_and_derived<MeshBase,_> >
#if 1
    , parameter::optional<parameter::deduced<tag::bases_list>, mpl::or_<boost::is_base_and_derived<detail::bases_base,_>,
                                                                        mpl::or_<fusion::traits::is_sequence<_>,
                                                                                 mpl::is_sequence<_> > > >
#else
    , parameter::optional<parameter::deduced<tag::bases_list>, fusion::traits::is_sequence<_> >
#endif
    , parameter::optional<parameter::deduced<tag::value_type>, boost::is_floating_point<_> >
    , parameter::optional<parameter::deduced<tag::periodicity_type>, boost::is_base_and_derived<detail::periodicity_base,_> >
    > functionspace_signature;


/**
 * \class FunctionSpace
 * \ingroup SpaceTime
 * @author Function Space Class
 *
 * \c FunctionSpace is a representation of a functional space parametrized by
 * the type of the mesh (\c MeshType)
 *
 * @author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
    typedef typename parameter::binding<args, tag::periodicity_type, NoPeriodicity>::type periodicity_type;
    typedef typename parameter::binding<args, tag::bases_list, detail::bases<Lagrange<1,Scalar> > >::type bases_list;

    BOOST_MPL_ASSERT_NOT( ( boost::is_same<mpl::at<bases_list,mpl::int_<0> >, mpl::void_> ) );

private:

    template<typename MeshType,typename BasisType>
    struct ChangeMesh
    {
        typedef boost::shared_ptr<FunctionSpace<MeshType,detail::bases<BasisType>,value_type, periodicity_type> > type;
    };
    template<typename BasisType>
    struct ChangeBasis
    {
        typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                  mpl::identity<boost::shared_ptr<FunctionSpace<meshes_list,detail::bases<BasisType>,value_type, periodicity_type> > >,
                                  mpl::identity<typename mpl::transform<meshes_list, ChangeMesh<mpl::_1,BasisType>, mpl::back_inserter<fusion::vector<> > > > >::type::type type;
    };
    typedef typename mpl::transform<bases_list, ChangeBasis<mpl::_1>, mpl::back_inserter<fusion::vector<> > >::type functionspace_vector_type;

    template<typename BasisType>
    struct ChangeBasisToComponentBasis
    {
        typedef typename BasisType::component_basis_type component_basis_type;
        typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                  mpl::identity<boost::shared_ptr<FunctionSpace<meshes_list,detail::bases<component_basis_type>,value_type, periodicity_type> > >,
                                  mpl::identity<typename mpl::transform<meshes_list, ChangeMesh<mpl::_1,component_basis_type>, mpl::back_inserter<fusion::vector<> > > > >::type::type type;
    };

    typedef typename mpl::transform<bases_list,
                                    ChangeBasisToComponentBasis<mpl::_1>,
                                    mpl::back_inserter<fusion::vector<> > >::type component_functionspace_vector_type;

    template<typename BasisType>
    struct GetComponentBasis
    {
        typedef typename BasisType::component_basis_type type;
    };

    typedef typename mpl::transform<bases_list,
                                    GetComponentBasis<mpl::_1>,
                                    mpl::back_inserter<fusion::vector<> > >::type component_basis_vector_type;



public:

    /** @name Constants
     */
    //@{
    static const bool is_composite = ( mpl::size<bases_list>::type::value > 1 );

    template<int N>
    struct GetMesh
    {
        typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                  mpl::identity<mpl::identity<meshes_list> >,
                                  mpl::identity<mpl::at_c<meshes_list,N> > >::type::type::type type;
    };
    // mesh
    typedef meshes_list MeshesListType;
    typedef typename GetMesh<0>::type mesh_0_type;
    typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                              mpl::identity<meshes_list>,
                              mpl::identity<mesh_0_type> >::type::type mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef typename mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                              mpl::identity<typename mesh_type::element_type>,
                              mpl::identity<mpl::void_> >::type::type convex_type;

    template<typename BasisType>
    struct GetNComponents
    {
        //typedef mpl::int_<BasisType::template apply<mesh_type::nDim,value_type,typename mesh_type::element_type>::type::nComponents> type;
        typedef mpl::int_<BasisType::template apply<mesh_type::nDim,
                                                    mesh_type::nRealDim,
                                                    value_type,
                                                    typename mesh_type::element_type>::type::nComponents> type;
    };

    static const uint16_type nDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                             mpl::identity<mpl::int_<meshes_list::nDim> >,
                                             mpl::identity<mpl::int_<-1> > >::type::type::value;
    static const uint16_type nRealDim = mpl::if_<boost::is_base_of<MeshBase, meshes_list >,
                                             mpl::identity<mpl::int_<meshes_list::nRealDim> >,
                                             mpl::identity<mpl::int_<-1> > >::type::type::value;



    //typedef typename mpl::at_c<bases_list,0>::type::template apply<mesh_type::nDim,value_type,typename mesh_type::element_type>::type basis_0_type;
    typedef typename mpl::at_c<bases_list,0>::type::template apply<GetMesh<0>::type::nDim,
                                                                   GetMesh<0>::type::nRealDim,
                                                                   value_type,
                                                                   typename GetMesh<0>::type::element_type>::type basis_0_type;

    static const uint16_type rank = ( is_composite? invalid_uint16_type_value : basis_0_type::rank );
    static const bool is_scalar = ( is_composite? false : basis_0_type::is_scalar );
    static const bool is_vectorial = ( is_composite? false : basis_0_type::is_vectorial );
    static const bool is_tensor2 = ( is_composite? false : basis_0_type::is_tensor2 );
    static const bool is_continuous = ( is_composite? false : basis_0_type::isContinuous);
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

    static const bool is_periodic = periodicity_type::is_periodic;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef typename ublas::type_traits<value_type>::real_type real_type;
    typedef typename node<value_type>::type node_type;


    typedef FunctionSpace<A0,A1,A2,A3,A4> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef boost::shared_ptr<functionspace_type> pointer_type;

    typedef FunctionSpace<meshes_list, component_basis_vector_type, value_type, periodicity_type> component_functionspace_type;
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
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB, geoelement_type> >,
                              mpl::identity<mpl::void_> >::type::type gmc_type;
    typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_ptrtype>, mpl::identity<mpl::void_> >::type::type geopc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_type>, mpl::identity<mpl::void_> >::type::type geopc_type;

    // dof
    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<DofComposite>,
                              mpl::identity<DofTable<mesh_type, basis_type, periodicity_type> > >::type::type dof_type;

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
     * \class Element
     */
    template<typename T = double,  typename Cont = VectorUblas<T> >
    class Element
        :
        public Cont,boost::addable<Element<T,Cont> >, boost::subtractable<Element<T,Cont> >
    {
        template<typename BasisType>
        struct ChangeElement
        {
            BOOST_MPL_ASSERT_NOT( ( boost::is_same<BasisType,mpl::void_> ) );
            typedef FunctionSpace<mesh_type,detail::bases<BasisType>,value_type,periodicity_type> fs_type;
            typedef typename fs_type::template Element<value_type, typename VectorUblas<T>::range::type > element_type;
            typedef element_type type;
        };
        typedef typename mpl::transform<bases_list, ChangeElement<mpl::_1>, mpl::back_inserter<mpl::vector<> > >::type element_vector_type;
        typedef typename VectorUblas<T>::range::type ct_type;


    public:

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

        /** @name Typedefs
         */
        //@{
        typedef T value_type;

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

        void init( functionspace_ptrtype const& __functionspace,
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
            *this = project( this->functionspace(), elements(this->functionspace()->mesh()), __expr );
            return *this;
        }
#endif
        template<typename VectorExpr>
        Element& operator=( VectorExpr const& __v );


        /**
         * \return the container read-only
         */
        super const& container() const { return *this; }

        /**
         * \return the container read-write
         */
        super & container() { return *this; }

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
            auto s = ublas::slice( THECOMP, N_COMPONENTS, _M_functionspace->nDofPerComponent() );
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
            auto s = ublas::slice( THECOMP, N_COMPONENTS, _M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( THECOMP );
            return component_type( compSpace(),
                                   typename component_type::container_type( (VectorUblas<value_type>&)*this, s ),
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
            FEEL_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, _M_functionspace->nDofPerComponent() );
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
            FEEL_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, _M_functionspace->nDofPerComponent() );
            std::string __name = this->name() + "_" + componentToString( i );

            component_type c( compSpace(),
                              typename component_type::container_type( (VectorUblas<value_type>&)*this, s ),
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
            FEEL_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, _M_functionspace->nDofPerComponent() );

            std::string __name = this->name() + "_" + componentToString( i );
            component_type c( compSpace(),
                              typename component_type::container_type( this->vec().data().expression(), s ),
                              __name,
                              start()+i,
                              i);
            return c;
        }
        component_type
        comp( ComponentType i, mpl::bool_<false> )
        {
            FEEL_ASSERT( i >= X && i < N_COMPONENTS );
            auto s = ublas::slice( i, N_COMPONENTS, _M_functionspace->nDofPerComponent() );

            std::string __name = this->name() + "_" + componentToString( i );
            component_type c( compSpace(),
                              typename component_type::container_type( (VectorUblas<value_type>&)*this, s ),
                              __name,
                              start()+i,
                              i);
            return c;
        }

        value_type localToGlobal( size_type e, size_type l, int c ) const
        {
            size_type index=start()+boost::get<0>(_M_functionspace->dof()->localToGlobal( e, l, c ));
            return super::operator()( index );
        }
#if 0
        value_type& localToGlobal( size_type e, size_type l, int c )
        {
            size_type index=start()+boost::get<0>(_M_functionspace->dof()->localToGlobal( e, l, c ));
            return super::operator()( index );
        }
#endif
        value_type  operator()( size_t i ) const { return super::operator()( i ); }
        value_type& operator()( size_t i )       { return super::operator()( i ); }
        Element& operator+=( Element const& _e )
            {
                for( int i=0; i < _e.nLocalDof(); ++i )
                    this->operator()(i) += _e(i);
                return *this;
            }
        Element& operator-=( Element const& _e )
            {
                for( int i=0; i < _e.nLocalDof(); ++i )
                    this->operator()(i) -= _e(i);
                return *this;
            }
        /**
         * update element when mesh has been changed
         */
        void operator()( MESH_CHANGES mesh_changes )
        {
            Debug( 5010 ) << "Update element after a change in the mesh\n";
        }

        template<typename AE>
        container_type& assign( const ublas::vector_expression<AE> &ae )
        {
            return super::assign( ae );
        }
        void assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index=start()+boost::get<0>(_M_functionspace->dof()->localToGlobal( ie, il, c ));
            this->operator[]( index ) = __v;
        }
        void plus_assign( size_type ie, uint16_type il, uint16_type c, value_type const& __v )
        {
            size_type index=start()+boost::get<0>(_M_functionspace->dof()->localToGlobal( ie, il, c ));
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
        DataMap const& map() const { return _M_functionspace->map(); }

        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh()  { return _M_functionspace->mesh(); }

        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh() const { return _M_functionspace->mesh(); }

        /**
         * \return the element-wise square root
         */
        Element sqrt() const { Element _e( _M_functionspace ); for( int i=0; i < _e.nLocalDof(); ++i ) _e(i)  = math::sqrt(this->operator()(i)); return _e; }

        /**
         * \return the element-wise power to n
         */
        Element pow( int n ) const { Element _e( _M_functionspace ); for( int i=0; i < _e.nLocalDof(); ++i ) _e(i)  = math::pow(this->operator()(i),n); return _e; }

        // Only works for scalar fields
        template < typename p0_space_type >
        typename p0_space_type::element_type extremeValue( boost::shared_ptr<p0_space_type> const& P0h, std::string extreme )
        {
            // check if the mesh coming from P0h and the class elements is the same
            FEEL_ASSERT( P0h->mesh() == this->mesh() ).error("mesh is not the same");
            FEEL_ASSERT( is_scalar ).error("only works for scalar fields");

            typename p0_space_type::element_type p0Element( P0h );

            for ( auto elt_it = P0h->mesh()->beginElement(); elt_it != P0h->mesh()->endElement(); ++elt_it )
                {
                    size_type eid = elt_it->id();

                    size_type dofp0 = boost::get<0>(P0h->dof()->localToGlobal(eid, 0, 0));
                    std::vector<value_type> values ( functionspace_type::fe_type::nLocalDof );

                    size_type dofpn = 0;
                    for ( uint16_type local_id=0; local_id < functionspace_type::fe_type::nLocalDof; ++local_id )
                        {
                            dofpn = boost::get<0>(this->functionSpace()->dof()->localToGlobal(eid, local_id, 0));
                            values[local_id] = this->operator()(dofpn);
                        }

                    if ( extreme == "max" )
                        p0Element.assign(eid, 0, 0, *( std::max_element( values.begin(), values.end() ) ));
                    else
                        p0Element.assign(eid, 0, 0, *( std::min_element( values.begin(), values.end() ) ));
                }

            return p0Element;
        }

        template < typename p0_space_type >
        typename p0_space_type::element_type max( boost::shared_ptr<p0_space_type> const& P0h ) { return this->extremeValue(P0h, "max" ); }

        value_type max() { return this->container().max(); }

        template < typename p0_space_type >
        typename p0_space_type::element_type min( boost::shared_ptr<p0_space_type> const& P0h ) { return this->extremeValue(P0h, "min" ); }

        value_type min() { return this->container().min(); }

        //! Interpolation at a set of points
        //@{
        /**
         * data structure that stores the interpolated values of the
         * element at a set of points
         */
        typedef detail::ID<value_type,nComponents1,nComponents2> id_type;


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


        void
        idInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const;

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
        ptsInContext( Context_t const & context, mpl::int_<1>) const
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
            gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId()) );

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
        /**
         * interpolate the function at node (real coordinate) x
         *
         * @return the interpolated value of the function at the real point x
         */
        id_type operator()( node_type const& __x, bool extrapolate = false ) const;

        //@}

        //! gradient interpolation tool
        //@{

        typedef detail::DD<value_type, nComponents1, nRealDim> grad_type;
        typedef detail::D<value_type,0,nComponents1, 1> dx_type;
        typedef detail::D<value_type,1,nComponents1, 1> dy_type;
        typedef detail::D<value_type,2,nComponents1, 1> dz_type;

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
        gradInterpolate(matrix_node_type __ptsReal, grad_array_type& v ) const;

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
        dxInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const;

        void
        dyInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const;

        void
        dzInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const;


        //@}

        typedef detail::Div<value_type> div_type;

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
        divInterpolate( matrix_node_type __ptsReal, div_array_type& v ) const;

        typedef detail::Curl<value_type,-1,nRealDim> curl_type;
        typedef detail::Curl<value_type,0,1> curlx_type;
        typedef detail::Curl<value_type,1,1> curly_type;
        typedef detail::Curl<value_type,2,1> curlz_type;

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
        curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v ) const;

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
        curlxInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v ) const;

        void
        curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v ) const;

        void
        curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v ) const;

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
        hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v ) const;

        typedef detail::H<value_type,nRealDim,nRealDim> hess_type;

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
        functionspace_ptrtype const& functionSpace() const { return _M_functionspace; }

        /**
         * get the \c FunctionSpace vector
         */
        //functionspace_vector_type const&
        //functionSpaces() const { return _M_functionspace->functionSpaces(); }

        functionspace_vector_type const&
        functionSpaces() const { return _M_functionspace->functionSpaces(); }

        /**
         * get the \p i -th \c FunctionSpace out the list
         */
        template<int i>
        typename mpl::at_c<functionspace_vector_type,i>::type
        functionSpace() const { return _M_functionspace->template functionSpace<i>(); }

        template<int i>
        typename mpl::at_c<element_vector_type,i>::type
        element()
        {
            size_type nbdof_start =  fusion::accumulate( this->functionSpaces(),
                                                         size_type( 0 ),
                                                         detail::NLocalDof( 0, i ) );

            typename mpl::at_c<functionspace_vector_type,i>::type space( _M_functionspace->template functionSpace<i>() );
            Debug( 5010 ) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
            Debug( 5010 ) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
            Debug( 5010 ) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
            Debug( 5010 ) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";
            ct_type ct( *this, ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ) );

            Debug( 5010 ) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
            Debug( 5010 ) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";

            return typename mpl::at_c<element_vector_type,i>::type( space, ct );
        }

        template<int i>
        typename mpl::at_c<element_vector_type,i>::type
        element() const
        {
            size_type nbdof_start =  fusion::accumulate( _M_functionspace->functionSpaces(),
                                                         size_type( 0 ),
                                                         detail::NLocalDof( 0, i ) );
            typename mpl::at_c<functionspace_vector_type,i>::type space( _M_functionspace->template functionSpace<i>() );

            Debug( 5010 ) << "Element <" << i << ">::start :  "<< nbdof_start << "\n";
            Debug( 5010 ) << "Element <" << i << ">::size :  "<<  space->nDof()<< "\n";
            Debug( 5010 ) << "Element <" << i << ">::local size :  "<<  space->nLocalDof()<< "\n";
            Debug( 5010 ) << "Element <" << -1 << ">::size :  "<<  this->size() << "\n";

            ct_type ct( const_cast<VectorUblas<value_type>&>(*this),
                        ublas::range( nbdof_start, nbdof_start+space->nLocalDof() ) );

            Debug( 5010 ) << "Element <" << i << ">::range.size :  "<<  ct.size()<< "\n";
            Debug( 5010 ) << "Element <" << i << ">::range.start :  "<<  ct.start()<< "\n";
            return typename mpl::at_c<element_vector_type,i>::type( space, ct );
        }

        /**
         *
         * @return the finite element space associated with the n-th component
         */
        component_functionspace_ptrtype const& compSpace() const { return _M_functionspace->compSpace(); }


        /**
         * \return the number of dof
         */
        size_type nDof() const { return _M_functionspace->nDof(); }

        /**
         * \return the number of local dof (dof per processor)
         */
        size_type nLocalDof() const { return _M_functionspace->nLocalDof(); }

        /**
           \return the number of dof per component
        */
        size_type nDofPerComponent() const { return _M_functionspace->nDofPerComponent(); }

        /**
           \return the name of the element
        */
        std::string const& name() const { return _M_name; }

        size_type start() const { return _M_start; }

        bool isAComponent() const { return _M_ct >= X && _M_ct <= Z; }

        ComponentType component() const
        {
            if (  _M_ct < X || _M_ct > Z )
                {
                    std::ostringstream __str;
                    __str << "invalid component extraction (should be 0 or 1 or 2): " << _M_ct;
                    throw std::logic_error( __str.str() );
                }
            return _M_ct;
        }

        std::string componentToString( ComponentType ct ) const
        {
            if (  ct < X || ct > Z )
                {
                    std::ostringstream __str;
                    __str << "invalid component extraction (should be 0 or 1 or 2): " << ct;
                    throw std::logic_error( __str.str() );
                }
            switch(ct)
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
            _M_name = __name;
        }

        void setFunctionSpace( functionspace_ptrtype space )
        {
            _M_functionspace = space;
            super::init( _M_functionspace->nDof(),  _M_functionspace->nLocalDof() );
        }

        //@}
    private:

        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int /*version*/)
        {
            //ar & BOOST_SERIALIZATION_NVP( boost::serialization::base_object<super>(*this) );
            ar & boost::serialization::make_nvp( "name", _M_name );
            Debug( 5010 ) << "got name " << _M_name << "\n";
            if ( Archive::is_saving::value )
                {
                    size_type s = this->size();
                    ar & boost::serialization::make_nvp("size", s);

                    typename container_type::const_iterator it = this->begin();
                    typename container_type::const_iterator en = this->end();
                    for ( size_type i = 0; it != en; ++it, ++i )
                        {
                            T value = this->operator[](i);
                            std::ostringstream v_str;
                            v_str << "value_" << i;

                            //                     Debug( 5010 ) << "save value " << value << " at " << v_str.str() << "\n";

                            ar & boost::serialization::make_nvp( v_str.str().c_str(), value );
                        }
                }
            if ( Archive::is_loading::value )
                {
                    size_type s( 0 );
                    ar & boost::serialization::make_nvp("size", s);

                    Debug( 5010 ) << "loading ublas::vector of size " << s << "\n";

                    this->resize(  s );
                    Debug( 5010 ) << "resize done: " << this->size() << "\n";

                    for ( size_type i = 0; i < s ; ++i )
                        {
                            value_type value(  0 );
                            std::ostringstream v_str;
                            v_str << "value_" << i;
                            //                     Debug( 5010 ) << "load value at " << v_str.str() << "\n";
                            ar & boost::serialization::make_nvp( v_str.str().c_str(), value );
                            //                     Debug( 5010 ) << "got value " << value << " at index " << i << "\n";
                            this->operator[](i) = value;
                        }
                }


        }

    private:

        /**
           Finite Element space
        */
        functionspace_ptrtype _M_functionspace;

        std::string _M_name;

        size_type _M_start;

        //! type of the component
        ComponentType _M_ct;

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
                   periodicity_type  periodicity = periodicity_type() )
    {
        this->init( mesh, mesh_components, periodicity );
    }

    FunctionSpace( mesh_ptrtype const& mesh, std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices ) { this->init( mesh, 0, dofindices ); }

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

    static pointer_type New( mesh_ptrtype const& __m, std::vector<boost::tuple<size_type, uint16_type, size_type> > const dofindices )
    {
        return pointer_type( new functionspace_type( __m, dofindices ) );
    }

    static pointer_type New( mesh_ptrtype const& __m,
                             size_type mesh_components = MESH_RENUMBER | MESH_CHECK,
                             periodicity_type periodicity = periodicity_type() )
    {
        return pointer_type( new functionspace_type( __m, mesh_components, periodicity ) );
    }
    /**
     * initialize the function space
     */
    void init( mesh_ptrtype const& mesh, size_type mesh_components = MESH_RENUMBER | MESH_CHECK, periodicity_type periodicity = periodicity_type() )
    {

        Context ctx( mesh_components );
        Debug( 5010 ) << "component     MESH_RENUMBER: " <<  ctx.test( MESH_RENUMBER ) << "\n";
        Debug( 5010 ) << "component MESH_UPDATE_EDGES: " <<  ctx.test( MESH_UPDATE_EDGES ) << "\n";
        Debug( 5010 ) << "component MESH_UPDATE_FACES: " <<  ctx.test( MESH_UPDATE_FACES ) << "\n";
        Debug( 5010 ) << "component    MESH_PARTITION: " <<  ctx.test( MESH_PARTITION ) << "\n";
        this->init( mesh, mesh_components, periodicity, std::vector<boost::tuple<size_type, uint16_type, size_type> >(), mpl::bool_<is_composite>() );
        mesh->addObserver( *this );
    }

    void init( mesh_ptrtype const& mesh, size_type mesh_components, std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices )
    {

        Context ctx( mesh_components );
        Debug( 5010 ) << "component     MESH_RENUMBER: " <<  ctx.test( MESH_RENUMBER ) << "\n";
        Debug( 5010 ) << "component MESH_UPDATE_EDGES: " <<  ctx.test( MESH_UPDATE_EDGES ) << "\n";
        Debug( 5010 ) << "component MESH_UPDATE_FACES: " <<  ctx.test( MESH_UPDATE_FACES ) << "\n";
        Debug( 5010 ) << "component    MESH_PARTITION: " <<  ctx.test( MESH_PARTITION ) << "\n";
        this->init( mesh, mesh_components, periodicity_type(), dofindices, mpl::bool_<is_composite>() );
        mesh->addObserver( *this );
    }

    //! destructor: do nothing thanks to shared_ptr<>
    ~FunctionSpace() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * update element when mesh has been changed
     */
    void operator()( MESH_CHANGES mesh_changes )
    {
        Debug( 5010 ) << "Update function space after a change in the mesh\n";

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
        gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId()) );

        return __c_interp->xReal();
    }


    /**
     * return 1 if scalar space or the number of components if vectorial space
     */
    uint16_type qDim() const { return N_COMPONENTS; }

    /**
     * \return the number of degrees of freedom for each space overall the entire domain
     */
    size_type nDof() const { return this->nDof( mpl::bool_<is_composite>() ); }

    /**
     * \return the number of degrees of freedom for each space on the current subdomain
     */
    size_type nLocalDof() const { return this->nLocalDof( mpl::bool_<is_composite>() ); }

    /**
     * \return the distribution of the dofs among the processors
     */
    proc_dist_map_type getProcDistMap() const { return procDistMap; }

    /**
     * \return the starting value of the global dof numbering
     */
    size_type nDofStart( size_type i = invalid_size_type_value ) const
    {
        size_type start =  fusion::accumulate( this->functionSpaces(), size_type( 0 ), detail::NbDof( 0, i ) );
        return start;
    }

    /**
     * \return the number of degrees of freedom per dim
     */
    size_type nDofPerComponent() const { return this->nDof()/qDim(); }


    /**
       \return the mesh associated with the approximation
    */
    mesh_ptrtype& mesh() { return _M_mesh; }

    /**
       \return the mesh associated with the approximation
    */
    mesh_ptrtype const& mesh() const { return _M_mesh; }

    /**
       \return the reference finite element
    */
    basis_ptrtype const& basis() const { return _M_ref_fe; }

    /**
       \return the reference finite element
    */
    reference_element_ptrtype const& fe() const { return _M_ref_fe; }

    /**
       \return the geometric mapping
    */
    gm_ptrtype const& gm() const { return _M_mesh->gm(); }

    /**
       \return the geometric mapping
    */
    gm1_ptrtype const& gm1() const { return _M_mesh->gm1(); }

    /**
       \return the degrees of freedom
    */
    DataMap const& map() const { return *_M_dof; }


    /**
       \return the degrees of freedom
    */
    dof_ptrtype dof() { return _M_dof; }

    /**
       \return the degrees of freedom
    */
    dof_ptrtype const& dof() const { return _M_dof; }

    /**
     * get the \c FunctionSpace vector
     */
    //functionspace_vector_type const&
    //functionSpaces() const { return _M_functionspaces; }

    functionspace_vector_type const&
    functionSpaces() const { return _M_functionspaces; }

    std::vector<std::vector<int> > dofIndexSplit()
        {
            if ( nSpaces > 1 )
            {
                uint cptSpaces=0;
                uint start=0;
                std::vector<std::vector<int> > is;
                auto result = boost::make_tuple(cptSpaces,start,is);

                boost::fusion::fold( functionSpaces(), result,  detail::computeNDofForEachSpace() );

                return result.template get<2>();
            }
            std::vector<std::vector<int> > is;
            is.push_back( std::vector<int>( nLocalDof() ) );
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
     * get the \p i -th \c FunctionSpace out the list
     */
    template<int i>
    typename mpl::at_c<functionspace_vector_type,i>::type
    functionSpace() { return fusion::at_c<i>( _M_functionspaces ); }

    /**
     *
     * @return the finite element space associated with the n-th component
     */
    component_functionspace_ptrtype const& compSpace() const { return _M_comp_space; }

    /**
     * \return trace space
     */
    trace_functionspace_ptrtype
    trace()  const
        {
            //return trace( mpl::greater<mpl::int_<nDim>,mpl::int_<1> >::type() )
            return trace_functionspace_type::New( mesh()->trace(boundaryfaces(mesh())) );
        }
    template<typename RangeT>
    trace_functionspace_ptrtype
    trace( RangeT range  )  const
        {
            return trace_functionspace_type::New( mesh()->trace(range) );
        }


    /**
       \return true if Space has a region tree to localize points
    */
    bool hasRegionTree() const { return _M_rt; }

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
        Log() << " number of components : " << qDim() << "\n";
        Log() << "         n Global Dof : " << nDof() << "\n";
        Log() << "         n Local  Dof : " << nLocalDof() << "\n";
    }

    //@}


    FunctionSpace( FunctionSpace const& __fe )
        :
        _M_mesh( __fe._M_mesh ),
        _M_ref_fe( __fe._M_ref_fe ),
        _M_comp_space( __fe._M_comp_space ),
        _M_dof( __fe._M_dof ),
        _M_rt( __fe._M_rt )
    {
        Debug( 5010 ) << "copying FunctionSpace\n";
    }

protected:

private:

    void init( mesh_ptrtype const& mesh,
               size_type mesh_components,
               periodicity_type const& periodicity,
               std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices,
               mpl::bool_<false> );
    void init( mesh_ptrtype const& mesh,
               size_type mesh_components,
               periodicity_type const& periodicity,
               std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices,
               mpl::bool_<true> );


    size_type nDof( mpl::bool_<false> ) const;
    size_type nDof( mpl::bool_<true> ) const;

    size_type nLocalDof( mpl::bool_<false> ) const;
    size_type nLocalDof( mpl::bool_<true> ) const;

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
            _M_functionspace( __functionspace ),
            _M_mesh( __m )
        {}
        component_functionspace_ptrtype operator()()
        {
            return operator()( mpl::bool_<functionspace_type::is_scalar>() );
        }
        component_functionspace_ptrtype operator()( mpl::bool_<true> )
        {
            FEEL_ASSERT( 0 ).error( "invalid call for component space extraction" );
            //return _M_functionspace;
            return component_functionspace_ptrtype();
        }
        component_functionspace_ptrtype operator()( mpl::bool_<false> )
        {
            return component_functionspace_type::NewPtr( _M_mesh );
        }

    private:

        FunctionSpace<A0,A1,A2,A3,A4> * _M_functionspace;
        mesh_ptrtype _M_mesh;
    };

#if 0
    template<typename ElementRange, typename OnExpr>
    void
    on( ElementRange const& _range,
        OnExpr const& _expr )
    {
        M_constraints.push_back( vf::project( this->shared_from_this(), _range, _expr )  );
    }
#endif // 0

protected:

    //friend class FunctionSpace<mesh_type, typename bases_list::component_basis_type, value_type>;
    //friend class FunctionSpace<mesh_type, bases_list, value_type>;


    // finite element mesh
    mesh_ptrtype _M_mesh;

    //! finite element reference type
    reference_element_ptrtype _M_ref_fe;

    //! component fe space
    component_functionspace_ptrtype _M_comp_space;

    //! Degrees of freedom
    dof_ptrtype _M_dof;

    /** region tree associated with the mesh */
    mutable boost::optional<region_tree_ptrtype> _M_rt;

    //mutable boost::prof::basic_profiler<boost::prof::basic_profile_manager<std::string, real_type, boost::high_resolution_timer, boost::prof::empty_logging_policy, boost::prof::default_stats_policy<std::string, real_type> > > _M_prof_find_points;
private:

    //! disable default constructor
    FunctionSpace();

    functionspace_vector_type _M_functionspaces;

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
                                         std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices,
                                         mpl::bool_<false> )
{
    Debug( 5010 ) << "calling init(<space>) begin\n";
    Debug( 5010 ) << "calling init(<space>) is_periodic: " << is_periodic << "\n";
    _M_mesh = __m;



    if ( basis_type::nDofPerEdge || nDim >= 3 )
        mesh_components |= MESH_UPDATE_EDGES;
    /*
     * update faces info in mesh only if dofs exists on faces or the
     * expansion is continuous between elements. This case handles strong
     * Dirichlet imposition
     */
    if ( basis_type::nDofPerFace || is_continuous  || nDim >= 3 )
        mesh_components |= MESH_UPDATE_FACES;

    _M_mesh->components().set( mesh_components );

    _M_mesh->updateForUse();

    _M_ref_fe = basis_ptrtype( new basis_type );
    _M_dof = dof_ptrtype( new dof_type( _M_ref_fe, periodicity ) );
    Debug( 5010 ) << "[functionspace] Dof indices is empty ? " << dofindices.empty() << "\n";
    _M_dof->setDofIndices( dofindices );
    Debug( 5010 ) << "[functionspace] is_periodic = " << is_periodic << "\n";

    _M_dof->build( _M_mesh );

    if ( is_vectorial )
        _M_comp_space = component_functionspace_ptrtype( new component_functionspace_type( _M_mesh, MESH_COMPONENTS_DEFAULTS, periodicity ) );

    Debug( 5010 ) << "nb dim : " << qDim() << "\n";
    Debug( 5010 ) << "nb dof : " << nDof() << "\n";
    Debug( 5010 ) << "nb dof per component: " << nDofPerComponent() << "\n";

    if ( is_vectorial )
        {
            Debug( 5010 ) << "component space :: nb dim : " << _M_comp_space->qDim() << "\n";
            Debug( 5010 ) << "component space :: nb dof : " << _M_comp_space->nDof() << "\n";
            Debug( 5010 ) << "component space :: nb dof per component: " << _M_comp_space->nDofPerComponent() << "\n";
        }

    detail::searchIndicesBySpace<proc_dist_map_type>( this, procDistMap);

    Debug( 5010 ) << "calling init(<space>) end\n";

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::init( mesh_ptrtype const& __m,
                                         size_type mesh_components,
                                         periodicity_type const& periodicity,
                                         std::vector<boost::tuple<size_type, uint16_type, size_type> > const& dofindices,
                                         mpl::bool_<true> )
{
    Debug( 5010 ) << "calling init(<composite>) begin\n";
    _M_mesh = __m;

#if 0
    if ( basis_type::nDofPerEdge )
        _M_mesh->components().set( MESH_UPDATE_EDGES );
    /*
     * update faces info in mesh only if dofs exists on faces or the
     * expansion is continuous between elements. This case handles
     * strong Dirichlet imposition
     */
    if ( basis_type::nDofPerFace || is_continuous  )
        _M_mesh->components().set( MESH_UPDATE_FACES );
#else
    _M_mesh->components().set( mesh_components | MESH_UPDATE_FACES | MESH_UPDATE_EDGES | MESH_PARTITION );
#endif

    _M_mesh->updateForUse();

    fusion::for_each( _M_functionspaces, detail::InitializeSpace<mesh_type>( __m, dofindices ) );

    _M_dof = dof_ptrtype( new dof_type( this->nDof(), this->nLocalDof() ) );
    Debug( 5010 ) << "calling nDof(<composite>)" << this->nDof() << "\n";
    Debug( 5010 ) << "calling init(<composite>) end\n";

    proc_dist_map_type emptyMap;
    procDistMap = fusion::accumulate( _M_functionspaces,
                                      emptyMap,
                                      detail::searchIndicesBySpace<proc_dist_map_type>() );

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nDof( mpl::bool_<true> ) const
{
    Debug( 5010 ) << "calling nDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( _M_functionspaces, size_type( 0 ), detail::NbDof() );
    Debug( 5010 ) << "calling nDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nDof( mpl::bool_<false> ) const
{
    return _M_dof->nDof();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDof( mpl::bool_<true> ) const
{
    Debug( 5010 ) << "calling nLocalDof(<composite>) begin\n";
    size_type ndof =  fusion::accumulate( _M_functionspaces, size_type( 0 ), detail::NLocalDof() );
    Debug( 5010 ) << "calling nLocalDof(<composite>) end\n";
    return ndof;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
size_type
FunctionSpace<A0, A1, A2, A3, A4>::nLocalDof( mpl::bool_<false> ) const
{
    return _M_dof->nLocalDof();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
void
FunctionSpace<A0, A1, A2, A3, A4>::updateRegionTree() const
{
    scalar_type EPS=1E-13;

    region_tree_ptrtype __rt( new region_tree_type );

    __rt->clear();
    BoundingBox<> __bb( _M_mesh->gm()->isLinear() );

    typedef typename mesh_type::element_iterator mesh_element_iterator;
    mesh_element_iterator it = _M_mesh->beginElementWithProcessId( _M_mesh->comm().rank() );
    mesh_element_iterator en = _M_mesh->endElementWithProcessId( _M_mesh->comm().rank() );
    for ( size_type __i = 0; it != en; ++__i, ++it )
        {
            __bb.make( it->G() );
            for (unsigned k=0; k < __bb.min.size(); ++k)
                {
                    __bb.min[k]-=EPS;
                    __bb.max[k]+=EPS;
                }
            __rt->addBox( __bb.min, __bb.max, it->id() );
        }

    //__rt->dump();
    _M_rt = __rt;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
region_tree_ptrtype const&
FunctionSpace<A0, A1, A2, A3, A4>::regionTree() const
{
    if ( !_M_rt )
        {
            scalar_type EPS=1E-13;

            region_tree_ptrtype __rt( new region_tree_type );

            __rt->clear();
            BoundingBox<> __bb( _M_mesh->gm()->isLinear() );

            typedef typename mesh_type::element_iterator mesh_element_iterator;
            mesh_element_iterator it = _M_mesh->beginElementWithProcessId( _M_mesh->comm().rank() );
            mesh_element_iterator en = _M_mesh->endElementWithProcessId( _M_mesh->comm().rank() );
            for ( size_type __i = 0; it != en; ++__i, ++it )
                {
                    __bb.make( it->G() );
                    for (unsigned k=0; k < __bb.min.size(); ++k)
                        {
                            __bb.min[k]-=EPS;
                            __bb.max[k]+=EPS;
                        }
                    __rt->addBox( __bb.min, __bb.max, it->id() );
                }

            //__rt->dump();
            _M_rt = __rt;
        }
    return _M_rt.get();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
bool
FunctionSpace<A0, A1, A2, A3, A4>::findPoint(node_type const& pt,size_type &cv , node_type &ptr ) const
{
    if ( !hasRegionTree() )
        regionTree();

    //_M_prof_find_points.resume();

    region_tree_type* __rt = _M_rt.get().get();

    region_tree_type::pbox_set_type boxlst;
    __rt->findBoxesAtPoint(pt, boxlst);

    typedef typename gm_type::Inverse inv_trans_type;
    typename gm_type::reference_convex_type refelem;

    std::pair<size_type,value_type> closest = std::make_pair( invalid_size_type_value, -1e30 );

    region_tree_type::pbox_set_type::const_iterator it = boxlst.begin();
    region_tree_type::pbox_set_type::const_iterator ite = boxlst.end();
    for (; it != ite; ++it)
        {
            inv_trans_type __git( _M_mesh->gm(), _M_mesh->element( ( *it )->id ) );

            size_type cv_stored = (*it)->id;


            Debug( 5010 ) << "[FunctionSpace::findPoint] id : " << cv_stored << "\n";

            __git.setXReal( pt );
            ptr = __git.xRef();


            bool isin;
            value_type dmin;
            boost::tie( isin, dmin ) = refelem.isIn( ptr );
            Debug( 5010 ) << "[FunctionSpace::findPoint] isin: " << isin << " dmin: " << dmin << "\n";

            closest =  ( dmin > closest.second )?std::make_pair( cv_stored, dmin ):closest;
            if ( isin )
                {
                    Debug( 5010 ) << "[FunctionSpace::findPoint] id of the convex where " << pt << " belongs : " << cv_stored << "\n";
                    Debug( 5010 ) << "[FunctionSpace::findPoint] ref coordinate: " << ptr << "\n";
                    cv = (*it)->id;
                    //_M_prof_find_points.pause();
                    return true;
                }
        }
    cv=closest.first;
    //_M_prof_find_points.pause();
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
    _M_start( 0 ),
    _M_ct( NO_COMPONENT )
{}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( Element const& __e )
    :
    super( __e ),
    _M_functionspace( __e._M_functionspace ),
    _M_name( __e._M_name ),
    _M_start( __e._M_start ),
    _M_ct( __e._M_ct )
{
    Debug( 5010 ) << "Element<copy>::range::start = " << this->start() << "\n";
    Debug( 5010 ) << "Element<copy>::range::size = " << this->size() << "\n";

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::Element( functionspace_ptrtype const& __functionspace,
                                                             std::string const& __name,
                                                             size_type __start,
                                                             ComponentType __ct )
    :
    super( __functionspace->map() ),
    _M_functionspace( __functionspace ),
    _M_name( __name ),
    _M_start( __start ),
    _M_ct( __ct )
{
    Debug( 5010 ) << "Element::start = " << this->start() << "\n";
    Debug( 5010 ) << "Element::size = " << this->size() << "\n";
    Debug( 5010 ) << "Element::ndof = " << this->nDof() << "\n";
    Debug( 5010 ) << "Element::nlocaldof = " << this->nLocalDof() << "\n";
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
    _M_functionspace( __functionspace ),
    _M_name( __name ),
    _M_start( __start ),
    _M_ct( __ct )
{
    Debug( 5010 ) << "Element<range>::range::start = " << __c.start() << "\n";
    Debug( 5010 ) << "Element<range>::range::size = " << __c.size() << "\n";
    Debug( 5010 ) << "Element<range>::start = " << this->start() << "\n";
    Debug( 5010 ) << "Element<range>::size = " << this->size() << "\n";
    Debug( 5010 ) << "Element<range>::ndof = " << this->nDof() << "\n";
    Debug( 5010 ) << "Element<range>::nlocaldof = " << this->nLocalDof() << "\n";
    _M_start = __c.start();
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::~Element()
{}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::init( functionspace_ptrtype const& __functionspace,
                                                          container_type const& __c )
{
    _M_functionspace = __functionspace;
    (container_type)*this = __c;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>&
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::operator=( Element const& __e )
{
    if (  this != &__e )
        {
            _M_functionspace = __e._M_functionspace;
            if ( __e._M_name != "unknown" )
                _M_name = __e._M_name;
            _M_start = __e._M_start;
            _M_ct = __e._M_ct;
            this->resize( _M_functionspace->nLocalDof() );
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
        Debug( 5010 ) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";
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

#if defined(HAVE_MPI)
        if ( nprocs > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, detail::vector_plus<int>() );
        }
#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* HAVE_MPI */

        id_type __id( this->id( *fectx ) );

        //Debug(5010) << "[interpolation]  id = " << __id << "\n";
#if defined(HAVE_MPI)
        Debug( 5010 ) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";
        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), __id, functionSpace()->mesh()->comm().rank() );
        }
        //Debug(5010) << "[interpolation] after broadcast id = " << __id << "\n";
#endif /* HAVE_MPI */
        return __id;
    }
    else
    {
#if defined(HAVE_MPI)
        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, detail::vector_plus<int>() );
        }
#endif /* HAVE_MPI */
        bool found = false;
        size_type i = 0;
        for(;i < global_found_pt.size();++i )
            if ( global_found_pt[i] != 0 )
            {
                Debug(5010) << "processor " << i << " has the point " << __x << "\n";
                found = true;
                break;
            }
        id_type __id;
        if ( found )
        {
            Debug( 5010 ) << "receiving interpolation context from processor " << i << "\n";
#if defined(HAVE_MPI)
            if ( functionSpace()->mesh()->comm().size() > 1 )
                mpi::broadcast( functionSpace()->mesh()->comm(), __id, i );
#endif /* HAVE_MPI */

            Debug(5010) << "[interpolation] after broadcast id = " << __id << "\n";
        }
        else
        {
            Warning() << "no processor seems to have the point " << __x << "\n";
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

    const uint16_type nq = context.xRefs().size2();
    //array_type v( boost::extents[nComponents1][nComponents2][context.xRefs().size2()] );
        for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for( typename array_type::index c1 = 0; c1 < ncdof; ++c1 )
            {
                typename array_type::index ldof = basis_type::nDof*c1+l;
                size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                //std::cout << "ldof = " << ldof << "\n";
                //std::cout << "gdof = " << gdof << "\n";

                FEEL_ASSERT( gdof < this->size() )
                    ( context.eId() )
                    ( l )( c1 )( ldof)( gdof )
                    ( this->size() )( this->localSize() )
                    ( this->firstLocalIndex() )( this->lastLocalIndex() )
                    .error( "FunctionSpace::Element invalid access index" );

                value_type v_ = this->globalValue( gdof );
                //std::cout << "v_ =" << v_ << "\n";
                //for( typename array_type::index c2 = 0; c2 < nComponents2; ++c2 )
                for( uint16_type q = 0; q < nq; ++q )
                {
                    for( typename array_type::index i = 0; i < nComponents1; ++i )
                        //for( typename array_type::index j = 0; j < nComponents2; ++j )
                    {
                        v[q](i,0) += v_*context.id( ldof, i, 0, q );
                        //v[q](i,0) += v_*context.gmc()->J(*)*context.pc()->phi( ldof, i, 0, q );
                    }
                }
            }
        }
    //return v;
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );

            //update precompute of basis functions
            __pc->update( pts );

            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            id_type __id( this->id( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
            {
                for( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[boost::get<0>(*itL)](i,0) =  __id(i,0,k);
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
        Debug( 5010 ) << "Point " << __x << " is in element " << __cv_id << " pt_ref=" << __x_ref << "\n";
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

#if defined(HAVE_MPI)
        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, detail::vector_plus<int>() );
        }
#else
        global_found_pt[ 0 ] = found_pt[ 0 ];
#endif /* HAVE_MPI */

        grad_type g_( this->grad( *fectx ) );
        //Debug(5010) << "[interpolation]  id = " << v << "\n";
#if defined(HAVE_MPI)
        Debug( 5010 ) << "sending interpolation context to all processors from " << functionSpace()->mesh()->comm().rank() << "\n";
        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            mpi::broadcast( functionSpace()->mesh()->comm(), g_, functionSpace()->mesh()->comm().rank() );
        }
        //Debug(5010) << "[interpolation] after broadcast g_ = " << g_ << "\n";
#endif /* HAVE_MPI */
        return g_;
    }
    else
    {
#if defined(HAVE_MPI)
        if ( functionSpace()->mesh()->comm().size() > 1 )
        {
            //mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, std::plus<std::vector<int> >() );
            mpi::all_reduce( functionSpace()->mesh()->comm(), found_pt, global_found_pt, detail::vector_plus<int>() );
        }
#endif /* HAVE_MPI */
        bool found = false;
        size_type i = 0;
        for(;i < global_found_pt.size();++i )
            if ( global_found_pt[i] != 0 )
            {
                Debug(5010) << "processor " << i << " has the point " << __x << "\n";
                found = true;
                break;
            }
        grad_type g_;
        if ( found )
        {
            Debug( 5010 ) << "receiving interpolation context from processor " << i << "\n";
#if defined(HAVE_MPI)
            if ( functionSpace()->mesh()->comm().size() > 1 )
                mpi::broadcast( functionSpace()->mesh()->comm(), g_, i );
#endif /* HAVE_MPI */

            //Debug(5010) << "[interpolation] after broadcast id = " << v << "\n";
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

    //std::cout << "coeff=" << coeff << "\n";
    //array_type v( boost::extents[nComponents1][nRealDim][context.xRefs().size2()] );
    //std::fill( v.data(), v.data()+v.num_elements(), value_type( 0 ) );
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    int ldof = c1*basis_type::nDof+l;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( l )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );

                    //value_type v_ = (*this)( gdof );
                    value_type v_ = this->globalValue( gdof );
                    for( size_type q = 0; q < context.xRefs().size2(); ++q )
                    {
                        for( int k = 0; k < nComponents1; ++k )
                            for( int j = 0; j < nRealDim; ++j )
                            {
                                v[q](k,j) += v_*context.grad( ldof, k, j, q );
                            }
                    }
                }
        }
#if 0
    std::cout << "xref = " << context.xRefs() << "\n";
    std::cout << "xreal = " << context.xReal() << "\n";
    std::cout << "pts = " << context.G() << "\n";
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            for( typename array_type::index i = 0; i < nDim; ++i )
                {
                    for( size_type q = 0; q < context.xRefs().size2(); ++q )
                        std::cout << "Gref[ " << l << "," << i << "," << q << "] = " << pc.grad( l, 0, i, q ) << "\n";

                }
        }
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            for( int j = 0; j < nRealDim; ++j )
                for( typename array_type::index i = 0; i < nDim; ++i )
                    {
                        for( size_type q = 0; q < context.xRefs().size2(); ++q )
                            std::cout << "G[ " << l << "," << j << "," << i << "," << q << "] = " << context.B(q)( j, i )*pc.grad( l, 0, i, q ) << "\n";

                    }
        }
    for( int k = 0; k < nComponents1; ++k )
        for( int j = 0; j < nRealDim; ++j )
            {
                for( size_type q = 0; q < context.xRefs().size2(); ++q )
                    {
                        std::cout << "v[ " << k << "," << j << "," << q << "]= " << v[q](k,j)  << "\n";
                        std::cout << "B(" << q << ")=" << context.B(q) << "\n";
                        std::cout << "J(" << q << ")=" << context.J(q) << "\n";
                        std::cout << "K(" << q << ")=" << context.K(q) << "\n";
                    }
            }
#endif
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::gradInterpolate(  matrix_node_type __ptsReal, grad_array_type& v ) const
{
    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );


            //update precompute of basis functions
            __pc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            grad_type __grad( this->grad( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
            {
                for( typename array_type::index i = 0; i < nComponents1; ++i )
                    for( uint j = 0; j < nRealDim; ++j )
                    {
                        v[boost::get<0>(*itL)](i,j) = __grad(i,j,k);
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
    const size_type Q = context.xRefs().size2();
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    int ldof = c1*basis_type::nDof+l;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( l )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );
                    //value_type v_ = (*this)( gdof );
                    value_type v_ = this->globalValue( gdof );
                        for( size_type q = 0; q < Q; ++q )
                            {
#if 0
                                std::cout << "v(" << gdof << ")=" << v_ << "\n";
                                std::cout << "context.div(" << ldof << "," << q << ")="
                                          << context.div( ldof, 0, 0, q ) << "\n" ;
#endif
                                v[q](0,0) += v_*context.div( ldof, 0, 0, q );
                            }
                }
        }
#else
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    int ldof = c1*basis_type::nDof+l;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( l )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );
                    //value_type v_ = (*this)( gdof );
                    value_type v_ = this->globalValue( gdof );
                    for( int k = 0; k < nComponents1; ++k )
                        {
                            for( typename array_type::index i = 0; i < nDim; ++i )
                                {
                                    for( size_type q = 0; q < context.xRefs().size2(); ++q )
                                        {
                                            v[q](0,0) += v_*context.gmContext()->B(q)( k, i )*context.pc()->grad( ldof, k, i, q );
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
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::divInterpolate( matrix_node_type __ptsReal, div_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            div_type __div( this->div( *__ctx ));

            //update the output data
            //for( typename array_type::index i = 0; i < nComponents1; ++i )
            //   {
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
                v[boost::get<0>(*itL)](0,0) =  __div(0,0,k);
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
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    int ldof = c1*basis_type::nDof+l;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( l )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );

                    //value_type v_ = (*this)( gdof );
                    value_type v_ = this->globalValue( gdof );
                    const uint16_type nq = context.xRefs().size2();
                    for( uint16_type q = 0; q < nq ; ++q )
                    {
                        if ( nDim == 3 )
                        {
                            for( typename array_type::index i = 0; i < nDim; ++i )
                            {
                                v[q](i,0) += v_*context.curl( ldof, i, 0, q );
                            }
                        }
                        else if ( nDim == 2 )
                        {
                            v[q](0,0) += v_*context.curl( ldof, 0, 0, q );
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
    for( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    int ldof = c1*basis_type::nDof+l;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), l, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( l )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );

                    //value_type v_ = (*this)( gdof );
                    value_type v_ = this->globalValue( gdof );
                    const uint16_type nq = context.xRefs().size2();
                    for( uint16_type q = 0; q < nq ; ++q )
                    {
                        if ( nDim == 3 )
                        {
                            v[q](0,0) += v_*context.curl( ldof, comp, 0, q );
                        }
                        else if ( nDim == 2 )
                        {
                            v[q](0,0) += v_*context.curl( ldof, 2, 0, q );
                        }

                    }
                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlInterpolate( matrix_node_type __ptsReal, curl_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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


    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            curl_type __curl( this->curl( *__ctx ));

            //update the output data
            itL=it->second.begin();
            if ( nDim == 3 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        {
                            v[boost::get<0>(*itL)](0,0) =  __curl(0,0,k);
                            v[boost::get<0>(*itL)](1,0) =  __curl(1,0,k);
                            v[boost::get<0>(*itL)](2,0) =  __curl(2,0,k);
                        }
                }
            else if ( nDim == 2 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        v[boost::get<0>(*itL)](0,0) =  __curl(0,0,k);
                }
        }

}


template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlxInterpolate(matrix_node_type __ptsReal, comp_curl_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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


    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            curlx_type __curlx( this->curlx( *__ctx ));

            //update the output data
            itL=it->second.begin();
            if ( nDim == 3 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        {
                            v[boost::get<0>(*itL)](0,0) =  __curlx(0,0,k);
                            v[boost::get<0>(*itL)](1,0) =  __curlx(1,0,k);
                            v[boost::get<0>(*itL)](2,0) =  __curlx(2,0,k);
                        }
                }
            else if ( nDim == 2 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        v[boost::get<0>(*itL)](0,0) =  __curlx(0,0,k);
                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlyInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            curly_type __curly( this->curly( *__ctx ));

            //update the output data
            itL=it->second.begin();
            if ( nDim == 3 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        {
                            v[boost::get<0>(*itL)](0,0) =  __curly(0,0,k);
                            v[boost::get<0>(*itL)](1,0) =  __curly(1,0,k);
                            v[boost::get<0>(*itL)](2,0) =  __curly(2,0,k);
                        }
                }
            else if ( nDim == 2 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        v[boost::get<0>(*itL)](0,0) =  __curly(0,0,k);
                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::curlzInterpolate( matrix_node_type __ptsReal, comp_curl_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            curlz_type __curlz( this->curlz( *__ctx ));

            //update the output data
            itL=it->second.begin();
            if ( nDim == 3 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        {
                            v[boost::get<0>(*itL)](0,0) =  __curlz(0,0,k);
                            v[boost::get<0>(*itL)](1,0) =  __curlz(1,0,k);
                            v[boost::get<0>(*itL)](2,0) =  __curlz(2,0,k);
                        }
                }
            else if ( nDim == 2 )
                {
                    itL=it->second.begin();
                    for (uint k=0;k<nbPtsElt;++k,++itL)
                        v[boost::get<0>(*itL)](0,0) =  __curlz(0,0,k);
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
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    size_type ldof = basis_type::nDof*c1 + i;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( i )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );

                    value_type v_ = this->globalValue( gdof );
                    for( size_type q = 0; q < context.xRefs().size2(); ++q )
                    {
                        for( typename array_type::index i = 0; i < nComponents1; ++i )
                        {
                            v[q](i,0) += v_*context.d( ldof, i, N, q );
                        }
                    }

                }
        }
}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dxInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            dx_type __dx( this->dx( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
                for( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[boost::get<0>(*itL)](i,0) =  __dx(i,0,k);
                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dyInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            dy_type __dy( this->dy( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
                for( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[boost::get<0>(*itL)](i,0) =  __dy(i,0,k);
                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::dzInterpolate( matrix_node_type __ptsReal, id_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );
            //evaluate element for these points
            dz_type __dz( this->dz( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
                for( typename array_type::index i = 0; i < nComponents1; ++i )
                {
                    v[boost::get<0>(*itL)](i,0) =  __dz(i,0,k);
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
    //FEEL_ASSERT( comp < nRealDim )( comp )( nRealDim ).error( "[FunctionSpace::Element] grad: invalid component" );
    for ( int i = 0; i < basis_type::nDof; ++i )
        {
            const int ncdof = is_product?nComponents1:1;
            for(int c1 = 0; c1 < ncdof; ++c1 )
                {
                    size_type ldof = basis_type::nDof*c1 + i;
                    size_type gdof = boost::get<0>(_M_functionspace->dof()->localToGlobal( context.eId(), i, c1 ) );
                    FEEL_ASSERT( gdof >= this->firstLocalIndex() &&
                                 gdof < this->lastLocalIndex() )
                        ( context.eId() )
                        ( i )( c1 )( ldof)( gdof )
                        ( this->size() )( this->localSize() )
                        ( this->firstLocalIndex() )( this->lastLocalIndex() )
                        .error( "FunctionSpace::Element invalid access index" );

                    value_type v_ = this->globalValue( gdof );

                    for( size_type q = 0; q < context.xRefs().size2(); ++q )
                    {
                        for( int i = 0; i < nRealDim; ++i )
                            for( int j = 0; j < nRealDim; ++j )
                            {
                                v[q](i,j) += v_*context.hess( ldof, i, j, q );
                            } // i,j
                    } // q

                }
        }

}

template<typename A0, typename A1, typename A2, typename A3, typename A4>
template<typename Y,  typename Cont>
void
FunctionSpace<A0, A1, A2, A3, A4>::Element<Y,Cont>::hessInterpolate( matrix_node_type __ptsReal, hess_array_type& v ) const
{

    typedef typename mesh_type::Localization::localization_ptrtype localization_ptrtype;
    typedef typename mesh_type::Localization::container_search_iterator_type analysis_iterator_type;
    typedef typename mesh_type::Localization::container_output_iterator_type analysis_output_iterator_type;

    // create analysys map : id -> List of pt
    localization_ptrtype __loc = this->functionSpace()->mesh()->tool_localization();
    __loc->run_analysis(__ptsReal);
    analysis_iterator_type it = __loc->result_analysis_begin();
    analysis_iterator_type it_end = __loc->result_analysis_end();
    analysis_output_iterator_type itL,itL_end;

    //geomap
    gm_ptrtype __gm = this->functionSpace()->gm();

    //if analysis map is empty : no interpolation
    if (it==it_end) return;

    //alocate a point matrix
    size_type nbPtsElt=it->second.size();
    uint nbCoord=boost::get<1>(*(it->second.begin())).size();
    matrix_node_type pts( nbCoord, nbPtsElt );

    //init the geomap context and precompute basis function
    geopc_ptrtype __geopc( new geopc_type( __gm, pts ) );
    pc_ptrtype __pc( new pc_type( this->functionSpace()->fe(), pts ));
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

    for ( ;it!=it_end;++it)
        {
            nbPtsElt = it->second.size();

            //iterate in the list pt for a element
            itL=it->second.begin();
            itL_end=it->second.end();

            //compute a point matrix with the list of point
            pts= matrix_node_type( nbCoord, nbPtsElt );
            for (size_type i=0;i<nbPtsElt;++i,++itL)
                ublas::column( pts, i ) = boost::get<1>(*itL);

            //update geomap context
            __geopc->update( pts );
            __c->update( this->functionSpace()->mesh()->element( it->first ), __geopc );
            //update precompute of basis functions
            __pc->update( pts );
            __ctx->update( __c, __pc );

            //evaluate element for these points
            hess_type __hess( this->hess( *__ctx ));

            //update the output data
            itL=it->second.begin();
            for (uint k=0;k<nbPtsElt;++k,++itL)
                for( int i = 0; i < nRealDim; ++i )
                    for( int j = 0; j < nRealDim; ++j )
                        v[boost::get<0>(*itL)](i,j) =  __hess(i,j,k);
        }

}

template<typename T,int M,int N>
std::ostream&
operator<<( std::ostream& os, detail::ID<T,M,N> const& id )
{
    const size_type* shape =  id.M_id.shape();
    for( size_type i = 0;i < shape[0]; ++i )
        {
            for( size_type j = 0;j < shape[1]; ++j )
                {
                    for( size_type k = 0;k < shape[2]; ++k )
                        {
                            os << id( i, j, k ) << ' ';
                        }
                    os << std::endl;
                }
            os << std::endl;
        }
    return os;
}
template<typename T,int M, int N>
inline
DebugStream&
operator<<( DebugStream& __os, detail::ID<T,M,N> const& id )
{
    if ( __os.doPrint() )
        {
            std::ostringstream __str;

            __str << id;

            __os << __str.str() << "\n";
        }
    return __os;
}
template<typename T, int M, int N>
inline
NdebugStream&
operator<<( NdebugStream& os, detail::ID<T,M,N> const& )
{
    return os;
}

/**
 * Computes the element wise product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the element wise product of \p v1 and \p v2
 */
template <typename ElementType>
ElementType
element_product( ElementType const& v1, ElementType const& v2 )
{
    FEEL_ASSERT( v1.functionSpace() == v2.functionSpace() ).error( "incompatible function spaces");

    typedef typename type_traits<typename ElementType::value_type>::real_type real_type;

    ElementType _t( v1.functionSpace() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();
    for( size_type i = 0; i < s; ++i )
        _t.operator()(start+i) = v1.operator()( start + i )* v2.operator()( start + i );
    return _t;
}

/**
 * Computes the element wise division of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the element wise division product of \p v1 and \p v2
 */
template <typename ElementType>
ElementType
element_div( ElementType const& v1, ElementType const& v2 )
{
    FEEL_ASSERT( v1.functionSpace() == v2.functionSpace() ).error( "incompatible function spaces");

    typedef typename type_traits<typename ElementType::value_type>::real_type real_type;

    ElementType _t( v1.functionSpace() );
    size_type s = v1.localSize();
    size_type start = v1.firstLocalIndex();
    for( size_type i = 0; i < s; ++i )
        _t.operator()(start+i) = v1.operator()( start + i )/v2.operator()( start + i );
    return _t;
}

template<typename MeshType>
auto
measurePointElements( boost::shared_ptr<FunctionSpace<MeshType,bases<Lagrange<MeshType::nOrder,Scalar> > > >& _space ) -> decltype( _space->element() )
{
    auto _fn = _space->element( "measurePointElements" );
    _fn.setZero();
    std::vector<bool> ptdone( _space->mesh()->numPoints(), false );
    auto elit = _space->mesh()->beginElement();
    auto elen = _space->mesh()->endElement();
    for( ; elit != elen; ++ elit )
    {
        for( int p = 0; p < elit->numPoints; ++p )
        {
            if ( ptdone[elit->point(p).id()] == false )
            {
                BOOST_FOREACH( auto pt, elit->point(p).elements() )
                {
                    _fn.plus_assign( elit->id(), p, 0, elit->measure() );
                }
                ptdone[elit->point(p).id()] = true;
            }
        }
    }
    return _fn;
}

} // Feel
#endif /* __FunctionSpace_H */
