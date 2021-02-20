/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2013-04-07

   Copyright (C) 2013-2016 Feel++ Consortium

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
   \file ReducedBasisSpace.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-07
*/
#ifndef __ReducedBasisSpace_H
#define __ReducedBasisSpace_H 1

#include <feel/feelcrb/crb.hpp>


namespace Feel
{

namespace detail
{
template<typename RbSpaceType>
struct InitializeRbSubSpace
{
    InitializeRbSubSpace( RbSpaceType const& rbSpaceType )
        :
        M_rbSpaceType( rbSpaceType )
    {}

    template <typename T>
    void operator()( T & x ) const
    {
        typedef typename T::first_type key_type;
        typedef typename T::second_type::element_type rbsubspace_type;
        if ( x.second && x.first == key_type() )
        {
            auto rbSubSpace = x.second;
            rbSubSpace->setFunctionSpace( M_rbSpaceType.template subFeFunctionSpace<key_type::value>() );
        }
        else
        {
            auto rbSubSpace = std::make_shared<rbsubspace_type>( M_rbSpaceType.worldCommPtr() );
            rbSubSpace->setFunctionSpace( M_rbSpaceType.template subFeFunctionSpace<key_type::value>() );
            x = std::make_pair(key_type(), rbSubSpace);
        }
    }
private :
    RbSpaceType const& M_rbSpaceType;
};

/**
 * add basis in subRbSpace
 * Type=0 : primal basis
 * Type=1 : dual basis
 */
template<typename RbSpaceType,int Type>
struct AddBasisElementRbSubSpace
{
    AddBasisElementRbSubSpace( RbSpaceType const& rbSpaceType, typename RbSpaceType::space_element_type const & e )
        :
        M_rbSpaceType( rbSpaceType ),
        M_e( e )
        {}

    template <typename T>
    void operator()( T & x ) const
        {
            this->operator()( x,mpl::int_<Type>() );
        }
    template <typename T>
    void operator()( T & x, mpl::int_<0> ) const
        {
            typedef typename T::first_type key_type;
            // create new element
            auto basiselt = M_rbSpaceType.template subFeFunctionSpace<key_type::value>()->element();
            // copy view in element
            basiselt = M_e.template element<key_type::value>();
            // add basis in subspace
            x.second->addPrimalBasisElement( basiselt );
        }
    template <typename T>
    void operator()( T & x, mpl::int_<1> ) const
        {
            typedef typename T::first_type key_type;
            // create new element
            auto basiselt = M_rbSpaceType.template subFeFunctionSpace<key_type::value>()->element();
            // copy view in element
            basiselt = M_e.template element<key_type::value>();
            // add basis in subspace
            x.second->addDualBasisElement( basiselt );
        }
private :
    RbSpaceType const& M_rbSpaceType;
    typename RbSpaceType::space_element_type const & M_e;
};

/**
 * delete n-last basis in subRbSpace
 * Type=0 : primal basis
 * Type=1 : dual basis
 */
template<typename RbSpaceType,int Type>
struct DeleteLastBasisElementsRbSubSpace
{
    DeleteLastBasisElementsRbSubSpace( RbSpaceType const& rbSpaceType, int number )
        :
        M_rbSpaceType( rbSpaceType ),
        M_number( number )
        {}

    template <typename T>
    void operator()( T & x ) const
        {
            this->operator()( x,mpl::int_<Type>() );
        }
    template <typename T>
    void operator()( T & x, mpl::int_<0> ) const
        {
            x.second->deleteLastPrimalBasisElements( M_number );
        }
    template <typename T>
    void operator()( T & x, mpl::int_<1> ) const
        {
            x.second->deleteLastDualBasisElements( M_number );
        }
private :
    RbSpaceType const& M_rbSpaceType;
    int M_number;
};

/**
 * delete all basis in subRbSpace
 * Type=0 : primal basis
 * Type=1 : dual basis
 */
template<int Type>
struct ClearBasisElementsRbSubSpace
{
    template <typename T>
    void operator()( T & x ) const
        {
            this->operator()( x,mpl::int_<Type>() );
        }
    template <typename T>
    void operator()( T & x, mpl::int_<0> ) const
        {
            x.second->primalRB().clear();
        }
    template <typename T>
    void operator()( T & x, mpl::int_<1> ) const
        {
            x.second->dualRB().clear();
        }
};


struct SetDimensionRbSubSpace
{
    SetDimensionRbSubSpace( int dim )
        :
        M_dim( dim )
        {}

    template <typename T>
    void operator()( T & x ) const
        {
            x.second->primalRB().resize( M_dim );
            x.second->dualRB().resize( M_dim );
        }

    int M_dim;
};

}

/**
 * \class ReducedBasisSpace
 * \brief Reduced Basis Space class
 *
 * @author Christophe Prud'homme
 * @see
 */

template<typename SpaceType>
class ReducedBasisSpace : public SpaceType
{

    typedef SpaceType super;
    typedef std::shared_ptr<super> super_ptrtype;

public :

    typedef double value_type;

    typedef typename super::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef super fespace_type;
    typedef super_ptrtype fespace_ptrtype;
    typedef super functionspace_type;
    typedef super_ptrtype functionspace_ptrtype;

    //typedef typename super::element_type space_element_type;
    typedef typename fespace_type::element_type space_element_type;
    typedef std::shared_ptr<space_element_type> space_element_ptrtype;

    typedef std::vector< space_element_ptrtype > rb_basis_type;
    typedef std::shared_ptr<rb_basis_type> rb_basis_ptrtype;

    typedef ReducedBasisSpace<SpaceType> this_type;
    typedef std::shared_ptr<this_type> this_ptrtype;

    typedef Eigen::Matrix<value_type,Eigen::Dynamic,1> eigen_vector_type;
    typedef Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;

    typedef typename functionspace_type::bases_list bases_list;
    typedef typename functionspace_type::meshes_list meshes_list;


    //static const bool is_composite = ( mpl::size<bases_list>::type::value > 1 );
    static const bool is_composite = super::is_composite;
    static const bool is_product = super::is_product ;
    static const uint16_type nSpaces = super::nSpaces;


    struct nodim { static const int nDim = -1; static const int nRealDim = -1; };
    static const uint16_type nDim = mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                                             mpl::identity<meshes_list >,
                                             mpl::identity<nodim> >::type::type::nDim;
    static const uint16_type nRealDim = mpl::if_<boost::is_base_of<MeshBase<>, meshes_list >,
                                                 mpl::identity<meshes_list>,
                                                 mpl::identity<nodim> >::type::type::nRealDim;


    // geomap
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm_type>, mpl::identity<mpl::void_> >::type::type gm_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::gm1_type>, mpl::identity<mpl::void_> >::type::type gm1_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename mesh_type::element_type>, mpl::identity<mpl::void_> >::type::type geoelement_type;
    typedef std::shared_ptr<gm_type> gm_ptrtype;
    typedef std::shared_ptr<gm1_type> gm1_ptrtype;
    static const size_type pts_gmc_context_v = vm::POINT;
    static const size_type gmc_context_v = vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context</*vm::POINT,*/ geoelement_type> >,
                              mpl::identity<mpl::void_> >::type::type pts_gmc_type;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::template Context</*vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB,*/ geoelement_type> >,
                              mpl::identity<mpl::void_> >::type::type gmc_type;
    typedef std::shared_ptr<gmc_type> gmc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_ptrtype>, mpl::identity<mpl::void_> >::type::type geopc_ptrtype;
    typedef typename mpl::if_<mpl::greater<mpl::int_<nDim>, mpl::int_<0> >,mpl::identity<typename gm_type::precompute_type>, mpl::identity<mpl::void_> >::type::type geopc_type;


    static const uint16_type nComponents = super::nComponents;

#if 0
    template<typename MeshListType,int N>
    struct GetMesh
    {
        typedef typename mpl::if_<mpl::or_<boost::is_base_of<MeshBase, MeshListType >,
                                           is_shared_ptr<MeshListType> >,
                                  mpl::identity<mpl::identity<MeshListType> >,
                                  mpl::identity<mpl::at_c<MeshListType,N> > >::type::type::type type;
    };



    typedef typename mpl::at_c<bases_list,0>::type::template apply<GetMesh<meshes_list,0>::type::nDim,
                                                                   GetMesh<meshes_list,0>::type::nRealDim,
                                                                   value_type,
                                                                   typename GetMesh<meshes_list,0>::type::element_type>::type basis_0_type;


    typedef typename mpl::if_<mpl::bool_<is_composite>,
                              mpl::identity<bases_list>,
                              mpl::identity<basis_0_type> >::type::type basis_type;


    typedef typename super::periodicity_type periodicity_type;
    typedef std::shared_ptr<periodicity_type> periodicity_ptrtype;
#else
    typedef typename super::basis_0_type basis_0_type;
    typedef typename super::basis_type basis_type;
#endif

    template<int i>
    struct sub_rbfunctionspace
    {
        typedef ReducedBasisSpace<typename fespace_type::template sub_functionspace<i>::type > type;
        typedef std::shared_ptr<type> ptrtype;
    };

    template<typename keyType>
    struct ChangeRbSpace
    {
        typedef std::pair< keyType, typename sub_rbfunctionspace<keyType::value>::ptrtype > the_type;
        typedef typename mpl::if_<mpl::bool_<fespace_type::is_composite>,
                                  mpl::identity<the_type>,
                                  mpl::identity<boost::none_t> >::type::type type;
    };
    typedef mpl::range_c<int,0, fespace_type::nSpaces> rangeRbSpaceType;
    typedef typename mpl::transform< rangeRbSpaceType, ChangeRbSpace<mpl::_1>,
                                    mpl::back_inserter<fusion::vector<> > >::type rbfunctionspace_vector_type;


#if 0
    ReducedBasisSpace( fespace_ptrtype const& Xh )
        :
        super( Xh ),
        M_mesh( Xh->mesh() )
        {
            this->init();
        }
#else
    explicit ReducedBasisSpace( worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr() )
        :
        super( worldcomm )
        {
            this->init( mpl::bool_<fespace_type::is_composite>() );
        }

    ReducedBasisSpace( fespace_ptrtype const& Xh, worldcomm_ptr_t const& worldcomm = Environment::worldCommPtr() )
        :
        super( worldcomm ),
        M_mesh()
        {
            if ( Xh )
            {
                this->setFunctionSpace( Xh );
            }
            else
                this->init( mpl::bool_<fespace_type::is_composite>() );
        }
    template< typename TheSpaceType>
    void setFunctionSpace( TheSpaceType const& feSpace )
        {
            this->setFunctionSpace( feSpace, mpl::bool_< boost::is_same<TheSpaceType,functionspace_ptrtype>::value> () );
        }
    template< typename TheSpaceType>
    void setFunctionSpace( TheSpaceType/*functionspace_ptrtype*/ const& feSpace, mpl::false_ )
        {}
    void setFunctionSpace( functionspace_ptrtype const& feSpace, mpl::true_ )
        {
            if ( feSpace )
            {
                super::shallowCopy( feSpace );
                M_feSpace = feSpace;
                M_mesh = feSpace->mesh();
                this->init( mpl::bool_<fespace_type::is_composite>() );
            }
        }

#endif

    //copy constructor
    ReducedBasisSpace( ReducedBasisSpace const& rb )
        :
        super ( rb ),
        M_primal_rb_basis( rb.M_primal_rb_basis ),
        M_dual_rb_basis( rb.M_dual_rb_basis ),
        M_mesh( rb.M_mesh )
        {}

    ReducedBasisSpace( this_ptrtype const& rb )
        :
        super ( *rb ),
        M_primal_rb_basis( rb->M_primal_rb_basis ),
        M_dual_rb_basis( rb->M_dual_rb_basis ),
        M_mesh( rb->M_mesh )
        {}

    void init( mpl::false_ )
        {
            //M_rbbasis = rbbasis_ptrtype( new rb_basis_type() );
        }
    void init( mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::InitializeRbSubSpace<this_type>( *this ) );
        }

    void setup( boost::property_tree::ptree const& ptree, std::string const& dbDir )
        {
            tic();

            int rbdim = ptree.template get<size_type>( "dimension" );
            this->setDimension( rbdim );
            std::string meshFilename = ptree.template get<std::string>( "mesh-filename" );
            if ( !dbDir.empty() && !fs::path(meshFilename).is_absolute() )
                meshFilename = (fs::path(dbDir)/fs::path(meshFilename).filename()).string();
            size_type meshUpdateContext = size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
            auto meshCtxInPtree = ptree.template get_optional<size_type>("mesh-context");
            if ( meshCtxInPtree )
                meshUpdateContext = *meshCtxInPtree;
            auto mesh = loadMesh(_mesh=new mesh_type(this->worldCommPtr()),_filename=meshFilename,
                                 _update=meshUpdateContext );
                                 //_update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES));
                                 //_update=size_type(MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES));
                                 //_update=size_type(MESH_UPDATE_FACES|MESH_UPDATE_EDGES));
            toc("ReducedBasisSpace::setup : load mesh",FLAGS_v>0);
            tic();
            auto spaceMeshSupport = typename functionspace_type::mesh_support_vector_type();
            auto feSpace = functionspace_type::New( _mesh=mesh,
                                                    _worldscomm=this->worldsComm(),
                                                    _range=spaceMeshSupport
                                                    );
            this->setFunctionSpace( feSpace );
            toc("ReducedBasisSpace::setup : init spaces",FLAGS_v>0);
        }


    /*
     * Get the mesh
     */
    mesh_ptrtype const& mesh() const
        {
            return M_mesh;
        }


    /*
     * add a new basis
     */
    void addPrimalBasisElement( space_element_type const & e )
        {
            space_element_ptrtype basis = M_feSpace->elementPtr();
            *basis = e;
            this->addPrimalBasisElement( basis );
        }

    void addPrimalBasisElement( space_element_ptrtype const & e )
        {
            M_primal_rb_basis.push_back( e );
            this->addPrimalBasisElementInSubSpace( *e, mpl::bool_<fespace_type::is_composite>() );

        }
    void addPrimalBasisElement( vector_ptrtype const & vec )
        {
            space_element_ptrtype e = M_feSpace->elementPtr();
            *e = *vec;
            this->addPrimalBasisElement( e );
        }

    space_element_type& primalBasisElement( int index )
        {
            int size = M_primal_rb_basis.size();
            CHECK( index < size ) << "bad index value, size of the RB "<<size<<" and index given : "<<index;
            return *M_primal_rb_basis[index];
        }
    void addDualBasisElement( space_element_type const & e )
        {
            space_element_ptrtype basis = M_feSpace->elementPtr();
            *basis = e;
            this->addDualBasisElement( basis );

        }

    void addDualBasisElement( space_element_ptrtype const & e )
        {
            M_dual_rb_basis.push_back( e );
            this->addDualBasisElementInSubSpace( *e, mpl::bool_<fespace_type::is_composite>() );
        }
    void addDualBasisElement( vector_ptrtype const & vec )
        {
            space_element_ptrtype e = M_feSpace->elementPtr();
            *e = *vec;
            this->addDualBasisElement( e );
        }


    space_element_type& dualBasisElement( int index )
        {
            int size = M_dual_rb_basis.size();
            CHECK( index < size ) << "bad index value, size of the RB "<<size<<" and index given : "<<index;
            return *M_dual_rb_basis[index];
        }

    void deleteLastPrimalBasisElements( int number )
        {
            int size = M_primal_rb_basis.size();
            CHECK( number < size )<<" error you want to delete "<<number<<" elements in a basis that contains only "<<size<<" elements\n";
            for(int i=0; i<number; i++)
            {
                M_primal_rb_basis.pop_back();
            }
            this->deleteLastPrimalBasisElementsInSubSpace( number, mpl::bool_<fespace_type::is_composite>() );
        }

    void deleteLastDualBasisElements( int number )
        {
            int size = M_dual_rb_basis.size();
            CHECK( number < size )<<" error you want to delete "<<number<<" elements in a basis that contains only "<<size<<" elements\n";
            for(int i=0; i<number; i++)
            {
                M_dual_rb_basis.pop_back();
            }
            this->deleteLastDualBasisElementsInSubSpace( number, mpl::bool_<fespace_type::is_composite>() );
        }

    /*
     * Get basis of the reduced basis space
     */
    rb_basis_type & primalRB()
        {
            return M_primal_rb_basis;
        }
    rb_basis_type & dualRB()
        {
            return M_dual_rb_basis;
        }
    rb_basis_type const& primalRB() const
        {
            return M_primal_rb_basis;
        }
    rb_basis_type const& dualRB() const
        {
            return M_dual_rb_basis;
        }

    /*
     * Set basis of the reduced basis space
     */
    void setBasis( boost::tuple< rb_basis_type , rb_basis_type > const& tuple )
        {
            this->setPrimalBasis( boost::get<0>( tuple ) );
            this->setDualBasis( boost::get<1>( tuple ) );
            /*auto primal = tuple.template get<0>();
            auto dual = tuple.template get<1>();
            M_primal_rb_basis = primal;
             M_dual_rb_basis = dual;*/
        }

    void setPrimalBasis( rb_basis_type const& rb)
        {
            M_primal_rb_basis = rb;
            this->setPrimalBasisInSubSpace( M_primal_rb_basis, mpl::bool_<fespace_type::is_composite>() );
        }
    void setDualBasis( rb_basis_type const& rb)
        {
            M_dual_rb_basis = rb;
            this->setDualBasisInSubSpace( M_dual_rb_basis, mpl::bool_<fespace_type::is_composite>() );
        }

    /**
     * update primal basis in subrbspace (do nothing if not composite)
     */
    void updatePrimalBasisForUse()
        {
            this->setPrimalBasisInSubSpace( M_primal_rb_basis, mpl::bool_<fespace_type::is_composite>() );
        }
    /**
     * update dual basis in subrbspace (do nothing if not composite)
     */
    void updateDualBasisForUse()
        {
            this->setDualBasisInSubSpace( M_dual_rb_basis, mpl::bool_<fespace_type::is_composite>() );
        }

    /**
     * orthonormalize the reduced basis using Gram-Schmidt process.
     * the tolerance is used to reiterate the process, by default we do not reiterate
     * @param norm norm to use (L2, H1, dot)
     * @param dual orthonormalize dual basis
     * @param tol tolerance above which we reiterate the process
     */
    void orthonormalize(std::string const& norm = "L2", bool dual = false, double tol = 1000)
        {
            auto wn = dual ? this->dualRB() : this->primalRB();
            auto N = wn.size();
            if( N == 0 )
                return;

            auto mf = form2(_test=this->functionSpace(), _trial=this->functionSpace());
            if( norm == "L2" || norm == "H1" )
                mf = integrate(_range=elements(support(this->functionSpace())),
                               _expr=inner(id(wn[0]), idt(wn[0])));
            if( norm == "H1" )
                mf += integrate(_range=elements(support(this->functionSpace())),
                                _expr=inner(grad(wn[0]), gradt(wn[0])));
            auto m = mf.matrixPtr();

            double max;
            int iter = 0;
            do
            {
                for( size_type i = 0; i < N; ++i )
                {
                    auto & wni = unwrap_ptr( wn[i] );
                    for( size_type j = 0; j < i; ++j )
                    {
                        auto const& wnj = unwrap_ptr( wn[j] );
                        double pij;
                        if( norm == "L2" || norm == "H1" )
                            pij = m->energy( wni, wnj );
                        else
                            pij = dot( wni, wnj );
                        wni.add( -pij, wnj );
                    }
                }
                for( size_type i = 0; i < N; ++i )
                {
                    auto& wni = unwrap_ptr( wn[i] );
                    double pii;
                    if( norm == "L2" || norm == "H1" )
                        pii = math::sqrt( m->energy( wni, wni ) );
                    else
                        pii = math::sqrt( dot(wni, wni ) );
                    wni.scale( 1./pii );
                }
                max = 0;
                for ( size_type i = 0; i < N; ++i )
                {
                    auto const & wni = unwrap_ptr( wn[i] );
                    for ( size_type j = 0; j < i; ++j )
                    {
                        auto const& wnj = unwrap_ptr( wn[j] );
                        double pij;
                        if( norm == "L2" || norm == "H1" )
                            pij = std::abs(m->energy( wni, wnj ) );
                        else
                            pij = std::abs(dot( wni, wnj ) );
                        if( pij > max )
                            max = pij;
                    }
                }
                iter++;
            } while(max > tol && iter < N);
        }

    //basis of RB space are elements of FEM function space
    //return value of the N^th basis ( vector ) at index idx
    //idx is the global dof ( fem )
    double basisValue(int N, int idx) const
        {
            return M_primal_rb_basis[N]->globalValue( idx );
        }

    /*
     * return true if the function space is composite
     */
    bool isComposite() const
        {
            return super::is_composite;
        }
    /*
     * size of the reduced basis space
     */
    size_type size() const
        {
            return this->dimension();
        }

    /*
     * dimension of the reduced basis space
     */
    size_type dimension() const
        {
            return M_primal_rb_basis.size();
        }

    /*
     * set dimension of the reduced basis space
     */
    void setDimension( size_type dim )
        {
            M_primal_rb_basis.resize( dim );
            M_dual_rb_basis.resize( dim );
            this->setDimensionInSubSpace( dim, mpl::bool_<fespace_type::is_composite>() );
        }

    /*
     * visualize basis of rb space
     */
    void visualize ()
        {
        }

    BOOST_PARAMETER_MEMBER_FUNCTION( ( this_ptrtype ),
                                     static New,
                                     tag,
                                     ( required
                                       ( space, *) ) )
        {
            //LOG( INFO ) << "ReducedBasis NEW (new impl)";
            return NewImpl( space );
        }

    static this_ptrtype NewImpl( fespace_ptrtype const& Xh )
        {

            return this_ptrtype( new this_type( Xh ) );
        }

    /*
     * Get the FEM functionspace
     */
    super_ptrtype /*const&*/ functionSpace() const
        {
            if ( M_feSpace )
                return M_feSpace;
            return std::const_pointer_cast<super>( super::shared_from_this() );
        }

    template<int i>
    typename sub_rbfunctionspace<i>::type::fespace_ptrtype const&
    subFeFunctionSpace() const
    {
        return super::template functionSpace<i>();
    }

    /**
     * save mesh on disk
     */
    void saveMesh( std::string const& meshFilename ) const
        {
#if defined(FEELPP_HAS_HDF5)
            tic();
            M_mesh->saveHDF5( meshFilename );
            toc("ReducedBasisSpace::saveMesh",FLAGS_v>0);
#endif
        }

    class ContextRB
        :
        public fespace_type::basis_context_type,
        public std::enable_shared_from_this<ContextRB>
    {
    public:
        typedef this_type reducedbasisspace_type;
        typedef this_ptrtype reducedbasisspace_ptrtype;
        typedef typename fespace_type::basis_context_type  super;
        typedef std::pair<int,super> ctx_type;
        typedef std::pair<int,std::shared_ptr<super>> ctx_ptrtype;

        static const uint16_type nComponents = reducedbasisspace_type::nComponents;

        ContextRB( ctx_ptrtype const& ctx, reducedbasisspace_ptrtype W )
        :
        super( *ctx.second ),
        M_index( ctx.first ),
        M_rbspace( W )

            {
                CHECK( M_rbspace ) << "Invalid rb space";
                LOG( INFO ) << "size : "<< M_rbspace->size();;
                M_phi.resize( nComponents, M_rbspace->size() );
                M_grad.resize(nComponents);
                for(int c=0; c<nComponents; c++)
                {
                    M_grad[c].resize( nDim );
                    for(int d=0; d<nDim; d++)
                    {
                        M_grad[c][d].resize(M_rbspace->size());
                    }
                }//end of resize

            }

        ContextRB() = default;
        ContextRB( ContextRB const& rbCtx ) = default;
        ContextRB& operator=( ContextRB const& ) = default;

        int pointIndex() const
            {
                return M_index;
            }

        reducedbasisspace_ptrtype rbSpace() const { return M_rbspace; }

        void setRbSpace( reducedbasisspace_ptrtype const& rbspace ) { M_rbspace = rbspace; }

        void setMeshForRbContext( mesh_ptrtype const& meshForRbContext )  { M_meshForRbContext = meshForRbContext; }

        void update()
            {
                int npts = this->nPoints();
                int N = M_rbspace->size();

                //
                // id
                //
                M_phi = M_rbspace->evaluateBasis( std::dynamic_pointer_cast<super>( this->shared_from_this() ) );

                //
                // grad
                //
                for( int i=0; i<N; i++)
                {
                    //matrix containing grad at all points
                    auto evaluation = M_rbspace->evaluateGradBasis( i , std::dynamic_pointer_cast<super>( this->shared_from_this() ) );
                    for(int c=0; c<nComponents; c++)
                    {
                        for(int d=0; d<nDim; d++)
                        {
                            //evaluation is a matrix
                            M_grad[c][d]( i ) = evaluation( d, c );
                        }//dimension
                    }//components
                }//rb functions

                for(int c=0; c<nComponents; c++)
                {
                    for(int d=0; d<nDim; d++)
                    {
                        DVLOG( 2 ) << "matrix M_grad of component "<<c<<"and dim "<<d<<" : \n"<<M_grad[c][d];
                    }
                }//components

            } // update


        eigen_vector_type id( eigen_vector_type const& coeffs, mpl::bool_<false> ) const
            {
                // this should work for scalar and vector fields we expect
                // coeffs to be (d,N) and M_phi (d,vN) where d is the number of
                // components and N the number of basis functions
                // TODO: check with vectorial functions
                //return (coeffs*M_phi.transpose()).diagonal();
                return M_phi*coeffs;

            }
        eigen_vector_type id( eigen_vector_type const& coeffs, mpl::bool_<true> ) const
            {
                // this should work for scalar and vector fields we expect
                // coeffs to be (d,N) and M_phi (d,vN) where d is the number of
                // components and N the number of basis functions
                // TODO: check with vectorial functions
                //return (coeffs*M_phi.transpose()).diagonal();
                return M_phi * coeffs;
            }
        //evaluation at only one node
        eigen_vector_type id( eigen_vector_type const& coeffs ) const
            {
                return id( coeffs, mpl::bool_<(nComponents>1)>() );
            }

        /*
         * for given element coefficients, gives the gradient of the element at node given in context_fem
         */
        eigen_vector_type grad( eigen_vector_type const& coeffs , int node_index=-1 ) const
            {
#if 0
                int npts = super::nPoints();
                DCHECK(npts > node_index)<<"node_index "<<node_index<<" must be lower that npts "<<npts;
                if( node_index >=0 )
                    npts=1; //we study only the node at node_index
#endif
                eigen_vector_type result ;
                for(int d=0; d<nDim; d++)
                {
                    for(int c=0; c<nComponents; c++)
                    {
                        //eigen_vector_type result_comp_dim( npts );
                        auto prod = coeffs.transpose()*M_grad[c][d];
                        eigen_vector_type result_comp_dim=prod.transpose();
                        //concatenate
                        int new_size = result.size() + result_comp_dim.size();
                        eigen_vector_type tmp( result );
                        result.resize( new_size );
                        result << tmp,result_comp_dim;
                    }
                }
#if 0
                //we now reorganize datas
                eigen_vector_type tmp( result );
                for(int p=0; p<npts; p++)
                {
                    for(int d=0; d<nDim; d++)
                    {
                        for(int c=0; c<nComponents; c++)
                            result(p*nComponents*nDim+c+d) = tmp(p+npts*c*d);
                    }
                }
#endif
                return result;

            }//grad


        /*
         * for given element coefficients, gives the dx,dy or dy of the element at node given in context_fem
         */
        eigen_vector_type d(int N, eigen_vector_type const& coeffs , int node_index=-1 ) const
            {
                int npts = super::nPoints();
                DCHECK(npts > node_index)<<"node_index "<<node_index<<" must be lower that npts "<<npts;
                if( node_index >=0 )
                    npts=1; //we study only the node at node_index

                eigen_vector_type result ;
                for(int c=0; c<nComponents; c++)
                {
                    eigen_vector_type result_comp_dim( npts );

                    auto prod = coeffs.transpose()*M_grad[c][N];
                    result_comp_dim=prod.transpose();
                    //concatenate
                    int new_size = result.size() + result_comp_dim.size();
                    eigen_vector_type tmp( result );
                    result.resize( new_size );
                    result << tmp,result_comp_dim;
                }
#if 0
                //we now reorganize datas
                eigen_vector_type tmp( result );
                for(int p=0; p<npts; p++)
                {
                    for(int c=0; c<nComponents; c++)
                        result(p*nComponents*nDim+c+N) = tmp(p+npts*c*N);
                }
#endif
                return result;
            }//d

        auto id( int ldof, int q ) const -> decltype( super::id( ldof, q ) )
            {
                return super::id( ldof, q );
            }
        value_type id( int ldof, int c1, int c2, int q ) const
            {
                return super::id( ldof, c1, c2, q );
            }
        auto grad( int ldof, int q ) const -> decltype( super::grad( ldof, q ) )
            {
                return super::grad( ldof, q );
            }
        value_type grad( int ldof, int c1, int c2, int q ) const
            {
                return super::grad( ldof, c1, c2, q );
            }

        value_type d( int ldof, int c1, int c2, int q ) const
            {
                return super::d( ldof, c1, c2, q );
            }

    private :

        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                ar & BOOST_SERIALIZATION_NVP( M_index );
                ar & BOOST_SERIALIZATION_NVP( M_phi );
                ar & BOOST_SERIALIZATION_NVP( M_grad );
                typename super::geometric_mapping_context_ptrtype/*auto*/ gmContext = this->gmContext();
                auto const& meshEltCtx = gmContext->element();
                ar & BOOST_SERIALIZATION_NVP( meshEltCtx );
                //ar & BOOST_SERIALIZATION_NVP( *gmContext );
                ar & boost::serialization::make_nvp( "gmContext", *gmContext );

            }

        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                // std::cout <<" load crbctx archive\n";
                ar & BOOST_SERIALIZATION_NVP( M_index );
                ar & BOOST_SERIALIZATION_NVP( M_phi );
                ar & BOOST_SERIALIZATION_NVP( M_grad );
                // std::cout << "M_phi="<<M_phi<<"\n";
                geoelement_type meshEltCtx;
                ar & BOOST_SERIALIZATION_NVP( meshEltCtx );
                if ( M_rbspace )
                {
                    // std::cout << "has rbspace\n";
                    if ( M_meshForRbContext )
                    {
                        if ( !M_meshForRbContext->hasElement( meshEltCtx.id() ) )
                        {
                            M_meshForRbContext->addElement( meshEltCtx, false );
                        }
                        auto const& meshEltCtxRegister = M_meshForRbContext->element( meshEltCtx.id() );
                        //typename super::geometric_mapping_context_ptrtype gmContext( new typename super::geometric_mapping_context_type( M_meshForRbContext->gm(),meshEltCtxRegister ) );
                        typename super::geometric_mapping_context_ptrtype gmContext = M_meshForRbContext->gm()->template context<0>( meshEltCtxRegister, typename super::geometric_mapping_context_type::precompute_ptrtype{} );
                        ar & boost::serialization::make_nvp( "gmContext", *gmContext );
                        this->setGmContext( gmContext );
                    }
                    else if ( M_rbspace->mesh() )
                    {
                        // std::cout << "has mesh in rbspace\n";
                        CHECK ( M_rbspace->mesh()->hasElement( meshEltCtx.id() ) ) << "fails because mesh doesnt have the element reloaded for gmc";
                        auto const& meshEltCtxRegister = M_rbspace->mesh()->element( meshEltCtx.id() );
                        //typename super::geometric_mapping_context_ptrtype gmContext( new typename super::geometric_mapping_context_type( M_rbspace->mesh()->gm(),meshEltCtxRegister ) );
                        typename super::geometric_mapping_context_ptrtype gmContext = M_rbspace->mesh()->gm()->template context<0>( meshEltCtxRegister, typename super::geometric_mapping_context_type::precompute_ptrtype{} );
                        ar & boost::serialization::make_nvp( "gmContext", *gmContext );
                        this->setGmContext( gmContext );
                    }
                    else
                        CHECK( false ) << "fails because no mesh defined";
                }
                //typename super::geometric_mapping_context_ptrtype gmc( new typename super::geometric_mapping_context_type( M_rbspace->mesh()->gm(),meshEltCtx ) );

            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

    private :
        int M_index;
        reducedbasisspace_ptrtype M_rbspace;
        mesh_ptrtype M_meshForRbContext; // only init when reload serialisation and no mesh defined in rbspace
        eigen_matrix_type M_phi;
        std::vector<std::vector< eigen_vector_type> > M_grad;

    };

    typedef ContextRB ctxrb_type;
    typedef std::shared_ptr<ctxrb_type> ctxrb_ptrtype;



    /*
     * store the evaluation of basis functions in given points
     */
    class ContextRBSet
        : public std::map<int,std::shared_ptr<ContextRB> >
    {
        typedef typename std::map<int,std::shared_ptr<ContextRB>> super;

    public :
        static const bool is_rb_context = true;


        typedef ReducedBasisSpace<SpaceType> rbspace_type;
        typedef std::shared_ptr<rbspace_type> rbspace_ptrtype;

        typedef typename rbspace_type::super_ptrtype functionspace_ptrtype;
        typedef typename rbspace_type::super functionspace_type;
        typedef typename functionspace_type::Context fe_context_type;
        typedef std::shared_ptr<fe_context_type> fe_context_ptrtype;

        typedef typename super::iterator iterator;

        typedef rbspace_type::eigen_vector_type eigen_vector_type;

        typedef rbspace_type::rb_basis_type rb_basis_type;
        typedef std::shared_ptr<rb_basis_type> basis_ptrtype;

        typedef Eigen::MatrixXd eigen_matrix_type;

        typedef rbspace_type::mesh_type mesh_type;

        //static const bool is_function_space_scalar = functionspace_type::is_scalar;
        //static const bool is_function_space_vectorial = functionspace_type::is_vectorial;
        //static const uint16_type nComponents = functionspace_type::nComponents;
        //static const uint16_type nDim = rbspace_type::nDim;

        ~ContextRBSet() {}

        ContextRBSet( super const& s )
            :
            super( s ),
            M_ctx_fespace(s.begin()->second->rbSpace()->functionSpace()),
            M_rbspace(s.begin()->second->rbSpace()),
            M_ctxHaveBeenMpiBroadcasted( false )
            {
                for( auto& m : s )
                {
                    M_ctx_fespace.addCtx( m.second, M_rbspace->worldComm().globalRank() );
                }
                M_nPoints = M_ctx_fespace.nPoints();
            }

        ContextRBSet( rbspace_ptrtype const& rbspace )
            :
            super(),
            M_ctx_fespace( rbspace->functionSpace() ),
            M_rbspace( rbspace ),
            M_nPoints( M_ctx_fespace.nPoints() ),
            M_ctxHaveBeenMpiBroadcasted( false )
            {
                update();
            }

        //only because the function functionSpace()
        //returns an object : rb_functionspaceptrtype
        ContextRBSet( functionspace_ptrtype const& functionspace )
            :
            super(),
            M_ctx_fespace( functionspace ),
            M_nPoints( M_ctx_fespace.nPoints() ),
            M_ctxHaveBeenMpiBroadcasted( false )
            {
                update();
            }

        ContextRBSet()
            :
            M_nPoints( 0 ),
            M_ctxHaveBeenMpiBroadcasted( false )
            {}
        ContextRBSet( ContextRBSet const& ctxRbSet ) = default;
        ContextRBSet& operator=( ContextRBSet const& ) = default;


        functionspace_ptrtype functionSpace() const
            {
                return M_rbspace->functionSpace();
            }

        void update()
            {
                for( auto& c : *this )
                {
                    rank_type procId = M_ctx_fespace.processorHavingPoint(c.first);
                    if ( procId == M_rbspace->worldComm().globalRank() )
                        c.second->update();
                }
                // update context for each process + define context mesh
                if ( this->nPoints() > 0 )
                {
                    for ( int k=0;k<this->nPoints();++k )
                    {
                        this->syncCtx( k );
                    }
                    M_ctxHaveBeenMpiBroadcasted = true;
                }
            }

        bool ctxHaveBeenMpiBroadcasted() const { return M_ctxHaveBeenMpiBroadcasted; }

        rbspace_ptrtype const& rbFunctionSpace() const
            {
                return M_rbspace;
            }

        rbspace_type* ptrFunctionSpace() const
            {
                return M_rbspace.get();
            }

        void setRbFunctionSpace( rbspace_ptrtype const& rbSpace )
            {
                M_rbspace = rbSpace;
            }

        std::pair<iterator,bool>
        add( node_type const& t )
        {
            return add( t, mpl::bool_<is_composite>() );
        }
        std::pair<iterator,bool>
        add( node_type const& t, mpl::bool_<true> )
        {
            return std::make_pair( this->end(), false );
        }
        std::pair<iterator,bool>
        add( node_type const& t, mpl::bool_<false> )
        {
            // we suppose here that the point will be added which means that it
            // has been located in the underlying mesh
            auto ret = M_ctx_fespace.add( t );
            M_nPoints = M_ctx_fespace.nPoints();
            // if the current processor handles the point previously located
            // add it to the ContextRB set
            if ( ret.second )
            {
                return this->insert( std::make_pair( ret.first->first, std::make_shared<ContextRB>( *ret.first, M_rbspace ) ) );
            }
            return std::make_pair( this->end(), false );
        }
        int nPoints() const
        {
            return M_nPoints;//M_ctx_fespace.nPoints();
        }

        fe_context_type const& feContext() const { return M_ctx_fespace; }

    private :

        void syncCtx( int ptId )
        {
            rank_type procId = M_ctx_fespace.processorHavingPoint(ptId);
            rank_type myrank = M_rbspace->worldComm().rank();

            if ( !M_meshForRbContext )
                M_meshForRbContext = std::make_shared<mesh_type>( M_rbspace->worldCommPtr() );

            ctxrb_ptrtype rbCtxReload;
            if ( myrank == procId )
            {
                // update additional mesh used in rb context
                CHECK( this->find( ptId ) != this->end() ) << "point id is not saved on this process";
                rbCtxReload = this->operator[]( ptId );
                CHECK( rbCtxReload ) << "rbCtxReload not init";
                rbCtxReload->setMeshForRbContext( M_meshForRbContext );
                auto const& modelMeshEltCtx = rbCtxReload->gmContext()->element();
                if ( !M_meshForRbContext->hasElement( modelMeshEltCtx.id(), modelMeshEltCtx.processId() ) )
                {
                    geoelement_type meshEltCtx = modelMeshEltCtx;;
                    M_meshForRbContext->addElement( meshEltCtx, false );
                }
            }
            else
            {
                // create a new rb context
                rbCtxReload.reset( new ctxrb_type );
                rbCtxReload->setRbSpace( M_rbspace );
                rbCtxReload->setMeshForRbContext( M_meshForRbContext );
            }
            // recv rb context from process which have localized the point
            mpi::broadcast( M_rbspace->worldComm().globalComm(), *rbCtxReload, procId );

            // save rb context with process which doesnt have localized the point
            if ( myrank != procId )
                this->operator[]( ptId ) = rbCtxReload;

            // std::cout << "["<<M_rbspace->worldComm().rank()<<"] M_meshForRbContext->numElements() : " << M_meshForRbContext->numElements() << "\n";
        }


        friend class boost::serialization::access;

        template<class Archive>
        void save( Archive & ar, const unsigned int version ) const
            {
                ar & BOOST_SERIALIZATION_NVP( M_nPoints );
                std::vector<int> rbCtxKeys;
                for ( auto const& rbCtxPair : *this )
                    rbCtxKeys.push_back( rbCtxPair.first );
                ar & BOOST_SERIALIZATION_NVP( rbCtxKeys );
                for ( int rbCtxKey : rbCtxKeys )
                {
                    std::string rbctxNameInSerialization = (boost::format("rbSpaceContext_%1%")%rbCtxKey).str();
                    ar & boost::serialization::make_nvp( rbctxNameInSerialization.c_str(), *(this->find( rbCtxKey )->second) );
                }

            }

        template<class Archive>
        void load( Archive & ar, const unsigned int version )
            {
                // std::cout << "reload rbctxset\n";
                if ( M_rbspace )
                {
                    if ( !M_meshForRbContext )
                        M_meshForRbContext = std::make_shared<mesh_type>( M_rbspace->worldCommPtr() );
                }

                ar & BOOST_SERIALIZATION_NVP( M_nPoints );
                std::vector<int> rbCtxKeys;
                ar & BOOST_SERIALIZATION_NVP( rbCtxKeys );
                for ( int rbCtxKey : rbCtxKeys )
                {
                    std::string rbctxNameInSerialization = (boost::format("rbSpaceContext_%1%")%rbCtxKey).str();
                    std::shared_ptr<ContextRB> rbCtxReload( new ContextRB );
                    if ( M_rbspace )
                    {
                        rbCtxReload->setRbSpace( M_rbspace );
                        rbCtxReload->setMeshForRbContext( M_meshForRbContext );
                    }
                    ar & boost::serialization::make_nvp( rbctxNameInSerialization.c_str(), *rbCtxReload );
                    this->operator[]( rbCtxKey ) = rbCtxReload;
                }
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()

    private :
        fe_context_type M_ctx_fespace;
        rbspace_ptrtype M_rbspace;
        int M_nPoints;
        mesh_ptrtype M_meshForRbContext; // only init when reload serialisation and no mesh defined in rbspace
        bool M_ctxHaveBeenMpiBroadcasted;
    };
    typedef ContextRBSet ctxrbset_type;
    typedef std::shared_ptr<ctxrbset_type> ctxrbset_ptrtype;

    /**
     * \return function space context
     */
    ContextRBSet context() { return ContextRBSet( std::dynamic_pointer_cast< ReducedBasisSpace<SpaceType> >( this->shared_from_this() ) ); }

    ctxrb_ptrtype
    contextBasis( std::pair<int,std::shared_ptr<ContextRB>> const& p, ContextRBSet const& c )
        {
            return p.second;
        }

    /**
     * evaluate all basis functions at (only one) node given by ctx
     * returns matrix of size nComponents x N if N is RB size
     */
    template <typename ContextBasisFem>
    eigen_matrix_type evaluateBasis( std::shared_ptr<ContextBasisFem> const& ctx )
        {
            eigen_matrix_type m( nComponents, M_primal_rb_basis.size() );

            auto s = this->functionSpace()->context();
            s.addCtx( ctx, this->functionSpace()->worldComm().localRank() );

            //index of the column to be filled
            int c = 0;

            //loop over each basis functions
            for( auto& b : M_primal_rb_basis )
            {
                //in a parallel context the point given by ctx is located only on one proc
                //so we need to indicate to evaluateFromContext not do "all_reduce" at the end
                //that is why we use a boolean mpi_communications
                m.col(c++) = evaluateFromContext( _context=s , _expr=idv(b) , _mpi_communications=false );
            }
            return m;
        }

    template <typename ContextBasisFem>
    eigen_matrix_type evaluateGradBasis( int i , std::shared_ptr<ContextBasisFem> const& ctxs )
        {
            auto ctx = this->functionSpace()->context();
            ctx.addCtx( ctxs, this->functionSpace()->worldComm().localRank() );

            int npts = ctx.nPoints();
            eigen_matrix_type evaluation;
            evaluation.resize( nDim , npts*nComponents );

            //in a parallel context the point given by ctx is located only on one proc
            //so we need to indicate to evaluateFromContext not do "all_reduce" at the end
            //that is why we use a boolean _mpi_communications

            if( nDim >= 1 )
            {
                auto evalx = evaluateFromContext( _context=ctx , _expr= dxv( M_primal_rb_basis[i] ) , _mpi_communications=false );
                for(int p=0; p<npts; p++)
                {
                    for(int c=0; c<nComponents;c++)
                        evaluation( 0 , p*nComponents+c ) = evalx( p*nComponents+c );
                }
            }
            if( nDim >= 2 )
            {
                auto evaly = evaluateFromContext( _context=ctx , _expr= dyv( M_primal_rb_basis[i] ), _mpi_communications=false );
                for(int p=0; p<npts; p++)
                {
                    for(int c=0; c<nComponents;c++)
                        evaluation( 1 , p*nComponents+c ) = evaly( p*nComponents+c );
                }
            }
            if( nDim == 3 )
            {
                auto evalz = evaluateFromContext( _context=ctx , _expr= dzv( M_primal_rb_basis[i] ), _mpi_communications=false );
                for(int p=0; p<npts; p++)
                {
                    for(int c=0; c<nComponents;c++)
                        evaluation( 2 , p*nComponents+c ) = evalz( p*nComponents+c );
                }
            }
            return evaluation;
        }

    /*
     *  contains coefficients in the RB expression
     *  i.e. contains coefficients u_i in
     *  u_N ( \mu ) = \sum_i^N u_i ( \mu ) \xi_i( x )
     *  where \xi_i( x ) are basis function
     */
    template<typename T = double,  typename Cont = Eigen::Matrix<T,Eigen::Dynamic,1> >
    class Element
        :
        public Cont,boost::addable<Element<T,Cont> >, boost::subtractable<Element<T,Cont> >, public FunctionSpaceBase::ElementBase
    {
    public:
        typedef ReducedBasisSpace<SpaceType> rbspace_type;
        typedef std::shared_ptr<rbspace_type> rbspace_ptrtype;

        typedef rbspace_type::super functionspace_type;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

        typedef Eigen::MatrixXd eigen_matrix_type;
        typedef rbspace_type::eigen_vector_type eigen_vector_type;

        //element of the FEM function space
        typedef rbspace_type::space_element_type space_element_type;

        //RB context
        typedef rbspace_type::ContextRBSet ctxrbset_type;
        typedef rbspace_type::ContextRB ctxrb_type;
        typedef std::shared_ptr<rbspace_type::ContextRB> ctxrb_ptrtype;

        typedef T value_type;

        typedef Cont super;
        typedef Cont container_type;


        typedef typename functionspace_type::gm_type gm_type;
        typedef std::shared_ptr<gm_type> gm_ptrtype;

        /**
         * geometry typedef
         */
        typedef typename mesh_type::element_type geoelement_type;
        //typedef typename gm_type::template Context<vm::POINT|vm::JACOBIAN|vm::HESSIAN|vm::KB, geoelement_type> gmc_type;
        //typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        using gmc_type = typename rbspace_type::gmc_type;
        using gmc_ptrtype = typename rbspace_type::gmc_ptrtype;
        typedef typename gm_type::precompute_ptrtype geopc_ptrtype;
        typedef typename gm_type::precompute_type geopc_type;


        static const uint16_type nRealDim = mesh_type::nRealDim;
        static const bool is_composite = functionspace_type::is_composite;
        static const bool is_mortar = functionspace_type::is_mortar;

        typedef typename mpl::if_<mpl::bool_<is_composite>,
                                  mpl::identity<boost::none_t>,
                                  mpl::identity<typename basis_0_type::polyset_type> >::type::type polyset_type;

        typedef typename mpl::if_<mpl::bool_<is_composite>,
                                  mpl::identity<boost::none_t>,
                                  mpl::identity<typename basis_0_type::PreCompute> >::type::type pc_type;
        typedef std::shared_ptr<pc_type> pc_ptrtype;

        static const uint16_type nComponents1 = functionspace_type::nComponents1;
        static const uint16_type nComponents2 = functionspace_type::nComponents2;

        typedef boost::multi_array<value_type,3> array_type;
#if 0
        typedef Eigen::Matrix<value_type,nComponents1,1> _id_type;
        typedef Eigen::Matrix<value_type,nComponents1,nRealDim> _grad_type;
#else
        using _id_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nComponents2>>;
        using _grad_type = Eigen::TensorFixedSize<value_type,Eigen::Sizes<nComponents1,nRealDim>>;
        using eigen_matrix_to_tensor_map = Eigen::TensorMap<Eigen::Matrix<value_type, nComponents1, 1> >;
#endif
        typedef boost::multi_array<_id_type,1> id_array_type;
        typedef boost::multi_array<_grad_type,1> grad_array_type;

        typedef typename matrix_node<value_type>::type matrix_node_type;

        friend class ReducedBasisSpace<SpaceType>;

        Element()
            {}

        Element( rbspace_ptrtype const& rbspace , std::string const& name="u")
            :
            M_femfunctionspace( rbspace->functionSpace() ),
            M_rbspace( rbspace ),
            M_name( name )
            {
                this->resize( M_rbspace.dimension() );
            }

        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh()
            {
                return M_femfunctionspace->mesh();
            }
        /**
         * \return the mesh associated to the function
         */
        mesh_ptrtype mesh() const
            {
                return M_femfunctionspace->mesh();
            }

        /*
         * return the reduced basis space associated
         */
        void setReducedBasisSpace( rbspace_ptrtype const& rbspace )
            {
                M_rbspace = rbspace;
            }

        /*
         * return the size of the element
         */
        int size() const
            {
                return super::size();
            }

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

        void setCoefficient( int index , double value )
            {
                int size=super::size();
                FEELPP_ASSERT( index < size )(index)(size).error("invalid index");
                this->operator()( index ) = value;
            }

        value_type operator()( size_t i ) const
            {
                int size=super::size();
                FEELPP_ASSERT( i < size )(i)(size).error("invalid index");
                return super::operator()( i );
            }


        int localSize() const
            {
                return super::size();
            }

        value_type& operator()( size_t i )
            {
                int size = super::size();
                FEELPP_ASSERT( i < size )(i)(size).error("invalid index");
                return super::operator()( i );
            }

        Element& operator+=( Element const& _e )
            {
                int size1=super::size();
                int size2=_e.size();
                FEELPP_ASSERT( size1 == size2 )(size1)(size2).error("invalid size");
                for ( int i=0; i <size1; ++i )
                    this->operator()( i ) += _e( i );
                return *this;
            }

        Element& operator-=( Element const& _e )
            {
                int size1=super::size();
                int size2=_e.size();
                FEELPP_ASSERT( size1 == size2 )(size1)(size2).error("invalid size");
                for ( int i=0; i <size2; ++i )
                    this->operator()( i ) -= _e( i );
                return *this;
            }

        space_element_type expansion( int  N=-1)
            {
                return M_rbspace.expansion( *this, N );
            }

        //basis of RB space are elements of FEM function space
        //return value of the N^th basis ( vector ) at index idx
        double basisValue( int N, int idx) const
            {
                return M_rbspace.basisValue( N , idx );
            }

        eigen_vector_type evaluate(  ctxrb_type const& context_rb )
            {
                return context_rb.id( *this );
            }

        /*
         * evaluate the element to nodes in contextRBSet
         */
        eigen_vector_type evaluate(  ctxrbset_type const& context_rb, bool do_comm=true )
            {
                bool do_communications = do_comm && M_rbspace.worldComm().globalSize() > 1 && !context_rb.ctxHaveBeenMpiBroadcasted();

                int npts = context_rb.nPoints();
                //from now we manipulate local datas
                //one proc can have zero, one or more context_rb
                eigen_vector_type local_result( npts*nComponents );
                eigen_vector_type global_result( npts*nComponents );
                local_result.setZero();
                global_result.setZero();
                auto it = context_rb.begin();
                auto en = context_rb.end();

                //loop on local points
                for ( int p = 0; it!=en ; ++it, ++p )
                {
                    int global_position=it->first;
                    auto ctx=it->second;
                    local_result.segment( nComponents*global_position, nComponents )  = this->evaluate( *ctx );
                }

                if( do_communications )
                    mpi::all_reduce( M_rbspace.worldComm() , local_result, global_result, std::plus< eigen_vector_type >() );
                else
                    global_result = local_result;

                return global_result;

                //int p = 0;
                //for( auto const& ctx : context_rb )
                //{
                    //result.segment( nComponents*p, nComponents )  = this->evaluate( *ctx.second );
                    //p++;
                //}

            }


        void updateGlobalValues() const
            {
                //nothing to do
            }
        bool areGlobalValuesUpdated() const
            {
                return true;
            }


        template <typename ... CTX>
        decltype(auto) //ctxrb_ptrtype
        selectContext( CTX const& ... ctx ) const
            {
                typedef boost::fusion::vector<CTX...> my_vector_ctx_type;
                typedef typename boost::fusion::result_of::distance<typename boost::fusion::result_of::begin<my_vector_ctx_type>::type,
                                                                    typename boost::fusion::result_of::find<my_vector_ctx_type,ctxrb_ptrtype>::type>::type pos_ctx_type;
                typedef typename boost::fusion::result_of::distance<typename boost::fusion::result_of::begin<my_vector_ctx_type>::type,
                                                                    typename boost::fusion::result_of::find<my_vector_ctx_type,
                                                                                                            std::shared_ptr<typename ctxrb_type::super> >::type>::type pos_fectx_type;

                static const int myNumberOfCtx = boost::mpl::size<my_vector_ctx_type>::type::value;
                CHECK( pos_ctx_type::value < myNumberOfCtx || pos_fectx_type::value < myNumberOfCtx ) << "no compatible context";
                static const int ctxPosition = (pos_ctx_type::value < myNumberOfCtx)? pos_ctx_type::value : ( (pos_fectx_type::value < myNumberOfCtx)? pos_fectx_type::value : 0 );
                my_vector_ctx_type ctxvec( ctx... );
                return boost::fusion::at_c<ctxPosition>( ctxvec );
            }


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
                return shape;
            }

        template<typename Context_t>
        id_type
        id( Context_t const & context ) const
            {
                return id_type( *this, context );
            }

        //the function id_ (called by idv, idt or id) will do not the same things
        //if it has a FEM context ( mpl::bool_<false> )
        //or if it has a RB context ( mpl::bool<true> )
        //with FEM context we use pre-computations of FEM basis functions whereas with RB context we use precomputations of RB basis functions
        template<typename Context_t> void  id_( Context_t const & context, id_array_type& v ) const;
        template<typename Context_t> void  id_( Context_t const & context, id_array_type& v , mpl::bool_<true> ) const;
        template<typename Context_t> void  id_( Context_t const & context, id_array_type& v , mpl::bool_<false>) const;

        template<typename Context_t>
        void
        id( Context_t const & context, id_array_type& v ) const
            {
                id_( context, v );
            }


        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        gradExtents( ContextType const & context ) const
            {
                boost::array<typename array_type::index, 1> shape;
                shape[0] = context.xRefs().size2();
                return shape;
            }

        template<typename Context_t> void  grad_( Context_t const & context, grad_array_type& v ) const;
        template<typename Context_t> void  grad_( Context_t const & context, grad_array_type& v , mpl::bool_<true> ) const;
        template<typename Context_t> void  grad_( Context_t const & context, grad_array_type& v , mpl::bool_<false>) const;

        template<typename Context_t>
        void
        grad( Context_t const & context, grad_array_type& v ) const
            {
                grad_( context , v);
            }


        template<typename ContextType>
        boost::array<typename array_type::index, 1>
        dxExtents( ContextType const & context ) const
            {
                boost::array<typename array_type::index, 1> shape;
                shape[0] = context.xRefs().size2();
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

        template<typename Context_t> void  d_( int N, Context_t const & context, id_array_type& v ) const;
        template<typename Context_t> void  d_( int N, Context_t const & context, id_array_type& v , mpl::bool_<true> ) const;
        template<typename Context_t> void  d_( int N, Context_t const & context, id_array_type& v , mpl::bool_<false>) const;

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
        idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const;

        /**
         * these functione are here for the compilation.
         * evaluateFromContext has the kywork _projection=true
         * and in that case we project the expression on
         * the function space associated to the context
         * and it may happens that user wants to do
         * evaluateFromContext( _context=ctx, _expr=dxv(...) , _projection=false )
         * and so even we don't project dxv(...), in order to compile
         * we need to have these functions
         */
        void
        gradInterpolate( matrix_node_type __ptsReal, grad_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
        {
            CHECK( false ) << "The function gradInterpolate is not yet implemented \n";
        }
        void
        dxInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
        {
            bool go = false;
            CHECK( go ) << "The function dxInterpolate is not yet implemented \n";
        }
        void
        dyInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
        {
            bool go = false;
            CHECK( go ) << "The function dyInterpolate is not yet implemented \n";
        }
        void
        dzInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
        {
            bool go = false;
            CHECK( go ) << "The function dzInterpolate is not yet implemented \n";
        }


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
                // typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
                // typedef std::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

                // gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(),  context.pc() ) );
                auto  __c_interp = context.geometricMapping()->template context<vm::POINT>( context.element_c(),  context.pc() );
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
                typedef typename Context_t::gm_type::template Context< Context_t::context|vm::POINT, typename Context_t::element_type> gmc_interp_type;
                typedef std::shared_ptr<gmc_interp_type> gmc_interp_ptrtype;

                typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::permutation_type permutation_type;
                typedef typename Context_t::gm_type::template Context<Context_t::context,typename Context_t::element_type>::precompute_ptrtype precompute_ptrtype;

                std::vector<std::map<permutation_type, precompute_ptrtype> > __geo_pcfaces = context.pcFaces();
                gmc_interp_ptrtype __c_interp( new gmc_interp_type( context.geometricMapping(), context.element_c(), __geo_pcfaces , context.faceId() ) );
#endif
                auto __geo_pcfaces = context.pcFaces();
                auto __c_interp = context.geometricMapping()->template context<vm::POINT>( context.element_c(), __geo_pcfaces , context.faceId() );

                return __c_interp->xReal();
            }


        /**
           \return the finite element space ( not reduced basis space )
        */
        functionspace_ptrtype const& functionSpace() const
            {
                return M_femfunctionspace;
            }

    private:
        functionspace_ptrtype M_femfunctionspace;
        rbspace_type M_rbspace;
        std::string M_name;
    };//Element of the reduced basis space

    typedef Element<value_type> element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    //Access to an element of the reduced basis space
    element_type element( std::string const& name="u" )
        {
            element_type u(std::dynamic_pointer_cast<ReducedBasisSpace<SpaceType>>(this->shared_from_this() ), name );
            u.setZero();
            return u;
        }

    element_ptrtype elementPtr( std::string const& name="u" )
        {
            element_ptrtype u( new element_type( std::dynamic_pointer_cast<ReducedBasisSpace<SpaceType>>(this->shared_from_this() ) , name ) );
            u->setZero();
            return u;
        }

    space_element_type expansion( element_type const& uRB, int  N=-1) const
        {
            int basis_size = M_primal_rb_basis.size();
            int number_of_coeff = ( N == -1 )? basis_size : N;
            CHECK( number_of_coeff <= basis_size ) << "invalid size : " << number_of_coeff << " must be less or equal than " << basis_size;
            return Feel::expansion( M_primal_rb_basis, uRB , number_of_coeff );
        }

    void expansion( element_type const& uRB, Vector<value_type> & uFE, int  N=-1) const
        {
            int basis_size = M_primal_rb_basis.size();
            int number_of_coeff = ( N == -1 )? basis_size : N;
            CHECK( number_of_coeff <= basis_size ) << "invalid size : " << number_of_coeff << " must be less or equal than " << basis_size;
            Feel::expansion( M_primal_rb_basis, uRB , uFE, number_of_coeff );
        }


    rbfunctionspace_vector_type const&
    rbfunctionspaces() const
    {
        return M_rbfunctionspaces;
    }

    /**
     * get the \p i -th \c FunctionSpace out the list
     */
    template<int i>
    typename mpl::at_c<rbfunctionspace_vector_type,i>::type::second_type const&
    rbFunctionSpace() const
    {
        return fusion::at_c<i>( M_rbfunctionspaces ).second;
    }

private :
    void addPrimalBasisElementInSubSpace( space_element_type const & e, mpl::false_ ) {}
    void addPrimalBasisElementInSubSpace( space_element_type const & e, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::AddBasisElementRbSubSpace<this_type,0>( *this, e ) );
        }
    void deleteLastPrimalBasisElementsInSubSpace( int number, mpl::false_ ) {}
    void deleteLastPrimalBasisElementsInSubSpace( int number, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::DeleteLastBasisElementsRbSubSpace<this_type,0>( *this, number ) );
        }
    void setPrimalBasisInSubSpace( rb_basis_type const& rb, mpl::false_ ) {}
    void setPrimalBasisInSubSpace( rb_basis_type const& rb, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::ClearBasisElementsRbSubSpace<0>() );
            for (auto const& thebasis : rb )
                this->addPrimalBasisElementInSubSpace( *thebasis, mpl::true_() );
        }

    void addDualBasisElementInSubSpace( space_element_type const & e, mpl::false_ ) {}
    void addDualBasisElementInSubSpace( space_element_type const & e, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::AddBasisElementRbSubSpace<this_type,1>( *this, e ) );
        }
    void deleteLastDualBasisElementsInSubSpace( int number, mpl::false_ ) {}
    void deleteLastDualBasisElementsInSubSpace( int number, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::DeleteLastBasisElementsRbSubSpace<this_type,1>( *this, number ) );
        }
    void setDualBasisInSubSpace( rb_basis_type const& rb, mpl::false_ ) {}
    void setDualBasisInSubSpace( rb_basis_type const& rb, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::ClearBasisElementsRbSubSpace<1>() );
            for (auto const& thebasis : rb )
                this->addDualBasisElementInSubSpace( *thebasis, mpl::true_() );
        }
    void setDimensionInSubSpace( size_type rbdim, mpl::false_ ) {}
    void setDimensionInSubSpace( size_type rbdim, mpl::true_ )
        {
            fusion::for_each( M_rbfunctionspaces,
                              Feel::detail::SetDimensionRbSubSpace( rbdim ) );
        }

private :

    fespace_ptrtype M_feSpace;
    rbfunctionspace_vector_type M_rbfunctionspaces;
    rb_basis_type M_primal_rb_basis;
    rb_basis_type M_dual_rb_basis;
    mesh_ptrtype M_mesh;

};//ReducedBasisSpace


template<typename SpaceType>
const uint16_type ReducedBasisSpace<SpaceType>::nComponents;

template<typename SpaceType>
const uint16_type ReducedBasisSpace<SpaceType>::ContextRB::nComponents;



template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v ) const
{
    return id_( context, v, boost::is_same<Context_t,ContextRB>() );
}


template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::grad_( Context_t const & context, grad_array_type& v ) const
{
    return grad_( context, v, boost::is_same<Context_t,ContextRB>() );
}

template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::d_( int N, Context_t const & context, id_array_type& v ) const
{
    return d_( N, context, v, boost::is_same<Context_t,ContextRB>() );
}


//warning :
//if we need to evaluate element at nodes in context
//use u.evaluate( rb_context ) should go faster
//because here we do the same thing ( context.id(*this) )
//but we add a loop over points
template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v , mpl::bool_<true> ) const
{
    //eigen_matrix_to_tensor_map m( context.id( *this ) );
    auto idEigenMatrix = context.id( *this );
    Eigen::TensorMap<_id_type> m( idEigenMatrix.data(), nComponents1, 1 );
    v[0] = m;
}

template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::grad_( Context_t const & context, grad_array_type& v , mpl::bool_<true> ) const
{
    //LOG( INFO ) << " grad_ with a RB context";

    //int index = rb_context.pointIndex();
    //evaluate the gradient at a specific point (called by evaluateFromContext)
    auto evaluation = context.grad( *this );
    for(int d=0;d<nDim;d++)
    {
        for(int c=0; c<nComponents; c++)
        {
            v[0]( c,d ) += evaluation(c+d*nComponents);
        }
    }
}

template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::d_( int N, Context_t const & context, id_array_type& v , mpl::bool_<true> ) const
{
    //LOG( INFO ) << " d_ with a RB context";

    auto evaluation = context.d(N, *this );

    for(int c=0; c<nComponents; c++)
    {
        v[0]( c,0 ) = evaluation(c);
    }

}

template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::id_( Context_t const & context, id_array_type& v , mpl::bool_<false> ) const
{
    //LOG( INFO ) << "id_ with a FEM context";
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<size_type> )
        return;

    const uint16_type nq = context.xRefs().size2();

    int rb_size = this->size();
    //loop on RB dof
    for(int N=0; N<rb_size; N++)
    {

        // the RB unknown can be written as
        // u^N = \sum_i^N u_i^N \PHI_i where \PHI_i are RB basis functions (i.e. elements of FEM function space)
        // u^N = \sum_i^N u_i^N \sum_j \PHI_ij \phi_j where PHI_ij is a scalar and phi_j is the j^th fem basis function
        value_type u_i = this->operator()( N );

        for ( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;

            for ( typename array_type::index c1 = 0; c1 < ncdof; ++c1 )
            {
                typename array_type::index ldof = basis_type::nDof*c1+l;
                size_type gdof = M_femfunctionspace->dof()->localToGlobal( elt_id, l, c1 ).index();

                //N is the index of the RB basis function (i.e fem element)
                //FEM coefficient associated to the global dof "gdof" of the N^th RB element in the basis
                value_type rb_basisij = this->basisValue( N , gdof );

                //coefficient u_i^N * \PHI_ij
                value_type coefficient = u_i*rb_basisij;

                for ( uint16_type q = 0; q < nq; ++q )
                {
                    for ( typename array_type::index i = 0; i < nComponents1; ++i )
                    {
                        //context contains evaluation of FEM basis functions
                        //context.id give access to M_phi in PolynomialSet
                        v[q]( i,0 ) += coefficient*context.id( ldof, i, 0, q );
                    }
                }
            }
        }
    }//end of loop on RB dof

}


template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::grad_( Context_t const & context, grad_array_type& v , mpl::bool_<false> ) const
{

    //LOG( INFO ) << "grad_ with a FEM context";
    if ( !this->areGlobalValuesUpdated() )
        this->updateGlobalValues();

    size_type elt_id = context.eId();
    if ( context.gmContext()->element().mesh()->isSubMeshFrom( this->mesh() ) )
        elt_id = context.gmContext()->element().mesh()->subMeshToMesh( context.eId() );
    if ( context.gmContext()->element().mesh()->isParentMeshOf( this->mesh() ) )
        elt_id = this->mesh()->meshToSubMesh( context.eId() );
    if ( elt_id == invalid_v<size_type> )
        return;

    int rb_size = this->size();

    //loop on RB dof
    for(int N=0; N<rb_size; N++)
    {

        value_type u_i = this->operator()( N );

        for ( int l = 0; l < basis_type::nDof; ++l )
        {
            const int ncdof = is_product?nComponents1:1;

            for ( int c1 = 0; c1 < ncdof; ++c1 )
            {
                int ldof = c1*basis_type::nDof+l;
                size_type gdof = M_femfunctionspace->dof()->localToGlobal( elt_id, l, c1 ).index();

                //N is the index of the RB basis function (i.e fem element)
                value_type rb_basisij = this->basisValue( N , gdof );

                //coefficient u_i^N * \PHI_ij
                value_type coefficient = u_i*rb_basisij;

                for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                {
                    for ( int k = 0; k < nComponents1; ++k )
                    {
                        for ( int j = 0; j < nRealDim; ++j )
                        {
                            v[q]( k,j ) += coefficient*context.grad( ldof, k, j, q );
                        }//j
                    }//k
                }//q
            }//c1
        }//l

    }//N


}


template<typename SpaceType>
template<typename Y,  typename Cont>
template<typename Context_t>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::d_( int N, Context_t const & context, id_array_type& v , mpl::bool_<false> ) const
{

    int rb_size = this->size();

    for(int rbN=0; rbN<rb_size; rbN++)
    {
        value_type u_i = this->operator()( rbN );

        for ( int i = 0; i < basis_type::nDof; ++i )
        {
            const int ncdof = is_product?nComponents1:1;

            for ( int c1 = 0; c1 < ncdof; ++c1 )
            {
                size_type ldof = basis_type::nDof*c1 + i;
                size_type gdof = M_femfunctionspace->dof()->localToGlobal( context.eId(), i, c1 ).index();

                //N is the index of the RB basis function (i.e fem element)
                value_type rb_basisij = this->basisValue( rbN , gdof );

                //coefficient u_i^N * \PHI_ij
                value_type coefficient = u_i*rb_basisij;

                for ( size_type q = 0; q < context.xRefs().size2(); ++q )
                {
                    for ( typename array_type::index i = 0; i < nComponents1; ++i )
                    {
                        v[q]( i,0 ) += coefficient*context.d( ldof, i, N, q );
                    }
                }//q
            }//c1
        }//i

    }//rbN

}




template<typename SpaceType>
template<typename Y,  typename Cont>
void
ReducedBasisSpace<SpaceType>::Element<Y,Cont>::idInterpolate( matrix_node_type __ptsReal, id_array_type& v, bool conformalEval, matrix_node_type const& setPointsConf ) const
{
    bool go = false;
    CHECK( go ) << "The function idInterpolate is not yet implemented \n";
    LOG( INFO ) << "not yet implemented";
}

template<typename SpaceType>
inline
std::shared_ptr<ReducedBasisSpace<SpaceType> >
RbSpacePch(  std::shared_ptr<SpaceType> const& Xh )
{
    return ReducedBasisSpace<SpaceType>::New( Xh );
}

template<int Order, typename SpaceType>
inline
std::shared_ptr<ReducedBasisSpace<SpaceType>>
RbSpacePchv(  std::shared_ptr<SpaceType> const& Xh)
{
    return ReducedBasisSpace<SpaceType>::New( Xh );
}


}//namespace Feel

#endif
