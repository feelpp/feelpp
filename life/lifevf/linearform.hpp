/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-01-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006,2007 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file linearform.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-18
 */
#ifndef __LinearForm_H
#define __LinearForm_H 1

#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/multi_array.hpp>
#include <life/lifevf/block.hpp>
#include <life/lifealg/vectorvalue.hpp>

namespace Life
{
namespace fusion = boost::fusion;
namespace vf
{
/// \cond detail
namespace detail
{
/**
 * \class LinearForm
 * \brief Class that represents Linear forms
 *
 * \p LinearForm is the glue between the variational language and its
 * representation. In particular \p LinearForm is parametrized by:
 *
 * -# the function space \p SpaceType that describe the type of
 *    functions the \p LinearForm is applied to.
 *
 * -# the type of representation the use want to use for the
 *    linearform, by default it is a ublas vector
 *
 * -# the element container that is used to represent the linear form
 *    in a mixed form context, see \p MixLinearForm
 *
 * @author Christophe Prud'homme
 * @see \p BilinearForm
 */
template<typename SpaceType,typename VectorType,typename ElemContType>
class LinearForm
{
public:


    /** @name Typedefs
     */
    //@{
    enum { nDim = SpaceType::nDim };

    typedef LinearForm<SpaceType, VectorType, ElemContType> self_type;

    typedef SpaceType space_type;
    typedef boost::shared_ptr<SpaceType> space_ptrtype;

    typedef typename space_type::value_type value_type;
    typedef typename space_type::real_type real_type;
    typedef VectorType vector_type;
    typedef typename space_type::template Element<value_type, ElemContType> element_type;

#if 0
    typedef typename space_type::component_fespace_type component_fespace_type;
    typedef typename space_type::element_type::component_type component_type;
#endif
    typedef typename space_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type::permutation_type permutation_type;
    typedef typename space_type::gm_type gm_type;
    typedef typename space_type::gm_ptrtype gm_ptrtype;

    typedef typename space_type::fe_type fe_type;
    typedef typename space_type::basis_0_type::precompute_type test_precompute_type;

    typedef boost::shared_ptr<test_precompute_type> test_precompute_ptrtype;
    //@}

    /**
     * @name Inner Structures
     */
    //@{

    /**
     * \class Context
     * \brief element-wise representation of the bilinear form
     *
     * Local represents the linear form on a element (or a face of
     * an element) The algebraic representation of the local form is a
     * dense vector. This local representation shall eventually be
     * used later to fill global vectors, or in matrix-vector product
     * depending on the global linear form representation.
     *
     *
     */
    template<typename GeomapContext,
             typename ExprT,
             typename IM
             >
    class Context
    {
    public:

        typedef LinearForm<SpaceType,VectorType, ElemContType> form_type;
        typedef typename SpaceType::dof_type dof_type;
        typedef typename form_type::value_type value_type;



        typedef GeomapContext map_geometric_mapping_context_type;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type geometric_mapping_context_ptrtype;
        typedef typename geometric_mapping_context_ptrtype::element_type geometric_mapping_context_type;
        typedef typename geometric_mapping_context_type::gm_type geometric_mapping_type;


        static const uint16_type nDim = geometric_mapping_type::nDim;

        typedef ExprT expression_type;

        typedef typename space_type::mesh_type mesh_type;
        typedef typename mesh_type::element_type mesh_element_type;
        typedef typename mesh_element_type::permutation_type permutation_type;
        typedef typename space_type::fe_type test_basis_type;
        typedef boost::shared_ptr<test_basis_type> test_basis_ptrtype;
        typedef typename test_basis_type::template Context< geometric_mapping_context_type::context,
                                                            test_basis_type,
                                                            geometric_mapping_type,
                                                            mesh_element_type> test_basiscontext_type;

        typedef mpl::int_<fusion::result_of::template size<GeomapContext>::type::value> map_size;

        typedef boost::shared_ptr<test_basiscontext_type> test_basiscontext_ptrtype;


        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
                                  gmc<1>,
                                  gmc<0> >::type gmc1;

        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type left_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type::element_type left_gmc_type;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type right_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type::element_type right_gmc_type;

        typedef fusion::map<fusion::pair<gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
        typedef fusion::map<fusion::pair<gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;

        typedef typename mpl::if_<mpl::equal_to<map_size, mpl::int_<1> >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_basiscontext_ptrtype> > >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_basiscontext_ptrtype>,
                                                            fusion::pair<gmc<1>, test_basiscontext_ptrtype> > > >::type::type map_test_basiscontext_type;

        typedef fusion::map<fusion::pair<gmc<0>, test_basiscontext_ptrtype> > map_left_test_basiscontext_type;
        typedef fusion::map<fusion::pair<gmc<1>, test_basiscontext_ptrtype> > map_right_test_basiscontext_type;


        typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map_left_test_basiscontext_type> eval0_expr_type;
        typedef boost::shared_ptr<eval0_expr_type> eval0_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map_right_test_basiscontext_type> eval1_expr_type;
        typedef boost::shared_ptr<eval1_expr_type> eval1_expr_ptrtype;
        //typedef typename ExprT::template tensor<map_right_gmc_type, map_test_basiscontext_type> eval1_expr_type;


        typedef typename test_basis_type::template Context< geometric_mapping_context_type::context,
                                                            test_basis_type,
                                                            geometric_mapping_type,
                                                            mesh_element_type>::template Index<> test_index_type;


        //typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map0_test_basiscontext_type> eval0_expr_type;
        //typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map1_test_basiscontext_type> eval1_expr_type;

        //typedef ublas::vector<value_type> local_vector_type;
        static const int rep_shape = 2;//1+(eval_expr_type::shape::M-1>0)+(eval_expr_type::shape::N-1>0);
        typedef boost::multi_array<value_type, rep_shape> local_vector_type;

    private:

        struct InitTestFec
        {
            template<typename Sig>
            struct result;

            template<typename T>
            struct result<InitTestFec(T)>
            {
                typedef fusion::pair<typename boost::remove_reference<T>::type::first_type,test_basiscontext_ptrtype> type;
            };


            InitTestFec( test_basis_ptrtype const& fe, form_type const& form )
                :
                _M_fe( fe ),
                _M_form( form )
                {}
            template<typename T>
            fusion::pair<typename T::first_type,test_basiscontext_ptrtype>
            operator()( T const& t) const
            {
                geometric_mapping_context_ptrtype gmcptr( t.second );

                return fusion::make_pair<typename boost::remove_reference<T>::type::first_type>( test_basiscontext_ptrtype( new test_basiscontext_type( _M_fe,
                                                                                                                                                        gmcptr,
                                                                                                                                                        _M_form.testPc( gmcptr->faceId(), gmcptr->permutation() ) ) ) );
            };

            test_basis_ptrtype const& _M_fe;
            form_type const& _M_form;
        };

        struct UpdateTestFec
        {
            UpdateTestFec( map_geometric_mapping_context_type const& mapgmc,
                           form_type const& form )
                :
                _M_mapgmc( mapgmc ),
                _M_form( form )
                {}

            template<typename T>
            void operator()( T& t) const
            {
                geometric_mapping_context_ptrtype gmcptr( fusion::at_key<typename boost::remove_reference<T>::type::first_type>( _M_mapgmc ) );
                t.second->update( gmcptr, _M_form.testPc( gmcptr->faceId(), gmcptr->permutation()  ) );
            };

            map_geometric_mapping_context_type const& _M_mapgmc;
            form_type const& _M_form;
        };

        struct __integrate
        {
            __integrate( IM const& im )
                :
                _M_np( im.nPoints() ),
                _M_prod( im.nPoints() ),
                _M_weight( im.weights() ),
                _M_is_face_im( IM::is_face_im )
            {
            }
            __integrate&
            operator=( IM const& im )

            {
                _M_weight = im.weights();
                return *this;
            }

            void update( geometric_mapping_context_type const& _gmc )
            {
                if ( _M_is_face_im )
                    for( uint16_type q = 0; q < _M_np; ++q )
                        {
                            _M_prod[q] = _M_weight( q )*_gmc.J(q)*_gmc.normalNorm(q);
                        }
                else
                    for( uint16_type q = 0; q < _M_np; ++q )
                        {
                            _M_prod[q] = _M_weight( q )*_gmc.J(q);
                        }
            }
            template<typename ExprType>
            value_type operator()( ExprType const& expr,
                                   test_index_type  const& indi,
                                   uint16_type c1,
                                   uint16_type c2 ) const
            {
                value_type res = value_type(0);

                for( uint16_type q = 0; q < _M_np; ++q )
                    {
                        const value_type val_expr = expr.evaliq( indi, c1, c2, q );
                        res += _M_prod[q]*val_expr;

                    }
                return res;
            }
        private:
            uint16_type _M_np;
            std::vector<value_type> _M_prod;
            ublas::vector<value_type> _M_weight;
            const bool _M_is_face_im;
        };

    public:

        Context( form_type& __form,
                 map_geometric_mapping_context_type const& _gmc,
                 ExprT const& expr,
                 IM const& im );

        template<typename IM2>
        Context( form_type& __form,
                 map_geometric_mapping_context_type const& _gmc,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2 );

        template<typename IM2>
        Context( form_type& __form,
                 map_geometric_mapping_context_type const& _gmc,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2,
                 mpl::int_<2> );

        void update( map_geometric_mapping_context_type const& _gmc );

        void update( map_geometric_mapping_context_type const& _gmc, mpl::int_<2> );


        void update( map_geometric_mapping_context_type const& _gmc, IM const& im );

        void update( map_geometric_mapping_context_type const& _gmc, IM const& im, mpl::int_<2> );


        void integrate()
        {
            integrate(mpl::int_<fusion::result_of::size<GeomapContext>::type::value>() );
        }


        void assemble( size_type elt_0 );

        void assemble( size_type elt_0, size_type elt_1 );

    private:


        void integrate( mpl::int_<1> );

        void integrate( mpl::int_<2> );
    private:

        form_type& _M_form;
        dof_type* _M_test_dof;
        const list_block_type& _M_lb;
        map_geometric_mapping_context_type _M_gmc;
        left_gmc_ptrtype _M_gmc_left;
        right_gmc_ptrtype _M_gmc_right;
        map_left_gmc_type _M_left_map;
        map_right_gmc_type _M_right_map;
        map_test_basiscontext_type _M_test_fec;
        map_left_test_basiscontext_type _M_test_fec0;
        map_right_test_basiscontext_type _M_test_fec1;

        local_vector_type _M_rep;
        eval0_expr_ptrtype _M_eval0_expr;
        eval1_expr_ptrtype _M_eval1_expr;

        __integrate M_integrator;

    }; // Context


    //@}

    /** @name Constructors, destructor
     */
    //@{

    LinearForm( LinearForm const & __vf );

    LinearForm( space_ptrtype const& __X,
                vector_type& __F,
                bool init = true,
                bool do_threshold = false,
                value_type threshold = type_traits<value_type>::epsilon() );
    LinearForm( space_ptrtype const& __X,
                vector_type& __F,
                list_block_type const& __lb,
                bool init = true,
                bool do_threshold = false,
                value_type threshold = type_traits<value_type>::epsilon()  );

    ~LinearForm()
        {}



    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * Construct the linear form given by the expression \p expr
     *
     * This function typically fills the representation that was used
     * for the \p LinearForm.
     */
    template <class ExprT>
    LinearForm& operator=( Expr<ExprT> const& expr);

    /**
     * Construct the linear form given by the expression \p expr
     *
     * This function typically fills the representation that was used
     * for the \p LinearForm by adding to it the evaluation of the
     * expression \p expr.
     */
    template <class ExprT>
    LinearForm& operator+=( Expr<ExprT> const& expr);

#if 0
    /**
     * \todo write documentation
     * \see MixedLinearForm
     */
    LinearForm<component_fespace_type,
               vector_type,
               typename component_type::container_type>
    operator()( component_type& __u )
        {
            typedef typename component_type::container_type container_type;

            typedef LinearForm<component_fespace_type, vector_type,
                               container_type> vflf_type;

            list_block_type __list_block;
            __list_block.push_back( Block( 0, 0,
                                           __u.start(), 0 ) );
            return vflf_type( __u, _M_F, __list_block, false );

        }
#endif

    /**
     * \return the entry \f$M_{i,j}\f$
     */
    value_type operator()( size_type i ) const { return _M_F( i ); }

    //@}

    /** @name Accessors
     */
    //@{

    /**
    * \return the test function space
     */
    space_ptrtype const& functionSpace() const { return _M_X; }

    /**
     * Geometric transformation
     */
    gm_ptrtype const& gm() const { return _M_X->gm(); }

   /**
     * Return the structure that holds the test basis functions
     * evaluated at a previously given set of points on a face of the
     * reference element
     * \see precomputeBasisAtPoints()
     */
    test_precompute_ptrtype const& testPc( uint16_type __f,
                                           permutation_type __p = permutation_type( permutation_type::NO_PERMUTATION ) ) const
    {
        if ( __f == invalid_uint16_type_value )
            return  _M_test_pc;
        return _M_test_pc_face.find( __f )->second.find( __p )->second;
    }

    vector_type& representation() const { return _M_F; }

    list_block_type const& blockList() const { return _M_lb; }

    /**
     * \return \c true if threshold applies, false otherwise
     */
    bool doThreshold( value_type const& v ) const { return ( math::abs( v ) > _M_threshold ); }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * Set the function space from which the linear form takes its
     * value from.
     */
    void setFunctionSpace( space_ptrtype const& __X )
    {
        _M_X = __X;
    }

    /**
     * set a threshold value for the matrix entries associated with the
     * bilinear form
     */
    void setThreshold( value_type eps ) { _M_threshold = eps; }

    /**
     * set the threshold strategy, true threshold the matrix entries,
     * false do not threshold
     */
    void setDoThreshold( bool do_threshold ) { _M_do_threshold = do_threshold; }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * precompute the basis function associated with the test and
     * trial space at a set of points
     */
    template<typename Pts>
    void precomputeBasisAtPoints( Pts const& pts )
    {
        _M_test_pc = test_precompute_ptrtype( new test_precompute_type( functionSpace()->fe(), pts ) );
    }

   /**
     * precompute the basis function associated with the test and
     * trial space at a set of points on a face of the reference
     * element
     */
    template<typename Pts>
    void precomputeBasisAtPoints( uint16_type __f, permutation_type __p, Pts const& pts )
    {
        _M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( functionSpace()->fe(), pts ) );
    }

    /**
     * add value \p v at position (\p i) of the vector
     * associated with the linear form
     */
    void add( size_type i,  value_type const& v )
    {
        if ( _M_do_threshold )
            {
                if ( doThreshold( v ) )
                    _M_F.add( i, v );
            }
        else
            _M_F.add( i, v );
    }

    /**
     * set value \p v at position (\p i) of the vector
     * associated with the linear form
     */
    void set( size_type i,  value_type const& v ) { _M_F.set( i, v ); }


    //@}

protected:

    LinearForm();

private:
    template <class ExprT> void assign( Expr<ExprT> const& expr, bool init, mpl::bool_<false> );
    template <class ExprT> void assign( Expr<ExprT> const& expr, bool init, mpl::bool_<true> );
private:

    space_ptrtype _M_X;

    vector_type& _M_F;

    list_block_type _M_lb;

    test_precompute_ptrtype _M_test_pc;

    std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > _M_test_pc_face;

    bool _M_do_threshold;
    value_type _M_threshold;
};

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( LinearForm const & __vf )
    :
    _M_X( __vf._M_X ),
    _M_F( __vf._M_F ),
    _M_lb( __vf._M_lb ),
    _M_test_pc(),
    _M_test_pc_face(),
    _M_do_threshold( __vf._M_do_threshold ),
    _M_threshold( __vf._M_threshold )

{
    Debug( 5060 ) << "LinearForm copy constructor\n";
    Debug( 5060 ) << "     n Dof : " << _M_X->nDof() << "\n";
    Debug( 5060 ) << "    F size : " << _M_F.size() << "\n";
    Debug( 5060 ) << "block size : " << _M_lb.size() << "\n";
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( space_ptrtype const& __X,
                                                             vector_type& __F,
                                                             bool init,
                                                             bool do_threshold,
                                                             value_type threshold  )
    :
    _M_X( __X ),
    _M_F( __F ),
    _M_lb(),
    _M_do_threshold( do_threshold ),
    _M_threshold( threshold )
{

    for( uint16_type __i = 0;__i < _M_X->qDim(); ++__i)
    {
        _M_lb.push_back( Block( __i, 0,
                                __i*_M_X->nDofPerComponent(),
                                0 ));
        Debug( 5050 ) << "[linearform::linearform] block: "
                      << Block( __i, 0, __i*_M_X->nDofPerComponent(), 0 )  << "\n";
    }

    if (  init )
        _M_F.zero();
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( space_ptrtype const& __X,
                                                             vector_type& __F,
                                                             list_block_type const& __lb,
                                                             bool init,
                                                             bool do_threshold,
                                                             value_type threshold )
    :
    _M_X( __X ),
    _M_F( __F ),
    _M_lb( __lb ),
    _M_do_threshold( do_threshold ),
    _M_threshold( threshold )
{
    if ( init )
        _M_F.zero();
}

template<typename LFType, typename TheSpaceType, typename ExprType>
struct LFAssign
{
    LFAssign( LFAssign const& lfa )
        :
        _M_lf( lfa._M_lf ),
        _M_Xh( lfa._M_Xh ),
        _M_expr( lfa._M_expr ),
        _M_index( lfa._M_index ),
        _M_init( lfa._M_init )
    {}
    LFAssign( LFType& lf, TheSpaceType const& Xh, ExprType const& expr, bool init )
        :
        _M_lf( lf ),
        _M_Xh( Xh ),
        _M_expr( expr ),
        _M_index( 0 ),
        _M_init( init)
    {}
    template<typename SpaceType>
    void operator()( boost::shared_ptr<SpaceType> const& X ) const
    {
        Debug(5050) << "expression has test functions ? :"
                    << ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result
                    << "\n";
        if ( !ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result )
            {
                ++_M_index;
                return;
            }

        list_block_type __list_block;
        __list_block.push_back( Block( 0, 0, _M_Xh->nDofStart( _M_index ), 0 ) );
        LinearForm<SpaceType,typename LFType::vector_type, typename LFType::element_type> lf( X,
                                                                                              _M_lf.representation(),
                                                                                              __list_block,
                                                                                              false );
        //
        // in composite integration, make sure that if _M_init is \p
        // true for the first space, it is set to \p false for the
        // next spaces otherwise it will erase/clear to 0 the previous
        // assemblies
        //
        if ( _M_init )
            {
                // assembly
                lf = _M_expr;

                // make sure we won't erase the assembly we just did
                _M_init = false;
            }
        else
            {
                lf += _M_expr;
            }
        ++_M_index;
    }
private:
    LFType& _M_lf;
    TheSpaceType const& _M_Xh;
    ExprType const& _M_expr;
    mutable size_type _M_index;
    mutable bool _M_init;
};
template<typename LFType, typename TheSpaceType, typename ExprType>
LFAssign<LFType,TheSpaceType,ExprType>
make_lfassign( LFType& lf, TheSpaceType const& Xh, ExprType const& expr, bool init )
{ return LFAssign<LFType,TheSpaceType,ExprType>( lf, Xh, expr, init ); }



template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
void
LinearForm<SpaceType, VectorType, ElemContType>::assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<true> )
{
    fusion::for_each( _M_X->functionSpaces(), make_lfassign( *this, _M_X, __expr, init ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
void
LinearForm<SpaceType, VectorType, ElemContType>::assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<false> )
{
    if ( _M_lb.empty() )
        {
            _M_lb.push_back( Block( 0, 0, 0, 0 ) );
        }
    if ( init )
        {
            typename list_block_type::const_iterator __bit = _M_lb.begin();
            typename list_block_type::const_iterator __ben = _M_lb.end();
            for ( ; __bit != __ben; ++__bit )
                {
                    Debug( 5050 ) << "LinearForm:: block: " << *__bit << "\n";
                    size_type g_ic_start = __bit->globalRowStart();
                    Debug( 5050 ) << "LinearForm:: g_ic_start: " << g_ic_start << "\n";

                    _M_F.zero( g_ic_start,g_ic_start + _M_X->nDof() );
                }
        }
    __expr.assemble( _M_X, *this );
    //vector_range_type r( _M_F, ublas::range( _M_v.start(),_M_v.start()+ _M_v.size() ) );
    //std::cout << "r = " << r << "\n";
    //std::cout << "after F= " << _M_F << "\n";
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
LinearForm<SpaceType, VectorType, ElemContType>&
LinearForm<SpaceType, VectorType, ElemContType>::operator=( Expr<ExprT> const& __expr )
{
    // loop(fusion::for_each) over sub-functionspaces in SpaceType
    // pass expression and initialize
    this->assign( __expr, true, mpl::bool_<(SpaceType::nSpaces > 1)>() );
    return *this;
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
LinearForm<SpaceType, VectorType, ElemContType>&
LinearForm<SpaceType, VectorType, ElemContType>::operator+=( Expr<ExprT> const& __expr )
{
    this->assign( __expr, false, mpl::bool_<(SpaceType::nSpaces > 1)>() );
    return *this;
}




//
// Context class for linear forms
//
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im )
    :
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_test_fec( fusion::transform( _M_gmc, InitTestFec(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_basiscontext_type::nDof]
           [test_basiscontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr(),
    M_integrator( im )
{
    _M_eval0_expr->init( im );
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
template<typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im,
                                                                                           IM2 const& im2 )
    :
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_test_fec( fusion::transform( _M_gmc, InitTestFec(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_basiscontext_type::nDof]
           [test_basiscontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr(),
    M_integrator( im )
{
    _M_eval0_expr->init( im2 );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT, typename IM>
template<typename IM2>
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                           map_geometric_mapping_context_type const& _gmc,
                                                                                           ExprT const& expr,
                                                                                           IM const& im,
                                                                                           IM2 const& im2,
                                                                                           mpl::int_<2> )
    :
    _M_form( __form ),
    _M_test_dof( __form.functionSpace()->dof().get() ),
    _M_lb( __form.blockList() ),
    _M_gmc( _gmc ),
    _M_gmc_left( fusion::at_key<gmc<0> >( _gmc ) ),
    _M_gmc_right( fusion::at_key<gmc1 >( _gmc ) ),
    _M_left_map( fusion::make_map<gmc<0> >( _M_gmc_left ) ),
    _M_right_map( fusion::make_map<gmc<0> >( _M_gmc_right ) ),
    _M_test_fec( fusion::transform( _M_gmc, InitTestFec(__form.functionSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_test_fec1( fusion::make_pair<gmc<1> >( fusion::at_key<gmc<1> >( _M_test_fec ) ) ),
    _M_rep(boost::extents[fusion::size(_gmc)*test_basiscontext_type::nDof]
           [test_basiscontext_type::nComponents1]),
    _M_eval0_expr( new eval0_expr_type( expr, _gmc, _M_test_fec0 ) ),
    _M_eval1_expr( new eval1_expr_type( expr, _gmc, _M_test_fec1 ) ),
    M_integrator( im )
{
    _M_eval0_expr->init( im2 );
    _M_eval1_expr->init( im2 );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc )
{
    _M_gmc = _gmc;
    _M_gmc_left = fusion::at_key<gmc<0> >( _gmc );
    _M_left_map = fusion::make_map<gmc<0> >( _M_gmc_left );
    fusion::for_each( _M_test_fec, UpdateTestFec( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval0_expr->update( _gmc, _M_test_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::int_<2> )
{
    typedef mpl::int_<fusion::result_of::template size<map_geometric_mapping_context_type>::type::value> map_size;
    BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,map_geometric_mapping_context_type ));
    _M_gmc = _gmc;
#if 0
    _M_gmc_left = fusion::at_key<gmc<0> >( _gmc );
    _M_gmc_right =  fusion::at_key<gmc1 >( _gmc );
    _M_left_map = fusion::make_map<gmc<0> >( _M_gmc_left );
    _M_right_map = fusion::make_map<gmc<0> >( _M_gmc_right );
#endif
    fusion::for_each( _M_test_fec, UpdateTestFec( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_test_fec1 = fusion::make_map<gmc<1> >( fusion::at_key<gmc<1> >( _M_test_fec ) );
    _M_eval0_expr->update( _gmc, _M_test_fec0 );
    _M_eval1_expr->update( _gmc, _M_test_fec1 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, IM const& im )
{
    M_integrator = im;
    this->update( _gmc );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, IM const& im, mpl::int_<2> )
{
    M_integrator = im;
    this->update( _gmc, mpl::int_<2>() );
}


template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<1> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    typedef geometric_mapping_context_type gmc_type;

    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( (shape::M == 1 && shape::N == 1),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

    test_index_type indi;
    int i, c1;
    value_type res;
    //#pragma omp parallel
    //{
        //#pragma omp single
        //Debug() << "[linearform::integrate] num threads: " << OMP_GET_NUM_THREADS << "\n";

    //#pragma omp for private(i,c1,indi, res)
        for ( i = 0; i < test_basiscontext_type::nDof; ++i )
            for( c1 = 0;c1 < test_basiscontext_type::nComponents1; ++ c1 )
                {
                    indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                    res = M_integrator( *_M_eval0_expr, indi, 0, 0 );
                    _M_rep[i][c1] = res;
                }
        //}
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<2> )
{
    //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    typedef geometric_mapping_context_type gmc_type;

    typedef mpl::int_<fusion::result_of::template size<map_geometric_mapping_context_type>::type::value> map_size;
    BOOST_MPL_ASSERT_MSG( map_size::value == 2, INVALID_GEOMAP, (map_size,map_geometric_mapping_context_type ));

    typedef typename eval0_expr_type::shape shape;
    BOOST_MPL_ASSERT_MSG( (shape::M == 1 && shape::N == 1),
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );
    test_index_type indi;

    for ( uint16_type i = 0; i < test_basiscontext_type::nDof; ++i )
        for( uint16_type c1 = 0;c1 < test_basiscontext_type::nComponents1; ++ c1 )
            {
                indi.setIndex( boost::make_tuple( i, c1, 0 ) );

                BOOST_MPL_ASSERT( IM::is_face_im );
                _M_rep[i][c1] = M_integrator( *_M_eval0_expr, indi, 0, 0 );
                uint16_type ii = i + test_basiscontext_type::nDof;
                _M_rep[ii][c1] = M_integrator( *_M_eval1_expr, indi, 0, 0 );
            }
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0 )
{
    test_index_type indi;
    size_type ig;
    int isign;
    size_type row_start = _M_lb.front().globalRowStart();
    for ( uint16_type k = 0 ; k < test_basiscontext_type::nDof; k++ )
        for( uint16_type c1 = 0;c1 < test_basiscontext_type::nComponents1; ++ c1 )
            {
                //uint16_type c = test_basiscontext_type::nComponents2*c1+c2;
                boost::tie( ig, isign, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, k, c1 );
                ig += row_start;
                _M_form.add( ig, value_type(isign)*_M_rep[k][c1] );
            }
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
LinearForm<SpaceType, VectorType, ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0, size_type elt_1 )
{
    test_index_type indi;
    size_type ig0,ig1;
    int isign0,isign1;
    size_type row_start = _M_lb.front().globalRowStart();
    for ( uint16_type k = 0 ; k < test_basiscontext_type::nDof; k++ )
        for( uint16_type c1 = 0;c1 < test_basiscontext_type::nComponents1; ++ c1 )
            //for( uint16_type c2 = 0;c2 < test_basiscontext_type::nComponents2; ++ c2 )
            {
                //uint16_type c = test_basiscontext_type::nComponents2*c1+c2;
                boost::tie( ig0, isign0, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, k, c1 );
                ig0 += row_start;
                _M_form.add( ig0, value_type(isign0)*_M_rep[k][c1] );

                boost::tie( ig1, isign1, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_1, k, c1 );
                ig1 += row_start;
                size_type l = test_basiscontext_type::nDof + k;
                _M_form.add( ig1, value_type(isign1)*_M_rep[l][c1] );
            }
}

template<typename X1, typename IntElts, typename Im, typename ExprT>
typename Expr<Integrator<IntElts, Im, ExprT> >::value_type
integrate( boost::shared_ptr<X1> const& __X1,
           IntElts const& elts,
           Im const& __im,
           ExprT const& __expr )
{
    VectorValue<typename Expr<Integrator<IntElts, Im, ExprT> >::value_type> M;
    form( __X1, M ) = integrate( elts, __im, __expr );
    return M.vec();
}

} // detail
/// \endcond
} // vf
} // life

#endif /* __LinearForm_H */

