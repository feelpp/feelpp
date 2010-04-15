/* -*- mode: c++ -*-

   This file is part of the Life library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2005-01-18

   Copyright (C) 2005,2006 EPFL
   Copyright (C) 2006,2007,2008 Université Joseph Fourier (Grenoble I)

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
 *  \file bilinearform.hpp
 *  \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
 *  \date 2005-01-18
 */
#ifndef __BilinearForm_H
#define __BilinearForm_H 1

#include <set>

#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/spirit/home/phoenix.hpp>
#include <boost/spirit/home/phoenix/core/argument.hpp>
#include <life/lifecore/context.hpp>
#include <life/lifealg/matrixvalue.hpp>
#include <life/lifealg/vectorublas.hpp>
#include <life/lifealg/graphcsr.hpp>
#include <life/lifevf/block.hpp>
#include <life/lifevf/fec.hpp>


namespace Life
{
namespace parameter = boost::parameter;
namespace fusion = boost::fusion;
namespace vf
{
enum DofGraph
    {
        DOF_PATTERN_COUPLED  = 1 << 0,
        DOF_PATTERN_NEIGHBOR = 1 << 1
    };

/// \cond detail
template<typename FE1,typename FE2,typename ElemContType> class BilinearForm;
namespace detail
{
template<typename BFType, typename Space1Type>
struct compute_graph2
{
    compute_graph2( BFType& bf, boost::shared_ptr<Space1Type> const& space1, size_type test_index, size_type hints  )
        :
        M_bf( bf ),
        M_space1( space1 ),
        M_test_index( test_index ),
        M_trial_index( 0 ),
        M_hints( hints )
    {}

    template <typename Space2>
    void operator()( boost::shared_ptr<Space2> const& space2 ) const;

    mutable BFType& M_bf;
    boost::shared_ptr<Space1Type> const& M_space1;
    size_type M_test_index;
    mutable size_type M_trial_index;
    size_type M_hints;
};


template<typename BFType>
struct compute_graph1
{
    compute_graph1( BFType& bf, size_type hints )
        :
        M_bf( bf ),
        M_test_index( 0 ),
        M_hints( hints )
    {}

    template <typename Space1>
    void operator()( boost::shared_ptr<Space1> const& space1 ) const
    {
        fusion::for_each( M_bf.trialSpace()->functionSpaces(),
                          compute_graph2<BFType,Space1>( M_bf, space1, M_test_index, M_hints ) );
        ++M_test_index;
    }
    mutable BFType& M_bf;
    mutable size_type M_test_index;
    size_type M_hints;
};


template<typename BFType, typename ExprType>
struct BFAssign2
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign2( BFAssign2 const& lfa )
        :
        _M_bf( lfa._M_bf ),
        _M_expr( lfa._M_expr ),
        _M_test_index( lfa._M_test_index )
    {}
    BFAssign2( BFType& lf, ExprType const& expr )
        :
        _M_bf( lf ),
        _M_expr( expr ),
        _M_test_index( 0 )
    {}
    template<typename SpaceType>
    void operator()( boost::shared_ptr<SpaceType> const& X ) const
    {
        Debug( 5050 ) << "[BFAssign2::operator()] start loop on trial spaces against test space index: " << _M_test_index << "\n";
        fusion::for_each( _M_bf.trialSpace()->functionSpaces(),
                          make_bfassign1( _M_bf, _M_expr, X, _M_test_index ) );
        Debug( 5050 ) << "[BFAssign2::operator()] stop loop on trial spaces against test space index: " << _M_test_index << "\n";
        ++_M_test_index;

    }
private:
    BFType& _M_bf;
    ExprType const& _M_expr;
    mutable size_type _M_test_index;
};
template<typename BFType, typename ExprType>
BFAssign2<BFType,ExprType>
make_bfassign2( BFType& lf, ExprType const& expr )
{ return BFAssign2<BFType,ExprType>( lf, expr ); }

template<typename BFType, typename ExprType, typename TestSpaceType>
struct BFAssign1
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign1( BFAssign1 const& lfa )
        :
        _M_bf( lfa._M_bf ),
        _M_test( lfa._M_test ),
        _M_expr( lfa._M_expr ),
        _M_trial_index( lfa._M_trial_index )
    {}
    BFAssign1( BFType& lf,
               ExprType const& expr,
               boost::shared_ptr<TestSpaceType> const& Testh,
               size_type test_index )
        :
        _M_bf( lf ),
        _M_test( Testh ),
        _M_expr( expr ),
        _M_trial_index( 0 ),
        _M_test_index( test_index )
    {}
    template<typename SpaceType>
    void operator()( boost::shared_ptr<SpaceType> const& trial ) const;

private:
    BFType& _M_bf;
    boost::shared_ptr<TestSpaceType> _M_test;
    ExprType const& _M_expr;
    mutable size_type _M_trial_index;
    size_type _M_test_index;
};
template<typename BFType, typename ExprType, typename TestSpaceType>
BFAssign1<BFType,ExprType,TestSpaceType>
make_bfassign1( BFType& lf,
                ExprType const& expr,
                boost::shared_ptr<TestSpaceType> const& test_space,
                size_type test_index )
{ return BFAssign1<BFType,ExprType,TestSpaceType>( lf, expr, test_space, test_index ); }


/*!
  \class BilinearForm
  \brief brief description

  @author Christophe Prud'homme
  @see
*/
template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename FE1::value_type> >
class BilinearForm
{
public:


    /** @name Typedefs
     */
    //@{
    enum { nDim = FE1::nDim };

    typedef FE1 space_1_type;
    typedef boost::shared_ptr<FE1> space_1_ptrtype;
    typedef space_1_type test_space_type;
    typedef boost::shared_ptr<space_1_type> test_space_ptrtype;

    typedef FE2 space_2_type;
    typedef boost::shared_ptr<FE2> space_2_ptrtype;
    typedef space_2_type trial_space_type;
    typedef boost::shared_ptr<space_2_type> trial_space_ptrtype;

    typedef typename FE1::value_type value_type;
    typedef typename space_1_type::template Element<value_type,ElemContType> element_1_type;

    typedef typename space_2_type::template Element<value_type,ElemContType> element_2_type;
    typedef BilinearForm<FE1, FE2, ElemContType> self_type;

#if 0
    typedef typename space_1_type::component_fespace_type component_space_1_type;
    typedef typename element_1_type::component_type component_1_type;
    typedef typename space_2_type::component_fespace_type component_space_2_type;
    typedef typename element_2_type::component_type component_2_type;
#endif // 0

    typedef typename space_1_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type mesh_element_type;
    typedef typename mesh_element_type::permutation_type permutation_type;
    typedef typename space_1_type::fe_type fe_1_type;
    typedef typename space_2_type::fe_type fe_2_type;

    typedef typename space_1_type::gm_type gm_type;
    typedef typename space_1_type::gm_ptrtype gm_ptrtype;

    //typedef ublas::compressed_matrix<value_type, ublas::row_major> csr_matrix_type;
    typedef MatrixSparse<value_type> matrix_type;
    static const bool is_row_major = true;//matrix_type::is_row_major;

    typedef typename mpl::if_<mpl::equal_to<mpl::bool_<is_row_major>, mpl::bool_<true> >,
                              mpl::identity<ublas::row_major>,
                              mpl::identity<ublas::column_major> >::type::type layout_type;

    typedef typename space_1_type::basis_0_type::precompute_type test_precompute_type;
    typedef boost::shared_ptr<test_precompute_type> test_precompute_ptrtype;

    typedef typename space_2_type::basis_0_type::precompute_type trial_precompute_type;
    typedef boost::shared_ptr<trial_precompute_type> trial_precompute_ptrtype;

    typedef GraphCSR graph_type;
    typedef boost::shared_ptr<graph_type> graph_ptrtype;
    //@}


    /**
     * @name Inner Structures
     */
    //@{

    /**
     * \class Context
     * \brief element-wise representation of the bilinear form
     *
     * Local represents the bilinear form on a element (or a face of
     * an element) The algebraic representation of the local form is a
     * dense matrix. This local representation shall eventually be
     * used later to fill global matrices, or in matrix-vector product
     * depending on the global bilinear form representation.
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
        typedef Context<GeomapContext,ExprT,IM> form_context_type;
        typedef BilinearForm<FE1, FE2, ElemContType> form_type;
        typedef typename FE1::dof_type dof_1_type;
        typedef typename FE2::dof_type dof_2_type;

        typedef typename form_type::value_type value_type;



        typedef GeomapContext map_geometric_mapping_context_type;

        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type geometric_mapping_context_ptrtype;

        typedef typename geometric_mapping_context_ptrtype::element_type geometric_mapping_context_type;
        typedef typename geometric_mapping_context_type::gm_type geometric_mapping_type;

        static const uint16_type nDim = geometric_mapping_type::nDim;

        typedef ExprT expression_type;


        typedef mpl::int_<fusion::result_of::template size<GeomapContext>::type::value> map_size;

        typedef typename FE2::fe_type trial_fe_type;
        typedef boost::shared_ptr<trial_fe_type> trial_fe_ptrtype;
        typedef typename trial_fe_type::template Context< geometric_mapping_context_type::context,
                                                          trial_fe_type,
                                                          geometric_mapping_type,
                                                          mesh_element_type> trial_fecontext_type;
        typedef boost::shared_ptr<trial_fecontext_type> trial_fecontext_ptrtype;
        typedef typename mpl::if_<mpl::equal_to<map_size, mpl::int_<1> >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype>,
                                                            fusion::pair<gmc<1>, trial_fecontext_ptrtype> > > >::type::type map_trial_fecontext_type;

        typedef fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > map_left_trial_fecontext_type;
        typedef fusion::map<fusion::pair<gmc<1>, trial_fecontext_ptrtype> > map_right_trial_fecontext_type;
        typedef typename trial_fe_type::template Context< geometric_mapping_context_type::context,
                                                          trial_fe_type,
                                                          geometric_mapping_type,
                                                          mesh_element_type>::template Index<> trial_index_type;
        typedef typename FE1::fe_type test_fe_type;
        typedef boost::shared_ptr<test_fe_type> test_fe_ptrtype;
        typedef typename test_fe_type::template Context< geometric_mapping_context_type::context,
                                                         test_fe_type,
                                                         geometric_mapping_type,
                                                         mesh_element_type> test_fecontext_type;
        typedef boost::shared_ptr<test_fecontext_type> test_fecontext_ptrtype;
        typedef typename mpl::if_<mpl::equal_to<map_size, mpl::int_<1> >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype>,
                                                            fusion::pair<gmc<1>, test_fecontext_ptrtype> > > >::type::type map_test_fecontext_type;
        typedef fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > map_left_test_fecontext_type;
        typedef fusion::map<fusion::pair<gmc<1>, test_fecontext_ptrtype> > map_right_test_fecontext_type;

        typedef typename test_fe_type::template Context< geometric_mapping_context_type::context,
                                                         test_fe_type,
                                                         geometric_mapping_type,
                                                         mesh_element_type>::template Index<> test_index_type;

        typedef typename ExprT::template tensor<map_geometric_mapping_context_type,
                                               map_left_test_fecontext_type,
                                               map_left_trial_fecontext_type> eval00_expr_type;
        typedef boost::shared_ptr<eval00_expr_type> eval00_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_context_type,
                                               map_left_test_fecontext_type,
                                               map_right_trial_fecontext_type> eval01_expr_type;
        typedef boost::shared_ptr<eval01_expr_type> eval01_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_context_type,
                                               map_right_test_fecontext_type,
                                               map_left_trial_fecontext_type> eval10_expr_type;
        typedef boost::shared_ptr<eval10_expr_type> eval10_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_context_type,
                                               map_right_test_fecontext_type,
                                               map_right_trial_fecontext_type> eval11_expr_type;
        typedef boost::shared_ptr<eval11_expr_type> eval11_expr_ptrtype;

        static const int rep_shape = 4;//2+(eval_expr_type::shape::M-1>0)+(eval_expr_type::shape::N-1>0);
        typedef boost::multi_array<value_type, rep_shape> local_matrix_type;

    public:

        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  left, Right right )
        { return getMap( left, right, boost::is_same<Left, map_trial_fecontext_type>() ); }
        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  left, Right /*right*/, mpl::bool_<true> ) { return left; }
        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  /*left*/, Right right, mpl::bool_<false> ) { return map_trial_fecontext_type(right); }

        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  left, Right right )
        { return getMap( left, right, boost::is_same<Left, map_left_trial_fecontext_type>() ); }
        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  left, Right /*right*/, mpl::bool_<true> ) { return left; }
        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  /*left*/, Right right, mpl::bool_<false> ) { return right; }

        Context( form_type& __form,
                 map_geometric_mapping_context_type const& gmc,
                 ExprT const& expr,
                 IM const& im );

        template<typename IM2>
        Context( form_type& __form,
                 map_geometric_mapping_context_type const& gmc,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2 );

        template<typename IM2>
        Context( form_type& __form,
                 map_geometric_mapping_context_type const& gmc,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2,
                 mpl::int_<2> );

        void update( map_geometric_mapping_context_type const& gmc );

        void update( map_geometric_mapping_context_type const& gmc, mpl::int_<2> );

        void update( map_geometric_mapping_context_type const& gmc, IM const& im )
        {
            M_integrator = im;
            update( gmc );
        }

        void update( map_geometric_mapping_context_type const& gmc, IM const& im, mpl::int_<2> )
        {
            M_integrator = im;
            update( gmc, mpl::int_<2>() );
        }

        void integrate()
        {
            integrate( mpl::int_<fusion::result_of::size<GeomapContext>::type::value>() );
        }



        void assemble( size_type elt_0 );

        void assemble( size_type elt_0, size_type elt_1  );

    private:
        void update( map_geometric_mapping_context_type const& gmc, mpl::bool_<false> );

        void update( map_geometric_mapping_context_type const& gmc, mpl::bool_<true> );

        void integrate( mpl::int_<1> );

        void integrate( mpl::int_<2> );
    private:

        form_type& _M_form;
        const list_block_type& _M_lb;
        dof_1_type* _M_test_dof;
        dof_2_type* _M_trial_dof;
        map_geometric_mapping_context_type _M_gmc;
        map_test_fecontext_type _M_test_fec;
        map_left_test_fecontext_type _M_test_fec0;
        map_right_test_fecontext_type _M_test_fec1;
        map_trial_fecontext_type _M_trial_fec;
        map_left_trial_fecontext_type _M_trial_fec0;
        map_right_trial_fecontext_type _M_trial_fec1;


        local_matrix_type _M_rep;
        eval00_expr_ptrtype _M_eval_expr00;
        eval01_expr_ptrtype _M_eval_expr01;
        eval10_expr_ptrtype _M_eval_expr10;
        eval11_expr_ptrtype _M_eval_expr11;

        IM M_integrator;

        test_index_type indi;
        trial_index_type indj;

        }; // Context

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**

     */
    BilinearForm( space_1_ptrtype const& __X1,
                  space_2_ptrtype const& __X2,
                  matrix_type& __M,
                  bool build = true,
                  bool do_threshold = false,
                  value_type threshold = type_traits<value_type>::epsilon(),
                  size_type graph_hints = DOF_PATTERN_COUPLED );

    BilinearForm( space_1_ptrtype const& __X1,
                  space_2_ptrtype const& __X2,
                  matrix_type& __M,
                  list_block_type const& __lb,
                  bool do_threshold = false,
                  value_type threshold = type_traits<value_type>::epsilon(),
                  size_type graph_hints = DOF_PATTERN_COUPLED );

    BilinearForm( BilinearForm const& __vf )
        :
        _M_X1( __vf._M_X1 ),
        _M_X2( __vf._M_X2 ),
        _M_matrix( __vf._M_matrix ),
        _M_lb( __vf._M_lb ),
        _M_test_pc( __vf._M_test_pc ),
        _M_test_pc_face( __vf._M_test_pc_face ),
        _M_trial_pc( __vf._M_trial_pc ),
        _M_trial_pc_face( __vf._M_trial_pc_face ),
        _M_do_threshold( __vf._M_do_threshold ),
        _M_threshold( __vf._M_threshold ),
        M_graph( __vf.M_graph )
    {}

    ~BilinearForm()
    {}



    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * copy operator
     */
    BilinearForm&
    operator=( BilinearForm const& form )
    {
        if ( this != &form )
        {
            _M_X1 = form._M_X1;
            _M_X2 = form._M_X2;
            _M_matrix = form._M_matrix;
            _M_lb = form._M_lb;
            M_graph = form.M_graph;
        }
        return *this;
    }

    /**
     * expression assignment
     */
    template <class ExprT>
    BilinearForm& operator=( Expr<ExprT> const& expr );

    /**
     * plus expression assignment
     */
    template <class ExprT>
    BilinearForm& operator+=( Expr<ExprT> const& expr );

    /**
     * Computes the energy norm associated with the bilinear form
     *
     * @param __v element of Space 1 (test space)
     * @param __u element of Space 2 (trial space)
     * @return the energy norm
     */
    value_type operator()( element_1_type const& __v,  element_2_type const& __u ) const
    {
        return _M_matrix.energy( __v, __u );
    }

    /**
     * \return the entry \f$M_{i,j}\f$
     */
    // value_type& operator()( size_type i, size_type j ) { return _M_matrix( i, j ); }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return the trial function space
     */
    space_2_ptrtype const& trialSpace() const { return _M_X2; }

    /**
     * \return the test function space
     */
    space_1_ptrtype const& testSpace() const { return _M_X1; }

    /**
     * Geometric transformation
     */
    gm_ptrtype const& gm() const { return _M_X1->gm(); }

    /**
     * Return the structure that holds the test basis functions
     * evaluated at a previously given set of points on the face of
     * the reference element
     *
     * \see precomputeBasisAtPoints()
     */
    test_precompute_ptrtype const& testPc( uint16_type __f,
                                           permutation_type __p = permutation_type( permutation_type::NO_PERMUTATION ) ) const
    {
        if ( __f == invalid_uint16_type_value )
            return  _M_test_pc;
        //LIFE_ASSERT( _M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );
        return _M_test_pc_face.find(__f )->second.find( __p )->second;
    }

    /**
     * Return the structure that holds the trial basis functions
     * evaluated at a previously given set of points on the face of
     * the reference element
     * \see precomputeBasisAtPoints()
     */
    trial_precompute_ptrtype const& trialPc( uint16_type __f,
                                             permutation_type __p = permutation_type( permutation_type::NO_PERMUTATION ) ) const
    {
        if ( __f == invalid_uint16_type_value )
            return  _M_trial_pc;
        //LIFE_ASSERT( _M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
        return _M_trial_pc_face.find(__f )->second.find( __p )->second;
    }

    /**
     * \return the matrix associated to the bilinear form
     */
    matrix_type const& matrix() const { return _M_matrix; }

    matrix_type& matrix() { return _M_matrix; }

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
        _M_test_pc = test_precompute_ptrtype( new test_precompute_type( testSpace()->fe(), pts ) );
        _M_trial_pc = trial_precompute_ptrtype( new trial_precompute_type( trialSpace()->fe(), pts ) );
    }

    /**
     * precompute the basis function associated with the test and
     * trial space at a set of points on the face of
     * the reference element
     */
    template<typename Pts>
    void precomputeBasisAtPoints( uint16_type __f, permutation_type const& __p, Pts const& pts )
    {
        _M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( testSpace()->fe(), pts ) );
        //LIFE_ASSERT( _M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );

        _M_trial_pc_face[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( trialSpace()->fe(), pts ) );
        //LIFE_ASSERT( _M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
    }



    /**
     * Diagonalize representation(matrix) associated to the \p
     * BilinearForm at selected dofs \p dofs by putting 0 on
     * the extra diagonal terms and 1 on the diagonal.
     *
     * If \p ON_ELIMINATION_KEEP_DIAGONAL is set in \p on_context then
     * the diagonal value of the matrix is kept and the right habd
     * side \p rhs is modified accordingly.
     */
    void zeroRows( std::vector<int> __dofs,
                   std::vector<value_type> __values,
                   Vector<value_type>& rhs,
                   Life::Context const& on_context );

    /**
     * add value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void add( size_type i,  size_type j,  value_type const& v )
    {
        if ( _M_do_threshold )
            {
                if ( doThreshold( v ) )
                    _M_matrix.add( i, j, v );
            }
        else
            _M_matrix.add( i, j, v );

    }



    /**
     * set value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void set( size_type i,  size_type j,  value_type const& v ) { _M_matrix.set( i, j, v ); }

    void addToNOz( size_type i, size_type n ) { M_n_oz[i] += n; }
    void addToNNz( size_type i, size_type n ) { M_n_nz[i] += n; }
    size_type nOz( size_type i ) const { return M_n_oz[i]; }
    size_type nNz( size_type i ) const { return M_n_nz[i]; }


    graph_ptrtype computeGraph( size_type hints, mpl::bool_<true> );
    graph_ptrtype computeGraph( size_type hints, mpl::bool_<false> );
    void mergeGraph( int row, int col, graph_ptrtype g );

    //@}

protected:

    BilinearForm();

private:


    template <class ExprT>
    void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<false> );

    template <class ExprT>
    void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<true> );
private:

    space_1_ptrtype _M_X1;
    space_2_ptrtype _M_X2;

    matrix_type& _M_matrix;

    bool _M_do_build;

    list_block_type _M_lb;

    test_precompute_ptrtype _M_test_pc;
    std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > _M_test_pc_face;
    trial_precompute_ptrtype _M_trial_pc;
    std::map<uint16_type, std::map<permutation_type, trial_precompute_ptrtype> > _M_trial_pc_face;

    bool _M_do_threshold;
    value_type _M_threshold;

    std::vector<size_type> M_n_nz;
    std::vector<size_type> M_n_oz;

    graph_ptrtype M_graph;
};
template<typename FE1,  typename FE2, typename ElemContType>
BilinearForm<FE1, FE2, ElemContType>::BilinearForm( space_1_ptrtype const& Xh,
                                                    space_2_ptrtype const& Yh,
                                                    matrix_type& __M,
                                                    bool build,
                                                    bool do_threshold,
                                                    value_type threshold,
                                                    size_type graph_hints
                                                    )
    :
    _M_X1( Xh ),
    _M_X2( Yh ),
    _M_matrix( __M ),
    _M_do_build( build ),
    _M_lb(),
    _M_do_threshold( do_threshold ),
    _M_threshold( threshold ),
    M_graph( new graph_type( Xh->nLocalDof(),
                             Xh->nDofStart(), Xh->nDofStart()+Xh->nLocalDof(),
                             Yh->nDofStart(), Yh->nDofStart()+Yh->nLocalDof() ) )
{
    Debug( 5050 ) << "begin constructor with default listblock\n";

    _M_lb.push_back( Block ( 0, 0, 0, 0 ) );

    if ( build )
        {
            const size_type n1_dof_on_proc = _M_X1->nLocalDof();
            M_n_nz.resize (n1_dof_on_proc);
            M_n_oz.resize (n1_dof_on_proc);

            boost::timer t;
            Debug( 5050 ) << "compute graph\n";
            graph_ptrtype graph = computeGraph( graph_hints, mpl::bool_<(FE1::nSpaces == 1)>() );
            Debug( 5050 ) << "computed graph in " << t.elapsed() << "s\n"; t.restart();
            //_M_X1.setState( space_1_type::SOLUTION );
            // _M_X2->setState( space_2_type::TRIAL );

            Debug( 5050 ) << "init matrix with graph\n";
            //_M_matrix.resize( _M_X2->nDof(), _M_X1->nDof(), false );
            _M_matrix.init( _M_X1->nDof(), _M_X2->nDof(), _M_X1->nLocalDof(), _M_X2->nLocalDof(),
                            graph );
            Debug( 5050 ) << "initialized matrix in " << t.elapsed() << "s\n"; t.restart();
            //_M_X2->dof()->nLocalDof(), _M_X1->dof()->nLocalDof() );


        }
    Debug( 5050 ) << "begin constructor with default listblock done\n";
}

template<typename FE1,  typename FE2, typename ElemContType>
BilinearForm<FE1, FE2, ElemContType>::BilinearForm( space_1_ptrtype const& Xh,
                                                    space_2_ptrtype const& Yh,
                                                    matrix_type& __M,
                                                    list_block_type const& __lb,
                                                    bool do_threshold ,
                                                    value_type threshold,
                                                    size_type graph_hints )
    :
    _M_X1( Xh ),
    _M_X2( Yh ),
    _M_matrix( __M ),
    _M_do_build( false ),
    _M_lb( __lb ),
    _M_do_threshold( do_threshold ),
    _M_threshold( threshold ),
    M_graph( new graph_type( Xh->nLocalDof(),
                             Xh->nDofStart(), Xh->nDofStart()+Xh->nLocalDof(),
                             Yh->nDofStart(), Yh->nDofStart()+Yh->nLocalDof() ) )
{

}

template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
                                              bool init,
                                              mpl::bool_<false> )
{
    if ( init )
        {
            Debug( 5055 ) << "[BilinearForm::assign<false>] start\n";
            typedef ublas::matrix_range<matrix_type> matrix_range_type;
            typename list_block_type::const_iterator __bit = _M_lb.begin();
            typename list_block_type::const_iterator __ben = _M_lb.end();
            for ( ; __bit != __ben; ++__bit )
                {
                    size_type g_ic_start = __bit->globalRowStart();
                    size_type g_jc_start = __bit->globalColumnStart();
                    _M_matrix.zero( g_ic_start, g_ic_start + _M_X1->nDof(),
                                    g_jc_start, g_jc_start + _M_X2->nDof() );
                }
        }
    __expr.assemble( _M_X1, _M_X2, *this );

    Debug( 5055 ) << "[BilinearForm::assign<false>] stop\n";

}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
                                              bool /*init*/,
                                              mpl::bool_<true> )
{
    Debug( 5050 ) << "BilinearForm::assign() start loop on test spaces\n";
    fusion::for_each( _M_X1->functionSpaces(), make_bfassign2( *this, __expr ) );
    Debug( 5050 ) << "BilinearForm::assign() stop loop on test spaces\n";
}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator=( Expr<ExprT> const& __expr )
{
    // loop(fusion::for_each) over sub-functionspaces in SpaceType
    // pass expression and initialize
    this->assign( __expr, true, mpl::bool_<(FE1::nSpaces > 1)>() );
    return *this;
}

template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator+=( Expr<ExprT> const& __expr )
{
    Debug( 5055 ) << "[BilinearForm::operator+=] start\n";
    this->assign( __expr, false, mpl::bool_<(FE1::nSpaces > 1)>() );
    Debug( 5055 ) << "[BilinearForm::operator+=] stop\n";
    return *this;
}


template<typename FE1,  typename FE2, typename ElemContType>
void
BilinearForm<FE1,FE2,ElemContType>::zeroRows( std::vector<int> __dofs,
                                              std::vector<value_type> __values,
                                              Vector<value_type>& rhs,
                                              Life::Context const& on_context )
{
    _M_matrix.zeroRows( __dofs, __values, rhs, on_context );
}

template<typename BidirectionalIterator>
inline
void
sortSparsityRow (const BidirectionalIterator begin,
                 BidirectionalIterator       middle,
                 const BidirectionalIterator end)
{
    if ((begin == middle) || (middle == end)) return;

    assert (std::distance (begin,  middle) > 0);
    assert (std::distance (middle, end)    > 0);
    LIFE_ASSERT( std::unique (begin,  middle) == middle )
        ( *begin )( *middle ).error( "duplicate dof(begin,middle)" );
    LIFE_ASSERT (std::unique (middle, end)    == end)
        (*begin)( *middle ).error( "duplicate dof (middle,end)" );

    while (middle != end)
        {
            BidirectionalIterator
                b = middle,
                a = b-1;

            // Bubble-sort the middle value downward
            while (!(*a < *b)) // *a & *b are less-than comparable, so use <
                {
                    std::swap (*a, *b);

                    if (a == begin) break;

                    b=a;
                    --a;
                }

            ++middle;
        }

    // Assure the algorithm worked if we are in DEBUG mode
#ifdef DEBUG
    {
        // SGI STL extension!
        // assert (std::is_sorted(begin,end));

        BidirectionalIterator
            prev  = begin,
            first = begin;

        for (++first; first != end; prev=first, ++first)
            if (*first < *prev)
                assert(false);
    }
#endif

    // Make sure the two ranges did not contain any common elements
    assert (std::unique (begin, end) == end);
} //

template<typename FE1,  typename FE2, typename ElemContType>
void
BilinearForm<FE1,FE2,ElemContType>::mergeGraph( int row, int col, graph_ptrtype g )
{
    Debug( 5050 ) << "[merge graph] for composite bilinear form\n";
    Debug( 5050 ) << "[mergeGraph] row = " << row << "\n";
    Debug( 5050 ) << "[mergeGraph] col = " << col << "\n";

    // nothing yet in store
    if ( !M_graph || M_graph->empty() )
        {
            Debug( 5050 ) << "[merge graph] nothing yet in store, copy graph\n";
            M_graph = g;

#if 0
            typename graph_type::const_iterator it = g->begin();
            typename graph_type::const_iterator en = g->end();
            for( ; it != en; ++it )
                {
                    std::vector<size_type> const& ivec = boost::get<2>( it->second );
                    for( int i = 0;i < ivec.size(); ++i )
                        {
                            Debug( 5050 ) << "[mergeGraph] ivec[" << i << "] = " << ivec[i] << "\n";
                        }
                }
#endif // 9
        }
    else
        {
            M_graph->setLastRowEntryOnProc( row + g->lastRowEntryOnProc() );
            M_graph->setLastColEntryOnProc( col + g->lastColEntryOnProc() );

            Debug( 5050 ) << "[merge graph] already something in store\n";
            typename graph_type::const_iterator it = g->begin();
            typename graph_type::const_iterator en = g->end();
            for( ; it != en; ++it )
                {
                    int theglobalrow = row+it->first;
                    int thelocalrow = row + boost::get<1>( it->second );



                    Debug( 5050 ) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
                    M_graph->row(theglobalrow).template get<1>() = thelocalrow;

                    // Get the row of the sparsity pattern
                    std::vector<size_type> ivec( boost::get<2>( it->second ).begin(), boost::get<2>( it->second ).end() );
                    std::for_each( ivec.begin(), ivec.end(), boost::phoenix::arg_names::arg1 += col );
                    //std::set<size_type> iout( ivec.size()+ M_graph->row(theglobalrow).template get<2>().size() );
                    std::set<size_type> iout( ivec.begin(), ivec.end() );

                    iout.insert( M_graph->row(theglobalrow).template get<2>().begin(),
                                 M_graph->row(theglobalrow).template get<2>().end() );
#if 0
                    size_type start1 = M_graph->row(theglobalrow).template get<2>().size();
                    std::copy( M_graph->row(theglobalrow).template get<2>().begin(),
                               M_graph->row(theglobalrow).template get<2>().end(),
                               iout.begin() );

                    for( size_type i = 0;i < ivec.size(); ++i )
                        {
                            //Debug( 5050 ) << "[mergeGraph] ivec[" << i << "] = " << ivec[i] << "\n";
                            iout[start1+i]=col+ivec[i];
                        }
#endif
                    M_graph->row(theglobalrow).template get<2>() = iout;
                }


        }



    Debug( 5050 ) << "merge graph for composite bilinear form done\n";
}
template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraph( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form\n";
    fusion::for_each( _M_X1->functionSpaces(), compute_graph1<self_type>( *this, hints ) );
    Debug( 5050 ) << "closing graph for composite bilinear form done in " << t.elapsed() << "s\n"; t.restart();

#if 0
    typename graph_type::const_iterator it = M_graph->begin();
    typename graph_type::const_iterator en = M_graph->end();
    for( ; it != en; ++it )
        {
            std::cout << "row " << it->first << ", " <<  boost::get<1>( it->second ) << ": ";
            std::copy( boost::get<2>( it->second ).begin(), boost::get<2>( it->second ).end(), std::ostream_iterator<int>( std::cout, " " ) );
            std::cout << "\n";
        }
#endif
    M_graph->close();
    Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";


    return M_graph;
}
#if 0
template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraph( size_type hints, mpl::bool_<true> )
{
    boost::timer t;
    // Compute the sparsity structure of the global matrix.  This can be
    // fed into a PetscMatrix to allocate exacly the number of nonzeros
    // necessary to store the matrix.  This algorithm should be linear
    // in the (# of elements)*(# nodes per element)
    const size_type nprocs           = _M_X1->mesh()->comm().size();
    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                                    first1_dof_on_proc, last1_dof_on_proc,
                                                    first2_dof_on_proc, last2_dof_on_proc ) );

    typedef typename mesh_type::element_const_iterator mesh_element_const_iterator;

    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    Life::Context graph( hints );
    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    Debug( 5050 ) << "[computeGraph] test : " << ( graph.test ( DOF_PATTERN_COUPLED ) || graph.test ( DOF_PATTERN_NEIGHBOR ) ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_COUPLED )=" <<  graph.test ( DOF_PATTERN_COUPLED ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_NEIGHBOR)=" <<  graph.test ( DOF_PATTERN_NEIGHBOR ) << "\n";
#if 0
    if ( graph.test ( DOF_PATTERN_COUPLED ) ||
         graph.test ( DOF_PATTERN_NEIGHBOR ) )
#else
        if (1)
#endif
        {
            Debug( 5050 ) << "[computeGraph] test (DOF_PATTERN_COUPLED || DOF_PATTERN_NEIGHBOR) ok\n";
            std::vector<size_type>
                element_dof1,
                element_dof2,
                neighbor_dof,
                dof_to_add;

            for ( ; elem_it != elem_en; ++elem_it)
                {
#if !defined(NDEBUG)
                    Debug( 5050 ) << "[BilinearForm::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
                    mesh_element_type const& elem = *elem_it;

                    // Get the global indices of the DOFs with support on this element
                    element_dof1 = _M_X1->dof()->getIndices( elem.id() );
                    element_dof2 = _M_X2->dof()->getIndices( elem.id() );

                    // We can be more efficient if we sort the element DOFs
                    // into increasing order
                    std::sort(element_dof1.begin(), element_dof1.end());
                    std::sort(element_dof2.begin(), element_dof2.end());

                    const uint16_type  n1_dof_on_element = element_dof1.size();
                    const uint16_type  n2_dof_on_element = element_dof2.size();

                    for (size_type i=0; i<n1_dof_on_element; i++)
                        {
                            const size_type ig1 = element_dof1[i];

                            // Only bother if this matrix row will be stored
                            // on this processor.
                            //if ((ig1 >= first1_dof_on_proc) &&
                            //(ig1 <= last1_dof_on_proc))
                                {
                                    // This is what I mean
                                    // assert ((ig - first_dof_on_proc) >= 0);
                                    // but do the test like this because ig and
                                    // first_dof_on_proc are size_types
#if 0
                                    LIFE_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                                    LIFE_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
                                        ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                                    graph_type::row_type& row = sparsity_graph->row(ig1);
                                    bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                                    row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                    row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                                    Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                                    // If the row is empty we will add *all* the element DOFs,
                                    // so just do that.
                                    if (row.get<2>().empty())
                                        {
                                            row.get<2>() = element_dof2;
                                        }
                                    else
                                        {
                                            // Build a list of the DOF indices not found in the
                                            // sparsity graph
                                            dof_to_add.clear();

                                            // Cache iterators.  Low will move forward, subsequent
                                            // searches will be on smaller ranges
                                            std::vector<size_type>::iterator
                                                low  = std::lower_bound (row.get<2>().begin(), row.get<2>().end(), element_dof2.front()),
                                                high = std::upper_bound (low,         row.get<2>().end(), element_dof2.back());

                                            for (size_type j=0; j<n2_dof_on_element; j++)
                                                {
                                                    const size_type jg = element_dof2[j];
                                                    //Debug() << "computeGraph : ig:" << ig1 << ", lig: " << ig1-first1_dof_on_proc  << ", jg = " << jg << "\n";
#if 0
                                                    // See if jg is in the sorted range
                                                    std::pair<std::vector<size_type>::iterator,
                                                        std::vector<size_type>::iterator>
                                                        pos = std::equal_range (low, high, jg);

                                                    // Must add jg if it wasn't found
                                                    if (pos.first == pos.second)
                                                        dof_to_add.push_back(jg);

                                                    // pos.first is now a valid lower bound for any
                                                    // remaining element Dof. (That's why we sorted them.)
                                                    // Use it for the next search
                                                    low = pos.first;
#else
                                                    // See if jg is in the sorted range
                                                    std::pair<std::vector<size_type>::iterator,
                                                        std::vector<size_type>::iterator>
                                                        pos = std::equal_range (row.get<2>().begin(), row.get<2>().end(), jg);

                                                    // Insert jg if it wasn't found
                                                    if (pos.first == pos.second)
                                                        dof_to_add.push_back(jg);
#endif
                                                }

                                            // Add to the sparsity graph
                                            if (!dof_to_add.empty())
                                                {
                                                    const size_type old_size = row.get<2>().size();

                                                    row.get<2>().insert (row.get<2>().end(),
                                                                dof_to_add.begin(),
                                                                dof_to_add.end());

                                                    //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                    sortSparsityRow (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                }

                                        }

                                    // Now (possibly) add dof from neighboring elements
                                    //if ( graph.test( DOF_PATTERN_NEIGHBOR ) )
                                        for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                                            {
                                                mesh_element_type const* neighbor = NULL;
                                                size_type neighbor_id = elem.neighbor(ms).first;
                                                size_type neighbor_process_id = elem.neighbor(ms).second;
                                                if ( neighbor_id != invalid_size_type_value )
                                                    //&& neighbor_process_id != proc_id )
                                                    {

#if 0
                                                        Debug() << "element id " << elem.id()
                                                                << ", element neighbor id " << neighbor_id
                                                                << " in proc " << neighbor_process_id << "\n";
#endif
                                                        neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                                                                             neighbor_process_id ) );

#if 0
                                                        Debug() << "found neighbor of element id " << elem.id()
                                                                << ", element neighbor id " << neighbor->id()
                                                                << " in proc " << neighbor->processId() << "\n";
#endif

                                                        if ( neighbor_id == neighbor->id()  )
                                                            {
                                                                neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );

                                                                const size_type n_dof_on_neighbor = neighbor_dof.size();
#if 0
                                                                for (size_type j=0; j<n_dof_on_neighbor; j++)
                                                                    {
                                                                        Debug( 5051 ) << "neighbor elem id: " << neighbor->id() << " dof " << neighbor_dof[j] << "\n";
                                                                    }
                                                                Debug( 5051 ) << "looking for dof " << ig1  << "\n";
#endif
#if 0
                                                                std::pair<std::vector<size_type>::iterator,
                                                                    std::vector<size_type>::iterator>
                                                                    posig = std::equal_range (neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#else
                                                                std::vector<size_type>::iterator posig = std::find( neighbor_dof.begin(), neighbor_dof.end(), ig1 );
#endif
                                                                // Insert jg if it wasn't found
                                                                //if (posig.first != posig.second)
                                                                if ( posig != neighbor_dof.end() ||
                                                                     graph.test ( DOF_PATTERN_NEIGHBOR ) )

                                                                    {
                                                                        //Debug() << "found element in proc " << neighbor_process_id << " that shares dof\n";
                                                                        for (size_type j=0; j<n_dof_on_neighbor; j++)
                                                                            {
                                                                                const size_type jg = neighbor_dof[j];

#if 0
                                                                                // See if jg is in the sorted range
                                                                                std::pair<std::vector<size_type>::iterator,
                                                                                    std::vector<size_type>::iterator>
                                                                                    pos = std::equal_range (row.get<2>().begin(), row.get<2>().end(), jg);
#else
                                                                                std::vector<size_type>::iterator pos = std::find( row.get<2>().begin(), row.get<2>().end(), jg );

#endif
                                                                                // Insert jg if it wasn't found
                                                                                if (pos == row.get<2>().end() )
                                                                                    {
                                                                                        const size_type old_size = row.get<2>().size();
                                                                                        row.get<2>().push_back (jg);
                                                                                        //std::inplace_merge (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                                                        sortSparsityRow (row.get<2>().begin(), row.get<2>().begin()+old_size, row.get<2>().end());
                                                                                    }
                                                                            }
                                                                    }
                                                            }
                                                    }

                                            } // neighbor graph
#if 0
                                        for( int k = 0; k < row.get<2>().size(); ++k )
                                            Debug() << "row[ " << ig1 - first1_dof_on_proc << ","<< k << " ]=" << row.get<2>()[k] << "\n";
#endif
                                } // only dof on proc

                        }// dof loop
                } // element iterator loop
        }
    else
        {}


    sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#else
template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraph( size_type hints, mpl::bool_<true> )
{
    boost::timer t;
    // Compute the sparsity structure of the global matrix.  This can be
    // fed into a PetscMatrix to allocate exacly the number of nonzeros
    // necessary to store the matrix.  This algorithm should be linear
    // in the (# of elements)*(# nodes per element)
    const size_type nprocs           = _M_X1->mesh()->comm().size();
    const size_type proc_id           = _M_X1->mesh()->comm().rank();
    const size_type n1_dof_on_proc    = _M_X1->nLocalDof();
    //const size_type n2_dof_on_proc    = _M_X2->nLocalDof();
    const size_type first1_dof_on_proc = _M_X1->dof()->firstDof( proc_id );
    const size_type last1_dof_on_proc = _M_X1->dof()->lastDof( proc_id );
    const size_type first2_dof_on_proc = _M_X2->dof()->firstDof( proc_id );
    const size_type last2_dof_on_proc = _M_X2->dof()->lastDof( proc_id );

    graph_ptrtype sparsity_graph( new graph_type( n1_dof_on_proc,
                                                    first1_dof_on_proc, last1_dof_on_proc,
                                                    first2_dof_on_proc, last2_dof_on_proc ) );

    typedef typename mesh_type::element_const_iterator mesh_element_const_iterator;

    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    Life::Context graph( hints );
    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    Debug( 5050 ) << "[computeGraph] test : " << ( graph.test ( DOF_PATTERN_COUPLED ) || graph.test ( DOF_PATTERN_NEIGHBOR ) ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_COUPLED )=" <<  graph.test ( DOF_PATTERN_COUPLED ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_NEIGHBOR)=" <<  graph.test ( DOF_PATTERN_NEIGHBOR ) << "\n";
#if 0
    if ( graph.test ( DOF_PATTERN_COUPLED ) ||
         graph.test ( DOF_PATTERN_NEIGHBOR ) )
#else
        if (1)
#endif
        {
            Debug( 5050 ) << "[computeGraph] test (DOF_PATTERN_COUPLED || DOF_PATTERN_NEIGHBOR) ok\n";
            std::vector<size_type>
                element_dof1,
                element_dof2,
                neighbor_dof;
            std::set<size_type> dof_to_add;

            for ( ; elem_it != elem_en; ++elem_it)
                {
#if !defined(NDEBUG)
                    Debug( 5050 ) << "[BilinearForm::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
                    mesh_element_type const& elem = *elem_it;

                    // Get the global indices of the DOFs with support on this element
                    element_dof1 = _M_X1->dof()->getIndices( elem.id() );
                    element_dof2 = _M_X2->dof()->getIndices( elem.id() );

                    // We can be more efficient if we sort the element DOFs
                    // into increasing order
                    std::sort(element_dof1.begin(), element_dof1.end());
                    std::sort(element_dof2.begin(), element_dof2.end());

                    const uint16_type  n1_dof_on_element = element_dof1.size();
                    const uint16_type  n2_dof_on_element = element_dof2.size();

                    for (size_type i=0; i<n1_dof_on_element; i++)
                        {
                            const size_type ig1 = element_dof1[i];

                            // Only bother if this matrix row will be stored
                            // on this processor.
                            //if ((ig1 >= first1_dof_on_proc) &&
                            //(ig1 <= last1_dof_on_proc))
                                {
                                    // This is what I mean
                                    // assert ((ig - first_dof_on_proc) >= 0);
                                    // but do the test like this because ig and
                                    // first_dof_on_proc are size_types
#if 0
                                    LIFE_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                                    LIFE_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
                                        ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                                    graph_type::row_type& row = sparsity_graph->row(ig1);
                                    bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                                    row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                    row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                                    Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                                    // Build a list of the DOF indices not found in the
                                    // sparsity graph
                                    //dof_to_add.clear();
                                    //dof_to_add.insert( row.get<2>().begin(), row.get<2>().end() );

                                    row.get<2>().insert( element_dof2.begin(), element_dof2.end() );
                                    // Now (possibly) add dof from neighboring elements
                                    if ( graph.test( DOF_PATTERN_NEIGHBOR ) )
                                    {
                                        for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                                        {
                                            mesh_element_type const* neighbor = NULL;
                                            size_type neighbor_id = elem.neighbor(ms).first;
                                            size_type neighbor_process_id = elem.neighbor(ms).second;
                                            if ( neighbor_id != invalid_size_type_value )
                                                //&& neighbor_process_id != proc_id )
                                            {

                                                neighbor = boost::addressof( _M_X1->mesh()->element( neighbor_id,
                                                                                                     neighbor_process_id ) );

                                                if ( neighbor_id == neighbor->id()  )
                                                {
                                                    neighbor_dof = _M_X2->dof()->getIndices( neighbor->id() );
                                                    row.get<2>().insert( neighbor_dof.begin(), neighbor_dof.end() );
                                                } // neighbor_id
                                            }

                                        } // neighbor graph
                                    }
                                } // only dof on proc

                        }// dof loop
                } // element iterator loop
        }
    else
        {}


    sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#endif

//
// Context
//
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                              map_geometric_mapping_context_type const& _gmc,
                                                                              ExprT const& expr,
                                                                              IM const& im )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc,
                                    detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(),
                                                                       _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmc,
                                                          detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(),
                                                                                              _M_form ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1] ),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),
    M_integrator( im )
{
    _M_eval_expr00->init( im );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                              map_geometric_mapping_context_type const& _gmc,
                                                                              ExprT const& expr,
                                                                              IM const& im,
                                                                              IM2 const& im2 )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_trial_fec( getMap( _M_test_fec, fusion::transform( _gmc, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), _M_form ) ) ) ),
    _M_trial_fec0( getMapL( _M_test_fec0, fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1] ),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01(),
    _M_eval_expr10(),
    _M_eval_expr11(),
    M_integrator( im )
{
    // faces
    _M_eval_expr00->init( im2 );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
template<typename IM2>
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::Context( form_type& __form,
                                                                                         map_geometric_mapping_context_type const& _gmc,
                                                                                         ExprT const& expr,
                                                                                         IM const& im,
                                                                                         IM2 const& im2,
                                                                                         mpl::int_<2> )
    :
    _M_form( __form ),
    _M_lb( __form.blockList() ),
    _M_test_dof( __form.testSpace()->dof().get() ),
    _M_trial_dof( __form.trialSpace()->dof().get() ),
    _M_gmc( _gmc ),
    _M_test_fec( fusion::transform( _gmc, detail::FEContextInit<0,form_context_type>(__form.testSpace()->fe(), _M_form ) ) ),
    _M_test_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) ) ),
    _M_test_fec1( fusion::make_map<gmc<1> >( fusion::at_key<gmc<1> >( _M_test_fec ) ) ),
    _M_trial_fec( fusion::transform( _gmc, detail::FEContextInit<1,form_context_type>( __form.trialSpace()->fe(), _M_form ) ) ),
    _M_trial_fec0( fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) ) ),
    _M_trial_fec1( fusion::make_map<gmc<1> >( fusion::at_key<gmc<1> >( _M_trial_fec ) ) ),
    _M_rep( boost::extents[fusion::size(_gmc)*test_fecontext_type::nDof][test_fecontext_type::nComponents1]
            [fusion::size(_gmc)*trial_fecontext_type::nDof][trial_fecontext_type::nComponents1]),
    _M_eval_expr00( new eval00_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec0 ) ),
    _M_eval_expr01( new eval01_expr_type( expr, _gmc, _M_test_fec0, _M_trial_fec1 ) ),
    _M_eval_expr10( new eval10_expr_type( expr, _gmc, _M_test_fec1, _M_trial_fec0 ) ),
    _M_eval_expr11( new eval11_expr_type( expr, _gmc, _M_test_fec1, _M_trial_fec1 ) ),
    M_integrator( im )
{
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 ).error( "invalid test_fec");
    LIFE_ASSERT( fusion::at_key<gmc<1> >( _M_test_fec1 ).get() != 0 ).error( "invalid test_fec");
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 ).error( "invalid trial_fec");
    LIFE_ASSERT( fusion::at_key<gmc<1> >( _M_trial_fec1 ).get() != 0 ).error( "invalid trial_fec");

    _M_eval_expr00->init( im2 );
    _M_eval_expr01->init( im2 );
    _M_eval_expr10->init( im2 );
    _M_eval_expr11->init( im2 );
}

template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc )
{
    update( _gmc,  boost::is_same<map_test_fecontext_type, map_trial_fecontext_type>() );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::bool_<false> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmc, _M_form ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_trial_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::bool_<true> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_test_fec0 );
    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
}
template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::update( map_geometric_mapping_context_type const& _gmc, mpl::int_<2> )
{
    fusion::for_each( _M_test_fec, detail::FEContextUpdate<0,form_context_type>( _gmc, _M_form ) );
    _M_test_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_test_fec ) );
    _M_test_fec1 = fusion::make_map<gmc<1> >( fusion::at_key<gmc<1> >( _M_test_fec ) );
    fusion::for_each( _M_trial_fec, detail::FEContextUpdate<1,form_context_type>( _gmc, _M_form ) );
    _M_trial_fec0 = fusion::make_map<gmc<0> >( fusion::at_key<gmc<0> >( _M_trial_fec ) );
    _M_trial_fec1 = fusion::make_map<gmc<1> >( fusion::at_key<gmc<1> >( _M_trial_fec ) );

    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_test_fec0 ).get() != 0 )
        ( 0 ).error( "invalid test_fec0" );
    LIFE_ASSERT( fusion::at_key<gmc<1> >( _M_test_fec1 ).get() != 0 )
        ( 1 ).error( "invalid test_fec1" );
    LIFE_ASSERT( fusion::at_key<gmc<0> >( _M_trial_fec0 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec0" );
    LIFE_ASSERT( fusion::at_key<gmc<1> >( _M_trial_fec1 ).get() != 0 )
        ( 0 ).error( "invalid trial_fec1" );

    _M_eval_expr00->update( _gmc, _M_test_fec0, _M_trial_fec0 );
    _M_eval_expr01->update( _gmc, _M_test_fec0, _M_trial_fec1 );
    _M_eval_expr10->update( _gmc, _M_test_fec1, _M_trial_fec0 );
    _M_eval_expr11->update( _gmc, _M_test_fec1, _M_trial_fec1 );

    std::fill( _M_rep.data(), _M_rep.data()+_M_rep.num_elements(), value_type(0) );
    M_integrator.update( *fusion::at_key<gmc<0> >( _gmc ) );
}



template<typename FE1,  typename FE2, typename ElemContType>
template<typename GeomapContext,typename ExprT,typename IM>
void
BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<1> )
{

    typedef geometric_mapping_context_type gmc_type;
    typedef typename eval00_expr_type::shape shape;
    static const bool cond = (shape::M == 1 && shape::N == 1);
    BOOST_MPL_ASSERT_MSG( cond,
                          INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                          (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

#if !defined(NDEBUG)
    geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
    Debug( 5050 ) << "[BilinearForm::integrate] local assembly in element " << _gmc.id() << "\n";
#endif /* NDEBUG */

    for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
        for ( uint16_type c1 = 0; c1 < test_fecontext_type::nComponents1; ++c1 )
            {
                indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                LIFE_ASSERT( size_type(indi.index()) ==
                             size_type(c1*test_fecontext_type::nDof+i) )
                    ( i )( c1 )
                    ( indi.index() )
                    ( c1*test_fecontext_type::nDof+i ).error( "invalid test index" );

                for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                    for ( uint16_type c2 = 0; c2 < trial_fecontext_type::nComponents1; ++c2 )
                        {
                            indj.setIndex( boost::make_tuple( j, c2, 0 ) );
                            LIFE_ASSERT( size_type(indj.index()) ==
                                         size_type(c2*trial_fecontext_type::nDof+j) )
                                ( j )( c2 )
                                ( indj.index() )
                                ( c2*trial_fecontext_type::nDof+j ).error( "invalid trial index" );
#if 0
                            if ( vm::has_symm<expression_type::context>::value )
                                {
                                    if ( indj.index() <= indi.index() )
                                        {
                                            _M_rep[i][c1][j][c2] = M_integrator( *_M_eval_expr00,
                                                                                 indi,
                                                                                 indj,
                                                                                 0, 0,
                                                                                 mpl::int_<vm::SYMM>() );
                                            _M_rep[j][c2][i][c1] = _M_rep[i][c1][j][c2];
                                        }
                                }
                            else
#endif
                                {
#if 0
                                    if ( (indj.index() <= indi.index()) &&
                                         (test_fecontext_type::nDof*test_fecontext_type::nComponents1 ==
                                          trial_fecontext_type::nDof*trial_fecontext_type::nComponents1) )
                                        {
                                            _M_rep[i][c1][j][c2] = M_integrator( *_M_eval_expr00,
                                                                                 indi,
                                                                                 indj,
                                                                                 0, 0 );
                                            _M_rep[j][c2][i][c1] = _M_rep[i][c1][j][c2];
                                        }
                                    else
#endif
                                        {
                                            _M_rep[i][c1][j][c2] = M_integrator( *_M_eval_expr00,
                                                                                 indi,
                                                                                 indj,
                                                                                 0, 0 );
                                        }
                                }
                        } // j, c2
            }
}

    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
    BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::integrate( mpl::int_<2> )
        {
            //geometric_mapping_context_type const& _gmc = *fusion::at_key<gmc<0> >( _M_gmc );
            typedef geometric_mapping_context_type gmc_type;
            typedef typename eval00_expr_type::shape shape;
            BOOST_MPL_ASSERT_MSG( (shape::M == 1 && shape::N == 1),
                                  INVALID_TENSOR_SHAPE_SHOULD_BE_RANK_0,
                                  (mpl::int_<shape::M>, mpl::int_<shape::N> ) );

            //BOOST_MPL_ASSERT_MSG( eval_expr_type::shape::M == 1, INVALID_TENSOR_EXTENT_M, mpl::int_<eval_expr_type::shape::M> );
            //BOOST_MPL_ASSERT_MSG( eval_expr_type::shape::N == 1, INVALID_TENSOR_EXTENT_N, mpl::int_<eval_expr_type::shape::N> );

            for ( uint16_type i = 0; i < test_fecontext_type::nDof; ++i )
                for ( uint16_type c1 = 0; c1 < test_fecontext_type::nComponents1; ++c1 )
                    {
                        indi.setIndex( boost::make_tuple( i, c1, 0 ) );
                        for ( uint16_type j = 0; j < trial_fecontext_type::nDof; ++j )
                            for ( uint16_type c2 = 0; c2 < trial_fecontext_type::nComponents1; ++c2 )
                                {
                                    indj.setIndex( boost::make_tuple( j, c2, 0 ) );
                                    _M_rep[i][c1][j][c2] += M_integrator( *_M_eval_expr00, indi, indj, 0, 0 );

                                    uint16_type ii = i;
                                    uint16_type jj = j + trial_fecontext_type::nDof;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr01, indi, indj, 0, 0 );

                                    ii = i+ test_fecontext_type::nDof;
                                    jj = j ;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr10, indi, indj, 0, 0 );

                                    ii = i+ test_fecontext_type::nDof;
                                    jj = j + trial_fecontext_type::nDof;

                                    _M_rep[ii][c1][jj][c2] += M_integrator( *_M_eval_expr11, indi, indj, 0, 0 );
                                }

                    }
        }
    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
        BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0 )
        {

            size_type ig,jg;
            size_type row_start = _M_lb.front().globalRowStart();
            size_type col_start = _M_lb.front().globalColumnStart();
            int isign, jsign;

#if !defined(NDEBUG)
            Debug( 5050 ) << "[BilinearForm::assemble] global assembly in element " << elt_0 << "\n";
            Debug( 5050 ) << "[BilinearForm::assemble] row start " << row_start << "\n";
            Debug( 5050 ) << "[BilinearForm::assemble] col start " << col_start << "\n";
#endif /* NDEBUG */

            for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                for ( uint16_type c1 = 0 ; c1 < test_fecontext_type::nComponents1; ++c1 )
                    {
                        boost::tie( ig, isign, boost::tuples::ignore) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        ig += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                            for ( uint16_type c2 = 0 ; c2 < trial_fecontext_type::nComponents1; ++c2 )
                                {
                                    boost::tie( jg, jsign, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c2 );
                                    jg += col_start;
                                    _M_form.add( ig, jg, value_type(isign*jsign)*_M_rep[i][c1][j][c2] );
                                }
                    }
        }

    template<typename FE1,  typename FE2, typename ElemContType>
        template<typename GeomapContext,typename ExprT,typename IM>
        void
        BilinearForm<FE1,FE2,ElemContType>::Context<GeomapContext,ExprT,IM>::assemble( size_type elt_0, size_type elt_1  )
        {
            test_index_type indi;
            trial_index_type indj;
            size_type ig0,ig1,jg0,jg1;
            size_type row_start = _M_lb.front().globalRowStart();
            size_type col_start = _M_lb.front().globalColumnStart();
            int isign0, isign1, jsign0, jsign1;
            const uint16_type test_ndof = test_fecontext_type::nDof;
            const uint16_type trial_ndof = trial_fecontext_type::nDof;
            //std::cout << "rep_" << elt_0 << "=" << _M_rep << "\n";
            for ( uint16_type i = 0 ; i < test_fecontext_type::nDof; ++i )
                for ( uint16_type c1 = 0 ; c1 < test_fecontext_type::nComponents1; ++c1 )
                    {
                        boost::tie( ig0, isign0, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_0, i, c1 );
                        boost::tie( ig1, isign1, boost::tuples::ignore ) = _M_test_dof->localToGlobal( elt_1, i, c1 );
                        ig0 += row_start;
                        ig1 += row_start;
                        for ( uint16_type j = 0 ; j < trial_fecontext_type::nDof; j++ )
                            for ( uint16_type c2 = 0 ; c2 < trial_fecontext_type::nComponents1; ++c2 )
                                {
                                    boost::tie( jg0, jsign0, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_0, j, c2 );
                                    boost::tie( jg1, jsign1, boost::tuples::ignore ) = _M_trial_dof->localToGlobal( elt_1, j, c2 );
                                    jg0 += col_start;
                                    jg1 += col_start;

                                    _M_form.add( ig0, jg0, value_type(isign0*jsign0)*_M_rep[i][c1][j][c2] );
                                    _M_form.add( ig0, jg1, value_type(isign0*jsign1)*_M_rep[i][c1][j+trial_ndof][c2] );
                                    _M_form.add( ig1, jg0, value_type(isign0*jsign0)*_M_rep[i+test_ndof][c1][j][c2] );
                                    _M_form.add( ig1, jg1, value_type(isign0*jsign1)*_M_rep[i+test_ndof][c1][j+trial_ndof][c2] );
                                }
                    }
        }
template<typename BFType, typename ExprType, typename TestSpaceType>
template<typename SpaceType>
void BFAssign1<BFType,ExprType,TestSpaceType>::operator()( boost::shared_ptr<SpaceType> const& trial ) const
{
    Debug(5050) << "[BFAssign1::operator()] expression has test functions index "
                << _M_test_index << " : "
                << ExprType::template HasTestFunction<typename TestSpaceType::reference_element_type>::result << " (0:no, 1:yes)\n";
    Debug(5050) << "[BFAssign1::operator()] expression has trial functions index "
                << _M_trial_index << " :"
                << ExprType::template HasTrialFunction<typename SpaceType::reference_element_type>::result << " (0:no, 1:yes)\n";

    if ( !ExprType::template HasTestFunction<typename TestSpaceType::reference_element_type>::result ||
         !ExprType::template HasTrialFunction<typename SpaceType::reference_element_type>::result )
        {
            ++_M_trial_index;
            return;
        }

    Debug( 5050 ) << "[BFAssign1::operator()] terms found with "
                  << "testindex: " << _M_test_index << " trialindex: " << _M_trial_index << "\n";
    typedef SpaceType trial_space_type;
    typedef TestSpaceType test_space_type;

    Life::vf::Block block ( 0, 0,
                            _M_bf.testSpace()->nDofStart( _M_test_index ),
                            _M_bf.trialSpace()->nDofStart( _M_trial_index ) );
    Debug( 5050 ) << "[BFAssign1::operator()] block: " << block << "\n";
    Life::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Life::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;

    bf_type bf( _M_test,trial, _M_bf.matrix(), list_block );

    bf += _M_expr;

    ++_M_trial_index;
}

template<typename BFType, typename TestSpaceType>
template<typename SpaceType>
void
compute_graph2<BFType,TestSpaceType>::operator()( boost::shared_ptr<SpaceType> const& trial ) const
{
    boost::timer t;
    Debug( 5050 ) << "[compute_graph2::operator()] begin ==================================================\n";
    Debug( 5050 ) << "[compute_graph2::operator()] testindex: " << M_test_index << "\n";
    Debug( 5050 ) << "[compute_graph2::operator()] trialindex: " << M_trial_index << "\n";
    typedef SpaceType trial_space_type;
    typedef TestSpaceType test_space_type;

    Life::vf::Block block ( 0, 0,
                            M_bf.testSpace()->nDofStart( M_test_index ),
                            M_bf.trialSpace()->nDofStart( M_trial_index ) );
    Debug( 5050 ) << "[compute_graph2::operator()] block: " << block << "\n";
    Life::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Life::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;


    bf_type bf( M_space1, trial, M_bf.matrix(), list_block );
    typename bf_type::graph_ptrtype graph = bf.computeGraph( M_hints, mpl::bool_<true>() );
    Debug( 5050 ) << "[compute_graph2::operator()] testindex: " << M_test_index << " row space index starts at " << M_bf.testSpace()->nDofStart( M_test_index ) << "\n";
    Debug( 5050 ) << "[compute_graph2::operator()] trialindex: " << M_trial_index << " col space index starts at " << M_bf.trialSpace()->nDofStart( M_trial_index ) << "\n";
    M_bf.mergeGraph( M_bf.testSpace()->nDofStart( M_test_index ), M_bf.trialSpace()->nDofStart( M_trial_index ) , graph );


#if 0
    const size_type n1_dof_on_proc    = M_space1->nLocalDof();
    size_type start = M_bf.testSpace()->nDofStart( M_test_index );
    for (size_type i=0; i< n1_dof_on_proc; i++)
        {
            M_bf.addToNOz( start+i, bf.nOz( i ) );
            M_bf.addToNNz( start+i, bf.nNz( i ) );
        }
#endif
    Debug( 5050 ) << "[compute_graph2] trial index " << M_trial_index << " done in " << t.elapsed() << "s\n";
    ++M_trial_index;

    Debug( 5050 ) << "[compute_graph2::operator()] end ==================================================\n";
}

} // detail
} // vf
/// \endcond
} // life

#endif /* __BilinearForm_H */

