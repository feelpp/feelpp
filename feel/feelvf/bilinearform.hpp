/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2005-01-18

   Copyright (C) 2005,2006 EPFL
   Copyright (C) 2006-2010 Université Joseph Fourier (Grenoble I)

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

#include <Eigen/Eigen>
#include <Eigen/StdVector>


#include <set>

#include <boost/parameter.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/spirit/home/phoenix.hpp>
#include <boost/spirit/home/phoenix/core/argument.hpp>
#include <feel/feelcore/context.hpp>
#include <feel/feelalg/matrixvalue.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/graphcsr.hpp>
#include <feel/feelvf/block.hpp>
#include <feel/feelvf/fec.hpp>
#include <feel/feelvf/formcontextbase.hpp>


namespace Feel
{
namespace parameter = boost::parameter;
namespace fusion = boost::fusion;
namespace vf
{
enum DofGraph
    {
        DOF_PATTERN_DEFAULT  = 1 << 0,
        DOF_PATTERN_COUPLED  = 1 << 1,
        DOF_PATTERN_NEIGHBOR = 1 << 2
    };

/// \cond detail
template<typename FE1,typename FE2,typename ElemContType> class BilinearForm;
namespace detail
{

template<typename BFType, typename Space1Type>
struct compute_graph3
{
    compute_graph3( BFType& bf, boost::shared_ptr<Space1Type> const& space1, size_type trial_index, size_type hints  )
        :
        M_bf( bf ),
        M_space1( space1 ),
        M_test_index( 0 ),
        M_trial_index( trial_index ),
        M_hints( hints )
    {}

    template <typename Space2>
    void operator()( boost::shared_ptr<Space2> const& space2 ) const;

    mutable BFType& M_bf;
    boost::shared_ptr<Space1Type> const& M_space1;
    mutable size_type M_test_index;
    size_type M_trial_index;
    size_type M_hints;
};

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


template<typename BFType, typename ExprType, typename TrialSpaceType>
struct BFAssign3
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign3( BFAssign3 const& lfa )
        :
        _M_bf( lfa._M_bf ),
        _M_trial( lfa._M_trial ),
        _M_expr( lfa._M_expr ),
        _M_test_index( lfa._M_test_index )
    {}
    BFAssign3( BFType& lf,
               ExprType const& expr,
               boost::shared_ptr<TrialSpaceType> const& Trialh,
               size_type trial_index )
        :
        _M_bf( lf ),
        _M_trial( Trialh ),
        _M_expr( expr ),
        _M_trial_index( trial_index ),
        _M_test_index( 0 )
    {}
    template<typename SpaceType>
    void operator()( boost::shared_ptr<SpaceType> const& test ) const;

private:
    BFType& _M_bf;
    boost::shared_ptr<TrialSpaceType> _M_trial;
    ExprType const& _M_expr;
    size_type _M_trial_index;
    mutable size_type _M_test_index;
};
template<typename BFType, typename ExprType, typename TrialSpaceType>
BFAssign3<BFType,ExprType,TrialSpaceType>
make_bfassign3( BFType& lf,
                ExprType const& expr,
                boost::shared_ptr<TrialSpaceType> const& trial_space,
                size_type trial_index )
{ return BFAssign3<BFType,ExprType,TrialSpaceType>( lf, expr, trial_space, trial_index ); }


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

    typedef typename space_1_type::mesh_type mesh_1_type;
    typedef typename mesh_1_type::element_type mesh_element_1_type;
    typedef typename mesh_element_1_type::permutation_type permutation_1_type;

    typedef typename space_2_type::mesh_type mesh_2_type;
    typedef typename mesh_2_type::element_type mesh_element_2_type;
    typedef typename mesh_element_2_type::permutation_type permutation_2_type;

    typedef typename space_1_type::fe_type fe_1_type;
    typedef typename space_2_type::fe_type fe_2_type;

    typedef typename space_1_type::gm_type gm_1_type;
    typedef typename space_1_type::gm_ptrtype gm_1_ptrtype;

    typedef typename space_1_type::gm1_type gm1_1_type;
    typedef typename space_1_type::gm1_ptrtype gm1_1_ptrtype;

    typedef typename space_2_type::gm_type gm_2_type;
    typedef typename space_2_type::gm_ptrtype gm_2_ptrtype;

    typedef typename space_2_type::gm1_type gm1_2_type;
    typedef typename space_2_type::gm1_ptrtype gm1_2_ptrtype;

    //typedef ublas::compressed_matrix<value_type, ublas::row_major> csr_matrix_type;
    typedef MatrixSparse<value_type> matrix_type;
    static const bool is_row_major = true;//matrix_type::is_row_major;

    typedef typename mpl::if_<mpl::equal_to<mpl::bool_<is_row_major>, mpl::bool_<true> >,
                              mpl::identity<ublas::row_major>,
                              mpl::identity<ublas::column_major> >::type::type layout_type;

    template<int _N = 0>
    struct test_precompute
    {
        //typedef typename space_1_type::basis_0_type::template precompute<_N>::type type;
        typedef typename space_1_type::basis_0_type::PreCompute type;
        typedef boost::shared_ptr<type> ptrtype;
    };
    template<int _N = 0>
    struct trial_precompute
    {
        //typedef typename space_2_type::basis_0_type::template precompute<_N>::type type;
        typedef typename space_2_type::basis_0_type::PreCompute type;
        typedef boost::shared_ptr<type> ptrtype;
    };

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
    template<typename GeomapTestContext,
             typename ExprT,
             typename IM,
             typename GeomapExprContext = GeomapTestContext,
             typename GeomapTrialContext = GeomapTestContext
             >
    class Context //: public FormContextBase<GeomapTestContext,IM,GeomapExprContext>
    {
        typedef FormContextBase<GeomapTestContext,IM,GeomapExprContext> super;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW


        typedef Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext> form_context_type;
        typedef BilinearForm<FE1, FE2, ElemContType> form_type;
        typedef typename FE1::dof_type dof_1_type;
        typedef typename FE2::dof_type dof_2_type;

        typedef typename form_type::value_type value_type;


        typedef typename super::map_geometric_mapping_context_type map_test_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_context_ptrtype test_geometric_mapping_context_ptrtype;
        typedef typename super::geometric_mapping_context_type test_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_type test_geometric_mapping_type;
        typedef typename super::map_size test_map_size;
        typedef typename super::gmc1 test_gmc1;
        typedef typename super::left_gmc_ptrtype test_left_gmc_ptrtype;
        typedef typename super::right_gmc_ptrtype test_right_gmc_ptrtype;
        typedef typename super::map_left_gmc_type test_map_left_gmc_type;
        typedef typename super::map_right_gmc_type test_map_right_gmc_type;


        typedef FormContextBase<GeomapTrialContext,IM,GeomapExprContext> super2;

        typedef typename super2::map_geometric_mapping_context_type map_trial_geometric_mapping_context_type;
        typedef typename super2::geometric_mapping_context_ptrtype trial_geometric_mapping_context_ptrtype;
        typedef typename super2::geometric_mapping_context_type trial_geometric_mapping_context_type;
        typedef typename super2::geometric_mapping_type trial_geometric_mapping_type;
        typedef typename super2::map_size trial_map_size;
        typedef typename super2::gmc1 trial_gmc1;
        typedef typename super2::left_gmc_ptrtype trial_left_gmc_ptrtype;
        typedef typename super2::right_gmc_ptrtype trial_right_gmc_ptrtype;
        typedef typename super2::map_left_gmc_type trial_map_left_gmc_type;
        typedef typename super2::map_right_gmc_type trial_map_right_gmc_type;





        static const uint16_type nDim = test_geometric_mapping_type::nDim;

        typedef ExprT expression_type;

        typedef typename test_precompute<>::type test_precompute_type;
        typedef typename test_precompute<>::ptrtype test_precompute_ptrtype;
        typedef typename trial_precompute<>::type trial_precompute_type;
        typedef typename trial_precompute<>::ptrtype trial_precompute_ptrtype;

        typedef typename FE2::fe_type trial_fe_type;
        typedef boost::shared_ptr<trial_fe_type> trial_fe_ptrtype;
        typedef typename trial_fe_type::template Context< trial_geometric_mapping_context_type::context,
                                                          trial_fe_type,
                                                          trial_geometric_mapping_type,
                                                          mesh_element_2_type> trial_fecontext_type;
        typedef boost::shared_ptr<trial_fecontext_type> trial_fecontext_ptrtype;
        typedef typename mpl::if_<mpl::equal_to<trial_map_size, mpl::int_<1> >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype>,
                                                            fusion::pair<gmc<1>, trial_fecontext_ptrtype> > > >::type::type map_trial_fecontext_type;

        typedef fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > map_left_trial_fecontext_type;
        typedef fusion::map<fusion::pair<trial_gmc1, trial_fecontext_ptrtype> > map_right_trial_fecontext_type;

        typedef typename FE1::fe_type test_fe_type;
        typedef boost::shared_ptr<test_fe_type> test_fe_ptrtype;
        typedef typename test_fe_type::template Context< test_geometric_mapping_context_type::context,
                                                         test_fe_type,
                                                         test_geometric_mapping_type,
                                                         mesh_element_1_type> test_fecontext_type;
        typedef boost::shared_ptr<test_fecontext_type> test_fecontext_ptrtype;
        typedef typename mpl::if_<mpl::equal_to<test_map_size, mpl::int_<1> >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > >,
                                  mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype>,
                                                            fusion::pair<gmc<1>, test_fecontext_ptrtype> > > >::type::type map_test_fecontext_type;
        typedef fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > map_left_test_fecontext_type;
        typedef fusion::map<fusion::pair<test_gmc1, test_fecontext_ptrtype> > map_right_test_fecontext_type;

        //--------------------------------------------------------------------------------------//

        typedef typename super::map_geometric_mapping_expr_context_type map_geometric_mapping_expr_context_type;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                                               map_left_test_fecontext_type,
                                               map_left_trial_fecontext_type> eval00_expr_type;
        typedef boost::shared_ptr<eval00_expr_type> eval00_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                                               map_left_test_fecontext_type,
                                               map_right_trial_fecontext_type> eval01_expr_type;
        typedef boost::shared_ptr<eval01_expr_type> eval01_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                                               map_right_test_fecontext_type,
                                               map_left_trial_fecontext_type> eval10_expr_type;
        typedef boost::shared_ptr<eval10_expr_type> eval10_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                                               map_right_test_fecontext_type,
                                               map_right_trial_fecontext_type> eval11_expr_type;
        typedef boost::shared_ptr<eval11_expr_type> eval11_expr_ptrtype;


        static const int rep_shape = 4;//2+(eval_expr_type::shape::M-1>0)+(eval_expr_type::shape::N-1>0);
        //typedef boost::multi_array<value_type, rep_shape> local_matrix_type;

        typedef typename FE1::dof_type test_dof_type;
        typedef typename FE2::dof_type trial_dof_type;
        static const int nDofPerElementTest = FE1::dof_type::nDofPerElement;
        static const int nDofPerElementTrial = FE2::dof_type::nDofPerElement;
        static const int nDofPerComponentTest = FE1::dof_type::fe_type::nLocalDof;
        static const int nDofPerComponentTrial = FE2::dof_type::fe_type::nLocalDof;
        static const int local_mat_traits = mpl::if_<mpl::equal_to<mpl::int_<nDofPerElementTrial>,mpl::int_<1> >,
                                                     mpl::int_<Eigen::ColMajor>,
                                                     mpl::int_<Eigen::RowMajor> >::type::value;
        typedef Eigen::Matrix<value_type, nDofPerElementTest, nDofPerElementTrial,local_mat_traits> local_matrix_type;
        typedef Eigen::Matrix<value_type, 2*nDofPerElementTest, 2*nDofPerElementTrial,Eigen::RowMajor> local2_matrix_type;
        typedef Eigen::Matrix<value_type, nDofPerComponentTest, nDofPerComponentTrial,local_mat_traits> c_local_matrix_type;
        typedef Eigen::Matrix<value_type, 2*nDofPerComponentTest, 2*nDofPerComponentTrial,Eigen::RowMajor> c_local2_matrix_type;
        typedef Eigen::Matrix<int, nDofPerElementTest, 1> local_row_type;
        typedef Eigen::Matrix<int, 2*nDofPerElementTest, 1> local2_row_type;
        typedef Eigen::Matrix<int, nDofPerElementTrial, 1> local_col_type;
        typedef Eigen::Matrix<int, 2*nDofPerElementTrial, 1> local2_col_type;

        typedef Eigen::Matrix<int, nDofPerComponentTest, 1> c_local_row_type;
        typedef Eigen::Matrix<int, 2*nDofPerComponentTest, 1> c_local2_row_type;
        typedef Eigen::Matrix<int, nDofPerComponentTrial, 1> c_local_col_type;
        typedef Eigen::Matrix<int, 2*nDofPerComponentTrial, 1> c_local2_col_type;


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
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2 );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2,
                 mpl::int_<2> );

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& _gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr );

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const& gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad );

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& _gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr,
                     mpl::int_<2> );

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr,
                     IM const& im )
        {
            M_integrator = im;
            update( gmcTest, gmcTrial, gmcExpr );
        }

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const& gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        IM const& im,
                                        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad )
        {
            M_integrator = im;
            updateInCaseOfInterpolate( gmcTest, gmcTrial, gmcExpr, indexLocalToQuad);
        }

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr,
                     IM const& im, mpl::int_<2> )
        {
            M_integrator = im;
            update( gmcTest, gmcTrial, gmcExpr, mpl::int_<2>() );
        }

        void integrate()
        {
            integrate( mpl::int_<fusion::result_of::size<GeomapTestContext>::type::value>() );
        }
        void integrateInCaseOfInterpolate(std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                          bool isFirstExperience )
        {
            integrateInCaseOfInterpolate( mpl::int_<fusion::result_of::size<GeomapTestContext>::type::value>(),
                                          indexLocalToQuad,
                                          isFirstExperience );
        }


        void assemble() { assemble( fusion::at_key<gmc<0> >( _M_test_gmc )->id());  }
        void assemble( size_type elt_0 );

        void assemble( mpl::int_<2> )
            {
                assemble( fusion::at_key<gmc<0> >( _M_test_gmc )->id(),
                          fusion::at_key<test_gmc1>( _M_test_gmc )->id() );
            }
        void assemble( size_type elt_0, size_type elt_1  );

        void assembleInCaseOfInterpolate();

        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points
         */
        template<typename Pts>
        void precomputeBasisAtPoints( Pts const& pts )
            {
                _M_test_pc = test_precompute_ptrtype( new test_precompute_type( _M_form.testSpace()->fe(), pts ) );
                _M_trial_pc = trial_precompute_ptrtype( new trial_precompute_type( _M_form.trialSpace()->fe(), pts ) );
            }

        template<typename PtsTest,typename PtsTrial>
        void precomputeBasisAtPoints( PtsTest const& pts1,PtsTrial const& pts2  )
            {
                _M_test_pc = test_precompute_ptrtype( new test_precompute_type( _M_form.testSpace()->fe(), pts1 ) );
                _M_trial_pc = trial_precompute_ptrtype( new trial_precompute_type( _M_form.trialSpace()->fe(), pts2 ) );
            }


        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points on the face of
         * the reference element
         */
        template<typename Pts>
        void precomputeBasisAtPoints( uint16_type __f, permutation_1_type const& __p, Pts const& pts )
            {
                _M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( _M_form.testSpace()->fe(), pts ) );
                //FEEL_ASSERT( _M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );

                _M_trial_pc_face[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( _M_form.trialSpace()->fe(), pts ) );
                //FEEL_ASSERT( _M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
            }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts )
        {
            typedef typename boost::is_same< permutation_1_type, typename QuadMapped<PtsSet>::permutation_type>::type is_same_permuation_type;
            return precomputeTestBasisAtPoints(pts, mpl::bool_<is_same_permuation_type::value>() );
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts, mpl::bool_<false> )
        {
            std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> > testpc;
            return testpc;
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts, mpl::bool_<true> )
            {
                QuadMapped<PtsSet> qm;
                typedef typename QuadMapped<PtsSet>::permutation_type permutation_type;
                typename QuadMapped<PtsSet>::permutation_points_type ppts( qm( pts ) );

                std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > testpc;

                for ( uint16_type __f = 0; __f < pts.nFaces(); ++__f )
                {
                    for( permutation_type __p( permutation_type::IDENTITY );
                         __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                    {
                        testpc[__f][__p] = test_precompute_ptrtype( new test_precompute_type( _M_form.testSpace()->fe(), ppts[__f].find(__p)->second ) );
                    }
                }
                return testpc;
            }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> >
        precomputeTrialBasisAtPoints( PtsSet const& pts )
        {
            typedef typename boost::is_same< permutation_2_type, typename QuadMapped<PtsSet>::permutation_type>::type is_same_permuation_type;
            return precomputeTrialBasisAtPoints(pts, mpl::bool_<is_same_permuation_type::value>() );
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> >
        precomputeTrialBasisAtPoints( PtsSet const& pts,  mpl::bool_<false> )
        {
            std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> > trialpc;
            return trialpc;
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> >
        precomputeTrialBasisAtPoints( PtsSet const& pts, mpl::bool_<true> )
            {
                QuadMapped<PtsSet> qm;
                typedef typename QuadMapped<PtsSet>::permutation_type permutation_type;
                typename QuadMapped<PtsSet>::permutation_points_type ppts( qm( pts ) );

                std::map<uint16_type, std::map<permutation_type,trial_precompute_ptrtype> > trialpc;

                for ( uint16_type __f = 0; __f < pts.nFaces(); ++__f )
                {
                    for( permutation_type __p( permutation_type::IDENTITY );
                         __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                    {
                        trialpc[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( _M_form.trialSpace()->fe(), ppts[__f].find(__p)->second ) );
                    }
                }
                return trialpc;
            }

        /**
         * Return the structure that holds the test basis functions
         * evaluated at a previously given set of points on the face of
         * the reference element
         *
         * \see precomputeBasisAtPoints()
         */
        test_precompute_ptrtype const& testPc( uint16_type __f,
                                               permutation_1_type __p = permutation_1_type( permutation_1_type::NO_PERMUTATION ) ) const
            {
                if ( __f == invalid_uint16_type_value )
                    return  _M_test_pc;
                //FEEL_ASSERT( _M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );
                return _M_test_pc_face.find(__f )->second.find( __p )->second;
            }

        /**
         * Return the structure that holds the trial basis functions
         * evaluated at a previously given set of points on the face of
         * the reference element
         * \see precomputeBasisAtPoints()
         */
        trial_precompute_ptrtype const& trialPc( uint16_type __f,
                                                 permutation_2_type __p = permutation_2_type( permutation_2_type::NO_PERMUTATION ) ) const
            {
                if ( __f == invalid_uint16_type_value )
                    return  _M_trial_pc;
                //FEEL_ASSERT( _M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
                return _M_trial_pc_face.find(__f )->second.find( __p )->second;
            }


    private:

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr,
                     mpl::bool_<false> );

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr,
                     mpl::bool_<true> );

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const& gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        mpl::bool_<false> );

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const& gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        mpl::bool_<true> );


        void integrate( mpl::int_<1> );

        void integrate( mpl::int_<2> );

        void integrateInCaseOfInterpolate( mpl::int_<1>,
                                           std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                           bool isFirstExperience );
    private:

        form_type& _M_form;
        const list_block_type& _M_lb;
        dof_1_type* _M_test_dof;
        dof_2_type* _M_trial_dof;

        test_precompute_ptrtype _M_test_pc;
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> > _M_test_pc_face;
        trial_precompute_ptrtype _M_trial_pc;
        std::map<uint16_type, std::map<permutation_2_type, trial_precompute_ptrtype> > _M_trial_pc_face;


        map_test_geometric_mapping_context_type _M_test_gmc;
        map_trial_geometric_mapping_context_type _M_trial_gmc;

        map_test_fecontext_type _M_test_fec;
        map_left_test_fecontext_type _M_test_fec0;
        map_right_test_fecontext_type _M_test_fec1;
        map_trial_fecontext_type _M_trial_fec;
        map_left_trial_fecontext_type _M_trial_fec0;
        map_right_trial_fecontext_type _M_trial_fec1;


        local_matrix_type _M_rep;
        local2_matrix_type _M_rep_2;
        local_row_type M_local_rows;
        local2_row_type M_local_rows_2;
        local_col_type M_local_cols;
        local2_col_type M_local_cols_2;
        local_row_type M_local_rowsigns;
        local2_row_type M_local_rowsigns_2;
        local_col_type M_local_colsigns;
        local2_col_type M_local_colsigns_2;

        c_local_matrix_type M_c_rep;
        c_local2_matrix_type M_c_rep_2;
        c_local_row_type M_c_local_rows;
        c_local2_row_type M_c_local_rows_2;
        c_local_col_type M_c_local_cols;
        c_local2_col_type M_c_local_cols_2;
        c_local_row_type M_c_local_rowsigns;
        c_local2_row_type M_c_local_rowsigns_2;
        c_local_col_type M_c_local_colsigns;
        c_local2_col_type M_c_local_colsigns_2;

        eval00_expr_ptrtype _M_eval_expr00;
        eval01_expr_ptrtype _M_eval_expr01;
        eval10_expr_ptrtype _M_eval_expr10;
        eval11_expr_ptrtype _M_eval_expr11;

        IM M_integrator;

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
        M_pattern( __vf.M_pattern ),
        _M_X1( __vf._M_X1 ),
        _M_X2( __vf._M_X2 ),
        _M_matrix( __vf._M_matrix ),
        _M_lb( __vf._M_lb ),
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
            M_pattern = form.M_pattern;
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
     * return the pattern
     */
    size_type pattern() const { return M_pattern; }

    /**
     * \return true if the pattern is coupled with respect to the components,
     * false otherwise
     */
    bool isPatternCoupled() const { Feel::Context ctx( M_pattern ); return ctx.test( DOF_PATTERN_COUPLED ); }

    /**
     * \return true if the pattern is the default one, false otherwise
     */
    bool isPatternDefault() const { Feel::Context ctx( M_pattern ); return ctx.test( DOF_PATTERN_DEFAULT ); }

    /**
     * \return true if the pattern adds the neighboring elements, false otherwise
     */
    bool isPatternNeighbor() const { Feel::Context ctx( M_pattern ); return ctx.test( DOF_PATTERN_NEIGHBOR ); }

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
    gm_1_ptrtype const& gm() const { return _M_X1->gm(); }

    /**
     * Geometric transformation
     */
    gm1_1_ptrtype const& gm1() const { return _M_X1->gm1(); }

    /**
     * Geometric transformation
     */
    gm_2_ptrtype const& gmTrial() const { return _M_X2->gm(); }

    /**
     * Geometric transformation
     */
    gm1_2_ptrtype const& gm1Trial() const { return _M_X2->gm1(); }


    /**
     * \return the matrix associated to the bilinear form
     */
    matrix_type const& matrix() const { return _M_matrix; }

    matrix_type& matrix() { return _M_matrix; }

    list_block_type const& blockList() const { return _M_lb; }

    /**
     * \return the threshold
     */
    value_type threshold() const { return _M_threshold; }

    /**
     * \return \c true if threshold applies, false otherwise
     */
    bool doThreshold( value_type const& v ) const { return ( math::abs( v ) > _M_threshold ); }

    /**
     * return true if do threshold. false otherwise
     */
    bool doThreshold() const { return _M_do_threshold; }

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
                   Feel::Context const& on_context );

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
     * add value \p v at position (\p i, \p j) of the matrix
     * associated with the bilinear form
     */
    void addMatrix( int* rows, int nrows,
                    int* cols, int ncols,
                    value_type* data )
        {
            _M_matrix.addMatrix( rows, nrows, cols, ncols, data );
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

    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> );
    graph_ptrtype computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> );

    //@}

protected:

    BilinearForm();

private:


    template <class ExprT>
    void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<false> );

    template <class ExprT>
    void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<true> );
private:

    size_type M_pattern;

    space_1_ptrtype _M_X1;
    space_2_ptrtype _M_X2;

    matrix_type& _M_matrix;

    bool _M_do_build;

    list_block_type _M_lb;


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
    M_pattern( graph_hints ),
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
    boost::timer tim;
    Debug( 5050 ) << "begin constructor with default listblock\n";

    _M_lb.push_back( Block ( 0, 0, 0, 0 ) );

    if ( build )
        {
            const size_type n1_dof_on_proc = _M_X1->nLocalDof();
            M_n_nz.resize (n1_dof_on_proc);
            M_n_oz.resize (n1_dof_on_proc);

            boost::timer t;
            Debug( 5050 ) << "compute graph\n";

            graph_ptrtype graph;
            if ( dynamic_cast<void*>( _M_X1->mesh().get()) == dynamic_cast<void*>( _M_X2->mesh().get()) )
                graph = computeGraph( graph_hints, mpl::bool_<mpl::and_< mpl::bool_< (FE1::nSpaces == 1)>,
                                                                         mpl::bool_< (FE2::nSpaces == 1)> >::type::value >() );
            else
                graph = computeGraphInCaseOfInterpolate( graph_hints, mpl::bool_<mpl::and_< mpl::bool_< (FE1::nSpaces == 1)>,
                                                                                            mpl::bool_< (FE2::nSpaces == 1)> >::type::value >() );

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
    Debug( 5050 ) << " - form init in " << tim.elapsed() << "\n";
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
    M_pattern( graph_hints ),
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
    if (FE1::nSpaces > 1 && FE2::nSpaces > 1 )
        fusion::for_each( _M_X1->functionSpaces(), make_bfassign2( *this, __expr ) );
    else if (FE1::nSpaces > 1 )
        fusion::for_each( _M_X1->functionSpaces(),
                          make_bfassign3( *this, __expr, _M_X2, 0 ) );
    else
        fusion::for_each( _M_X2->functionSpaces(),
                          make_bfassign1( *this, __expr, _M_X1, 0 ) );

    Debug( 5050 ) << "BilinearForm::assign() stop loop on test spaces\n";
}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator=( Expr<ExprT> const& __expr )
{
    // loop(fusion::for_each) over sub-functionspaces in SpaceType
    // pass expression and initialize
    this->assign( __expr, true, mpl::bool_<mpl::or_< mpl::bool_< (FE1::nSpaces > 1)>,
                                           mpl::bool_< (FE2::nSpaces > 1)> >::type::value >() );
    return *this;
}

template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator+=( Expr<ExprT> const& __expr )
{
    Debug( 5055 ) << "[BilinearForm::operator+=] start\n";
    this->assign( __expr, false, mpl::bool_<mpl::or_< mpl::bool_< (FE1::nSpaces > 1)>,
                                            mpl::bool_< (FE2::nSpaces > 1)> >::type::value >() );
    Debug( 5055 ) << "[BilinearForm::operator+=] stop\n";
    return *this;
}


template<typename FE1,  typename FE2, typename ElemContType>
void
BilinearForm<FE1,FE2,ElemContType>::zeroRows( std::vector<int> __dofs,
                                              std::vector<value_type> __values,
                                              Vector<value_type>& rhs,
                                              Feel::Context const& on_context )
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
    FEEL_ASSERT( std::unique (begin,  middle) == middle )
        ( *begin )( *middle ).error( "duplicate dof(begin,middle)" );
    FEEL_ASSERT (std::unique (middle, end)    == end)
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
    boost::timer tim;
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
                //auto row1_entries = boost::unwrap_ref( boost::ref( M_graph->row(theglobalrow).template get<2>() ) );
                std::set<size_type>& row1_entries = M_graph->row(theglobalrow).template get<2>();
                std::set<size_type> const& row2_entries = boost::get<2>( it->second );

                Debug( 5050 ) << "[mergeGraph] adding information to global row [" << theglobalrow << "], localrow=" << thelocalrow << "\n";
                M_graph->row(theglobalrow).template get<1>() = thelocalrow;

                if ( row1_entries.empty() )
                {
                    // if row is empty then no need to shift the dof in
                    // composite case since the merge in done block-row-wise
                    row1_entries = row2_entries;
                }
                else
                {
                    // ensure unique sorted ids
                    auto itg = boost::prior(row1_entries.end());
                    // shift dofs in case of composite spaces
                    std::for_each( row2_entries.begin(), row2_entries.end(),[&]( size_type o ){ itg = row1_entries.insert( itg, o+col); });
                }
            }
        }
    Debug( 5050 ) << " -- merge_graph (" << row << "," << col << ") in " << tim.elapsed() << "\n";
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

    Feel::Context graph( hints );
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
                                    FEEL_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                                    FEEL_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
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

    Debug( 5050 ) << "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
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

    typedef typename mesh_1_type::element_const_iterator mesh_element_const_iterator;

    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    Feel::Context graph( hints );
    // If the user did not explicitly specify the DOF coupling
    // then all the DOFS are coupled to each other.  Furthermore,
    // we can take a shortcut and do this more quickly here.  So
    // we use an if-test.
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_DEFAULT )=" <<  graph.test ( DOF_PATTERN_DEFAULT ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_COUPLED )=" <<  graph.test ( DOF_PATTERN_COUPLED ) << "\n";
    Debug( 5050 ) << "[computeGraph]  : graph.test ( DOF_PATTERN_NEIGHBOR)=" <<  graph.test ( DOF_PATTERN_NEIGHBOR ) << "\n";
    bool do_less =  ( ( graph.test( DOF_PATTERN_DEFAULT ) &&
                        ( _M_X1->dof()->nComponents ==
                          _M_X2->dof()->nComponents ) ) &&
                      !graph.test( DOF_PATTERN_COUPLED ) );
    std::vector<size_type>
        element_dof2(_M_X2->dof()->getIndicesSize()),
        neighbor_dof;

    for ( ; elem_it != elem_en; ++elem_it)
    {
#if !defined(NDEBUG)
        Debug( 5050 ) << "[BilinearForm::computePatter] element " << elem_it->id() << " on proc " << elem_it->processId() << "\n";
#endif /* NDEBUG */
        mesh_element_1_type const& elem = *elem_it;

        // Get the global indices of the DOFs with support on this element
        //element_dof1 = _M_X1->dof()->getIndices( elem.id() );
        _M_X2->dof()->getIndicesSet( elem.id(), element_dof2 );

        // We can be more efficient if we sort the element DOFs
        // into increasing order
        //std::sort(element_dof1.begin(), element_dof1.end());
        std::sort(element_dof2.begin(), element_dof2.end());

        //const uint16_type  n1_dof_on_element = element_dof1.size();
        const uint16_type  n1_dof_on_element = _M_X1->dof()->getIndicesSize();
        const uint16_type  n2_dof_on_element = element_dof2.size();

        for (size_type i=0; i<n1_dof_on_element; i++)
        //BOOST_FOREACH( auto ig1, _M_X1->dof()->getIndices( elem.id() ) )
        {
            const size_type ig1 = _M_X1->dof()->localToGlobalId( elem.id(), i );
            //const size_type ig1 = element_dof1[i];
            const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
            const int ncomp1 = i / ndofpercomponent1;
            const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

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
                FEEL_ASSERT (ig1 >= first1_dof_on_proc )( ig1 )( first1_dof_on_proc ).error ("invalid dof index");
                FEEL_ASSERT ((ig1 - first1_dof_on_proc) < sparsity_graph->size() )
                    ( ig1 )( first1_dof_on_proc )( sparsity_graph->size() ).error( "invalid dof index" );
#endif
                graph_type::row_type& row = sparsity_graph->row(ig1);
                bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                Debug( 5051 ) << "work with row " << ig1 << " local index " << ig1 - first1_dof_on_proc << "\n";

                if ( do_less )
                {
                    if ( ncomp1 == (_M_X2->dof()->nComponents-1) )
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.end() );
                    else
                        row.get<2>().insert( element_dof2.begin()+ncomp1*ndofpercomponent2,
                                             element_dof2.begin()+(ncomp1+1)*ndofpercomponent2 );
                }
                else
                {
                    row.get<2>().insert( element_dof2.begin(), element_dof2.end() );
                }
                // Now (possibly) add dof from neighboring elements
                if ( graph.test( DOF_PATTERN_NEIGHBOR ) )
                {
                    for (uint16_type ms=0; ms < elem.nNeighbors(); ms++)
                    {
                        mesh_element_1_type const* neighbor = NULL;
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
                                if ( do_less )
                                {
                                    if ( ncomp1 == (_M_X2->dof()->nComponents-1) )
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.end() );
                                    else
                                        row.get<2>().insert( neighbor_dof.begin()+ncomp1*ndofpercomponent2,
                                                             neighbor_dof.begin()+(ncomp1+1)*ndofpercomponent2 );

                                }
                                else
                                {
                                    row.get<2>().insert( neighbor_dof.begin(), neighbor_dof.end() );
                                }

                            } // neighbor_id
                        }

                    } // neighbor graph
                }
            } // only dof on proc

        }// dof loop
    } // element iterator loop
    Debug( 5050 )<< "[computeGraph<true>] before calling close in " << t.elapsed() << "s\n";
    sparsity_graph->close();
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    Debug( 5050 ) << "[computeGraph<true>] done in " << t.elapsed() << "s\n";
    return sparsity_graph;
}
#endif

template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false> )
{
    boost::timer t;
    Debug( 5050 ) << "compute graph for composite bilinear form with interpolation\n";

    auto graph = computeGraphInCaseOfInterpolate( hints, mpl::bool_< ( FE1::nSpaces > 1)>(), mpl::bool_< ( FE2::nSpaces > 1)>() );

    Debug( 5050 ) << "closing graph for composite bilinear form with interpolation done in " << t.elapsed() << "s\n"; t.restart();
    graph->close();
    Debug( 5050 ) << "compute graph for composite bilinear form done in " << t.elapsed() << "s\n";

    return graph;
}

template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<true> )
{
    fusion::for_each( _M_X1->functionSpaces(), compute_graph1<self_type>( *this, hints ) );

    return M_graph;
}

template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true>, mpl::bool_<false> )
{
    fusion::for_each( _M_X1->functionSpaces(),
                      compute_graph3<self_type,space_2_type>( *this, _M_X2, 0, hints ) );
    return M_graph;
}

template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<false>, mpl::bool_<true> )
{
    fusion::for_each( _M_X2->functionSpaces(),
                      compute_graph2<self_type,space_1_type>( *this, _M_X1, 0, hints ) );
    return M_graph;
}


template<typename FE1,  typename FE2, typename ElemContType>
typename BilinearForm<FE1,FE2,ElemContType>::graph_ptrtype
BilinearForm<FE1,FE2,ElemContType>::computeGraphInCaseOfInterpolate( size_type hints, mpl::bool_<true> )
{

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

    typedef typename mesh_1_type::element_const_iterator mesh_element_const_iterator;
    mesh_element_const_iterator       elem_it  = _M_X1->mesh()->beginElementWithProcessId( proc_id );
    const mesh_element_const_iterator elem_en  = _M_X1->mesh()->endElementWithProcessId( proc_id );

    auto locTool = _M_X2->mesh()->tool_localization();
    locTool->updateForUse();
    locTool->setExtrapolation(false);

    std::vector<size_type>
        element_dof1,
        element_dof2;

    for ( ; elem_it != elem_en; ++elem_it)
    {
        mesh_element_1_type const& elem = *elem_it;

        // Get the global indices of the DOFs with support on this element
        element_dof1 = _M_X1->dof()->getIndices( elem.id() );

        const uint16_type  n1_dof_on_element = element_dof1.size();
        //const uint16_type  n2_dof_on_element = element_dof2.size();

        for (size_type i=0; i<n1_dof_on_element; i++)
            {
                const size_type ig1 = element_dof1[i];
                const int ndofpercomponent1 = n1_dof_on_element / _M_X1->dof()->nComponents;
                const int ncomp1 = i / ndofpercomponent1;
                //const int ndofpercomponent2 = n2_dof_on_element / _M_X2->dof()->nComponents;

                auto ptRealDof = boost::get<0>(_M_X1->dof()->dofPoint(ig1));

                //std::cout << "\n Pt dof " << ptRealDof;

                auto res = locTool->searchElements(ptRealDof);
                auto hasFind = res.get<0>();

                if ( hasFind )
                    {
                        auto res_it = res.get<1>().begin();
                        auto res_en = res.get<1>().end();
                        for ( ; res_it != res_en ; ++res_it)
                            {
                                element_dof2 = _M_X2->dof()->getIndices( res_it->get<0>() );

                                graph_type::row_type& row = sparsity_graph->row(ig1);
                                bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                                row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                                row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;

                                //if ( do_less ) {}
                                //else
                                row.get<2>().insert( element_dof2.begin(), element_dof2.end() );


                            }

                    } // if (hasFind)
                else
                    {
                        // row empty
                        graph_type::row_type& row = sparsity_graph->row(ig1);
                        bool is_on_proc = ( ig1 >= first1_dof_on_proc) && (ig1 <= last1_dof_on_proc);
                        row.get<0>() = is_on_proc?proc_id:invalid_size_type_value;
                        row.get<1>() = is_on_proc?ig1 - first1_dof_on_proc:invalid_size_type_value;
                    }

            } // for (size_type i=0; i<n1_dof_on_element; i++)



    } // for ( ; elem_it ... )

    locTool->setExtrapolation(true);

    sparsity_graph->close();

    return sparsity_graph;
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

    Feel::vf::Block block ( 0, 0,
                            _M_bf.testSpace()->nDofStart( _M_test_index ),
                            _M_bf.trialSpace()->nDofStart( _M_trial_index ) );
    Debug( 5050 ) << "[BFAssign1::operator()] block: " << block << "\n";
    Feel::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Feel::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;

    bf_type bf( _M_test,trial, _M_bf.matrix(), list_block, _M_bf.doThreshold(), _M_bf.threshold(), _M_bf.pattern()  );

    bf += _M_expr;

    ++_M_trial_index;
}

template<typename BFType, typename ExprType, typename TrialSpaceType>
template<typename SpaceType>
void BFAssign3<BFType,ExprType,TrialSpaceType>::operator()( boost::shared_ptr<SpaceType> const& test ) const
{
    Debug(5050) << "[BFAssign3::operator()] expression has trial functions index "
                << _M_test_index << " : "
                << ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result << " (0:no, 1:yes)\n";
    Debug(5050) << "[BFAssign3::operator()] expression has test functions index "
                << _M_trial_index << " :"
                << ExprType::template HasTrialFunction<typename TrialSpaceType::reference_element_type>::result << " (0:no, 1:yes)\n";

    if ( !ExprType::template HasTrialFunction<typename TrialSpaceType::reference_element_type>::result ||
         !ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result )
        {
            ++_M_test_index;
            return;
        }

    Debug( 5050 ) << "[BFAssign3::operator()] terms found with "
                  << "testindex: " << _M_test_index << " trialindex: " << _M_trial_index << "\n";
    typedef SpaceType test_space_type;
    typedef TrialSpaceType trial_space_type;

    Feel::vf::Block block ( 0, 0,
                            _M_bf.testSpace()->nDofStart( _M_test_index ),
                            _M_bf.trialSpace()->nDofStart( _M_trial_index ) );
    Debug( 5050 ) << "[BFAssign1::operator()] block: " << block << "\n";
    Feel::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Feel::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;

    bf_type bf( test, _M_trial, _M_bf.matrix(), list_block, _M_bf.doThreshold(), _M_bf.threshold(), _M_bf.pattern() );

    bf += _M_expr;

    ++_M_test_index;
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

    Feel::vf::Block block ( 0, 0,
                            M_bf.testSpace()->nDofStart( M_test_index ),
                            M_bf.trialSpace()->nDofStart( M_trial_index ) );
    Debug( 5050 ) << "[compute_graph2::operator()] block: " << block << "\n";
    Feel::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Feel::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;


    bf_type bf( M_space1, trial, M_bf.matrix(), list_block );

    typename bf_type::graph_ptrtype graph;

    if ( dynamic_cast<void*>( bf.testSpace()->mesh().get()) == dynamic_cast<void*>( bf.trialSpace()->mesh().get()) )
        graph = bf.computeGraph( M_hints, mpl::bool_<true>() );
    else
        graph = bf.computeGraphInCaseOfInterpolate( M_hints, mpl::bool_<true>() );


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

template<typename BFType, typename TrialSpaceType>
template<typename SpaceType>
void
compute_graph3<BFType,TrialSpaceType>::operator()( boost::shared_ptr<SpaceType> const& test ) const
{
    boost::timer t;
    Debug( 5050 ) << "[compute_graph2::operator()] begin ==================================================\n";
    Debug( 5050 ) << "[compute_graph2::operator()] testindex: " << M_test_index << "\n";
    Debug( 5050 ) << "[compute_graph2::operator()] trialindex: " << M_trial_index << "\n";
    typedef SpaceType test_space_type;
    typedef TrialSpaceType trial_space_type;

    Feel::vf::Block block ( 0, 0,
                            M_bf.testSpace()->nDofStart( M_test_index ),
                            M_bf.trialSpace()->nDofStart( M_trial_index ) );
    Debug( 5050 ) << "[compute_graph2::operator()] block: " << block << "\n";
    Feel::vf::list_block_type list_block;
    list_block.push_back( block );

    typedef typename BFType::matrix_type matrix_type;
    typedef Feel::vf::detail::BilinearForm<test_space_type,
        trial_space_type,
        ublas::vector_range<ublas::vector<double> > > bf_type;


    bf_type bf( test, M_space1, M_bf.matrix(), list_block );

    typename bf_type::graph_ptrtype graph;

    if ( dynamic_cast<void*>( bf.testSpace()->mesh().get()) == dynamic_cast<void*>( bf.trialSpace()->mesh().get()) )
        graph = bf.computeGraph( M_hints, mpl::bool_<true>() );
    else
        graph = bf.computeGraphInCaseOfInterpolate( M_hints, mpl::bool_<true>() );


    Debug( 5050 ) << "[compute_graph2::operator()] testindex: " << M_test_index << " row space index starts at " << M_bf.testSpace()->nDofStart( M_test_index ) << "\n";
    Debug( 5050 ) << "[compute_graph2::operator()] trialindex: " << M_trial_index << " col space index starts at " << M_bf.trialSpace()->nDofStart( M_trial_index ) << "\n";
    M_bf.mergeGraph( M_bf.testSpace()->nDofStart( M_test_index ), M_bf.trialSpace()->nDofStart( M_trial_index ) , graph );


    Debug( 5050 ) << "[compute_graph2] test index " << M_test_index << " done in " << t.elapsed() << "s\n";
    ++M_test_index;

    Debug( 5050 ) << "[compute_graph2::operator()] end ==================================================\n";
}

} // detail
} // vf
/// \endcond
} // feel

#include <feel/feelvf/bilinearformcontext.hpp>

#endif /* __BilinearForm_H */

