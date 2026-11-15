/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2006-2010 Université Joseph Fourier (Grenoble I)
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
 *  \file bilinearform.hpp
 *  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *  \date 2005-01-18
 */
#ifndef __BilinearForm_H
#define __BilinearForm_H 1

#include <Eigen/Eigen>
#include <Eigen/StdVector>


#include <set>

#include <boost/smart_ptr/make_shared.hpp>
#include <boost/fusion/algorithm.hpp>

#include <feel/feelcore/context.hpp>
#include <feel/feelalg/matrixvalue.hpp>
#include <feel/feelalg/vectorublas.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/stencil.hpp>

#include <feel/feelvf/assemblystorage.hpp>
#include <feel/feelvf/bilinearformbase.hpp>
#include <feel/feelvf/linearform.hpp>
#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/block.hpp>
#include <feel/feelvf/fec.hpp>
#include <feel/feelvf/formcontextbase.hpp>


namespace Feel
{
namespace parameter = boost::parameter;
namespace fusion = boost::fusion;
namespace vf
{

/// \cond detail
template<typename FE1,typename FE2,typename ElemContType> class BilinearForm;
namespace detail
{


template<typename BFType, typename ExprType>
struct BFAssign2
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign2( BFAssign2 const& lfa )
        :
        M_bf( lfa.M_bf ),
        M_expr( lfa.M_expr ),
        M_test_index( lfa.M_test_index )
    {}
    BFAssign2( BFType& lf, ExprType const& expr )
        :
        M_bf( lf ),
        M_expr( expr ),
        M_test_index( 0 )
    {}
    template<typename SpaceType>
    void operator()( std::shared_ptr<SpaceType> const& X ) const
    {
        DVLOG(2) << "[BFAssign2::operator()] start loop on trial spaces against test space index: " << M_test_index << "\n";

        if ( M_bf.testSpace()->worldsComm()[M_test_index]->isActive() )
        {
            fusion::for_each( M_bf.trialSpace()->functionSpaces(),
                              make_bfassign1( M_bf, M_expr, X, M_test_index ) );
        }

        DVLOG(2) << "[BFAssign2::operator()] stop loop on trial spaces against test space index: " << M_test_index << "\n";
        ++M_test_index;

    }
private:
    BFType& M_bf;
    ExprType const& M_expr;
    mutable size_type M_test_index;
};
template<typename BFType, typename ExprType>
BFAssign2<BFType,ExprType>
make_bfassign2( BFType& lf, ExprType const& expr )
{
    return BFAssign2<BFType,ExprType>( lf, expr );
}

template<typename BFType, typename ExprType, typename TestSpaceType>
struct BFAssign1
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign1( BFAssign1 const& lfa )
        :
        M_bf( lfa.M_bf ),
        M_test( lfa.M_test ),
        M_expr( lfa.M_expr ),
        M_trial_index( lfa.M_trial_index ),
        M_test_index( lfa.M_test_index )
    {}
    BFAssign1( BFType& lf,
               ExprType const& expr,
               std::shared_ptr<TestSpaceType> const& Testh,
               size_type test_index )
        :
        M_bf( lf ),
        M_test( Testh ),
        M_expr( expr ),
        M_trial_index( 0 ),
        M_test_index( test_index )
    {}
    template<typename SpaceType>
    void operator()( std::shared_ptr<SpaceType> const& trial ) const;

private:
    BFType& M_bf;
    std::shared_ptr<TestSpaceType> M_test;
    ExprType const& M_expr;
    mutable size_type M_trial_index;
    size_type M_test_index;
};
template<typename BFType, typename ExprType, typename TestSpaceType>
BFAssign1<BFType,ExprType,TestSpaceType>
make_bfassign1( BFType& lf,
                ExprType const& expr,
                std::shared_ptr<TestSpaceType> const& test_space,
                size_type test_index )
{
    return BFAssign1<BFType,ExprType,TestSpaceType>( lf, expr, test_space, test_index );
}


template<typename BFType, typename ExprType, typename TrialSpaceType>
struct BFAssign3
{
    typedef typename BFType::matrix_type matrix_type;
    BFAssign3( BFAssign3 const& lfa )
        :
        M_bf( lfa.M_bf ),
        M_trial( lfa.M_trial ),
        M_expr( lfa.M_expr ),
        M_trial_index( lfa.M_trial_index ),
        M_test_index( lfa.M_test_index )
    {}
    BFAssign3( BFType& lf,
               ExprType const& expr,
               std::shared_ptr<TrialSpaceType> const& Trialh,
               size_type trial_index )
        :
        M_bf( lf ),
        M_trial( Trialh ),
        M_expr( expr ),
        M_trial_index( trial_index ),
        M_test_index( 0 )
    {}
    template<typename SpaceType>
    void operator()( std::shared_ptr<SpaceType> const& test ) const;

private:
    BFType& M_bf;
    std::shared_ptr<TrialSpaceType> M_trial;
    ExprType const& M_expr;
    size_type M_trial_index;
    mutable size_type M_test_index;
};
template<typename BFType, typename ExprType, typename TrialSpaceType>
BFAssign3<BFType,ExprType,TrialSpaceType>
make_bfassign3( BFType& lf,
                ExprType const& expr,
                std::shared_ptr<TrialSpaceType> const& trial_space,
                size_type trial_index )
{
    return BFAssign3<BFType,ExprType,TrialSpaceType>( lf, expr, trial_space, trial_index );
}




//!
//! BilinearForm class
//!
template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename FE1::value_type> >
class BilinearForm : public BilinearFormBase<typename FE1::value_type>
{
public:

    /** @name Typedefs
     */
    //@{
    enum { nDim = FE1::nDim };

    using value_type = typename FE1::value_type;
    using super = BilinearFormBase<value_type>;

    using space_1_type = functionspace_type<FE1>;
    typedef std::shared_ptr<space_1_type> space_1_ptrtype;
    typedef space_1_type test_space_type;
    typedef std::shared_ptr<space_1_type> test_space_ptrtype;

    using space_2_type = functionspace_type<FE2>;
    typedef std::shared_ptr<space_2_type> space_2_ptrtype;
    typedef space_2_type trial_space_type;
    typedef std::shared_ptr<space_2_type> trial_space_ptrtype;


    typedef typename space_1_type::template Element<value_type,ElemContType> element_1_type;

    typedef typename space_2_type::template Element<value_type,ElemContType> element_2_type;

    using self_type = BilinearForm<FE1, FE2, ElemContType>;


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

    using index_type = typename mesh_1_type::index_type;
    using size_type = typename mesh_1_type::size_type;

    using matrix_ptrtype = typename super::matrix_ptrtype;

    template<typename SpaceType, bool UseMortar = false>
    struct finite_element
    {
        typedef  typename mpl::if_<mpl::bool_<UseMortar&&SpaceType::is_mortar>,
                                   mpl::identity<typename SpaceType::mortar_fe_type>,
                                   mpl::identity<typename SpaceType::fe_type> >::type::type type;
        typedef std::shared_ptr<type> ptrtype;
    };

    template<int _N = 0, bool UseMortar = false>
    struct test_precompute
    {
        //typedef typename space_1_type::basis_0_type::template precompute<_N>::type type;
        typedef typename finite_element<space_1_type,UseMortar>::type::PreCompute type;
        typedef std::shared_ptr<type> ptrtype;
    };
    template<int _N = 0, bool UseMortar = false>
    struct trial_precompute
    {
        //typedef typename space_2_type::basis_0_type::template precompute<_N>::type type;
        //typedef typename space_2_type::basis_0_type::PreCompute type;
        typedef typename finite_element<space_2_type,UseMortar>::type::PreCompute type;
        typedef std::shared_ptr<type> ptrtype;
    };


    // return test finite element
    template<bool UseMortar=false>
    typename finite_element<space_1_type,UseMortar>::ptrtype
    testFiniteElement() const
        {
            if constexpr ( UseMortar && space_1_type::is_mortar )
                return std::make_shared<typename finite_element<space_1_type,UseMortar>::type>();
            else
                return M_X1->fe();
        }
    // return trial finite element
    template<bool UseMortar=false>
    typename finite_element<space_2_type,UseMortar>::ptrtype
    trialFiniteElement() const
        {
            if constexpr ( UseMortar && space_2_type::is_mortar )
                return std::make_shared<typename finite_element<space_2_type,UseMortar>::type>();
            else
                return M_X2->fe();
        }
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
             typename GeomapTrialContext = GeomapTestContext,
             int UseMortarType = 0
             >
    class Context //: public FormContextBase<GeomapTestContext,IM,GeomapExprContext>
    {
        typedef FormContextBase<GeomapTestContext,IM,GeomapExprContext> super;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef Context<GeomapTestContext,ExprT,IM,GeomapExprContext,GeomapTrialContext,UseMortarType> form_context_type;
        static inline const bool UseMortar = UseMortarType > 0;
        static inline const bool UseMortarTest = (UseMortarType == 1) || (UseMortarType == 3);
        static inline const bool UseMortarTrial = (UseMortarType == 2) || (UseMortarType == 3);
        using form_type = self_type;
        typedef typename space_1_type::dof_type dof_1_type;
        typedef typename space_2_type::dof_type dof_2_type;

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





        static inline const uint16_type nDim = test_geometric_mapping_type::nDim;

        static inline const uint16_type nDimTest = test_space_type::mesh_type::nDim;
        static inline const uint16_type nDimTrial = trial_space_type::mesh_type::nDim;
        static inline const uint16_type nDimDiffBetweenTestTrial = ( nDimTest > nDimTrial )? nDimTest-nDimTrial : nDimTrial-nDimTest;

        typedef ExprT expression_type;

        typedef typename test_precompute<0,UseMortar>::type test_precompute_type;
        typedef typename test_precompute<0,UseMortar>::ptrtype test_precompute_ptrtype;
        typedef typename trial_precompute<0,UseMortar>::type trial_precompute_type;
        typedef typename trial_precompute<0,UseMortar>::ptrtype trial_precompute_ptrtype;

        typedef typename mpl::if_<mpl::bool_<UseMortar&&space_2_type::is_mortar>,
                                  mpl::identity<typename space_2_type::mortar_fe_type>,
                                  mpl::identity<typename space_2_type::fe_type> >::type::type trial_fe_type;
        typedef std::shared_ptr<trial_fe_type> trial_fe_ptrtype;
        typedef typename trial_fe_type::template Context< expression_type::context, // trial_geometric_mapping_context_type::context,
                                                          trial_fe_type,
                                                          trial_geometric_mapping_type,
                                                          mesh_element_2_type,
                                                          0,
                                                          trial_geometric_mapping_context_type::subEntityCoDim> trial_fecontext_type;
        typedef std::shared_ptr<trial_fecontext_type> trial_fecontext_ptrtype;
        typedef typename mpl::if_<mpl::equal_to<trial_map_size, mpl::int_<1> >,
                mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > >,
                mpl::identity<fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype>,
                fusion::pair<gmc<1>, trial_fecontext_ptrtype> > > >::type::type map_trial_fecontext_type;

        typedef fusion::map<fusion::pair<gmc<0>, trial_fecontext_ptrtype> > map_left_trial_fecontext_type;
        typedef fusion::map<fusion::pair<trial_gmc1, trial_fecontext_ptrtype> > map_right_trial_fecontext_type;

        typedef typename mpl::if_<mpl::bool_<UseMortar&&space_1_type::is_mortar>,
                                  mpl::identity<typename space_1_type::mortar_fe_type>,
                                  mpl::identity<typename space_1_type::fe_type> >::type::type test_fe_type;
        typedef std::shared_ptr<test_fe_type> test_fe_ptrtype;
        typedef typename test_fe_type::template Context< expression_type::context, // test_geometric_mapping_context_type::context,
                                                         test_fe_type,
                                                         test_geometric_mapping_type,
                                                         mesh_element_1_type,
                                                         0,
                                                         test_geometric_mapping_context_type::subEntityCoDim> test_fecontext_type;
        typedef std::shared_ptr<test_fecontext_type> test_fecontext_ptrtype;
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
        typedef std::shared_ptr<eval00_expr_type> eval00_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                map_left_test_fecontext_type,
                map_right_trial_fecontext_type> eval01_expr_type;
        typedef std::shared_ptr<eval01_expr_type> eval01_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                map_right_test_fecontext_type,
                map_left_trial_fecontext_type> eval10_expr_type;
        typedef std::shared_ptr<eval10_expr_type> eval10_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type,
                map_right_test_fecontext_type,
                map_right_trial_fecontext_type> eval11_expr_type;
        typedef std::shared_ptr<eval11_expr_type> eval11_expr_ptrtype;


        static const int rep_shape = 4;//2+(eval_expr_type::shape::M-1>0)+(eval_expr_type::shape::N-1>0);
        //typedef boost::multi_array<value_type, rep_shape> local_matrix_type;

        typedef typename space_1_type::dof_type test_dof_type;
        typedef typename space_2_type::dof_type trial_dof_type;
        static constexpr bool hasDynamicOrder = test_dof_type::is_order_dynamic || trial_dof_type::is_order_dynamic;
        static const int nDofPerElementTest = space_1_type::dof_type::nDofPerElement;
        static const int nDofPerElementTrial = space_2_type::dof_type::nDofPerElement;
        static const int nDofPerComponentTest =
            finiteElementFixedLocalDofCountOrZero<test_fe_type>( 1 );
        static const int nDofPerComponentTrial =
            finiteElementFixedLocalDofCountOrZero<trial_fe_type>( 1 );
        // Keep the test and trial extent contracts independent in mixed-order forms.
        using ElementStorage = ::Feel::detail::AssemblyStorage<value_type,
                                                               test_dof_type::is_order_dynamic ? Dynamic : nDofPerElementTest,
                                                               trial_dof_type::is_order_dynamic ? Dynamic : nDofPerElementTrial>;
        using ComponentStorage = ::Feel::detail::AssemblyStorage<value_type,
                                                                 test_dof_type::is_order_dynamic ? Dynamic : nDofPerComponentTest,
                                                                 trial_dof_type::is_order_dynamic ? Dynamic : nDofPerComponentTrial>;

        /** Ordinary local element matrix. */
        using local_matrix_type = typename ElementStorage::Ordinary;
        /** Local matrix after eliminating one test mortar constraint. */
        using mortar_test_local_matrix_type = typename ElementStorage::MortarTest;
        /** Local matrix after eliminating one trial mortar constraint. */
        using mortar_trial_local_matrix_type = typename ElementStorage::MortarTrial;
        /** Two-element face contribution. */
        using local2_matrix_type = typename ElementStorage::TwoCell;
        /** One component block, with its own fixed or dynamic cardinalities. */
        using c_local_matrix_type = typename ComponentStorage::Ordinary;
        /** Test-mortar component block. */
        using c_mortar_test_local_matrix_type = typename ComponentStorage::MortarTest;
        /** Trial-mortar component block. */
        using c_mortar_trial_local_matrix_type = typename ComponentStorage::MortarTrial;
        /** Two-element component block. */
        using c_local2_matrix_type = typename ComponentStorage::TwoCell;

        // Orientation signs and index vectors share the same local extent policy.
        using local_row_sign_type = typename ElementStorage::RowSigns;
        using local_col_sign_type = typename ElementStorage::ColSigns;
        using local2_row_sign_type = typename ElementStorage::TwoCellRowSigns;
        using local2_col_sign_type = typename ElementStorage::TwoCellColSigns;
        // DofTable's output-buffer API owns dynamic index vectors; retain its capacity reuse.
        using local_row_type = Eigen::VectorXi;
        using local_col_type = Eigen::VectorXi;
        using local2_row_type = Eigen::VectorXi;
        using local2_col_type = Eigen::VectorXi;
        using mortar_test_local_row_type = Eigen::VectorXi;
        using mortar_trial_local_col_type = Eigen::VectorXi;
        using c_local_row_type = Eigen::VectorXi;
        using c_local_col_type = Eigen::VectorXi;
        using c_local2_row_type = Eigen::VectorXi;
        using c_local2_col_type = Eigen::VectorXi;
        using c_mortar_test_local_row_type = Eigen::VectorXi;
        /** Component orientation signs preserve their own static axis extents. */
        using ComponentRowSigns = typename ComponentStorage::RowSigns;
        /** Trial component orientation signs are independent of test cardinality. */
        using ComponentColSigns = typename ComponentStorage::ColSigns;
        /** Orientation signs for the test components on two elements. */
        using ComponentPairRowSigns = typename ComponentStorage::TwoCellRowSigns;
        /** Orientation signs for the trial components on two elements. */
        using ComponentPairColSigns = typename ComponentStorage::TwoCellColSigns;

    public:

        bool canReuseTestContext() const
        {
            if constexpr ( !std::is_same_v<map_test_fecontext_type, map_trial_fecontext_type> )
                return false;
            else if constexpr ( hasDynamicOrder )
            {
                // Matching C++ types do not imply matching runtime bases.
                // Reuse only the same FE instance on the same mesh.
                return M_form.testSpace()->mesh()->isSameMesh( M_form.trialSpace()->mesh() ) &&
                       M_form.template testFiniteElement<UseMortarTest>().get() ==
                           M_form.template trialFiniteElement<UseMortarTrial>().get();
            }
            else
                return M_form.testSpace()->mesh()->isRelatedTo( M_form.trialSpace()->mesh() );
        }

        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  left, Right right, bool same_mesh = true )
        {
            if ( same_mesh && this->canReuseTestContext() )
                return getMap( left, right, std::is_same<Left, map_trial_fecontext_type>() );
            return getMap( left, right, std::integral_constant<bool, false>() );
        }
        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  left, Right /*right*/, std::true_type )
        {
            return left;
        }
        template<typename Left, typename Right>
        map_trial_fecontext_type getMap( Left  /*left*/, Right right, std::false_type )
        {
            return map_trial_fecontext_type( right );
        }

        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  left, Right right, bool same_mesh = true )
        {
            if ( same_mesh && this->canReuseTestContext() )
                return getMapL( left, right, std::is_same<Left, map_left_trial_fecontext_type>() );
            return getMapL( left, right, std::integral_constant<bool, false>() );
        }
        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  left, Right /*right*/, std::true_type )
        {
            return left;
        }
        template<typename Left, typename Right>
        map_left_trial_fecontext_type getMapL( Left  /*left*/, Right right, std::false_type )
        {
            return right;
        }

        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im );

        template<typename IMTest,typename IMTrial>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im, IMTest const& imTest, IMTrial const& imTrial );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2 );

        template<typename IM2,typename IMTest,typename IMTrial>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2, IMTest const& imTest, IMTrial const& imTrial );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2,
                 mpl::int_<2> );

        template<typename IM2,typename IMTest,typename IMTrial>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& gmcTest,
                 map_trial_geometric_mapping_context_type const& _gmcTrial,
                 map_geometric_mapping_expr_context_type const & gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2, IMTest const& imTest, IMTrial const& imTrial,
                 mpl::int_<2> );


        void initDynamicEigenMatrix()
        {
            using ::Feel::detail::assemblyDofCount;
            using ::Feel::detail::assemblyExtent;
            using ::Feel::detail::resizeAssemblyMatrix;

            // Obtain effective cardinalities before narrowing or adjusting either side.
            int const rows = assemblyDofCount( M_test_dof->nLocalDof() );
            int const cols = assemblyDofCount( M_trial_dof->nLocalDof() );
            int const componentRows = [&]
            {
                if constexpr ( test_dof_type::is_order_dynamic )
                    return assemblyDofCount( M_form.template testFiniteElement<UseMortarTest>()->localDof() );
                else
                    return nDofPerComponentTest;
            }();
            int const componentCols = [&]
            {
                if constexpr ( trial_dof_type::is_order_dynamic )
                    return assemblyDofCount( M_form.template trialFiniteElement<UseMortarTrial>()->localDof() );
                else
                    return nDofPerComponentTrial;
            }();

            // Reuse an unchanged shape; fixed axes remain checked in Release builds.
            auto bind = [&]( auto& ordinary, auto& mortar, auto& pair,
                             auto& rowIndices, auto& colIndices, auto& rowSigns, auto& colSigns,
                             auto& pairRows, auto& pairCols, auto& pairRowSigns, auto& pairColSigns,
                             int r, int c )
            {
                resizeAssemblyMatrix( ordinary, r, c );
                resizeAssemblyMatrix( mortar, assemblyExtent( r, 1 ), c );
                resizeAssemblyMatrix( pair, assemblyExtent( r, 0, 2 ), assemblyExtent( c, 0, 2 ) );
                resizeAssemblyMatrix( rowIndices, r );
                resizeAssemblyMatrix( colIndices, c );
                resizeAssemblyMatrix( rowSigns, r );
                resizeAssemblyMatrix( colSigns, c );
                resizeAssemblyMatrix( pairRows, assemblyExtent( r, 0, 2 ) );
                resizeAssemblyMatrix( pairCols, assemblyExtent( c, 0, 2 ) );
                resizeAssemblyMatrix( pairRowSigns, assemblyExtent( r, 0, 2 ) );
                resizeAssemblyMatrix( pairColSigns, assemblyExtent( c, 0, 2 ) );
            };
            bind( M_rep, M_mortarTest_rep, M_rep_2, M_local_rows, M_local_cols,
                  M_local_rowsigns, M_local_colsigns, M_local_rows_2, M_local_cols_2,
                  M_local_rowsigns_2, M_local_colsigns_2, rows, cols );
            bind( M_c_rep, M_c_mortarTest_rep, M_c_rep_2, M_c_local_rows, M_c_local_cols,
                  M_c_local_rowsigns, M_c_local_colsigns, M_c_local_rows_2, M_c_local_cols_2,
                  M_c_local_rowsigns_2, M_c_local_colsigns_2, componentRows, componentCols );

            // The trial-mortar owner is distinct from the test-mortar component path.
            resizeAssemblyMatrix( M_mortarTrial_rep, rows, assemblyExtent( cols, 1 ) );
            resizeAssemblyMatrix( M_mortarTest_local_rows, assemblyExtent( rows, 1 ) );
            resizeAssemblyMatrix( M_mortarTrial_local_cols, assemblyExtent( cols, 1 ) );
            resizeAssemblyMatrix( M_c_mortarTest_local_rows, assemblyExtent( componentRows, 1 ) );
        }

        bool trialElementIsOnBoundary( index_type test_eid ) const { return M_form.trialSpace()->mesh()->element( this->trialElementId( test_eid ) ).isOnBoundary(); }

        index_type trialElementId( index_type trial_eid ) const
            {
                if constexpr ( nDimDiffBetweenTestTrial == 0 )
                {
                    index_type idElem = trial_eid;
                    index_type domain_eid = idElem;
                    const bool test_related_to_trial = M_form.testSpace()->mesh()->isSubMeshFrom( M_form.trialSpace()->mesh() );
                    const bool trial_related_to_test = M_form.trialSpace()->mesh()->isSubMeshFrom( M_form.testSpace()->mesh() );
                    const bool test_sibling_of_trial = M_form.testSpace()->mesh()->isSiblingOf( M_form.trialSpace()->mesh() );
                    if ( test_related_to_trial )
                    {
                        domain_eid = M_form.testSpace()->mesh()->subMeshToMesh( idElem );
                        DVLOG(2) << "[test_related_to_trial] test element id: "  << idElem << " trial element id : " << domain_eid << "\n";
                    }
                    if( trial_related_to_test )
                    {
                        domain_eid = M_form.trialSpace()->mesh()->meshToSubMesh( idElem );
                        DVLOG(2) << "[trial_related_to_test] test element id: "  << idElem << " trial element id : " << domain_eid << "\n";
                    }
                    if ( test_sibling_of_trial )
                    {
                        domain_eid = M_form.testSpace()->mesh()->meshToSubMesh( M_form.trialSpace()->mesh(), trial_eid );
                        DVLOG(2) << "[trial_sibling_of_test] test element id: "  << idElem << " trial element id : " << domain_eid << "\n";
                    }
                    return domain_eid;
                }
                else // nDimDiffBetweenTestTrial == 1
                {
                    index_type idElem = trial_eid;
                    index_type domain_eid = idElem;
                    const bool test_related_to_trial = M_form.testSpace()->mesh()->isSubMeshFrom( M_form.trialSpace()->mesh() );
                    const bool trial_related_to_test = M_form.trialSpace()->mesh()->isSubMeshFrom( M_form.testSpace()->mesh() );
                    if ( test_related_to_trial )
                    {
                        domain_eid = M_form.trialSpace()->mesh()->face( M_form.testSpace()->mesh()->subMeshToMesh( idElem )).element0().id();
                        DVLOG(2) << "[test_related_to_trial] test element id: "  << idElem << " trial element id : " << domain_eid << "\n";
                    }
                    if( trial_related_to_test )
                    {
                        auto const& eltTest = M_form.testSpace()->mesh()->element(idElem);
                        std::set<index_type> idsFind;
                        for (uint16_type f=0;f< M_form.testSpace()->mesh()->numLocalFaces();++f)
                            {
                                const index_type idFind = M_form.trialSpace()->mesh()->meshToSubMesh( eltTest.face(f).id() );
                                if ( idFind != invalid_v<index_type> ) idsFind.insert( idFind );
                            }
                        if ( idsFind.size()>1 ) std::cout << " TODO trialElementId " << std::endl;

                        if ( idsFind.size()>0 )
                            domain_eid = *idsFind.begin();
                        else
                            domain_eid = invalid_v<index_type>;

                        DVLOG(2) << "[trial_related_to_test] test element id: "  << idElem << " trial element id : " << domain_eid << "\n";
                    }
                    return domain_eid;
                }
            }
        index_type trialElementId( typename mesh_1_type::element_iterator it ) const
            {
                return trialElementId( it->id() );
            }
        index_type trialElementId( typename mesh_1_type::element_type const& e ) const
            {
                return trialElementId( e.id() );
            }
        bool isZero( index_type i ) const
            {
                index_type domain_eid = trialElementId( i );
                if ( domain_eid == invalid_v<index_type> )
                    return true;
                return false;
            }
        bool isZero( typename mesh_1_type::element_iterator it ) const
            {
                return this->isZero( it->id() );
            }
        bool isZero( typename mesh_1_type::element_type const& e ) const
            {
                return this->isZero( e.id() );
            }

        void update( map_test_geometric_mapping_context_type const& gmcTest,
                     map_trial_geometric_mapping_context_type const& _gmcTrial,
                     map_geometric_mapping_expr_context_type const& gmcExpr );

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const& gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad );

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
                                        std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad )
        {
            M_integrator = im;
            updateInCaseOfInterpolate( gmcTest, gmcTrial, gmcExpr, indexLocalToQuad );
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
        void integrateInCaseOfInterpolate( std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad,
                                           bool isFirstExperience )
        {
            integrateInCaseOfInterpolate( mpl::int_<fusion::result_of::size<GeomapTestContext>::type::value>(),
                                          indexLocalToQuad,
                                          isFirstExperience );
        }


        void assemble()
        {
            assemble( fusion::at_key<gmc<0> >( M_test_gmc )->id() );
        }
        void assemble( index_type elt_0 );

        void assemble( mpl::int_<2> )
        {
            assemble( fusion::at_key<gmc<0> >( M_test_gmc )->id(),
                      fusion::at_key<test_gmc1>( M_test_gmc )->id() );
        }
        void assemble( index_type elt_0, index_type elt_1  );

        /**
         * local to global assembly given a pair of test and trial element ids
         *  - \p elt.first provides the test element
         *  - \p elt.second provides the trial element
         */
        void assemble( std::pair<index_type,index_type> const& elt );

        /**
         * local to global assembly associated to internal faces given a pair of
         * test and trial element ids
         *
         *  - \p elt.first provides the test element
         *  - \p elt.second provides the trial element
         */
        void assemble( std::pair<index_type,index_type> const& elt_0,
                       std::pair<index_type,index_type> const& elt_1 );

        void assembleInCaseOfInterpolate();


        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points
         */
        template<typename Pts>
        void precomputeBasisAtPoints( Pts const& pts )
        {
            M_test_pc = test_precompute_ptrtype( new test_precompute_type( M_form.testFiniteElement<UseMortarTest>(), pts ) );
            M_trial_pc = trial_precompute_ptrtype( new trial_precompute_type( M_form.trialFiniteElement<UseMortarTrial>(), pts ) );
        }

        template<typename PtsTest,typename PtsTrial>
        void precomputeBasisAtPoints( PtsTest const& pts1,PtsTrial const& pts2  )
        {
            M_test_pc = test_precompute_ptrtype( new test_precompute_type( M_form.testFiniteElement<UseMortarTest>(), pts1 ) );
            M_trial_pc = trial_precompute_ptrtype( new trial_precompute_type( M_form.trialFiniteElement<UseMortarTrial>(), pts2 ) );
        }


        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points on the face of
         * the reference element
         */
        template<typename Pts>
        void precomputeBasisAtPoints( uint16_type __f, permutation_1_type const& __p, Pts const& pts )
        {
            M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testFiniteElement<UseMortarTest>(), pts ) );
            //FEELPP_ASSERT( M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );

            M_trial_pc_face[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( M_form.trialFiniteElement<UseMortarTrial>(), pts ) );
            //FEELPP_ASSERT( M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts )
        {
            std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> > testpc;
            
            if constexpr ( std::is_same_v<permutation_1_type, typename QuadMapped<PtsSet>::permutation_type> )
            {
                //QuadMapped<PtsSet> qm;
                typedef typename QuadMapped<PtsSet>::permutation_type permutation_type;
                //typename QuadMapped<PtsSet>::permutation_points_type ppts( qm( pts ) );

                for ( uint16_type __f = 0; __f < pts.nFaces(); ++__f )
                {
                    for ( permutation_type __p( permutation_type::IDENTITY );
                            __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                    {
                        //testpc[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testFiniteElement<UseMortar>(), ppts[__f].find( __p )->second ) );
                        testpc[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testFiniteElement<UseMortarTest>(), pts.fpoints( __f,__p.value() ) ) );
                    }
                }
            }
            
            return testpc;
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> >
        precomputeTrialBasisAtPoints( PtsSet const& pts )
        {
            std::map<uint16_type, std::map<permutation_2_type,trial_precompute_ptrtype> > trialpc;
            
            if constexpr ( std::is_same_v<permutation_2_type, typename QuadMapped<PtsSet>::permutation_type> )
            {
                //QuadMapped<PtsSet> qm;
                typedef typename QuadMapped<PtsSet>::permutation_type permutation_type;
                //typename QuadMapped<PtsSet>::permutation_points_type ppts( qm( pts ) );

                for ( uint16_type __f = 0; __f < pts.nFaces(); ++__f )
                {
                    for ( permutation_type __p( permutation_type::IDENTITY );
                            __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                    {
                        //trialpc[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( M_form.trialFiniteElement<UseMortar>(), ppts[__f].find( __p )->second ) );
                        trialpc[__f][__p] = trial_precompute_ptrtype( new trial_precompute_type( M_form.trialFiniteElement<UseMortarTrial>(), pts.fpoints(__f, __p.value() ) ) );
                    }
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
                return  M_test_pc;

            //FEELPP_ASSERT( M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );
            return M_test_pc_face.find( __f )->second.find( __p )->second;
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
                return  M_trial_pc;

            //FEELPP_ASSERT( M_trial_pc_face.find(__f )->second )( __f ).error( "invalid trial precompute type" );
            return M_trial_pc_face.find( __f )->second.find( __p )->second;
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
                                           std::vector<boost::tuple<index_type,index_type> > const& indexLocalToQuad,
                                           bool isFirstExperience );
    private:

        form_type& M_form;
        const list_block_type& M_lb;
        dof_1_type* M_test_dof;
        dof_2_type* M_trial_dof;

        test_precompute_ptrtype M_test_pc;
        std::map<uint16_type, std::map<permutation_1_type,test_precompute_ptrtype> > M_test_pc_face;
        trial_precompute_ptrtype M_trial_pc;
        std::map<uint16_type, std::map<permutation_2_type, trial_precompute_ptrtype> > M_trial_pc_face;


        map_test_geometric_mapping_context_type M_test_gmc;
        map_trial_geometric_mapping_context_type M_trial_gmc;

        map_test_fecontext_type M_test_fec;
        map_left_test_fecontext_type M_test_fec0;
        map_right_test_fecontext_type M_test_fec1;
        map_trial_fecontext_type M_trial_fec;
        map_left_trial_fecontext_type M_trial_fec0;
        map_right_trial_fecontext_type M_trial_fec1;


        local_matrix_type M_rep;
        mortar_test_local_matrix_type M_mortarTest_rep;
        mortar_trial_local_matrix_type M_mortarTrial_rep;
        local2_matrix_type M_rep_2;
        local_row_type M_local_rows;
        mortar_test_local_row_type M_mortarTest_local_rows;
        mortar_trial_local_col_type M_mortarTrial_local_cols;
        local2_row_type M_local_rows_2;
        local_col_type M_local_cols;
        local2_col_type M_local_cols_2;

        local_row_sign_type/*local_row_type*/ M_local_rowsigns;
        local2_row_sign_type/*local2_row_type*/ M_local_rowsigns_2;
        local_col_sign_type/*local_col_type*/ M_local_colsigns;
        local2_col_sign_type/*local2_col_type*/ M_local_colsigns_2;

        c_local_matrix_type M_c_rep;
        c_mortar_test_local_matrix_type M_c_mortarTest_rep;
        c_local2_matrix_type M_c_rep_2;
        c_local_row_type M_c_local_rows;
        c_mortar_test_local_row_type M_c_mortarTest_local_rows;
        c_local2_row_type M_c_local_rows_2;
        c_local_col_type M_c_local_cols;
        c_local2_col_type M_c_local_cols_2;
        ComponentRowSigns M_c_local_rowsigns;
        ComponentPairRowSigns M_c_local_rowsigns_2;
        ComponentColSigns M_c_local_colsigns;
        ComponentPairColSigns M_c_local_colsigns_2;

        eval00_expr_ptrtype M_eval_expr00;
        eval01_expr_ptrtype M_eval_expr01;
        eval10_expr_ptrtype M_eval_expr10;
        eval11_expr_ptrtype M_eval_expr11;

        IM M_integrator;

    }; // Context

    //@}

    /** @name Constructors, destructor
     */
    //@{


    BilinearForm() = default;
    /**


     */
    BilinearForm( std::string name,
                  space_1_ptrtype const& __X1,
                  space_2_ptrtype const& __X2,
                  matrix_ptrtype& __M,
                  size_type rowstart = 0,
                  size_type colstart = 0,
                  bool build = true,
                  bool do_threshold = false,
                  value_type threshold = type_traits<value_type>::epsilon(),
                  size_type graph_hints = Pattern::COUPLED );

    BilinearForm( std::string name,
                  space_1_ptrtype const& __X1,
                  space_2_ptrtype const& __X2,
                  matrix_ptrtype& __M,
                  list_block_type const& __lb,
                  size_type rowstart = 0,
                  size_type colstart = 0,
                  bool do_threshold = false,
                  value_type threshold = type_traits<value_type>::epsilon(),
                  size_type graph_hints = Pattern::COUPLED );

    /**
     * @brief Construct a new Bilinear Form object
     *
     * @param f bilinear form to copy
     */
    BilinearForm( BilinearForm const& f )
    {
        bool same_spaces = (M_X1 == f.M_X1) && ( M_X2 == f.M_X2 );
        //tic(M_name, Environment::logVerbosityLevel() > 0 );
        M_X1 = f.M_X1;
        M_X2 = f.M_X2;
        if ( !this->isMatrixAllocated() || !same_spaces )
            this->allocateMatrix( M_X1, M_X2 );
        super::operator=( f );
        //toc(M_name, Environment::logVerbosityLevel() > 0 );
    }
    /**
     * @brief Construct a new Bilinear Form object
     *
     * @param __vf
     */
    BilinearForm( BilinearForm && __vf ) = default;

    ~BilinearForm() override
    {
        //toc(M_name, Environment::logVerbosityLevel() > 0 );
    }



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
            bool same_spaces = (M_X1 == form.M_X1) && ( M_X2 == form.M_X2 );
            M_X1 = form.M_X1;
            M_X2 = form.M_X2;
            if ( !this->isMatrixAllocated() || !same_spaces )
                this->allocateMatrix( M_X1, M_X2 );
            super::operator=( form );

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
        return this->M_matrix->energy( __v, __u );
    }

    /**
     * @brief Computes the linear form associated with the bilinear form for a given trial function
     *
     * @param __u
     * @return form1_t<test_space_type>
     */
    form1_t<test_space_type> operator()( element_2_type const& __u ) const
        {
            auto l = form1_t<test_space_type>( "linearform.l", M_X1, backend()->newVector( M_X1 ) );
            this->M_matrix->multVector( __u, l.vector() );
            return l;
        }

    BilinearForm& operator+=( BilinearForm const& a )
    {
        static_cast<super&>(*this) += static_cast<const super&>(a);
        return *this;
    }
    BilinearForm& operator-=( BilinearForm const& a )
    {
         static_cast<super&>(*this) -= static_cast<const super&>(a);
        return *this;
    }
    BilinearForm& operator*=( value_type const& a )
    {
        static_cast<super&>(*this) *= a;
        return *this;
    }
    BilinearForm& operator/=( value_type const& a )
    {
        static_cast<super&>(*this) /= a;
        return *this;
    }
    /**
     * @brief unary minus operator
     *
     * @return new BilinearForm
     */
    BilinearForm operator-() const
    {
        return BilinearForm( *this ) *= -1;
    }

    /**
     * \return the entry \f$M_{i,j}\f$
     */
    // value_type& operator()( size_type i, size_type j ) { return M_matrix( i, j ); }
    //@}

    /** @name Accessors
     */
    //@{


    /**
     * \return the trial function space
     */
    space_2_ptrtype const& trialSpace() const
    {
        return M_X2;
    }

    /**
     * \return the test function space
     */
    space_1_ptrtype const& testSpace() const
    {
        return M_X1;
    }

    /**
     * Geometric transformation
     */
    gm_1_ptrtype const& gm() const
    {
        return M_X1->gm();
    }

    /**
     * Geometric transformation
     */
    gm1_1_ptrtype const& gm1() const
    {
        return M_X1->gm1();
    }

    /**
     * Geometric transformation
     */
    gm_2_ptrtype const& gmTrial() const
    {
        return M_X2->gm();
    }

    /**
     * Geometric transformation
     */
    gm1_2_ptrtype const& gm1Trial() const
    {
        return M_X2->gm1();
    }



    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    //@}

private:


    template <class ExprT> void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<false> );
    template <class ExprT> void assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<true> );
    template <class ExprT> void assign( Expr<ExprT> const& __expr, mpl::bool_<true>, mpl::bool_<true> );
    template <class ExprT> void assign( Expr<ExprT> const& __expr, mpl::bool_<true>, mpl::bool_<false> );
    template <class ExprT> void assign( Expr<ExprT> const& __expr, mpl::bool_<true>, mpl::bool_<false>, mpl::bool_<true> );
    template <class ExprT> void assign( Expr<ExprT> const& __expr, mpl::bool_<true>, mpl::bool_<false>, mpl::bool_<false> );

private:

    space_1_ptrtype M_X1;
    space_2_ptrtype M_X2;

};
template<typename FE1,  typename FE2, typename ElemContType>
BilinearForm<FE1, FE2, ElemContType>::BilinearForm( std::string name,
                                                    space_1_ptrtype const& Xh,
                                                    space_2_ptrtype const& Yh,
                                                    matrix_ptrtype& __M,
                                                    size_type rowstart,
                                                    size_type colstart,
                                                    bool build,
                                                    bool do_threshold,
                                                    value_type threshold,
                                                    size_type graph_hints )
:
    super( name, Xh, Yh, __M, rowstart, colstart, build, do_threshold, threshold, graph_hints ),
    M_X1( Xh ),
    M_X2( Yh )
{}

template<typename FE1,  typename FE2, typename ElemContType>
BilinearForm<FE1, FE2, ElemContType>::BilinearForm( std::string name,
                                                    space_1_ptrtype const& Xh,
                                                    space_2_ptrtype const& Yh,
                                                    matrix_ptrtype& __M,
                                                    list_block_type const& __lb,
                                                    size_type rowstart,
                                                    size_type colstart,
                                                    bool do_threshold ,
                                                    value_type threshold,
                                                    size_type graph_hints )
:
    super( name, Xh, Yh, __M, __lb, rowstart, colstart, do_threshold, threshold, graph_hints ),
    M_X1( Xh ),
    M_X2( Yh )
{}


template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
                                              bool init,
                                              mpl::bool_<false> )
{
    // do nothing if process not active for this space
    if ( !this->M_X1->worldComm().isActive() )
        return;

    if ( init )
    {
        this->M_matrix->zero();
        this->resetDeferredDirichletAfterOperatorReset();
    }

    __expr.assemble( M_X1, M_X2, *this );

    DVLOG(2) << "[BilinearForm::assign<false>] stop\n";

}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
                                              bool init,
                                              mpl::bool_<true> )
{
    // do nothing if process not active for this space
    if ( !this->M_X1->worldComm().isActive() )
        return;

    DVLOG(2) << "BilinearForm::assign() start loop on test spaces\n";

    if ( init )
    {
        this->M_matrix->zero();
        this->resetDeferredDirichletAfterOperatorReset();
    }

    assign( __expr, mpl::bool_<true>(), mpl::bool_<( space_1_type::nSpaces > 1 && space_2_type::nSpaces > 1 )>() );
    DVLOG(2) << "BilinearForm::assign() stop loop on test spaces\n";

}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
        mpl::bool_<true>,
        mpl::bool_<true> )
{
    fusion::for_each( M_X1->functionSpaces(), make_bfassign2( *this, __expr ) );
}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
        mpl::bool_<true>,
        mpl::bool_<false> )
{
    assign( __expr, mpl::bool_<true>(), mpl::bool_<false>(), mpl::bool_<( space_1_type::nSpaces > 1 )>() );
}

template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
        mpl::bool_<true>,
        mpl::bool_<false>,
        mpl::bool_<true> )
{
    fusion::for_each( M_X1->functionSpaces(), make_bfassign3( *this, __expr, M_X2, 0 ) );
}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
void
BilinearForm<FE1, FE2, ElemContType>::assign( Expr<ExprT> const& __expr,
        mpl::bool_<true>,
        mpl::bool_<false>,
        mpl::bool_<false> )
{
    fusion::for_each( M_X2->functionSpaces(), make_bfassign1( *this, __expr, M_X1, 0 ) );

}
template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator=( Expr<ExprT> const& __expr )
{
    // loop(fusion::for_each) over sub-functionspaces in SpaceType
    // pass expression and initialize
#if 1
    this->assign( __expr, true, mpl::bool_<mpl::or_< mpl::bool_< ( space_1_type::nSpaces > 1 )>,
                  mpl::bool_< ( space_2_type::nSpaces > 1 )> >::type::value >() );
#else
    using has_multiple_spaces_t = mpl::bool_<mpl::or_< mpl::bool_< ( space_1_type::nSpaces > 1 )>,
                                             mpl::bool_< ( space_2_type::nSpaces > 1 )> >::type::value >;

    //M_futs_assign.push_back( std::async( std::launch::async, &BilinearForm<FE1, FE2, ElemContType>::template assign<ExprT>, this, __expr, true, has_multiple_spaces_t() ) );
    this->push_back( std::async( std::launch::async, [this,__expr](){ this->assign( __expr, true, has_multiple_spaces_t() ); } ) );
    std::cout << "futs: " << this->M_fut_assign.size() << std::endl;
#endif
    return *this;
}

template<typename FE1,  typename FE2,  typename ElemContType>
template<typename ExprT>
BilinearForm<FE1, FE2, ElemContType>&
BilinearForm<FE1, FE2, ElemContType>::operator+=( Expr<ExprT> const& __expr )
{
    DVLOG(2) << "[BilinearForm::operator+=] start\n";
    this->assign( __expr, false, mpl::bool_<mpl::or_< mpl::bool_< ( space_1_type::nSpaces > 1 )>,
                  mpl::bool_< ( space_2_type::nSpaces > 1 )> >::type::value >() );
    DVLOG(2) << "[BilinearForm::operator+=] stop\n";
    return *this;
}


template<typename BFType, typename ExprType, typename TestSpaceType>
template<typename SpaceType>
void BFAssign1<BFType,ExprType,TestSpaceType>::operator()( std::shared_ptr<SpaceType> const& trial ) const
{
    if ( M_bf.testSpace()->worldsComm()[M_test_index]->isActive() )
    {
        DVLOG(2) << "[BFAssign1::operator()] expression has test functions index "
                      << M_test_index << " : "
                 << ExprType::template HasTestFunction<typename TestSpaceType::reference_element_type>::result << " (0:no, 1:yes)\n"
                 << ExprType::template HasTestFunction<typename TestSpaceType::component_basis_type>::result << " (0:no, 1:yes)\n";
        DVLOG(2) << "[BFAssign1::operator()] expression has trial functions index "
                 << M_trial_index << " :"
                 << ExprType::template HasTrialFunction<typename SpaceType::reference_element_type>::result << " (0:no, 1:yes)\n"
                 << ExprType::template HasTrialFunction<typename SpaceType::component_basis_type>::result << " (0:no, 1:yes)\n";

        if ( !ExprType::template HasTestFunction<typename TestSpaceType::reference_element_type>::result ||
             //!ExprType::template HasTestFunction<typename TestSpaceType::component_basis_type>::result ||
             !ExprType::template HasTrialFunction<typename SpaceType::reference_element_type>::result //||
             //!ExprType::template HasTrialFunction<typename SpaceType::component_basis_type>::result
             )
        {
            ++M_trial_index;
            return;
        }

        DVLOG(2) << "[BFAssign1::operator()] terms found with "
                      << "testindex: " << M_test_index << " trialindex: " << M_trial_index << "\n";
        typedef SpaceType trial_space_type;
        typedef TestSpaceType test_space_type;
#if 0
        Feel::vf::list_block_type list_block;

        if ( M_bf.testSpace()->worldsComm()[M_test_index].globalSize()>1 )
        {
            //DVLOG(2) << "[BFAssign1::operator()] block: " << block << "\n";
            if (M_bf.testSpace()->hasEntriesForAllSpaces())
                list_block.push_back( Feel::vf::Block( 0, 0,
                                                       0/*M_bf.testSpace()->nLocalDofStart( M_test_index )*/,
                                                       0/*M_bf.trialSpace()->nLocalDofStart( M_trial_index )*/ ) );
            else
                list_block.push_back( Feel::vf::Block( 0, 0, 0, M_bf.trialSpace()->nLocalDofStart( M_trial_index ) ) );

        }

        else
        {
            //DVLOG(2) << "[BFAssign1::operator()] block: " << block << "\n";
            list_block.push_back( Feel::vf::Block( 0,0,
                                                   0/*M_bf.testSpace()->nDofStart( M_test_index )*/,
                                                   0/*M_bf.trialSpace()->nDofStart( M_trial_index )*/ ) );
        }
#endif
        typedef typename BFType::matrix_type matrix_type;
        typedef Feel::vf::detail::BilinearForm<test_space_type,
                trial_space_type,
                ublas::vector_range<ublas::vector<double> > > bf_type;
#if 0
        bf_type bf( M_test,trial, M_bf.matrixPtr(), list_block,  M_bf.rowStartInMatrix(), M_bf.colStartInMatrix(), M_bf.doThreshold(), M_bf.threshold(), M_bf.pattern()  );
        datamap_ptrtype dmTest = M_bf.matrixPtr()->mapRowPtr();//M_bf.testSpace()->dof()
        datamap_ptrtype dmTrial = M_bf.matrixPtr()->mapColPtr();//M_bf.trialSpace()->dof()
        bf.setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_bf.rowStartInMatrix() + M_test_index ) );
        bf.setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_bf.colStartInMatrix() + M_trial_index ) );
#else
        bf_type bf( M_bf.name()+"["+std::to_string(M_test_index)+","+std::to_string(M_trial_index)+"]", M_test, trial, M_bf.matrixPtr(),
                    M_bf.rowStartInMatrix() + M_test_index, M_bf.colStartInMatrix() + M_trial_index,
                    false, M_bf.doThreshold(), M_bf.threshold(), M_bf.pattern() );
#endif
        bf.shareDeferredDirichletState( M_bf.deferredDirichletStatePtr() );
        bf.setDirichletPolicy( M_bf.dirichletPolicy() );

        bf += M_expr;
    }

    ++M_trial_index;
}

template<typename BFType, typename ExprType, typename TrialSpaceType>
template<typename SpaceType>
void BFAssign3<BFType,ExprType,TrialSpaceType>::operator()( std::shared_ptr<SpaceType> const& test ) const
{
    if ( M_bf.testSpace()->worldsComm()[M_test_index]->isActive() )
    {

        DVLOG(2) << "[BFAssign3::operator()] expression has trial functions index "
                 << M_test_index << " : "
                 << ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result << " (0:no, 1:yes)\n"
                 << ExprType::template HasTestFunction<typename SpaceType::component_basis_type>::result << " (0:no, 1:yes)\n";
        DVLOG(2) << "[BFAssign3::operator()] expression has test functions index "
                 << M_trial_index << " :"
                 << ExprType::template HasTrialFunction<typename TrialSpaceType::reference_element_type>::result << " (0:no, 1:yes)\n"
                 << ExprType::template HasTrialFunction<typename TrialSpaceType::component_basis_type>::result << " (0:no, 1:yes)\n";


        if ( !ExprType::template HasTrialFunction<typename TrialSpaceType::reference_element_type>::result ||
             //!ExprType::template HasTrialFunction<typename TrialSpaceType::component_basis_type>::result ||
             !ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result //||
             //!ExprType::template HasTestFunction<typename SpaceType::component_basis_type>::result
             )
        {
            ++M_test_index;
            return;
        }

        DVLOG(2) << "[BFAssign3::operator()] terms found with "
                      << "testindex: " << M_test_index << " trialindex: " << M_trial_index << "\n";
        typedef SpaceType test_space_type;
        typedef TrialSpaceType trial_space_type;

        //Feel::vf::Block block ( 0, 0,
        //                        M_bf.testSpace()->nDofStart( M_test_index ),
        //                        M_bf.trialSpace()->nDofStart( M_trial_index ) );
        Feel::vf::list_block_type list_block;

        // with mpi, dof start to 0 (thanks to the LocalToGlobal mapping).
        if ( M_bf.testSpace()->worldsComm()[M_test_index]->globalSize()>1 )
        {
            //DVLOG(2) << "[BFAssign1::operator()] block: " << block << "\n";
            if (M_bf.testSpace()->hasEntriesForAllSpaces())
                list_block.push_back( Feel::vf::Block( 0, 0,
                                                       0/*M_bf.testSpace()->nLocalDofStart( M_test_index )*/,
                                                       0/*M_bf.trialSpace()->nLocalDofStart( M_trial_index )*/ ) );
            else
                list_block.push_back( Feel::vf::Block( 0, 0, 0, M_bf.trialSpace()->nLocalDofStart( M_trial_index ) ) );

        }

        else
        {
            //DVLOG(2) << "[BFAssign1::operator()] block: " << block << "\n";
            list_block.push_back( Feel::vf::Block( 0,0,
                                                   0/*M_bf.testSpace()->nDofStart( M_test_index )*/,
                                                   0/*M_bf.trialSpace()->nDofStart( M_trial_index )*/ ) );
        }

        typedef typename BFType::matrix_type matrix_type;
        typedef Feel::vf::detail::BilinearForm<test_space_type,
                trial_space_type,
                ublas::vector_range<ublas::vector<double> > > bf_type;
#if 0
        bf_type bf( test, M_trial, M_bf.matrixPtr(), list_block, M_bf.rowStartInMatrix(), M_bf.colStartInMatrix(), M_bf.doThreshold(), M_bf.threshold(), M_bf.pattern() );
        datamap_ptrtype dmTest = M_bf.matrixPtr()->mapRowPtr(); //M_bf.testSpace()->dof()
        datamap_ptrtype dmTrial = M_bf.matrixPtr()->mapColPtr(); //M_bf.trialSpace()->dof()
        bf.setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_bf.rowStartInMatrix() + M_test_index ) );
        bf.setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_bf.colStartInMatrix() + M_trial_index ) );
#else
        bf_type bf( M_bf.name()+"["+std::to_string(M_test_index)+","+std::to_string(M_trial_index)+"]",
                    test, M_trial, M_bf.matrixPtr(),
                    M_bf.rowStartInMatrix() + M_test_index, M_bf.colStartInMatrix()+ M_trial_index,
                    false, M_bf.doThreshold(), M_bf.threshold(), M_bf.pattern() );
#endif
        bf.shareDeferredDirichletState( M_bf.deferredDirichletStatePtr() );
        bf.setDirichletPolicy( M_bf.dirichletPolicy() );

        bf += M_expr;
    }

    ++M_test_index;
}


} // detail
} // vf

namespace meta
{
template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
struct BilinearForm
{
    typedef Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> type;
};


}

/// \endcond

template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> operator+( Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> const& __a,
                                                                Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> const& __b )
{
    BilinearForm<FE1,FE2,ElemContType> __c( __a );
    __c += __b;
    return __c;
}
template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> operator-( Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> const& __a,
                                                                Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> const& __b )
{
    BilinearForm<FE1,FE2,ElemContType> __c( __a );
    __c -= __b;
    return __c;
}

/**
 * @brief sum of bilinear forms
 *
 * @tparam FE1
 * @tparam FE2
 * @tparam BinaryOperation
 * @tparam VectorUblas<typename functionspace_type<FE1>::value_type>
 * @param v vector of bilinear forms
 * @param op is the binary operation to apply
 * @param init is the initial value
 * @return BilinearForm<FE1,FE2,ElemContType>
 */
template<typename FE1,
         typename FE2,
         class BinaryOperation,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> sum( std::vector<Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> > const& v,
                                                          BinaryOperation op,
                                                          Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType> init = Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType>() )
{
    for(auto const& a : v)
        init = op( std::move(init), a );
    return init;
}

template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
using form2_type = Feel::vf::detail::BilinearForm<FE1,FE2,ElemContType>;

template<typename FE1,
         typename FE2,
         typename ElemContType = VectorUblas<typename functionspace_type<FE1>::value_type> >
using form2_t = form2_type<FE1,FE2,ElemContType>;


template<typename T>
struct is_bilinear_form : std::is_base_of<BilinearFormBase<typename std::decay_t<T>::value_type>, std::decay_t<T>>::type {};
template <typename T>
inline constexpr bool is_bilinear_form_v = is_bilinear_form<T>::value;


template <typename T>
struct form_test_mesh : mp11::mp_identity<typename std::decay_t<T>::mesh_type> {};
template <typename T>
    requires is_bilinear_form_v<T>
struct form_test_mesh< T > : mp11::mp_identity<typename std::decay_t<T>::mesh_1_type> {};
template <typename T>
using form_test_mesh_t = typename form_test_mesh<T>::type;

} // feel

#include <feel/feelvf/bilinearformcontext.hpp>

#endif /* __BilinearForm_H */
