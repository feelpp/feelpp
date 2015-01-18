/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2011 Université Joseph Fourier (Grenoble I)

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
   \file linearform.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-18
 */
#ifndef __LinearForm_H
#define __LinearForm_H 1

#include <Eigen/Eigen>
#include <Eigen/StdVector>

#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/multi_array.hpp>
#include <feel/feelvf/block.hpp>
#include <feel/feelalg/vectorvalue.hpp>
#include <feel/feelvf/fec.hpp>
#include <feel/feelvf/formcontextbase.hpp>

namespace Feel
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
    typedef space_type test_space_type;
    typedef space_type trial_space_type;

    typedef typename space_type::value_type value_type;
    typedef typename space_type::real_type real_type;
    typedef VectorType vector_type;
    typedef boost::shared_ptr<VectorType> vector_ptrtype;
    //typedef typename space_type::template Element<value_type, ElemContType> element_type;
    typedef typename space_type::template Element<value_type> element_type;

#if 0
    typedef typename space_type::component_fespace_type component_fespace_type;
    typedef typename space_type::element_type::component_type component_type;
#endif
    typedef typename space_type::mesh_type mesh_type;
    typedef typename mesh_type::element_type mesh_test_element_type;
    typedef typename mesh_type::element_type::permutation_type permutation_type;
    typedef typename space_type::gm_type gm_type;
    typedef typename space_type::gm_ptrtype gm_ptrtype;
    typedef typename space_type::gm1_type gm1_type;
    typedef typename space_type::gm1_ptrtype gm1_ptrtype;

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
             typename IM,
             typename GeomapExprContext = GeomapContext
             >
    class Context //: public FormContextBase<GeomapContext, IM>
    {
        typedef FormContextBase<GeomapContext, IM, GeomapExprContext> super;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef Context<GeomapContext,ExprT,IM,GeomapExprContext> form_context_type;
        typedef LinearForm<SpaceType,VectorType, ElemContType> form_type;
        typedef typename SpaceType::dof_type dof_type;
        typedef typename form_type::value_type value_type;


#if 0
        typedef GeomapContext map_geometric_mapping_context_type;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type geometric_mapping_context_ptrtype;
        typedef typename geometric_mapping_context_ptrtype::element_type geometric_mapping_context_type;
        typedef typename geometric_mapping_context_type::gm_type geometric_mapping_type;
#else
        typedef typename super::map_geometric_mapping_context_type map_test_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_context_ptrtype test_geometric_mapping_context_ptrtype;
        typedef typename super::geometric_mapping_context_type test_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_type test_geometric_mapping_type;

        typedef typename super::map_geometric_mapping_context_type map_trial_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_context_ptrtype trial_geometric_mapping_context_ptrtype;
        typedef typename super::geometric_mapping_context_type trial_geometric_mapping_context_type;
        typedef typename super::geometric_mapping_type trial_geometric_mapping_type;

#endif

        static const uint16_type nDim = test_geometric_mapping_type::nDim;

        typedef ExprT expression_type;

        typedef typename space_type::mesh_type mesh_type;
        typedef typename mesh_type::element_type mesh_element_type;
        typedef typename mesh_element_type::permutation_type permutation_type;
        typedef typename space_type::fe_type test_fe_type;
        typedef typename space_type::fe_type trial_fe_type;
        typedef boost::shared_ptr<test_fe_type> test_fe_ptrtype;
        typedef typename test_fe_type::template Context< test_geometric_mapping_context_type::context,
                test_fe_type,
                test_geometric_mapping_type,
                mesh_element_type> test_fecontext_type;
        typedef test_fecontext_type trial_fecontext_type;

#if 0
        typedef mpl::int_<fusion::result_of::template size<GeomapContext>::type::value> map_size;
#else
        typedef typename super::map_size map_size;
#endif

        typedef boost::shared_ptr<test_fecontext_type> test_fecontext_ptrtype;

#if 0
        typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >,
                gmc<1>,
                gmc<0> >::type gmc1;
#else
        typedef typename super::gmc1 gmc1;
#endif

#if 0
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type left_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type::element_type left_gmc_type;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type right_gmc_ptrtype;
        typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type::element_type right_gmc_type;

        typedef fusion::map<fusion::pair<gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
        typedef fusion::map<fusion::pair<gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;
#else
        typedef typename super::left_gmc_ptrtype left_gmc_ptrtype;
        typedef typename super::left_gmc_type left_gmc_type;
        typedef typename super::right_gmc_ptrtype right_gmc_ptrtype;
        typedef typename super::right_gmc_type right_gmc_type;

        typedef typename super::map_left_gmc_type map_left_gmc_type;
        typedef typename super::map_right_gmc_type map_right_gmc_type;
#endif

        typedef typename mpl::if_<mpl::equal_to<map_size, mpl::int_<1> >,
                mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > >,
                mpl::identity<fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype>,
                fusion::pair<gmc<1>, test_fecontext_ptrtype> > > >::type::type map_test_fecontext_type;

        typedef fusion::map<fusion::pair<gmc<0>, test_fecontext_ptrtype> > map_left_test_fecontext_type;
        typedef fusion::map<fusion::pair<gmc<1>, test_fecontext_ptrtype> > map_right_test_fecontext_type;

        typedef typename super::map_geometric_mapping_expr_context_type map_geometric_mapping_expr_context_type;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type, map_left_test_fecontext_type> eval0_expr_type;
        typedef boost::shared_ptr<eval0_expr_type> eval0_expr_ptrtype;

        typedef typename ExprT::template tensor<map_geometric_mapping_expr_context_type, map_right_test_fecontext_type> eval1_expr_type;
        typedef boost::shared_ptr<eval1_expr_type> eval1_expr_ptrtype;
        //typedef typename ExprT::template tensor<map_right_gmc_type, map_test_fecontext_type> eval1_expr_type;


        typedef typename test_fe_type::template Context< test_geometric_mapping_context_type::context,
                test_fe_type,
                test_geometric_mapping_type,
                mesh_element_type>::template Index<> test_index_type;


        //typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map0_test_fecontext_type> eval0_expr_type;
        //typedef typename ExprT::template tensor<map_geometric_mapping_context_type, map1_test_fecontext_type> eval1_expr_type;

        //typedef ublas::vector<value_type> local_vector_type;
        static const int rep_shape = 2;//1+(eval_expr_type::shape::M-1>0)+(eval_expr_type::shape::N-1>0);
        //typedef boost::multi_array<value_type, rep_shape> local_vector_type;
        typedef typename space_type::dof_type test_dof_type;
        static const int nDofPerElementTest = space_type::dof_type::nDofPerElement;
        typedef Eigen::Matrix<value_type, nDofPerElementTest, 1> local_vector_type;
        typedef Eigen::Matrix<value_type, 2*nDofPerElementTest, 1> local2_vector_type;
        typedef Eigen::Matrix<int, nDofPerElementTest, 1> local_row_type;
        typedef Eigen::Matrix<int, 2*nDofPerElementTest, 1> local2_row_type;

    public:

        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& _gmcTest,
                 map_trial_geometric_mapping_context_type const & _gmcTrial, //useless(fix compilation)
                 map_geometric_mapping_expr_context_type const& _gmcExpr,
                 ExprT const& expr,
                 IM const& im );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& _gmcTest,
                 map_trial_geometric_mapping_context_type const & _gmcTrial, //useless(fix compilation)
                 map_geometric_mapping_expr_context_type const& _gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2 );

        template<typename IMExpr, typename IMTest>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& _gmcTest,
                 map_geometric_mapping_expr_context_type const& _gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IMExpr const& imExpr, IMTest const& imTest );

        template<typename IM2>
        Context( form_type& __form,
                 map_test_geometric_mapping_context_type const& _gmcTest,
                 map_trial_geometric_mapping_context_type const & _gmcTrial, //useless(fix compilation)
                 map_geometric_mapping_expr_context_type const& _gmcExpr,
                 ExprT const& expr,
                 IM const& im,
                 IM2 const& im2,
                 mpl::int_<2> );

        bool isZero( size_type i ) const
            {
                return false;
            }
        bool isZero( typename mesh_type::element_iterator it ) const
            {
                return this->isZero( it->id() );
            }
        bool isZero( typename mesh_type::element_type const& e ) const
            {
                return this->isZero( e.id() );
            }

        void update( map_test_geometric_mapping_context_type const& _gmcTest,
                     map_trial_geometric_mapping_context_type const & gmcTrial,
                     map_geometric_mapping_expr_context_type const& _gmcExpr );


        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const & gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad );

        void update( map_test_geometric_mapping_context_type const& _gmcTest,
                     map_trial_geometric_mapping_context_type const & gmcTrial,
                     map_geometric_mapping_expr_context_type const& _gmcExpr,
                     mpl::int_<2> );


        void update( map_test_geometric_mapping_context_type const& _gmcTest,
                     map_trial_geometric_mapping_context_type const & gmcTrial,
                     map_geometric_mapping_expr_context_type const& _gmcExpr,
                     IM const& im );

        void updateInCaseOfInterpolate( map_test_geometric_mapping_context_type const& gmcTest,
                                        map_trial_geometric_mapping_context_type const & gmcTrial,
                                        map_geometric_mapping_expr_context_type const& gmcExpr,
                                        IM const& im,
                                        std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad )
        {
            M_integrator = im;
            this->updateInCaseOfInterpolate( gmcTest, gmcTrial, gmcExpr, indexLocalToQuad );
        }


        void update( map_test_geometric_mapping_context_type const& _gmcTest,
                     map_trial_geometric_mapping_context_type const & gmcTrial,
                     map_geometric_mapping_expr_context_type const& _gmcExpr,
                     IM const& im, mpl::int_<2> );


        void integrate()
        {
            integrate( mpl::int_<fusion::result_of::size<GeomapContext>::type::value>() );
        }

        void integrateInCaseOfInterpolate( std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                           bool isFirstExperience )
        {
            integrateInCaseOfInterpolate( mpl::int_<fusion::result_of::size<GeomapContext>::type::value>(),
                                          indexLocalToQuad, isFirstExperience );
        }


        void assemble()
        {
            assemble( M_gmc_left->id() );
        }
        void assemble( size_type elt_0 );

        void assemble( mpl::int_<2> )
        {
            assemble( M_gmc_left->id(), M_gmc_right->id() );
        }
        void assemble( size_type elt_0, size_type elt_1 );

        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points
         */
        template<typename Pts>
        void precomputeBasisAtPoints( Pts const& pts )
        {
            M_test_pc = test_precompute_ptrtype( new test_precompute_type( M_form.testSpace()->fe(), pts ) );
        }

        /**
         * precompute the basis function associated with the test and
         * trial space at a set of points on the face of
         * the reference element
         */
        template<typename Pts>
        void precomputeBasisAtPoints( uint16_type __f, permutation_type const& __p, Pts const& pts )
        {
            M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testSpace()->fe(), pts ) );
            //FEELPP_ASSERT( M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );
        }
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
                return  M_test_pc;

            //FEELPP_ASSERT( M_test_pc_face.find(__f )->second )( __f ).error( "invalid test precompute type" );
            return M_test_pc_face.find( __f )->second.find( __p )->second;
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts )
        {
            typedef typename boost::is_same< permutation_type, typename QuadMapped<PtsSet>::permutation_type>::type is_same_permuation_type;
            return precomputeTestBasisAtPoints( pts, mpl::bool_<is_same_permuation_type::value>() );
        }


        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts, mpl::bool_<false> )
        {
            std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > testpc;
            return testpc;
        }

        template<typename PtsSet>
        std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> >
        precomputeTestBasisAtPoints( PtsSet const& pts, mpl::bool_<true> )
        {
            //QuadMapped<PtsSet> qm;
            typedef typename QuadMapped<PtsSet>::permutation_type permutation_type;
            //typename QuadMapped<PtsSet>::permutation_points_type ppts( qm( pts ) );

            std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > testpc;

            for ( uint16_type __f = 0; __f < pts.nFaces(); ++__f )
            {
                for ( permutation_type __p( permutation_type::IDENTITY );
                        __p < permutation_type( permutation_type::N_PERMUTATIONS ); ++__p )
                {
                    //testpc[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testSpace()->fe(), ppts[__f].find( __p )->second ) );
                    testpc[__f][__p] = test_precompute_ptrtype( new test_precompute_type( M_form.testSpace()->fe(), pts.fpoints( __f,__p.value() ) ) );
                }
            }

            return testpc;
        }

    private:


        void integrate( mpl::int_<1> );

        void integrate( mpl::int_<2> );

        void integrateInCaseOfInterpolate( mpl::int_<1>,
                                           std::vector<boost::tuple<size_type,size_type> > const& indexLocalToQuad,
                                           bool isFirstExperience );

    private:

        form_type& M_form;
        dof_type* M_test_dof;
        const list_block_type& M_lb;

        test_precompute_ptrtype M_test_pc;
        std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > M_test_pc_face;

        map_test_geometric_mapping_context_type M_gmc;
        left_gmc_ptrtype M_gmc_left;
        right_gmc_ptrtype M_gmc_right;
        map_left_gmc_type M_left_map;
        map_right_gmc_type M_right_map;
        map_test_fecontext_type M_test_fec;
        map_left_test_fecontext_type M_test_fec0;
        map_right_test_fecontext_type M_test_fec1;

        local_vector_type M_rep;
        local2_vector_type M_rep_2;
        local_row_type M_local_rows;
        local2_row_type M_local_rows_2;
        local_row_type M_local_rowsigns;
        local2_row_type M_local_rowsigns_2;

        eval0_expr_ptrtype M_eval0_expr;
        eval1_expr_ptrtype M_eval1_expr;

        IM M_integrator;

    }; // Context


    //@}

    /** @name Constructors, destructor
     */
    //@{

    LinearForm( ) {}
    LinearForm( LinearForm const & __vf );

    LinearForm( space_ptrtype const& __X,
                vector_ptrtype __F,
                size_type rowstart = 0,
                bool init = true,
                bool do_threshold = false,
                value_type threshold = type_traits<value_type>::epsilon() );
    LinearForm( space_ptrtype const& __X,
                vector_ptrtype __F,
                list_block_type const& __lb,
                size_type rowstart = 0,
                bool init = true,
                bool do_threshold = false,
                value_type threshold = type_traits<value_type>::epsilon()  );

    ~LinearForm()
    {}



    //@}

    /** @name Operator overloads
     */
    //@{

    LinearForm& operator=( LinearForm const& lf )
        {
            if ( this != &lf )
            {
                M_X = lf.M_X;
                M_F = lf.M_F;
                M_lb = lf.M_lb;
                M_row_startInVector = lf.M_row_startInVector;
                M_test_pc = lf.M_test_pc;
                M_test_pc_face = lf.M_test_pc_face;
                M_do_threshold = lf.M_do_threshold;
                M_threshold = lf.M_threshold;
            }
            return *this;
        }
    /**
     * Construct the linear form given by the expression \p expr
     *
     * This function typically fills the representation that was used
     * for the \p LinearForm.
     */
    template <class ExprT>
    LinearForm& operator=( Expr<ExprT> const& expr );

    /**
     * Construct the linear form given by the expression \p expr
     *
     * This function typically fills the representation that was used
     * for the \p LinearForm by adding to it the evaluation of the
     * expression \p expr.
     */
    template <class ExprT>
    LinearForm& operator+=( Expr<ExprT> const& expr );

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
        return vflf_type( __u, M_F, __list_block, false );

    }
#endif

    /**
     * \return the entry \f$M_{i,j}\f$
     */
    value_type operator()( size_type i ) const
    {
        return M_F( i );
    }

#if 0
    /**
     * Computes the application of the form on an element of the function space
     *
     * @param __v element of Space 1 (test space)
     * @return f(v)
     */
    value_type operator()( element_type const& __v ) const
    {
        return M_F->dot( __v );
    }
#endif
    /**
     * Computes the application of the form on an element of the function space
     *
     * @param __v element of Space 1 (test space)
     * @return f(v)
     */
    value_type operator()( typename space_type::element_type const& __v ) const
    {
        return M_F->dot( __v );
    }


    //@}

    /** @name Accessors
     */
    //@{

    /**
    * \return the test function space
     */
    space_ptrtype const& functionSpace() const
    {
        return M_X;
    }

    /**
    * \return the test function space
     */
    space_ptrtype const& testSpace() const
    {
        return M_X;
    }

    /**
     * Geometric transformation
     */
    gm_ptrtype const& gm() const
    {
        return M_X->gm();
    }

    /**
     * Geometric transformation
     */
    gm1_ptrtype const& gm1() const
    {
        return M_X->gm1();
    }

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
            return  M_test_pc;

        return M_test_pc_face.find( __f )->second.find( __p )->second;
    }

    vector_type& representation() const
    {
        return *M_F;
    }
    vector_ptrtype vectorPtr() const
    {
        return M_F;
    }
    vector_ptrtype vectorPtr()
    {
        return M_F;
    }

    vector_type& vector()
    {
        return *M_F;
    }

    vector_type const& vector() const
    {
        return *M_F;
    }

    list_block_type const& blockList() const
    {
        return M_lb;
    }

    size_type rowStartInVector() const
    {
        return M_row_startInVector;
    }


    /**
     * \return \c true if threshold applies, false otherwise
     */
    bool doThreshold( value_type const& v ) const
    {
        return ( math::abs( v ) > M_threshold );
    }

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
        M_X = __X;
    }

    /**
     * set a threshold value for the matrix entries associated with the
     * bilinear form
     */
    void setThreshold( value_type eps )
    {
        M_threshold = eps;
    }

    /**
     * set the threshold strategy, true threshold the matrix entries,
     * false do not threshold
     */
    void setDoThreshold( bool do_threshold )
    {
        M_do_threshold = do_threshold;
    }

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
        M_test_pc = test_precompute_ptrtype( new test_precompute_type( functionSpace()->fe(), pts ) );
    }

    /**
      * precompute the basis function associated with the test and
      * trial space at a set of points on a face of the reference
      * element
      */
    template<typename Pts>
    void precomputeBasisAtPoints( uint16_type __f, permutation_type __p, Pts const& pts )
    {
        M_test_pc_face[__f][__p] = test_precompute_ptrtype( new test_precompute_type( functionSpace()->fe(), pts ) );
    }

    /**
     * add value \p v at position (\p i) of the vector
     * associated with the linear form
     */
    void add( size_type i,  value_type const& v )
    {
        if ( M_do_threshold )
        {
            if ( doThreshold( v ) )
                M_F->add( i+this->rowStartInVector(), v );
        }

        else
            M_F->add( i+this->rowStartInVector(), v );
    }

    /**
     * add data \p v at indices \c i of the vector
     * associated with the linear form
     */
    void addVector( int* i, int n,  value_type* v )
    {
        if ( this->rowStartInVector()!=0 )
            for ( int k = 0; k< n ; ++k )
                i[k]+=this->rowStartInVector();

        M_F->addVector( i, n, v );
    }

    /**
     * set value \p v at position (\p i) of the vector
     * associated with the linear form
     */
    void set( size_type i,  value_type const& v )
    {
        M_F->set( i, v );
    }

    /**
     * set linear form to 0
     */
    void zero() { M_F->zero(); }

    /**
     * scale linear form by \p s
     */
    void scale( value_type s ) { M_F->scale( s ); }
    
    LinearForm& operator+=( LinearForm& f )
        {
            if ( this == &f )
                return *this;

            *M_F += *f.M_F;

            return *this;
        }
    //@}

private:
    template <class ExprT> void assign( Expr<ExprT> const& expr, bool init, mpl::bool_<false> );
    template <class ExprT> void assign( Expr<ExprT> const& expr, bool init, mpl::bool_<true> );
private:

    space_ptrtype M_X;

    vector_ptrtype M_F;

    list_block_type M_lb;

    size_type M_row_startInVector;

    test_precompute_ptrtype M_test_pc;

    std::map<uint16_type, std::map<permutation_type,test_precompute_ptrtype> > M_test_pc_face;

    bool M_do_threshold;
    value_type M_threshold;
};

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( LinearForm const & __vf )
    :
    M_X( __vf.M_X ),
    M_F( __vf.M_F ),
    M_lb( __vf.M_lb ),
    M_row_startInVector( __vf.M_row_startInVector ),
    M_test_pc(),
    M_test_pc_face(),
    M_do_threshold( __vf.M_do_threshold ),
    M_threshold( __vf.M_threshold )

{
    DVLOG(2) << "LinearForm copy constructor\n";
    DVLOG(2) << "     n Dof : " << M_X->nDof() << "\n";
    DVLOG(2) << "    F size : " << M_F->size() << "\n";
    DVLOG(2) << "block size : " << M_lb.size() << "\n";
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( space_ptrtype const& __X,
        vector_ptrtype __F,
        size_type rowstart,
        bool init,
        bool do_threshold,
        value_type threshold  )
    :
    M_X( __X ),
    M_F( __F ),
    M_lb(),
    M_row_startInVector( rowstart ),
    M_do_threshold( do_threshold ),
    M_threshold( threshold )
{

    if ( !this->M_X->worldComm().isActive() ) return;

    for ( uint16_type __i = 0; __i < M_X->qDim(); ++__i )
    {
        M_lb.push_back( Block( __i, 0,
                                __i*M_X->nDofPerComponent(),
                                0 ) );
        DVLOG(2) << "[linearform::linearform] block: "
                      << Block( __i, 0, __i*M_X->nDofPerComponent(), 0 )  << "\n";
    }

    if (  init )
        M_F->zero();
}

template<typename SpaceType, typename VectorType,  typename ElemContType>
LinearForm<SpaceType, VectorType, ElemContType>::LinearForm( space_ptrtype const& __X,
        vector_ptrtype __F,
        list_block_type const& __lb,
        size_type rowstart,
        bool init,
        bool do_threshold,
        value_type threshold )
    :
    M_X( __X ),
    M_F( __F ),
    M_lb( __lb ),
    M_row_startInVector( rowstart ),
    M_do_threshold( do_threshold ),
    M_threshold( threshold )
{
    if ( !this->M_X->worldComm().isActive() ) return;

    if ( init )
        M_F->zero();
}

template<typename LFType, typename TheSpaceType, typename ExprType>
struct LFAssign
{
    LFAssign( LFAssign const& lfa )
        :
        M_lf( lfa.M_lf ),
        M_Xh( lfa.M_Xh ),
        M_expr( lfa.M_expr ),
        M_index( lfa.M_index ),
        M_init( lfa.M_init )
    {}
    LFAssign( LFType& lf, TheSpaceType const& Xh, ExprType const& expr, bool init )
        :
        M_lf( lf ),
        M_Xh( Xh ),
        M_expr( expr ),
        M_index( 0 ),
        M_init( init )
    {}
    template<typename SpaceType>
    void operator()( boost::shared_ptr<SpaceType> const& X ) const
    {
        if ( M_lf.testSpace()->worldsComm()[M_index].isActive() )
        {
            DVLOG(2) << "expression has test functions ? :"
                          << ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result
                    << "\n";

            if ( !ExprType::template HasTestFunction<typename SpaceType::reference_element_type>::result )
            {
                ++M_index;
                return;
            }

            list_block_type __list_block;

            if ( M_lf.testSpace()->worldsComm()[M_index].globalSize()>1 )
                {
                    if (M_lf.testSpace()->hasEntriesForAllSpaces())
                        __list_block.push_back( Block( 0, 0, M_Xh->nLocalDofStart( M_index ), 0 ) );
                    else
                        __list_block.push_back( Block( 0, 0, 0, 0 ) );
                }
            else
                __list_block.push_back( Block( 0, 0, M_Xh->nDofStart( M_index ), 0 ) );

            LinearForm<SpaceType,typename LFType::vector_type, typename LFType::element_type> lf( X,
                    M_lf.vectorPtr(),
                    __list_block,
                    M_lf.rowStartInVector(),
                    false );

            //
            // in composite integration, make sure that if M_init is \p
            // true for the first space, it is set to \p false for the
            // next spaces otherwise it will erase/clear to 0 the previous
            // assemblies
            //
            if ( M_init )
            {
                // assembly
                lf = M_expr;

                // make sure we won't erase the assembly we just did
                M_init = false;
            }

            else
            {
                lf += M_expr;
            }
        }

        else // not active : there is the init case with a close in zero
        {
            if ( M_init ) M_lf.representation().zero();
        }

        ++M_index;
    }
private:
    LFType& M_lf;
    TheSpaceType const& M_Xh;
    ExprType const& M_expr;
    mutable size_type M_index;
    mutable bool M_init;

};

// implementation
template<typename LFType, typename TheSpaceType, typename ExprType>
LFAssign<LFType,TheSpaceType,ExprType>
make_lfassign( LFType& lf, TheSpaceType const& Xh, ExprType const& expr, bool init )
{
    return LFAssign<LFType,TheSpaceType,ExprType>( lf, Xh, expr, init );
}



template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
void
LinearForm<SpaceType, VectorType, ElemContType>::assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<true> )
{
    fusion::for_each( M_X->functionSpaces(), make_lfassign( *this, M_X, __expr, init ) );
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
void
LinearForm<SpaceType, VectorType, ElemContType>::assign( Expr<ExprT> const& __expr, bool init, mpl::bool_<false> )
{
    if ( M_lb.empty() )
    {
        M_lb.push_back( Block( 0, 0, 0, 0 ) );
    }

    if ( init )
    {
        typename list_block_type::const_iterator __bit = M_lb.begin();
        typename list_block_type::const_iterator __ben = M_lb.end();

        for ( ; __bit != __ben; ++__bit )
        {
            DVLOG(2) << "LinearForm:: block: " << *__bit << "\n";
            size_type g_ic_start = __bit->globalRowStart();
            DVLOG(2) << "LinearForm:: g_ic_start: " << g_ic_start << "\n";

            M_F->zero( g_ic_start,g_ic_start + M_X->nDof() );
        }
    }

    __expr.assemble( M_X, *this );
    //vector_range_type r( M_F, ublas::range( M_v.start(),M_v.start()+ M_v.size() ) );
    //std::cout << "r = " << r << "\n";
    //std::cout << "after F= " << M_F << "\n";
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
LinearForm<SpaceType, VectorType, ElemContType>&
LinearForm<SpaceType, VectorType, ElemContType>::operator=( Expr<ExprT> const& __expr )
{
    // loop(fusion::for_each) over sub-functionspaces in SpaceType
    // pass expression and initialize
    this->assign( __expr, true, mpl::bool_<( SpaceType::nSpaces > 1 )>() );
    return *this;
}
template<typename SpaceType, typename VectorType,  typename ElemContType>
template<typename ExprT>
LinearForm<SpaceType, VectorType, ElemContType>&
LinearForm<SpaceType, VectorType, ElemContType>::operator+=( Expr<ExprT> const& __expr )
{
    this->assign( __expr, false, mpl::bool_<( SpaceType::nSpaces > 1 )>() );
    return *this;
}


} // detail
/// \endcond
} // vf

namespace meta
{
template<typename SpaceType,typename VectorType=typename Backend<typename SpaceType::value_type>::vector_type,typename ElemContType=typename Backend<typename SpaceType::value_type>::vector_type>
struct LinearForm
{
    typedef Feel::vf::detail::LinearForm<SpaceType,VectorType,ElemContType> type;
};
}
} // feel

#include <feel/feelvf/linearformcontext.hpp>
#endif /* __LinearForm_H */
