/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 2013-11-19

  Copyright (C) 2012 Universite Joseph Fourier (Grenoble I)

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
   \file fluidmecstresstensor.hpp
   \author Vincent Chabannes
   \date 2013-11-19
 */
#ifndef FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H
#define FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{

enum FMSTExprApplyType { FM_ST_EVAL=0,FM_ST_JACOBIAN=1,FM_VISCOSITY_EVAL=2 };


    namespace detail
    {
        template < int nDim, typename SpecificExprType >
        struct GetShapeExpr
        {
            typedef Shape<nDim, Tensor2, false, false> shape_grad_type;
            typedef Shape<nDim, Scalar, false, false> shape_scalar_type;
            typedef typename mpl::if_<mpl::or_< boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                                boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> > >,
                                      shape_grad_type,
                                      shape_scalar_type >::type shape_type;
        };
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorBase : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                                           typename detail::GetShapeExpr< ExprType::expr_grad_velocity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape::nDim, typename ExprType::specific_expr_type>::shape_type,
                                                           typename ExprType::value_type >
    {
        typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                           typename detail::GetShapeExpr<ExprType::expr_grad_velocity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t>::shape::nDim, typename ExprType::specific_expr_type>::shape_type,
                           typename ExprType::value_type> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::key_type key_type;

        //typedef typename expr_type::geoelement_type geoelement_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::gmc_type gmc_type;

        typedef typename expr_type::expr_grad_velocity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_grad_velocity_type;
        typedef typename expr_type::expr_id_pressure_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_id_pressure_type;

        typedef typename expr_type::value_type value_type;

        typedef typename super_type::loc_tensor2_type loc_grad_type;
        typedef typename super_type::new_array_tensor2_type array_grad_type;
        typedef typename super_type::loc_scalar_type loc_scalar_type;
        typedef typename super_type::new_array_scalar_type array_scalar_type;
        typedef typename super_type::array_scalar_type old_array_scalar_type;

        typedef typename super_type::shape_type shape;
        typedef typename super_type::matrix_shape_type matrix_shape_type;

        typedef typename super_type::new_array_shape_type array_shape_type;
        using ret_type = Eigen::Map<Eigen::Matrix<value_type, shape::M, shape::N> const>;

        tensorFluidStressTensorBase( expr_type const& expr,
                                     Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom,fev,feu ),
            M_tensorExprIdPressure( expr.expressionIdPressure(),geom,fev,feu ),
            M_locRes( 0 ),
            M_locGradVelocity( 0 ),
            M_locPressure( 0 )
        {}
        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom,fev ),
            M_tensorExprIdPressure( expr.expressionIdPressure(),geom,fev ),
            M_locRes( 0 ),
            M_locGradVelocity( 0 ),
            M_locPressure( 0 )
        {}

        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom ),
            M_tensorExprIdPressure( expr.expressionIdPressure(),geom ),
            M_locRes( 0 ),
            M_locGradVelocity( 0 ),
            M_locPressure( 0 )
        {}


        expr_type const& expr() const { return M_expr; }

        array_grad_type const& locGradVelocity() const { return M_locGradVelocity; }
        loc_grad_type const& locGradVelocity( uint16_type q ) const { return M_locGradVelocity[q]; }
        array_scalar_type const& locPressure() const { return M_locPressure; }
        loc_scalar_type const& locPressure( uint16_type q ) const { return M_locPressure[q]; }

        array_shape_type & locRes() { return M_locRes; }
        array_shape_type const& locRes() const { return M_locRes; }
        matrix_shape_type & locRes( uint16_type q ) { return M_locRes[q]; }
        matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }
    protected :
        void updateCommon( Geo_t const& geom, bool upGradVelocity = true )
        {
            this->setGmc( geom );

            if ( upGradVelocity )
                M_tensorExprGradVelocity.update( geom );

            bool upIdPressure = this->expr().withPressureTerm();
            if ( upIdPressure )
                M_tensorExprIdPressure.update( geom );

            this->updateDataContainer( upGradVelocity, upIdPressure );
        }

        void updateCommon( Geo_t const& geom, uint16_type face, bool upGradVelocity = true  )
        {
            this->setGmc( geom );

            if ( upGradVelocity )
                M_tensorExprGradVelocity.update( geom, face );

            bool upIdPressure = this->expr().withPressureTerm();
            if ( upIdPressure )
                M_tensorExprIdPressure.update( geom, face );

            this->updateDataContainer( upGradVelocity, upIdPressure );
        }
    private :
        void updateDataContainer( bool upGradVelocity, bool upIdPressure )
            {
                uint16_type nPoints = this->gmc()->nPoints();
                if ( upGradVelocity )
                {
                    if ( M_locGradVelocity.size() != nPoints )
                        M_locGradVelocity.setConstant( nPoints, this->M_zeroLocTensor2 );
                    for ( uint16_type q=0;q< nPoints;++q )
                    {
                        if constexpr ( false /*expr_type::expr_grad_velocity_type::is_terminal*/ )
                                     {
                                         M_locGradVelocity[q] = M_tensorExprGradVelocity.evalq( q ); //not compile, need to investigate
                                     }
                        else
                        {
                            typename super_type::loc_tensor2_type& locData = M_locGradVelocity[q];
                            for (uint16_type c1=0;c1<super_type::shape_tensor2::M;++c1 )
                                for (uint16_type c2=0;c2<super_type::shape_tensor2::N;++c2 )
                                    locData(c1,c2) = M_tensorExprGradVelocity.evalq( c1,c2,q );
                        }
                    }
                }
                if ( upIdPressure )
                {
                    if ( M_locPressure.size() != nPoints )
                        M_locPressure.setConstant( nPoints, this->M_zeroLocScalar );
                    for ( uint16_type q=0;q< nPoints;++q )
                        M_locPressure[q](0,0) = M_tensorExprIdPressure.evalq( 0,0,q );
                }

                if ( M_locRes.size() != nPoints )
                    M_locRes.setConstant( nPoints, super_type::matrix_shape_type::Zero() );
            }
    public :

        /*virtual*/
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return ret_type(this->M_locMatrixShape.data());//evalijq( i,j,q,mpl::int_<gmc_type::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if constexpr ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL ||
                 expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                return evalq( c1,c2,q );
            else
            {
                //CHECK( false ) << "not allow";
                //LOG(WARNING) << "evalijq non optimized";
                return this->evalijq( i,j,q )( c1,c2 );
                //return value_type(0);
            }
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return ret_type(M_locRes[q].data());
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1, c2, q, mpl::int_<expr_type::specific_expr_type::value>() );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return ret_type(M_locRes[q].data());
        }

    private :
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<FMSTExprApplyType::FM_ST_EVAL> /**/ ) const
            {
                return M_locRes[q]( c1,c2 );
            }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> /**/ ) const
            {
                CHECK( false ) << "not allow";
                return M_locRes[q]( c1,c2 );
            }
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL> /**/ ) const
            {
                return M_locRes[q]( 0,0 );
            }

        expr_type const& M_expr;

        tensor_expr_grad_velocity_type M_tensorExprGradVelocity;
        tensor_expr_id_pressure_type M_tensorExprIdPressure;
        array_shape_type M_locRes;
        array_grad_type M_locGradVelocity;
        array_scalar_type M_locPressure;
    };


/**
     * Newtonian :
     *   sigma = -pI + 2 \mu D(u) with D(u) = 0.5*(\grad u + (\grad u)^T)
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorNewtonian : public tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
    {
        typedef tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    public :
        typedef ExprType expr_type;
        //typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        using ret_type = typename super_type::ret_type;

        tensorFluidStressTensorNewtonian( expr_type const& expr,
                                          Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            //M_muExprTensor( this->expr().materialProperties().property("dynamic-viscosity").template expr<1,1>().evaluator( geom ) )
            M_muExprTensor( this->expr().template materialPropertyExpr<1,1>("dynamic-viscosity").evaluator( geom ) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_muExprTensor( this->expr().template materialPropertyExpr<1,1>("dynamic-viscosity").evaluator( geom ) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_muExprTensor( this->expr().template materialPropertyExpr<1,1>("dynamic-viscosity").evaluator( geom ) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            M_muExprTensor.update( geom );
            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            M_muExprTensor.update( geom, face );
            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type muEval = M_muExprTensor.evalq(0,0,q);
                auto & locResAtQ = this->locRes(q);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                    locResAtQ(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    locResAtQ(1,0) = muEval*( du2vdx + du1vdy );
                    locResAtQ(0,1) = this->locRes(q)(1,0);
                    locResAtQ(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        locResAtQ(0,0) -= pres;
                        locResAtQ(1,1) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL ||
                          expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    locResAtQ(0,0) = muEval;
                }

            }
        }
        void updateImpl( mpl::int_<3> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type muEval = M_muExprTensor.evalq(0,0,q);
                auto & locResAtQ = this->locRes(q);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                    locResAtQ(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    locResAtQ(1,0) = muEval*( du2vdx + du1vdy );
                    locResAtQ(2,0) = muEval*( du3vdx + du1vdz );
                    locResAtQ(0,1) = this->locRes(q)(1,0);
                    locResAtQ(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    locResAtQ(1,2) = muEval*( du2vdz + du3vdy );
                    locResAtQ(2,1) = this->locRes(q)(1,2);
                    locResAtQ(0,2) = this->locRes(q)(2,0);
                    locResAtQ(2,2) = /*-pres +*/ 2*muEval*du3vdz;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        locResAtQ(0,0) -= pres;
                        locResAtQ(1,1) -= pres;
                        locResAtQ(2,2) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL ||
                          expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    locResAtQ(0,0) = muEval;
                }
            }
        }
        using super_type::evalijq; // fix clang warning

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            const value_type muEval =  this->locRes(q)(0,0);//M_muExprTensor.evalq(0,0,q);
            //M_locMatrixShape( matrix_shape_type::Zero() )
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;

            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            const value_type muEval =  this->locRes(q)(0,0);//M_muExprTensor.evalq(0,0,q);
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(2,0) = muEval*( du3tdx + du1tdz );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;
            this->locMatrixShape()(2,1) = muEval*( du2tdz + du3tdy );
            this->locMatrixShape()(0,2) = this->locMatrixShape()(2,0);
            this->locMatrixShape()(1,2) = this->locMatrixShape()(2,1);
            this->locMatrixShape()(2,2) = 2*muEval*du3tdz;

            return ret_type(this->locMatrixShape().data());
        }

    private :

        typename expr_type::template material_property_expr_type<1,1>::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/>  M_muExprTensor;
    };

    /**
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     * Power Law :
     *   sigma = -pI + 2 \mu(\gammap) D(u) with D(u) = 0.5*(\grad u + (\grad u)^T)
     *   \mu(\gammap) = ...
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorPowerLaw : public tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
    {
        typedef tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::value_type value_type;
        typedef typename super_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename super_type::ret_type;

        tensorFluidStressTensorPowerLaw( expr_type const& expr,
                                         Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() ),
            M_kExprTensor( this->expr().template materialPropertyExpr<1,1>("consistency-index").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("power-law-index").evaluator( geom ) ),
            M_muMinExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-min").evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-max").evaluator( geom ) )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() ),
            M_kExprTensor( this->expr().template materialPropertyExpr<1,1>("consistency-index").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("power-law-index").evaluator( geom ) ),
            M_muMinExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-min").evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-max").evaluator( geom ) )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() ),
            M_kExprTensor( this->expr().template materialPropertyExpr<1,1>("consistency-index").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("power-law-index").evaluator( geom ) ),
            M_muMinExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-min").evaluator( geom ) ),
            M_muMaxExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-max").evaluator( geom ) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );
            M_kExprTensor.update( geom );
            M_nExprTensor.update( geom );
            M_muMinExprTensor.update( geom );
            M_muMaxExprTensor.update( geom );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            M_muIsInInterval.resize( this->locRes().size() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );
            M_kExprTensor.update( geom,face );
            M_nExprTensor.update( geom,face );
            M_muMinExprTensor.update( geom,face );
            M_muMaxExprTensor.update( geom,face );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            M_muIsInInterval.resize( this->locRes().size() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type power_k_generic = M_kExprTensor.evalq(0,0,q);
                const value_type power_n_generic = ( M_nExprTensor.evalq(0,0,q) - 1.)/2.;
                const value_type muMin = M_muMinExprTensor.evalq(0,0,q);
                const value_type muMax = M_muMaxExprTensor.evalq(0,0,q);

                //auto const mu_powerlaw = power_k_generic*pow( 2.0*inner(defv,defv) /*+chiInv*/ , cst( power_n_generic ) )/**chiSup*/;
                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                const value_type gammapoint2v = 2*DxD;
                value_type muEval = power_k_generic*math::pow( gammapoint2v , power_n_generic );
                if ( muEval < muMin )
                {
                    muEval = muMin;
                    M_muIsInInterval[q] = false;
                }
                else if ( muEval > muMax )
                {
                    muEval = muMax;
                    M_muIsInInterval[q] = false;
                }
                else
                    M_muIsInInterval[q] = true;
                //-idv(p)*Id + 2*mu_powerlaw*defv;
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    M_locEvalPrecompute[q] = power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        void updateImpl( mpl::int_<3> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type power_k_generic = M_kExprTensor.evalq(0,0,q);
                const value_type power_n_generic = ( M_nExprTensor.evalq(0,0,q) - 1.)/2.;
                const value_type muMin = M_muMinExprTensor.evalq(0,0,q);
                const value_type muMax = M_muMaxExprTensor.evalq(0,0,q);

                //auto const mu_powerlaw = power_k_generic*pow( 2.0*inner(defv,defv) /*+chiInv*/ , cst( power_n_generic ) )/**chiSup*/;
                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                const value_type DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                    0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                const value_type gammapoint2v = 2*DxD;
                value_type muEval = power_k_generic*math::pow( gammapoint2v , power_n_generic );
                if ( muEval < muMin )
                {
                    muEval = muMin;
                    M_muIsInInterval[q] = false;
                }
                else if ( muEval > muMax )
                {
                    muEval = muMax;
                    M_muIsInInterval[q] = false;
                }
                else
                    M_muIsInInterval[q] = true;

                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(2,0) = muEval*( du3vdx + du1vdz );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    this->locRes(q)(1,2) = muEval*( du2vdz + du3vdy );
                    this->locRes(q)(2,1) = this->locRes(q)(1,2);
                    this->locRes(q)(0,2) = this->locRes(q)(2,0);
                    this->locRes(q)(2,2) = /*-pres +*/ 2*muEval*du3vdz;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                        this->locRes(q)(2,2) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    M_locEvalPrecompute[q] = power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }

        using super_type::evalijq; // fix clang warning

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q];

            matrix_shape_type & thelocRes = this->locMatrixShape();
            if ( M_muIsInInterval[q] )
            {
                //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
                const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
                //const value_type mut = gammapoint2t*power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                const value_type mut = gammapoint2t*M_locEvalPrecompute[q];
                thelocRes(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
                thelocRes(1,0) = muEval*( du2tdx + du1tdy ) + mut*( du2vdx + du1vdy );
                thelocRes(0,1) = thelocRes(1,0);
                thelocRes(1,1) = 2*(muEval*du2tdy + mut*du2vdy );
            }
            else
            {
                thelocRes(0,0) = 2*(muEval*du1tdx);
                thelocRes(1,0) = muEval*( du2tdx + du1tdy );
                thelocRes(0,1) = thelocRes(1,0);
                thelocRes(1,1) = 2*(muEval*du2tdy );
            }

            return ret_type(thelocRes.data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q];
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);

            const value_type DxDt01 = du1tdy+du2tdx;
            const value_type DxDv01 = du1vdy+du2vdx;
            const value_type DxDt02 = du1tdz+du3tdx;
            const value_type DxDv02 = du1vdz+du3vdx;
            const value_type DxDt12 = du2tdz+du3tdy;
            const value_type DxDv12 = du2vdz+du3vdy;

            matrix_shape_type & thelocRes = this->locMatrixShape();
            if ( M_muIsInInterval[q] )
            {
                const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                    DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
                //const value_type mut = gammapoint2t*power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                const value_type mut = gammapoint2t*M_locEvalPrecompute[q];

                thelocRes(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
                thelocRes(1,0) = muEval*( DxDt01 ) + mut*( DxDv01 );
                thelocRes(2,0) = muEval*( DxDt02 ) + mut*( DxDv02 );
                thelocRes(0,1) = thelocRes(1,0);
                thelocRes(1,1) = 2*(muEval*du2tdy + mut*du2vdy);
                thelocRes(2,1) = muEval*( DxDt12 ) + mut*( DxDv12 );
                thelocRes(0,2) = thelocRes(2,0);
                thelocRes(1,2) = thelocRes(2,1);
                thelocRes(2,2) = 2*(muEval*du3tdz + mut*du3vdz);
            }
            else
            {
                thelocRes(0,0) = 2*(muEval*du1tdx);
                thelocRes(1,0) = muEval*( DxDt01 );
                thelocRes(2,0) = muEval*( DxDt02 );
                thelocRes(0,1) = thelocRes(1,0);
                thelocRes(1,1) = 2*(muEval*du2tdy);
                thelocRes(2,1) = muEval*( DxDt12 );
                thelocRes(0,2) = thelocRes(2,0);
                thelocRes(1,2) = thelocRes(2,1);
                thelocRes(2,2) = 2*(muEval*du3tdz);
            }

            return ret_type(thelocRes.data());
        }
    private :

        typename super_type::array_value_type M_locEvalPrecompute;
        typename super_type::array_value_type M_locEvalMu;
        std::vector<bool> M_muIsInInterval;
        typename expr_type::template material_property_expr_type<1,1>::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/>  M_kExprTensor, M_nExprTensor, M_muMinExprTensor, M_muMaxExprTensor;

    };

    /**
     *
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorCarreau : public tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
    {
        typedef tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::value_type value_type;
        using ret_type = typename super_type::ret_type;

        tensorFluidStressTensorCarreau( expr_type const& expr,
                                        Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-n").evaluator( geom ) )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-n").evaluator( geom ) )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-law-n").evaluator( geom ) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );

            M_lambdaExprTensor.update( geom );
            M_nExprTensor.update( geom );
            M_mu0ExprTensor.update( geom );
            M_muInfExprTensor.update( geom );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );

            M_lambdaExprTensor.update( geom, face );
            M_nExprTensor.update( geom, face );
            M_mu0ExprTensor.update( geom, face );
            M_muInfExprTensor.update( geom, face );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreau_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreau_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreau_lambda_pow2_time2 = math::pow(carreau_lambda,2.)*2.0;
                const value_type carreauValPower =  (carreau_n-1)/2.0;
                const value_type carreau_lambda2 = math::pow(carreau_lambda,2);
                const value_type part1_carreauLaw = ( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2;

                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreau_lambda_pow2_time2*DxD, carreauValPower );

                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    const value_type gammapoint2v = 2*DxD;
                    M_locEvalPrecompute[q] = part1_carreauLaw*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }
            }
        }
        void updateImpl( mpl::int_<3> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreau_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreau_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreau_lambda_pow2_time2 = math::pow(carreau_lambda,2.)*2.0;
                const value_type carreauValPower =  (carreau_n-1)/2.0;
                const value_type carreau_lambda2 = math::pow(carreau_lambda,2);
                const value_type part1_carreauLaw = ( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2;

                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                const value_type DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                    0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreau_lambda_pow2_time2*DxD, carreauValPower );

                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(2,0) = muEval*( du3vdx + du1vdz );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    this->locRes(q)(1,2) = muEval*( du2vdz + du3vdy );
                    this->locRes(q)(2,1) = this->locRes(q)(1,2);
                    this->locRes(q)(0,2) = this->locRes(q)(2,0);
                    this->locRes(q)(2,2) = /*-pres +*/ 2*muEval*du3vdz;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                        this->locRes(q)(2,2) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    //const value_type mut = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                    const value_type gammapoint2v = 2*DxD;
                    M_locEvalPrecompute[q] = part1_carreauLaw*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        using super_type::evalijq; // fix clang warning

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim>() );
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/ , mpl::int_<2> /**/ ) const
        {
            const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
            const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
            /*const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);*/
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q];
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
            //const value_type mut = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q];

            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy ) + mut*( du2vdx + du1vdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy );

            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
            const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
            /*const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);*/
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q];
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);

            const value_type DxDt01 = du1tdy+du2tdx;
            const value_type DxDv01 = du1vdy+du2vdx;
            const value_type DxDt02 = du1tdz+du3tdx;
            const value_type DxDv02 = du1vdz+du3vdx;
            const value_type DxDt12 = du2tdz+du3tdy;
            const value_type DxDv12 = du2vdz+du3vdy;
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
            //const value_type mut = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q];

            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( DxDt01 ) + mut*( DxDv01 );
            this->locMatrixShape()(2,0) = muEval*( DxDt02 ) + mut*( DxDv02 );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy);
            this->locMatrixShape()(2,1) = muEval*( DxDt12 ) + mut*( DxDv12 );
            this->locMatrixShape()(0,2) = this->locMatrixShape()(2,0);
            this->locMatrixShape()(1,2) = this->locMatrixShape()(2,1);
            this->locMatrixShape()(2,2) = 2*(muEval*du3tdz + mut*du3vdz);

            return ret_type(this->locMatrixShape().data());
        }
    private :
        typename super_type::array_value_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_value_type M_locEvalMu;
        typename expr_type::template material_property_expr_type<1,1>::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/> M_mu0ExprTensor, M_muInfExprTensor, M_lambdaExprTensor, M_nExprTensor;
    };


    /**
     *
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorCarreauYasuda : public tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
    {
        typedef tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::value_type value_type;
        using ret_type = typename super_type::ret_type;

        tensorFluidStressTensorCarreauYasuda( expr_type const& expr,
                                              Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-n").evaluator( geom ) ),
            M_aExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-a").evaluator( geom ) )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-n").evaluator( geom ) ),
            M_aExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-a").evaluator( geom ) )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_mu0ExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-zero-shear").evaluator( geom ) ),
            M_muInfExprTensor( this->expr().template materialPropertyExpr<1,1>("viscosity-infinite-shear").evaluator( geom ) ),
            M_lambdaExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-lambda").evaluator( geom ) ),
            M_nExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-n").evaluator( geom ) ),
            M_aExprTensor( this->expr().template materialPropertyExpr<1,1>("carreau-yasuda-law-a").evaluator( geom ) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );

            M_mu0ExprTensor.update( geom );
            M_muInfExprTensor.update( geom );
            M_lambdaExprTensor.update( geom );
            M_nExprTensor.update( geom );
            M_aExprTensor.update( geom );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );

            M_mu0ExprTensor.update( geom, face );
            M_muInfExprTensor.update( geom, face );
            M_lambdaExprTensor.update( geom, face );
            M_nExprTensor.update( geom, face );
            M_aExprTensor.update( geom, face );

            //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_a = M_aExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda_pow_a = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type carreauYasudaValPower = carreauYasuda_a/2.;
                const value_type carreauYasudaValPower2 = (carreauYasuda_n-1)/carreauYasuda_a;
                const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);

                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                const value_type gammapoint2v = 2*DxD;
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreauYasuda_lambda_pow_a*math::pow( gammapoint2v, carreauYasudaValPower) , carreauYasudaValPower2 );

                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
                    const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
                    //const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                    M_locEvalPrecompute[q] = part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        void updateImpl( mpl::int_<3> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
                const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda = M_lambdaExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_n = M_nExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_a = M_aExprTensor.evalq(0,0,q);
                const value_type carreauYasuda_lambda_pow_a = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type carreauYasudaValPower = carreauYasuda_a/2.;
                const value_type carreauYasudaValPower2 = (carreauYasuda_n-1)/carreauYasuda_a;

                const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
                const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);


                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
                const value_type DxD = math::pow(du1vdx,2) + math::pow(du2vdy,2) + math::pow(du3vdz,2) +
                    0.5*( math::pow(du2vdx+du1vdy,2) + math::pow(du1vdz+du3vdx,2) + math::pow(du2vdz+du3vdy,2) );
                const value_type gammapoint2v = 2*DxD;
                const value_type muEval = mu_inf + (mu_0 - mu_inf)*math::pow( 1. + carreauYasuda_lambda_pow_a*math::pow( gammapoint2v, carreauYasudaValPower) , carreauYasudaValPower2 );

                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    this->locRes(q)(0,0) = /*-pres +*/ 2*muEval*du1vdx;
                    this->locRes(q)(1,0) = muEval*( du2vdx + du1vdy );
                    this->locRes(q)(2,0) = muEval*( du3vdx + du1vdz );
                    this->locRes(q)(0,1) = this->locRes(q)(1,0);
                    this->locRes(q)(1,1) = /*-pres +*/ 2*muEval*du2vdy;
                    this->locRes(q)(1,2) = muEval*( du2vdz + du3vdy );
                    this->locRes(q)(2,1) = this->locRes(q)(1,2);
                    this->locRes(q)(0,2) = this->locRes(q)(2,0);
                    this->locRes(q)(2,2) = /*-pres +*/ 2*muEval*du3vdz;
                    if ( withPressure )
                    {
                        const value_type pres = this->locPressure(q)(0,0);
                        this->locRes(q)(0,0) -= pres;
                        this->locRes(q)(1,1) -= pres;
                        this->locRes(q)(2,2) -= pres;
                    }
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN )
                {
                    //M_locEvalDxD[q](0,0) = DxD;
                    M_locEvalMu[q] = muEval;
                    const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
                    const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
                    //const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                    M_locEvalPrecompute[q] = part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );

                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            //return evalijq( i,j,q,mpl::int_<expr_type::nRealDim>() );
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        using super_type::evalijq; // fix clang warning

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
            const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);

            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q];
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
#if 0
            const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);
            const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
            const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
            const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
#else
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q];
#endif
            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy ) + mut*( du2vdx + du1vdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy );

            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type mu_inf = M_muInfExprTensor.evalq(0,0,q);
            const value_type mu_0 = M_mu0ExprTensor.evalq(0,0,q);

            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q];//(0,0);

            const value_type DxDt01 = du1tdy+du2tdx;
            const value_type DxDv01 = du1vdy+du2vdx;
            const value_type DxDt02 = du1tdz+du3tdx;
            const value_type DxDv02 = du1vdz+du3vdx;
            const value_type DxDt12 = du2tdz+du3tdy;
            const value_type DxDv12 = du2vdz+du3vdy;
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
#if 0
            const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);
            const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
            const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
            const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
#else
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q];
#endif
            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( DxDt01 ) + mut*( DxDv01 );
            this->locMatrixShape()(2,0) = muEval*( DxDt02 ) + mut*( DxDv02 );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy);
            this->locMatrixShape()(2,1) = muEval*( DxDt12 ) + mut*( DxDv12 );
            this->locMatrixShape()(0,2) = this->locMatrixShape()(2,0);
            this->locMatrixShape()(1,2) = this->locMatrixShape()(2,1);
            this->locMatrixShape()(2,2) = 2*(muEval*du3tdz + mut*du3vdz);

            return ret_type(this->locMatrixShape().data());
        }
    private :
        typename super_type::array_value_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_value_type M_locEvalMu;
        typename expr_type::template material_property_expr_type<1,1>::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/> M_mu0ExprTensor, M_muInfExprTensor, M_lambdaExprTensor, M_nExprTensor, M_aExprTensor;
    };

/// \cond detail
/**
 * \class FluidMecStressTensorImpl
 * \brief det of a matrix
 *
 * @author Vincent Chabannes
 * @see
 */
template<typename ExprGradVelocityType, typename FiniteElementVelocityType, typename ExprIdPressureType, typename ModelPhysicFluidType, typename SymbolsExprType, typename SpecificExprType>
class FluidMecStressTensorImpl
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FluidMecStressTensorImpl<ExprGradVelocityType,FiniteElementVelocityType,ExprIdPressureType,ModelPhysicFluidType,SymbolsExprType,SpecificExprType> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_pressure = vm::JACOBIAN;
    static const size_type context_muP0 = vm::JACOBIAN;
    static const size_type context = context_velocity;

    typedef ExprGradVelocityType expr_grad_velocity_type;
    typedef ExprIdPressureType expr_id_pressure_type;
    using model_physic_fluid_type = ModelPhysicFluidType;
    using symbols_expr_type = SymbolsExprType;

    typedef SpecificExprType specific_expr_type;
    //------------------------------------------------------------------------------//
    typedef typename ExprGradVelocityType::value_type value_type;
    //------------------------------------------------------------------------------//
 
    static const bool is_terminal = true;

    //static const uint16_type orderpressure = functionspace_pressure_type::basis_type::nOrder;

    template<typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> >,
                                            typename mpl::if_<boost::is_same<Func,FiniteElementVelocityType>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type,
                                            mpl::bool_<false> >::type::value;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    FluidMecStressTensorImpl( expr_grad_velocity_type const& exprGradVelocity, expr_id_pressure_type const& exprIdPressure,
                              model_physic_fluid_type const& physicFluid,
                              MaterialProperties const& matProps,
                              uint16_type polyOrder,
                              bool withPressure,
                              SymbolsExprType const& se)
        :
        M_exprGradVelocity( exprGradVelocity ),
        M_exprIdPressure( exprIdPressure ),
        M_physicFluid( physicFluid ),
        M_matProps( matProps ),
        M_polynomialOrder( polyOrder ),
        M_withPressureTerm( withPressure ),
        M_se( se )
    {}

    FluidMecStressTensorImpl( FluidMecStressTensorImpl const & op ) = default;

    ~FluidMecStressTensorImpl()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! polynomial order
    uint16_type polynomialOrder() const
        {
            if ( M_polynomialOrder != invalid_uint16_type_value )
                return M_polynomialOrder;
            uint16_type orderGradVelocity = M_exprGradVelocity.polynomialOrder();
            if ( this->dynamicViscosity().isNewtonianLaw() )
            {
                if ( SpecificExprType::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                    return this->materialPropertyExpr<1,1>("dynamic-viscosity").polynomialOrder();
                else if ( SpecificExprType::value == FMSTExprApplyType::FM_ST_EVAL )
                    return std::max( orderGradVelocity, M_exprIdPressure.polynomialOrder() );
                else
                    return orderGradVelocity;
            }
            else
                return 2*(orderGradVelocity+1);
        }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }


    template <typename OtherSymbolsExprType>
    auto applySymbolsExpr( OtherSymbolsExprType const& se ) const
        {
            auto newse = Feel::vf::symbolsExpr( M_se, se );
            using new_se_type = std::decay_t<decltype(newse)>;
            using new_this_type = FluidMecStressTensorImpl<ExprGradVelocityType,FiniteElementVelocityType,ExprIdPressureType,ModelPhysicFluidType,new_se_type,SpecificExprType>;
            return new_this_type( M_exprGradVelocity,M_exprIdPressure,M_physicFluid,M_matProps,M_polynomialOrder,M_withPressureTerm,newse );
        }


    expr_grad_velocity_type const& expressionGradVelocity() const { return M_exprGradVelocity; }
    expr_id_pressure_type const& expressionIdPressure() const { return M_exprIdPressure; }
    auto const& dynamicViscosity() const { return M_physicFluid.dynamicViscosity(); }
    MaterialProperties const& materialProperties() const { return M_matProps; }
    bool withPressureTerm() const { return M_withPressureTerm; }


    template <int M,int N>
    using material_property_expr_type = std::decay_t<decltype( expr( ModelExpression{}.template expr<M,N>(),symbols_expr_type{} ) )>;

    template <int M,int N>
    material_property_expr_type<M,N> materialPropertyExpr( std::string const& prop ) const { return expr( this->materialProperties().property( prop ).template expr<M,N>(), M_se ); }

    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename this_type::value_type value_type;

        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        //----------------------------------------------------------------------------------------------------//

        typedef typename this_type::expr_grad_velocity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_grad_velocity_type;
        typedef typename detail::GetShapeExpr<tensor_expr_grad_velocity_type::shape::nDim, SpecificExprType>::shape_type shape;

        //----------------------------------------------------------------------------------------------------//

        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t,shape,value_type> tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename tensorbase_type::ret_type;
        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            this->initTensor( expr, geom, fev, feu );
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            this->initTensor( expr, geom, fev );
        }
        tensor( this_type const& expr, Geo_t const& geom )
        {
            this->initTensor( expr, geom );
        }

        template<typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update(geom);
        }
        void update( Geo_t const& geom )
        {
            M_tensorbase->update( geom );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            M_tensorbase->update( geom, face );
        }

        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i,c1,c2,q );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return M_tensorbase->evalq( q );
        }
        template<typename... TheArgsType>
        void initTensor( this_type const& expr, const TheArgsType&... theInitArgs )
            {
                auto const& dynamicViscosity = expr.dynamicViscosity();
                if ( dynamicViscosity.isNewtonianLaw() )
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,theInitArgs...) );
                else if ( dynamicViscosity.isPowerLaw() )
                    M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,theInitArgs...) );
                else if ( dynamicViscosity.isCarreauLaw() )
                    M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,theInitArgs...) );
                else if ( dynamicViscosity.isCarreauYasudaLaw() )
                    M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,theInitArgs...) );
                else
                    CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
            }
    private:
        tensorbase_ptrtype M_tensorbase;
    };

private:
    expr_grad_velocity_type M_exprGradVelocity;
    expr_id_pressure_type M_exprIdPressure;

    model_physic_fluid_type const& M_physicFluid;
    MaterialProperties const& M_matProps;
    uint16_type M_polynomialOrder;
    bool M_withPressureTerm;
    SymbolsExprType M_se;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */

template<class ExprGradVelocityType,class ExprIdPressureType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecNewtonianStressTensor( Expr<ExprGradVelocityType> const& grad_u, Expr<ExprIdPressureType> const& p,
                               ModelPhysicFluidType const& physicFluid,
                               MaterialProperties const& matProps,
                               bool withPressure,
                               SymbolsExprType const& se = symbols_expression_empty_t{},
                               uint16_type polyOrder = invalid_uint16_type_value )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,Expr<ExprIdPressureType>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,p,physicFluid,matProps,polyOrder,withPressure,se ) );
}

template<class ExprGradVelocityType,class ElementVelocityType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecNewtonianViscousStressTensorJacobian( Expr<ExprGradVelocityType> const& grad_u, ElementVelocityType const& /*u*/,
                                              ModelPhysicFluidType const& physicFluid,
                                              MaterialProperties const& matProps,
                                              SymbolsExprType const& se = symbols_expression_empty_t{},
                                              uint16_type polyOrder = invalid_uint16_type_value )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,typename unwrap_ptr_t<ElementVelocityType>::functionspace_type::reference_element_type,Expr<Cst<double>>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,/*p*/cst(0.),physicFluid,matProps,polyOrder,false,se ) );
}

template<class ExprGradVelocityType, class ModelPhysicFluidType,typename SymbolsExprType = symbols_expression_empty_t >
inline
auto
fluidMecViscosity( Expr<ExprGradVelocityType> const& grad_u,
                   ModelPhysicFluidType const& physicFluid,
                   MaterialProperties const& matProps,
                   SymbolsExprType const& se = symbols_expression_empty_t{},
                   uint16_type polyOrder = invalid_uint16_type_value )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,Expr<Cst<double>>,ModelPhysicFluidType,SymbolsExprType,mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,/*p*/cst(0.),physicFluid,matProps,polyOrder,false,se ) );
}


} // namespace FeelModels
} // namespace Feel
#endif /* FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H */
