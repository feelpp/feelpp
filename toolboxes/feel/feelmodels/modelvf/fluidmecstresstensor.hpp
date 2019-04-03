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
//#include <feel/feelmodels/fluid/viscositymodeldescription.hpp>
//template<class SpaceType> class DynamicViscosityModel;

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

        typedef typename expr_type::fe_pressure_type fe_pressure_type;
        typedef typename expr_type::geoelement_type geoelement_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::gmc_type gmc_type;

        typedef typename expr_type::expr_grad_velocity_type::template tensor<Geo_t, Basis_i_t, Basis_j_t> tensor_expr_grad_velocity_type;

        // fe pressure context
        typedef typename fe_pressure_type::PreCompute pc_pressure_type;
        typedef std::shared_ptr<pc_pressure_type> pc_pressure_ptrtype;
        typedef typename fe_pressure_type::template Context<expr_type::context_pressure, fe_pressure_type, gm_type,geoelement_type,gmc_type::context> ctx_pressure_type;
        typedef std::shared_ptr<ctx_pressure_type> ctx_pressure_ptrtype;

        typedef typename expr_type::value_type value_type;

        typedef typename super_type::loc_tensor2_type loc_grad_type;
        typedef typename super_type::new_array_tensor2_type array_grad_type;
        typedef typename super_type::loc_scalar_type loc_scalar_type;
        typedef typename super_type::array_scalar_type array_scalar_type;

        typedef typename super_type::shape_type shape;
        typedef typename super_type::matrix_shape_type matrix_shape_type;
        typedef typename super_type::new_array_shape_type array_shape_type;

        tensorFluidStressTensorBase( expr_type const& expr,
                                     Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom,fev,feu ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( this->gmc()->nPoints() ),
            M_locGradVelocity( this->gmc()->nPoints() ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
        {}
        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom,fev ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( this->gmc()->nPoints() ),
            M_locGradVelocity( this->gmc()->nPoints() ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
        {}

        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_tensorExprGradVelocity( expr.expressionGradVelocity(),geom ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( this->gmc()->nPoints() ),
            M_locGradVelocity( this->gmc()->nPoints() ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
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

        void updateCommon( Geo_t const& geom, bool upGradVelocity = true )
        {
            this->setGmc( geom );

            if ( upGradVelocity )
            {
                M_tensorExprGradVelocity.update( geom );
                M_locGradVelocity.resize( this->gmc()->nPoints() );
                for ( int q=0;q< this->gmc()->nPoints();++q )
                {
                    if constexpr ( false /*expr_type::expr_grad_velocity_type::is_terminal*/ )
                    {
                        M_locGradVelocity[q] = M_tensorExprGradVelocity.evalq( q ); //not compile, need to investigate
                    }
                    else
                    {
                        typename super_type::loc_tensor2_type& locData = M_locGradVelocity[q];
                        locData = this->M_zeroLocTensor2;
                        for (int c1=0;c1<super_type::shape_tensor2::M;++c1 )
                            for (int c2=0;c2<super_type::shape_tensor2::N;++c2 )
                                locData(c1,c2) = M_tensorExprGradVelocity.evalq( c1,c2,q );
                    }
                }
            }

            if ( this->expr().withPressureTerm() )
            {
                bool faceEval = this->gmc()->faceId() != invalid_uint16_type_value;
                if ( faceEval )
                    M_pcPressure->update( this->gmc()->pc()->nodes() );
                M_ctxPressure->update( fusion::at_key<key_type>( geom ),  (pc_pressure_ptrtype const&) M_pcPressure );
                std::fill( M_locPressure.data(), M_locPressure.data()+M_locPressure.num_elements(), this->M_zeroLocScalar );
                this->expr().pressure().id( *M_ctxPressure, M_locPressure );
            }

            //std::fill( M_locRes.data(), M_locRes.data()+M_locRes.num_elements(), loc_res_type::Zero() );
            //update(mpl::int_<gmc_type::nDim>(), mpl::int_<SpecificExprType::value>() );
        }
        void updateCommon( Geo_t const& geom, uint16_type face, bool upGradVelocity = true  )
        {
            //CHECK(false) << "not implemented";
            this->setGmc( geom );

            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            {
                if ( this->expr().withPressureTerm() )
                    M_pcPressure->update( this->gmc()->pc()->nodes() );
            }

            if ( upGradVelocity )
            {
                M_tensorExprGradVelocity.update( geom, face );
                M_locGradVelocity.resize( this->gmc()->nPoints() );
                for ( int q=0;q< this->gmc()->nPoints();++q )
                {
                    if constexpr ( false /*expr_type::expr_grad_velocity_type::is_terminal*/ )
                    {
                        M_locGradVelocity[q] = M_tensorExprGradVelocity.evalq( q );
                    }
                    else
                    {
                        typename super_type::loc_tensor2_type& locData = M_locGradVelocity[q];
                        locData = this->M_zeroLocTensor2;
                        for (int c1=0;c1<super_type::shape_tensor2::M;++c1 )
                            for (int c2=0;c2<super_type::shape_tensor2::N;++c2 )
                                locData(c1,c2) = M_tensorExprGradVelocity.evalq( c1,c2,q );
                    }
                }
            }

            if ( this->expr().withPressureTerm() )
            {
                M_ctxPressure->update( fusion::at_key<key_type>( geom ),  (pc_pressure_ptrtype const&) M_pcPressure );
                std::fill( M_locPressure.data(), M_locPressure.data()+M_locPressure.num_elements(), this->M_zeroLocScalar );
                this->expr().pressure().id( *M_ctxPressure, M_locPressure );
            }

        }


        /*virtual*/
        Eigen::Matrix<value_type, shape::M, shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return this->M_locMatrixShape;//evalijq( i,j,q,mpl::int_<gmc_type::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL ||
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
        matrix_shape_type const&
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_locRes[q];
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return this->evalq( c1, c2, q, mpl::int_<expr_type::specific_expr_type::value>() );
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_locRes[q];
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

        pc_pressure_ptrtype M_pcPressure;
        ctx_pressure_ptrtype M_ctxPressure;

        array_shape_type M_locRes;
        array_grad_type M_locGradVelocity;
        array_scalar_type M_locPressure;

        //mutable matrix_shape_type M_locMatrixShape;
    };

    /**
     * Newtonian :
     *   sigma = -pI + 2 \mu D(u) with D(u) = 0.5*(\grad u + (\grad u)^T)
     */
    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorNewtonianOLD : public tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
    {
        typedef tensorFluidStressTensorBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        // fe muP0 context
        typedef typename expr_type::fe_muP0_type fe_muP0_type;
        static const size_type context_muP0 = expr_type::context_muP0;
        typedef typename fe_muP0_type::PreCompute pc_muP0_type;
        typedef std::shared_ptr<pc_muP0_type> pc_muP0_ptrtype;
        typedef typename fe_muP0_type::template Context<context_muP0, fe_muP0_type, gm_type,geoelement_type,gmc_type::context> ctx_muP0_type;
        typedef std::shared_ptr<ctx_muP0_type> ctx_muP0_ptrtype;

        tensorFluidStressTensorNewtonianOLD( expr_type const& expr,
                                             Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}
        tensorFluidStressTensorNewtonianOLD( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}
        tensorFluidStressTensorNewtonianOLD( expr_type const& expr, Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            this->updateImpl();
        }

        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            this->updateImpl();
        }

        void updateImpl()
        {
            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
                M_pcMuP0->update( this->gmc()->pc()->nodes() );
            M_ctxMuP0->update( this->gmc(),  (pc_muP0_ptrtype const&) M_pcMuP0 );
            std::fill( M_locMuP0.data(), M_locMuP0.data()+M_locMuP0.num_elements(), this->M_zeroLocScalar );
            this->expr().viscosityModelDesc().fieldMu().id( *M_ctxMuP0, M_locMuP0 );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type muEval = M_locMuP0[q](0,0);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
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
                const value_type muEval = M_locMuP0[q](0,0);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
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
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }
            }
        }
        using super_type::evalijq; // fix clang warning

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim>() );
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            const value_type muEval = M_locMuP0[q](0,0);
            //M_locMatrixShape( matrix_shape_type::Zero() )
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            const value_type muEval = M_locMuP0[q](0,0);
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(2,0) = muEval*( du3tdx + du1tdz );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;
            this->locMatrixShape()(2,1) = muEval*( du2tdz + du3tdy );
            this->locMatrixShape()(0,2) = this->locMatrixShape()(2,0);
            this->locMatrixShape()(1,2) = this->locMatrixShape()(2,1);
            this->locMatrixShape()(2,2) = 2*muEval*du3tdz;

            return this->locMatrixShape();
        }

    private :

        pc_muP0_ptrtype M_pcMuP0;
        ctx_muP0_ptrtype M_ctxMuP0;
        typename super_type::array_scalar_type M_locMuP0;
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
        typedef typename super_type::geoelement_type geoelement_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::value_type value_type;

        // fe muP0 context
        typedef typename expr_type::fe_muP0_type fe_muP0_type;
        static const size_type context_muP0 = expr_type::context_muP0;
        typedef typename fe_muP0_type::PreCompute pc_muP0_type;
        typedef std::shared_ptr<pc_muP0_type> pc_muP0_ptrtype;
        typedef typename fe_muP0_type::template Context<context_muP0, fe_muP0_type, gm_type,geoelement_type,gmc_type::context> ctx_muP0_type;
        typedef std::shared_ptr<ctx_muP0_type> ctx_muP0_ptrtype;

        tensorFluidStressTensorNewtonian( expr_type const& expr,
                                          Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_muExprTensor( this->expr().dynamicViscosity().newtonian().expr().evaluator( geom ) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_muExprTensor( this->expr().dynamicViscosity().newtonian().expr().evaluator( geom ) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_muExprTensor( this->expr().dynamicViscosity().newtonian().expr().evaluator( geom ) )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            M_muExprTensor.update( geom );
            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face, expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL );
            M_muExprTensor.update( geom, face );
            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                const value_type muEval = M_muExprTensor.evalq(0,0,q);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
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
                const value_type muEval = M_muExprTensor.evalq(0,0,q);
                if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL )
                {
                    auto const& gradVelocityEval = this->locGradVelocity(q);
                    const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
                    const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
                    const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
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
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }
            }
        }
        using super_type::evalijq; // fix clang warning

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0);
            const value_type muEval = M_muExprTensor.evalq(0,0,q);
            //M_locMatrixShape( matrix_shape_type::Zero() )
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0,0), du1tdy = gradTrial(0,1,0), du1tdz = gradTrial(0,2,0);
            const value_type du2tdx = gradTrial(1,0,0), du2tdy = gradTrial(1,1,0), du2tdz = gradTrial(1,2,0);
            const value_type du3tdx = gradTrial(2,0,0), du3tdy = gradTrial(2,1,0), du3tdz = gradTrial(2,2,0);
            const value_type muEval = M_muExprTensor.evalq(0,0,q);
            this->locMatrixShape()(0,0) = 2*muEval*du1tdx;
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy );
            this->locMatrixShape()(2,0) = muEval*( du3tdx + du1tdz );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*muEval*du2tdy;
            this->locMatrixShape()(2,1) = muEval*( du2tdz + du3tdy );
            this->locMatrixShape()(0,2) = this->locMatrixShape()(2,0);
            this->locMatrixShape()(1,2) = this->locMatrixShape()(2,1);
            this->locMatrixShape()(2,2) = 2*muEval*du3tdz;

            return this->locMatrixShape();
        }

    private :

        typename ModelExpressionScalar::expr_scalar_type::template tensor<Geo_t/*, Basis_i_t, Basis_j_t*/>  M_muExprTensor;
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

        tensorFluidStressTensorPowerLaw( expr_type const& expr,
                                         Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu,
                                         Feel::FeelModels::DynamicViscosityPowerLaw const& powerLawModel )
            :
            super_type( expr,geom,fev,feu ),
            M_powerLawModel( powerLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev,
                                         Feel::FeelModels::DynamicViscosityPowerLaw const& powerLawModel )
            :
            super_type( expr,geom,fev ),
            M_powerLawModel( powerLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr,Geo_t const& geom,
                                         Feel::FeelModels::DynamicViscosityPowerLaw const& powerLawModel )
            :
            super_type( expr,geom ),
            M_powerLawModel( powerLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_muIsInInterval( this->gmc()->xRefs().size2() )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            M_muIsInInterval.resize( this->locRes().size() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            M_muIsInInterval.resize( this->locRes().size() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type power_k_generic = M_powerLawModel.k();
            const value_type power_n_generic = (M_powerLawModel.n()-1.)/2.;
            const value_type muMin = M_powerLawModel.muMin();
            const value_type muMax = M_powerLawModel.muMax();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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
            const value_type power_k_generic = M_powerLawModel.k();
            const value_type power_n_generic = (M_powerLawModel.n()-1.)/2.;
            const value_type muMin = M_powerLawModel.muMin();
            const value_type muMax = M_powerLawModel.muMax();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
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

            return thelocRes;
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
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

            return thelocRes;
        }
    private :
        Feel::FeelModels::DynamicViscosityPowerLaw const& M_powerLawModel;
        typename super_type::array_value_type M_locEvalPrecompute;
        typename super_type::array_value_type M_locEvalMu;
        std::vector<bool> M_muIsInInterval;
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

        tensorFluidStressTensorCarreau( expr_type const& expr,
                                        Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu,
                                        Feel::FeelModels::DynamicViscosityCarreauLaw const& carreauLawModel )
            :
            super_type( expr,geom,fev,feu ),
            M_carreauLawModel( carreauLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev,
                                        Feel::FeelModels::DynamicViscosityCarreauLaw const& carreauLawModel )
            :
            super_type( expr,geom,fev ),
            M_carreauLawModel( carreauLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr,Geo_t const& geom,
                                        Feel::FeelModels::DynamicViscosityCarreauLaw const& carreauLawModel )
            :
            super_type( expr,geom ),
            M_carreauLawModel( carreauLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }
        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type mu_inf = M_carreauLawModel.muInf();
            const value_type mu_0 = M_carreauLawModel.mu0();
            const value_type carreau_lambda = M_carreauLawModel.lambda();
            const value_type carreau_n = M_carreauLawModel.n();
            const value_type carreau_lambda_pow2_time2 = math::pow(carreau_lambda,2.)*2.0;
            const value_type carreauValPower =  (carreau_n-1)/2.0;
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);
            const value_type part1_carreauLaw = ( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2;

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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
            const value_type mu_inf = M_carreauLawModel.muInf();
            const value_type mu_0 = M_carreauLawModel.mu0();
            const value_type carreau_lambda = M_carreauLawModel.lambda();
            const value_type carreau_n = M_carreauLawModel.n();
            const value_type carreau_lambda_pow2_time2 = math::pow(carreau_lambda,2.)*2.0;
            const value_type carreauValPower =  (carreau_n-1)/2.0;
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);
            const value_type part1_carreauLaw = ( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2;

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim>() );
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/ , mpl::int_<2> /**/ ) const
        {
            const value_type mu_inf = M_carreauLawModel.muInf();
            const value_type mu_0 = M_carreauLawModel.mu0();
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

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type mu_inf = M_carreauLawModel.muInf();
            const value_type mu_0 = M_carreauLawModel.mu0();
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

            return this->locMatrixShape();
        }
    private :
        Feel::FeelModels::DynamicViscosityCarreauLaw const& M_carreauLawModel;
        typename super_type::array_value_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_value_type M_locEvalMu;
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

        tensorFluidStressTensorCarreauYasuda( expr_type const& expr,
                                              Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu,
                                              Feel::FeelModels::DynamicViscosityCarreauYasudaLaw const& carreauYasudaLawModel )
            :
            super_type( expr,geom,fev,feu ),
            M_carreauYasudaLawModel( carreauYasudaLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev,
                                              Feel::FeelModels::DynamicViscosityCarreauYasudaLaw const& carreauYasudaLawModel )
            :
            super_type( expr,geom,fev ),
            M_carreauYasudaLawModel( carreauYasudaLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr,Geo_t const& geom,
                                              Feel::FeelModels::DynamicViscosityCarreauYasudaLaw const& carreauYasudaLawModel )
            :
            super_type( expr,geom ),
            M_carreauYasudaLawModel( carreauYasudaLawModel ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

        void update( Geo_t const& geom ) override
        {
            super_type::updateCommon( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void update( Geo_t const& geom, uint16_type face ) override
        {
            super_type::updateCommon( geom, face );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().size(), super_type::matrix_shape_type::Zero() );
            updateImpl( mpl::int_<super_type::shape_type::nDim>() );
        }

        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type mu_inf = M_carreauYasudaLawModel.muInf();
            const value_type mu_0 = M_carreauYasudaLawModel.mu0();
            const value_type carreauYasuda_lambda = M_carreauYasudaLawModel.lambda();
            const value_type carreauYasuda_n = M_carreauYasudaLawModel.n();
            const value_type carreauYasuda_a = M_carreauYasudaLawModel.a();
            const value_type carreauYasuda_lambda_pow_a = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type carreauYasudaValPower = carreauYasuda_a/2.;
            const value_type carreauYasudaValPower2 = (carreauYasuda_n-1)/carreauYasuda_a;
            const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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
            const value_type mu_inf = M_carreauYasudaLawModel.muInf();
            const value_type mu_0 = M_carreauYasudaLawModel.mu0();
            const value_type carreauYasuda_lambda = M_carreauYasudaLawModel.lambda();
            const value_type carreauYasuda_n = M_carreauYasudaLawModel.n();
            const value_type carreauYasuda_a = M_carreauYasudaLawModel.a();
            const value_type carreauYasuda_lambda_pow_a = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type carreauYasudaValPower = carreauYasuda_a/2.;
            const value_type carreauYasudaValPower2 = (carreauYasuda_n-1)/carreauYasuda_a;

            const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);

            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
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
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const override
        {
            //return evalijq( i,j,q,mpl::int_<expr_type::nRealDim>() );
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<super_type::shape_type::nDim/*expr_type::nRealDim*/>() );
        }
        using super_type::evalijq; // fix clang warning

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/, mpl::int_<2> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::false_ /**/ , mpl::int_<3> /**/ ) const
        {
            CHECK( false) << "not allow\n";
            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<2> /**/ ) const
        {
            const value_type mu_inf = M_carreauYasudaLawModel.muInf();
            const value_type mu_0 = M_carreauYasudaLawModel.mu0();
            /*const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
             const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();*/
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

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type mu_inf = M_carreauYasudaLawModel.muInf();
            const value_type mu_0 = M_carreauYasudaLawModel.mu0();
            /*const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
             const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();*/
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

            return this->locMatrixShape();
        }
    private :
        Feel::FeelModels::DynamicViscosityCarreauYasudaLaw const& M_carreauYasudaLawModel;
        typename super_type::array_value_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_value_type M_locEvalMu;
    };

/// \cond detail
/**
 * \class FluidMecStressTensorImpl
 * \brief det of a matrix
 *
 * @author Vincent Chabannes
 * @see
 */
template<typename ExprGradVelocityType, typename FiniteElementVelocityType, typename ElementPressureType, typename ViscosityModelType, typename SpecificExprType>
class FluidMecStressTensorImpl
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FluidMecStressTensorImpl<ExprGradVelocityType,FiniteElementVelocityType,ElementPressureType,ViscosityModelType,SpecificExprType> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_pressure = vm::JACOBIAN;
    static const size_type context_muP0 = vm::JACOBIAN;
    static const size_type context = context_velocity;

    typedef ExprGradVelocityType expr_grad_velocity_type;
    typedef ElementPressureType element_pressure_type;
    typedef ViscosityModelType viscosity_model_type;
    typedef typename viscosity_model_type::element_type element_muP0_type;
    typedef SpecificExprType specific_expr_type;
    //------------------------------------------------------------------------------//
    // pressure functionspace
    typedef typename element_pressure_type::functionspace_type functionspace_pressure_type;
    typedef typename functionspace_pressure_type::reference_element_type fe_pressure_type;
    //------------------------------------------------------------------------------//
    // muP0 functionspace
    typedef typename element_muP0_type::functionspace_type functionspace_muP0_type;
    typedef typename functionspace_muP0_type::reference_element_type fe_muP0_type;
    //------------------------------------------------------------------------------//
    // expression desc
    typedef typename functionspace_pressure_type/*functionspace_velocity_type*/::geoelement_type geoelement_type;
    //typedef typename functionspace_velocity_type::gm_type gm_type;
    typedef typename functionspace_pressure_type/*functionspace_velocity_type*/::value_type value_type;
    typedef value_type evaluate_type;

    static const bool is_terminal = true;

    static const uint16_type orderpressure = functionspace_pressure_type::basis_type::nOrder;

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

    FluidMecStressTensorImpl( expr_grad_velocity_type const& exprGradVelocity, element_pressure_type const& p,
                              viscosity_model_type const& viscosityModelDesc, std::string const& matName,
                              uint16_type polyOrder,
                              bool withPressure, bool useNewNewtonianVersion )
        :
        M_exprGradVelocity( exprGradVelocity ),
        M_pressure( boost::cref(p) ),
        M_viscosityModelDesc( viscosityModelDesc ),
        M_matName( matName ),
        M_polynomialOrder( polyOrder ),
        M_withPressureTerm( withPressure ),
        M_useNewNewtonianVersion( useNewNewtonianVersion )
    {}

    FluidMecStressTensorImpl( FluidMecStressTensorImpl const & op )
        :
        M_exprGradVelocity( op.M_exprGradVelocity ),
        M_pressure( op.M_pressure ),
        M_viscosityModelDesc( op.M_viscosityModelDesc ),
        M_matName( op.M_matName ),
        M_polynomialOrder( op.M_polynomialOrder ),
        M_withPressureTerm( op.M_withPressureTerm ),
        M_useNewNewtonianVersion( op.M_useNewNewtonianVersion )
    {}

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
                    return this->dynamicViscosity().newtonian().expr().polynomialOrder();
                else if ( SpecificExprType::value == FMSTExprApplyType::FM_ST_EVAL )
                    return std::max( orderGradVelocity, orderpressure );
                else
                    return orderGradVelocity;
            }
            else
                return 2*(orderGradVelocity+1);
        }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    expr_grad_velocity_type const& expressionGradVelocity() const { return M_exprGradVelocity; }

    element_pressure_type const& pressure() const { return M_pressure; }
    viscosity_model_type const& viscosityModelDesc() const { return M_viscosityModelDesc; }
    auto const& dynamicViscosity() const { return M_viscosityModelDesc.dynamicViscosity( M_matName ); }
    bool withPressureTerm() const { return M_withPressureTerm; }
    bool useNewNewtonianVersion() const { return M_useNewNewtonianVersion; }
    //@}


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename element_pressure_type/*element_velocity_type*/::value_type value_type;

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

        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            auto const& dynamicViscosity = expr.dynamicViscosity();
            if ( dynamicViscosity.isNewtonianLaw() )
            {
                if ( expr.useNewNewtonianVersion() )
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
                else
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonianOLD<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
            }
            else if ( dynamicViscosity.isPowerLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu,dynamicViscosity.powerLaw()) );
            else if ( dynamicViscosity.isWalburnSchneckLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu,dynamicViscosity.walburnSchneckLaw().powerLaw()) );
            else if ( dynamicViscosity.isCarreauLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu,dynamicViscosity.carreauLaw()) );
            else if ( dynamicViscosity.isCarreauYasudaLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu,dynamicViscosity.carreauYasudaLaw()) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            auto const& dynamicViscosity = expr.dynamicViscosity();
            if ( dynamicViscosity.isNewtonianLaw() )
            {
                if ( expr.useNewNewtonianVersion() )
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
                else
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonianOLD<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
            }
            else if ( dynamicViscosity.isPowerLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,dynamicViscosity.powerLaw()) );
            else if ( dynamicViscosity.isWalburnSchneckLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,dynamicViscosity.walburnSchneckLaw().powerLaw()) );
            else if ( dynamicViscosity.isCarreauLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,dynamicViscosity.carreauLaw()) );
            else if ( dynamicViscosity.isCarreauYasudaLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,dynamicViscosity.carreauYasudaLaw()) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
        }
        tensor( this_type const& expr, Geo_t const& geom )
        {
            auto const& dynamicViscosity = expr.dynamicViscosity();
            if ( dynamicViscosity.isNewtonianLaw() )
            {
                if ( expr.useNewNewtonianVersion() )
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
                else
                    M_tensorbase.reset( new tensorFluidStressTensorNewtonianOLD<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
            }
            else if ( dynamicViscosity.isPowerLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,dynamicViscosity.powerLaw()) );
            else if ( dynamicViscosity.isWalburnSchneckLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,dynamicViscosity.walburnSchneckLaw().powerLaw()) );
            else if ( dynamicViscosity.isCarreauLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,dynamicViscosity.carreauLaw()) );
            else if ( dynamicViscosity.isCarreauYasudaLaw() )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,dynamicViscosity.carreauYasudaLaw()) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< dynamicViscosity.lawName() <<"\n";
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

        Eigen::Matrix<value_type, shape::M, shape::N> const&
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
        matrix_shape_type const&
        evaliq( uint16_type i, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }
        matrix_shape_type const&
        evalq( uint16_type q ) const
        {
            return M_tensorbase->evalq( q );
        }

    private:
        tensorbase_ptrtype M_tensorbase;
    };

private:
    expr_grad_velocity_type M_exprGradVelocity;
    boost::reference_wrapper<const element_pressure_type> M_pressure;
    viscosity_model_type const& M_viscosityModelDesc;
    std::string M_matName;
    uint16_type M_polynomialOrder;
    bool M_withPressureTerm;
    bool M_useNewNewtonianVersion;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */

template<class ExprGradVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> > >
fluidMecNewtonianStressTensor( Expr<ExprGradVelocityType> const& grad_u, ElementPressureType const& p,
                               ViscosityModelType const& viscosityModelDesc, std::string const& matName, bool withPressure, uint16_type polyOrder = invalid_uint16_type_value, bool useNewNewtonianVersion = false )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,p,viscosityModelDesc,matName,polyOrder,withPressure,useNewNewtonianVersion ) );
}

template<class ExprGradVelocityType,class ElementVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,typename ElementVelocityType::functionspace_type::reference_element_type,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> > >
fluidMecNewtonianStressTensorJacobian( Expr<ExprGradVelocityType> const& grad_u, ElementVelocityType const& /*u*/, ElementPressureType const& p,
                                       ViscosityModelType const& viscosityModelDesc, std::string const& matName, uint16_type polyOrder = invalid_uint16_type_value, bool useNewNewtonianVersion = false )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,typename ElementVelocityType::functionspace_type::reference_element_type,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,p,viscosityModelDesc,matName,polyOrder,false,useNewNewtonianVersion ) );
}

template<class ExprGradVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL> > >
fluidMecViscosity( Expr<ExprGradVelocityType> const& grad_u, ElementPressureType const& p,
                   ViscosityModelType const& viscosityModelDesc, std::string const& matName, uint16_type polyOrder = invalid_uint16_type_value, bool useNewNewtonianVersion = false )
{
    typedef FluidMecStressTensorImpl<Expr<ExprGradVelocityType>,std::nullptr_t,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL> > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( grad_u,p,viscosityModelDesc,matName,polyOrder,false,useNewNewtonianVersion ) );
}


} // namespace FeelModels
} // namespace Feel
#endif /* FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H */
