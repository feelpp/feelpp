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
namespace vf
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


#if 0
            typedef typename mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                      shape_grad_type,
                                      shape_grad_type/*shape_scalar_type*/ >::type shape_type;
#else
            typedef typename mpl::if_<mpl::or_< boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                                boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> > >,
                                      shape_grad_type,
                                      shape_scalar_type >::type shape_type;
#endif
        };
    }

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t,typename ExprType>
    struct tensorFluidStressTensorBase : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                                                           typename detail::GetShapeExpr<ExprType::nRealDim, typename ExprType::specific_expr_type>::shape_type,
                                                           typename ExprType::value_type >
    {
        typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,
                           typename detail::GetShapeExpr<ExprType::nRealDim, typename ExprType::specific_expr_type>::shape_type,
                           typename ExprType::value_type> super_type;
    public :
        typedef ExprType expr_type;
        typedef typename super_type::key_type key_type;

        typedef typename expr_type::fe_velocity_type fe_velocity_type;
        typedef typename expr_type::fe_pressure_type fe_pressure_type;
        typedef typename expr_type::geoelement_type geoelement_type;
        typedef typename super_type::gm_type gm_type;
        typedef typename super_type::gmc_type gmc_type;
        // fe velocity context
        typedef typename fe_velocity_type::PreCompute pc_velocity_type;
        typedef boost::shared_ptr<pc_velocity_type> pc_velocity_ptrtype;
        typedef typename fe_velocity_type::template Context<expr_type::context_velocity, fe_velocity_type, gm_type,geoelement_type,gmc_type::context> ctx_velocity_type;
        typedef boost::shared_ptr<ctx_velocity_type> ctx_velocity_ptrtype;
        // fe pressure context
        typedef typename fe_pressure_type::PreCompute pc_pressure_type;
        typedef boost::shared_ptr<pc_pressure_type> pc_pressure_ptrtype;
        typedef typename fe_pressure_type::template Context<expr_type::context_pressure, fe_pressure_type, gm_type,geoelement_type,gmc_type::context> ctx_pressure_type;
        typedef boost::shared_ptr<ctx_pressure_type> ctx_pressure_ptrtype;

        typedef typename expr_type::value_type value_type;
        typedef Shape<gmc_type::nDim, Tensor2, false, false> shape_grad;
        typedef Eigen::Matrix<value_type,shape_grad::M,shape_grad::N> loc_grad_type;
        typedef boost::multi_array<loc_grad_type,1> array_grad_type;

        typedef Shape<gmc_type::nDim, Scalar, false, false> shape_scalar;
        typedef Eigen::Matrix<value_type,shape_scalar::M,shape_scalar::N> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type,1> array_scalar_type;

#if 0
        typedef typename mpl::if_<boost::is_same<typename expr_type::specific_expr_type,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                  shape_grad,
                                  shape_grad/*shape_scalar*/ >::type shape;
        typedef Eigen::Matrix<value_type,shape::M,shape::N> matrix_shape_type;
        typedef boost::multi_array<matrix_shape_type,1> array_shape_type;
#else
        //typedef typename detail::GetShapeExpr<ExprType::nRealDim, typename ExprType::specific_expr_type>::shape_type shape;
        //typedef Eigen::Matrix<value_type,shape::M,shape::N> matrix_shape_type;
        //typedef boost::multi_array<matrix_shape_type,1> array_shape_type;
        typedef typename super_type::shape_type shape;
        typedef typename super_type::matrix_shape_type matrix_shape_type;
        typedef typename super_type::array_shape_type array_shape_type;
#endif

#if 0
        typedef loc_grad_type loc_res_type;
        typedef array_grad_type array_res_type;
#else

        typedef matrix_shape_type loc_res_type;
        typedef array_shape_type array_res_type;
#endif

        tensorFluidStressTensorBase( expr_type const& expr,
                                     Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom,fev,feu ),
            M_expr( expr ),
            M_pcVelocity( new pc_velocity_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxVelocity( new ctx_velocity_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_velocity_ptrtype const&)M_pcVelocity ) ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradVelocity( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_locMatrixShape( matrix_shape_type::Zero() )
        {}
        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom,fev ),
            M_expr( expr ),
            M_pcVelocity( new pc_velocity_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxVelocity( new ctx_velocity_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_velocity_ptrtype const&)M_pcVelocity ) ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradVelocity( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_locMatrixShape( matrix_shape_type::Zero() )
        {}

        tensorFluidStressTensorBase( expr_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_pcVelocity( new pc_velocity_type( expr.velocity().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxVelocity( new ctx_velocity_type( expr.velocity().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_velocity_ptrtype const&)M_pcVelocity ) ),
            M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), fusion::at_key<key_type>( geom )->xRefs() ) ),
            M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),fusion::at_key<key_type>( geom ),(pc_pressure_ptrtype const&)M_pcPressure ) ),
            M_locRes( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locGradVelocity( expr.velocity().gradExtents(*fusion::at_key<key_type>( geom )) ),
            M_locPressure( expr.pressure().idExtents(*fusion::at_key<key_type>( geom )) )
            //M_locMatrixShape( matrix_shape_type::Zero() )
        {}


        expr_type const& expr() const { return M_expr; }

        array_grad_type const& locGradVelocity() const { return M_locGradVelocity; }
        loc_grad_type const& locGradVelocity( uint16_type q ) const { return M_locGradVelocity[q]; }
        array_scalar_type const& locPressure() const { return M_locPressure; }
        loc_scalar_type const& locPressure( uint16_type q ) const { return M_locPressure[q]; }

        array_res_type & locRes() { return M_locRes; }
        array_res_type const& locRes() const { return M_locRes; }
        loc_res_type/*loc_grad_type*/ & locRes( uint16_type q ) { return M_locRes[q]; }
        loc_res_type/*loc_grad_type*/ const& locRes( uint16_type q ) const { return M_locRes[q]; }

        //matrix_shape_type /*const*/& locMatrixShape() const { return M_locMatrixShape; }

        /*virtual*/ void update( Geo_t const& geom )
        {
            M_ctxVelocity->update( fusion::at_key<key_type>( geom ),  (pc_velocity_ptrtype const&) M_pcVelocity );
            std::fill( M_locGradVelocity.data(), M_locGradVelocity.data()+M_locGradVelocity.num_elements(), loc_grad_type::Zero() );
            this->expr().velocity().grad( *M_ctxVelocity, M_locGradVelocity );

            if ( this->expr().withPressureTerm() )
            {
                M_ctxPressure->update( fusion::at_key<key_type>( geom ),  (pc_pressure_ptrtype const&) M_pcPressure );
                std::fill( M_locPressure.data(), M_locPressure.data()+M_locPressure.num_elements(), loc_scalar_type::Zero() );
                this->expr().pressure().id( *M_ctxPressure, M_locPressure );
            }

            //std::fill( M_locRes.data(), M_locRes.data()+M_locRes.num_elements(), loc_res_type::Zero() );
            //update(mpl::int_<gmc_type::nDim>(), mpl::int_<SpecificExprType::value>() );
        }
        void update( Geo_t const& geom, uint16_type face )
        {
            //CHECK(false) << "not implemented";
            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            {
                M_pcVelocity->update( this->gmc()->pc()->nodes() );
                if ( this->expr().withPressureTerm() )
                    M_pcPressure->update( this->gmc()->pc()->nodes() );
            }
            this->update( geom );
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
            return evalq( c1,c2,q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evalq( c1,c2,q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_locRes[q]( c1,c2 );
        }

    private :
        expr_type const& M_expr;

        pc_velocity_ptrtype M_pcVelocity;
        ctx_velocity_ptrtype M_ctxVelocity;
        pc_pressure_ptrtype M_pcPressure;
        ctx_pressure_ptrtype M_ctxPressure;

        array_res_type M_locRes;
        array_grad_type M_locGradVelocity;
        array_scalar_type M_locPressure;

        //mutable matrix_shape_type M_locMatrixShape;
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
        typedef boost::shared_ptr<pc_muP0_type> pc_muP0_ptrtype;
        typedef typename fe_muP0_type::template Context<context_muP0, fe_muP0_type, gm_type,geoelement_type,gmc_type::context> ctx_muP0_type;
        typedef boost::shared_ptr<ctx_muP0_type> ctx_muP0_ptrtype;

        tensorFluidStressTensorNewtonian( expr_type const& expr,
                                          Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}
        tensorFluidStressTensorNewtonian( expr_type const& expr, Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_pcMuP0( new pc_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxMuP0( new ctx_muP0_type( this->expr().viscosityModelDesc().fieldMu().functionSpace()->fe(),this->gmc(),(pc_muP0_ptrtype const&)M_pcMuP0 ) ),
            M_locMuP0( this->expr().viscosityModelDesc().fieldMu().idExtents(*this->gmc()) )
        {}

        void update( Geo_t const& geom )
        {
            if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_EVAL ||
                 expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                super_type::update( geom );

            M_ctxMuP0->update( this->gmc(),  (pc_muP0_ptrtype const&) M_pcMuP0 );
            std::fill( M_locMuP0.data(), M_locMuP0.data()+M_locMuP0.num_elements(), /*typename*/ super_type::loc_scalar_type::Zero() );
            this->expr().viscosityModelDesc().fieldMu().id( *M_ctxMuP0, M_locMuP0 );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().num_elements(), super_type::loc_res_type::Zero() );
            updateImpl( mpl::int_<gmc_type::nDim>()/*, expr_type::specific_expr_type()*/ );
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
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<expr_type::nRealDim>() );
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
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1);
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
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1), du1tdz = gradTrial(0,2);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1), du2tdz = gradTrial(1,2);
            const value_type du3tdx = gradTrial(2,0), du3tdy = gradTrial(2,1), du3tdz = gradTrial(2,2);
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

        tensorFluidStressTensorPowerLaw( expr_type const& expr,
                                         Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorPowerLaw( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

        void update( Geo_t const& geom )
        {
            super_type::update( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().num_elements(), super_type::loc_res_type::Zero() );
            updateImpl( mpl::int_<expr_type::nRealDim>() );
        }
        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type power_k_generic = this->expr().viscosityModelDesc().powerLaw_k_generic();
            const value_type power_n_generic = this->expr().viscosityModelDesc().powerLaw_n_generic();
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                //auto const mu_powerlaw = power_k_generic*pow( 2.0*inner(defv,defv) /*+chiInv*/ , cst( power_n_generic ) )/**chiSup*/;
                auto const& gradVelocityEval = this->locGradVelocity(q);
                const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
                const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
                const value_type DxD = math::pow(du1vdx,2) + 0.5*math::pow(du2vdx+du1vdy,2)  + math::pow(du2vdy,2);
                const value_type gammapoint2v = 2*DxD;
                const value_type muEval = power_k_generic*math::pow( gammapoint2v , power_n_generic );
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
                    M_locEvalMu[q](0,0) = muEval;
                    M_locEvalPrecompute[q](0,0) = power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
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
            const value_type power_k_generic = this->expr().viscosityModelDesc().powerLaw_k_generic();
            const value_type power_n_generic = this->expr().viscosityModelDesc().powerLaw_n_generic();
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
                const value_type muEval = power_k_generic*math::pow( gammapoint2v , power_n_generic );

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
                    M_locEvalMu[q](0,0) = muEval;
                    M_locEvalPrecompute[q](0,0) = power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }

        using super_type::evalijq; // fix clang warning

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<expr_type::nRealDim>() );
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
            const value_type power_k_generic = this->expr().viscosityModelDesc().powerLaw_k_generic();
            const value_type power_n_generic = this->expr().viscosityModelDesc().powerLaw_n_generic();
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q](0,0);
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
            //const value_type mut = gammapoint2t*power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);

            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy ) + mut*( du2vdx + du1vdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy );

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type power_k_generic = this->expr().viscosityModelDesc().powerLaw_k_generic();
            const value_type power_n_generic = this->expr().viscosityModelDesc().powerLaw_n_generic();
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1), du1tdz = gradTrial(0,2);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1), du2tdz = gradTrial(1,2);
            const value_type du3tdx = gradTrial(2,0), du3tdy = gradTrial(2,1), du3tdz = gradTrial(2,2);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q](0,0);
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);

            const value_type DxDt01 = du1tdy+du2tdx;
            const value_type DxDv01 = du1vdy+du2vdx;
            const value_type DxDt02 = du1tdz+du3tdx;
            const value_type DxDv02 = du1vdz+du3vdx;
            const value_type DxDt12 = du2tdz+du3tdy;
            const value_type DxDv12 = du2vdz+du3vdy;
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy + du3tdz*du3vdz) +
                DxDt01*DxDv01 + DxDt02*DxDv02 + DxDt12*DxDv12;
            //const value_type mut = gammapoint2t*power_k_generic*power_n_generic*math::pow( gammapoint2v , power_n_generic-1.0 );
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);

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
        typename super_type::array_scalar_type M_locEvalPrecompute;
        typename super_type::array_scalar_type M_locEvalMu;
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
                                         Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreau( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

        void update( Geo_t const& geom )
        {
            super_type::update( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().num_elements(), super_type::loc_res_type::Zero() );
            updateImpl( mpl::int_<expr_type::nRealDim>() );
        }
        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
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
                    M_locEvalMu[q](0,0) = muEval;
                    const value_type gammapoint2v = 2*DxD;
                    M_locEvalPrecompute[q](0,0) = part1_carreauLaw*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
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
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
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
                    M_locEvalMu[q](0,0) = muEval;
                    //const value_type mut = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                    const value_type gammapoint2v = 2*DxD;
                    M_locEvalPrecompute[q](0,0) = part1_carreauLaw*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        using super_type::evalijq; // fix clang warning

        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<expr_type::nRealDim>() );
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
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            /*const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);*/
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q](0,0);
            //const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
            //const value_type mut = gammapoint2t*( (carreau_n-1)/2.0 )*( mu_0 - mu_inf )*carreau_lambda2*math::pow( 1. + carreau_lambda2*gammapoint2v, (carreau_n-3)/2.0 );
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);

            this->locMatrixShape()(0,0) = 2*(muEval*du1tdx + mut*du1vdx);
            this->locMatrixShape()(1,0) = muEval*( du2tdx + du1tdy ) + mut*( du2vdx + du1vdy );
            this->locMatrixShape()(0,1) = this->locMatrixShape()(1,0);
            this->locMatrixShape()(1,1) = 2*(muEval*du2tdy + mut*du2vdy );

            return this->locMatrixShape();
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::true_ /**/, mpl::int_<3> /**/ ) const
        {
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            /*const value_type carreau_lambda = this->expr().viscosityModelDesc().carreau_lambda();
            const value_type carreau_n = this->expr().viscosityModelDesc().carreau_n();
            const value_type carreau_lambda2 = math::pow(carreau_lambda,2);*/
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1), du1tdz = gradTrial(0,2);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1), du2tdz = gradTrial(1,2);
            const value_type du3tdx = gradTrial(2,0), du3tdy = gradTrial(2,1), du3tdz = gradTrial(2,2);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q](0,0);
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
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);

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
        typename super_type::array_scalar_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_scalar_type M_locEvalMu;
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
                                         Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( expr,geom,fev,feu ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( expr,geom,fev ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}
        tensorFluidStressTensorCarreauYasuda( expr_type const& expr,Geo_t const& geom )
            :
            super_type( expr,geom ),
            M_locEvalPrecompute( boost::extents[ this->gmc()->xRefs().size2()] ),
            M_locEvalMu( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

        void update( Geo_t const& geom )
        {
            super_type::update( geom );

            std::fill( this->locRes().data(), this->locRes().data()+this->locRes().num_elements(), super_type::loc_res_type::Zero() );
            updateImpl( mpl::int_<expr_type::nRealDim>() );
        }
        void updateImpl( mpl::int_<2> /**/ )
        {
            const bool withPressure = this->expr().withPressureTerm();
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
            const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();
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
                    M_locEvalMu[q](0,0) = muEval;
                    const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
                    const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
                    //const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                    M_locEvalPrecompute[q](0,0) = part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
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
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
            const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();
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
                    M_locEvalMu[q](0,0) = muEval;
                    const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
                    const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
                    //const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
                    M_locEvalPrecompute[q](0,0) = part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );

                }
                else if ( expr_type::specific_expr_type::value == FMSTExprApplyType::FM_VISCOSITY_EVAL )
                {
                    this->locRes(q)(0,0) = muEval;
                }

            }
        }
        Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            //return evalijq( i,j,q,mpl::int_<expr_type::nRealDim>() );
            return evalijq( i,j,q,
                            mpl::bool_< expr_type::specific_expr_type::value == FMSTExprApplyType::FM_ST_JACOBIAN >(),
                            mpl::int_<expr_type::nRealDim>() );
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
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
            const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1);
            const value_type muEval = M_locEvalMu[q](0,0);
            const value_type gammapoint2t = 2*(du1tdx*du1vdx + du2tdy*du2vdy) + (du1tdy+du2tdx)*(du1vdy+du2vdx);
#if 0
            const value_type gammapoint2v = 2*M_locEvalDxD[q](0,0);
            const value_type carreauYasuda_lambdaA = math::pow(carreauYasuda_lambda,carreauYasuda_a);
            const value_type part1_carreauYasudaLaw = ( (carreauYasuda_n-1)/carreauYasuda_a )*( mu_0 - mu_inf);
            const value_type part2_carreauYasudaLaw = 1. + carreauYasuda_lambdaA*math::pow( gammapoint2v, carreauYasuda_a/2.);
            const value_type part3_carreauYasudaLaw = 0.5*carreauYasuda_a*carreauYasuda_lambdaA*math::pow( gammapoint2v, (carreauYasuda_a-2.)/2.);
            const value_type mut = gammapoint2t*part1_carreauYasudaLaw*math::pow( part2_carreauYasudaLaw, (carreauYasuda_n-3)/carreauYasuda_a );
#else
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);
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
            const value_type mu_inf = this->expr().viscosityModelDesc().mu_inf();
            const value_type mu_0 = this->expr().viscosityModelDesc().mu_0();
            const value_type carreauYasuda_lambda = this->expr().viscosityModelDesc().carreauYasuda_lambda();
            const value_type carreauYasuda_n = this->expr().viscosityModelDesc().carreauYasuda_n();
            const value_type carreauYasuda_a = this->expr().viscosityModelDesc().carreauYasuda_a();
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type du1tdx = gradTrial(0,0), du1tdy = gradTrial(0,1), du1tdz = gradTrial(0,2);
            const value_type du2tdx = gradTrial(1,0), du2tdy = gradTrial(1,1), du2tdz = gradTrial(1,2);
            const value_type du3tdx = gradTrial(2,0), du3tdy = gradTrial(2,1), du3tdz = gradTrial(2,2);
            auto const& gradVelocityEval = this->locGradVelocity(q);
            const value_type du1vdx = gradVelocityEval(0,0), du1vdy = gradVelocityEval(0,1), du1vdz = gradVelocityEval(0,2);
            const value_type du2vdx = gradVelocityEval(1,0), du2vdy = gradVelocityEval(1,1), du2vdz = gradVelocityEval(1,2);
            const value_type du3vdx = gradVelocityEval(2,0), du3vdy = gradVelocityEval(2,1), du3vdz = gradVelocityEval(2,2);
            const value_type muEval = M_locEvalMu[q](0,0);

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
            const value_type mut = gammapoint2t*M_locEvalPrecompute[q](0,0);
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
        typename super_type::array_scalar_type M_locEvalPrecompute;//M_locEvalDxD;
        typename super_type::array_scalar_type M_locEvalMu;
    };

/// \cond detail
/**
 * \class FluidMecStressTensorImpl
 * \brief det of a matrix
 *
 * @author Vincent Chabannes
 * @see
 */
template<typename ElementVelocityType, typename ElementPressureType, typename ViscosityModelType, typename SpecificExprType, int QuadOrder>
class FluidMecStressTensorImpl
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,SpecificExprType,QuadOrder> this_type;

    static const size_type context_velocity = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_pressure = vm::JACOBIAN;
    static const size_type context_muP0 = vm::JACOBIAN;
    static const size_type context = context_velocity;

    typedef ElementVelocityType element_velocity_type;
    typedef ElementPressureType element_pressure_type;
    typedef ViscosityModelType viscosity_model_type;
    typedef typename viscosity_model_type::element_type element_muP0_type;
    typedef SpecificExprType specific_expr_type;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_velocity_type::functionspace_type functionspace_velocity_type;
    typedef typename functionspace_velocity_type::reference_element_type* fe_velocity_ptrtype;
    typedef typename functionspace_velocity_type::reference_element_type fe_velocity_type;
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
    typedef typename functionspace_velocity_type::geoelement_type geoelement_type;
    //typedef typename functionspace_velocity_type::gm_type gm_type;
    typedef typename functionspace_velocity_type::value_type value_type;
    typedef value_type evaluate_type;
    static const uint16_type nDim = fe_velocity_type::nDim;
    static const uint16_type nRealDim = fe_velocity_type::nRealDim;
    static const uint16_type rank = fe_velocity_type::rank;
    static const uint16_type nComponents1 = fe_velocity_type::nComponents1;
    static const uint16_type nComponents2 = fe_velocity_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type ordervelocity = functionspace_velocity_type::basis_type::nOrder;
    static const uint16_type orderpressure = functionspace_pressure_type::basis_type::nOrder;

    //QuadOrder
    static const uint16_type imorderAuto = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                                typename mpl::max< mpl::int_<ordervelocity-1>, mpl::int_<orderpressure> >::type,
                                                mpl::int_<ordervelocity-1> >::type::value;
    static const uint16_type imorder = (QuadOrder>=0)?QuadOrder:imorderAuto;

    /*static const uint16_type imorder = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                                typename mpl::max< mpl::int_<ordervelocity-1>, mpl::int_<orderpressure> >::type,
                                                mpl::int_<ordervelocity-1> >::type::value*2;*/
    static const bool imIsPoly = true;



    template<typename Func>
    struct HasTestFunction
    {
#if 0
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                            mpl::bool_<false>,
                                            typename mpl::if_<boost::is_same<Func,fe_velocity_type>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false> >::type >::type::value;
#else
        static const bool result = false;
#endif
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN> >,
                                            typename mpl::if_<boost::is_same<Func,fe_velocity_type>,
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

    FluidMecStressTensorImpl( element_velocity_type const & u, element_pressure_type const& p,
                              viscosity_model_type const& viscosityModelDesc, bool withPressure )
        :
        M_velocity( boost::cref(u) ),
        M_pressure( boost::cref(p) ),
        M_viscosityModelDesc( viscosityModelDesc ),
        M_withPressureTerm( withPressure )
    {}

    FluidMecStressTensorImpl( FluidMecStressTensorImpl const & op )
        :
        M_velocity( op.M_velocity ),
        M_pressure( op.M_pressure ),
        M_viscosityModelDesc( op.M_viscosityModelDesc ),
        M_withPressureTerm( op.M_withPressureTerm )
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

    element_velocity_type const& velocity() const { return M_velocity; }
    element_pressure_type const& pressure() const { return M_pressure; }
    viscosity_model_type const& viscosityModelDesc() const { return M_viscosityModelDesc; }
    bool withPressureTerm() const { return M_withPressureTerm; }
    //@}





    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename element_velocity_type::value_type value_type;
#if 1
        // geomap context
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<0> >,
                mpl::identity<vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef boost::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;
#endif
        //----------------------------------------------------------------------------------------------------//

        /*typedef Shape<gmc_type::nDim, Vectorial, false, false> shape_id;
        typedef Eigen::Matrix<value_type,shape_id::M,shape_id::N> loc_id_type;
        typedef boost::multi_array<loc_id_type,1> array_id_type;*/

        typedef Shape<gmc_type::nDim, Tensor2, false, false> shape_grad;
        typedef Eigen::Matrix<value_type,shape_grad::M,shape_grad::N> loc_grad_type;
        typedef boost::multi_array<loc_grad_type,1> array_grad_type;

        typedef Shape<gmc_type::nDim, Scalar, false, false> shape_scalar;
        typedef Eigen::Matrix<value_type,shape_scalar::M,shape_scalar::N> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type,1> array_scalar_type;

#if 0
        typedef typename mpl::if_<boost::is_same<SpecificExprType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL> >,
                                  shape_grad,
                                  shape_grad/*shape_scalar*/ >::type shape;
#else
        typedef typename detail::GetShapeExpr<this_type::nRealDim, SpecificExprType>::shape_type shape;
#endif

        typedef Eigen::Matrix<value_type,shape::M,shape::N> matrix_shape_type;
        //typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> locAssembly_LinearForm_type;

        //----------------------------------------------------------------------------------------------------//

        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t,shape,value_type> tensorbase_type;
        typedef boost::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        //----------------------------------------------------------------------------------------------------//
        struct is_zero
        {
            static const bool value = false;
        };

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "newtonian" )
                M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "power_law" || expr.viscosityModelDesc().dynamicViscosityLaw() == "walburn-schneck_law" )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau-yasuda_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev,feu) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< expr.viscosityModelDesc().dynamicViscosityLaw() <<"\n";
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "newtonian" )
                M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "power_law" || expr.viscosityModelDesc().dynamicViscosityLaw() == "walburn-schneck_law" )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau-yasuda_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom,fev) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< expr.viscosityModelDesc().dynamicViscosityLaw() <<"\n";
        }
        tensor( this_type const& expr, Geo_t const& geom )
        {
            if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "newtonian" )
                M_tensorbase.reset( new tensorFluidStressTensorNewtonian<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "power_law" || expr.viscosityModelDesc().dynamicViscosityLaw() == "walburn-schneck_law" )
                M_tensorbase.reset( new tensorFluidStressTensorPowerLaw<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreau<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
            else if ( expr.viscosityModelDesc().dynamicViscosityLaw() == "carreau-yasuda_law" )
                M_tensorbase.reset( new tensorFluidStressTensorCarreauYasuda<Geo_t, Basis_i_t, Basis_j_t,this_type>(expr,geom) );
            else
                CHECK ( false ) << "invalid viscosity model : "<< expr.viscosityModelDesc().dynamicViscosityLaw() <<"\n";
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

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1,c2,q );
        }

    private:
        tensorbase_ptrtype M_tensorbase;
    };

private:
    boost::reference_wrapper<const element_velocity_type> M_velocity;
    boost::reference_wrapper<const element_pressure_type> M_pressure;
    viscosity_model_type const& M_viscosityModelDesc;
    bool M_withPressureTerm;
};
/// \endcond

/**
 * \brief det of the expression tensor
 */

template<int QuadOrder=-1,class ElementVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL>,QuadOrder > >
fluidMecNewtonianStressTensor( ElementVelocityType const& u, ElementPressureType const& p,
                               ViscosityModelType const& viscosityModelDesc, bool withPressure )
{
    typedef FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_EVAL>, QuadOrder > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( u,p,viscosityModelDesc, withPressure ) );
}

template<int QuadOrder=-1,class ElementVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN>,QuadOrder > >
fluidMecNewtonianStressTensorJacobian( ElementVelocityType const& u, ElementPressureType const& p,
                                       ViscosityModelType const& viscosityModelDesc, bool withPressure )
{
    typedef FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_ST_JACOBIAN>, QuadOrder > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( u,p,viscosityModelDesc, withPressure ) );
}

template<int QuadOrder=-1,class ElementVelocityType,class ElementPressureType, class ViscosityModelType >
inline
Expr< FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL>,QuadOrder > >
fluidMecViscosity( ElementVelocityType const& u, ElementPressureType const& p,
                   ViscosityModelType const& viscosityModelDesc )
{
    typedef FluidMecStressTensorImpl<ElementVelocityType,ElementPressureType,ViscosityModelType,mpl::int_<FMSTExprApplyType::FM_VISCOSITY_EVAL>, QuadOrder > fmstresstensor_t;
    return Expr< fmstresstensor_t >(  fmstresstensor_t( u,p,viscosityModelDesc, false ) );
}



} // namespace FeelModels
} // namespace vf
} // namespace Feel
#endif /* FEELPP_MODELS_VF_FLUIDMEC_STRESSTENSOR_H */
