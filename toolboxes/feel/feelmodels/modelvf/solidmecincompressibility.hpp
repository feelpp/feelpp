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
   \file solidmecresidual.hpp
   \author Vincent Chabannes
   \date 2013-11-19
 */
#ifndef __FEELPP_MODELS_VF_SOLIDMECINCOMPRESSIBILITY_H
#define __FEELPP_MODELS_VF_SOLIDMECINCOMPRESSIBILITY_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{
#if 0
struct ExprApplySolidMecPresFormType
{
    enum ExprApplyType { EVAL=0,JACOBIAN_TRIAL_DISP=1,JACOBIAN_TRIAL_PRES=2 };
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationMultiplierBase : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,typename ExprType::my_shape_type,typename ExprType::value_type >
{
    typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,typename ExprType::my_shape_type,typename ExprType::value_type > super_type;
    typedef ExprType expr_type;
    typedef typename expr_type::spec_expr_type SpecificExprType;
    typedef typename expr_type::value_type value_type;
    typedef typename expr_type::geoelement_type geoelement_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::shape_type shape;
    using ret_type = typename super_type::ret_type;
    
    // fe disp context
    typedef typename expr_type::fe_disp_type::PreCompute pc_disp_type;
    typedef std::shared_ptr<pc_disp_type> pc_disp_ptrtype;
    typedef typename expr_type::fe_disp_type::template Context<expr_type::context_disp, typename expr_type::fe_disp_type,
                                                               gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_disp_type;
    typedef std::shared_ptr<ctx_disp_type> ctx_disp_ptrtype;

    // fe pressure context
    typedef typename expr_type::fe_pressure_type::PreCompute pc_pressure_type;
    typedef std::shared_ptr<pc_pressure_type> pc_pressure_ptrtype;
    typedef typename expr_type::fe_pressure_type::template Context<expr_type::context_pressure, typename expr_type::fe_pressure_type,
                                                                   gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_pressure_type;
    typedef std::shared_ptr<ctx_pressure_type> ctx_pressure_ptrtype;

    tensorSolidMecPressureFormulationMultiplierBase( expr_type const& expr,
                                                     Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( geom, fev, feu ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(), this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(), this->gmc(),(pc_pressure_ptrtype const&)M_pcPressure ) ),
        M_locRes( boost::extents[ this->gmc()->xRefs().size2()] ),
        M_locGradDisplacement( boost::extents[ this->gmc()->xRefs().size2()] ),
        //M_locMatrixGradDisplacement( boost::extents[ this->gmc()->xRefs().size2()] ),
        M_locIdPressure( boost::extents[ this->gmc()->xRefs().size2()] )
        {}

    tensorSolidMecPressureFormulationMultiplierBase( expr_type const& expr,
                                                     Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( geom, fev ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),this->gmc(),(pc_pressure_ptrtype const&)M_pcPressure ) ),
        M_locRes( expr.disp().gradExtents(*this->gmc()) ),
        M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) ),
        //M_locMatrixGradDisplacement( boost::extents[ this->gmc()->xRefs().size2()] ),
        M_locIdPressure( expr.pressure().idExtents(*this->gmc()) )
        {}

    tensorSolidMecPressureFormulationMultiplierBase( expr_type const& expr, Geo_t const& geom )
        :
        super_type( geom ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_pcPressure( new pc_pressure_type( expr.pressure().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxPressure( new ctx_pressure_type( expr.pressure().functionSpace()->fe(),this->gmc(),(pc_pressure_ptrtype const&)M_pcPressure ) ),
        M_locRes( expr.disp().gradExtents(*this->gmc()) ),
        M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) ),
        //M_locMatrixGradDisplacement( boost::extents[ this->gmc()->xRefs().size2()] ),
        M_locIdPressure( expr.pressure().idExtents(*this->gmc()) )
        {}

    void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) { update(geom); }
    void update( Geo_t const& geom, Basis_i_t const& /*fev*/ ) { update(geom); }
    virtual void update( Geo_t const& geom ) = 0;

    using super_type::evalijq; // fix clang warning

    virtual
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const { return super_type::evalijq(i,j,q); }

    virtual
    value_type evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const = 0;

    //virtual
    //value_type evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const = 0;

    value_type evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        DCHECK( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL ) << "only for EVAL expression";
        return M_locRes[q]( c1,c2 );
    }
    ret_type
    evaliq( uint16_type i, uint16_type q ) const
    {
        DCHECK( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL ) << "only for EVAL expression";
        return ret_type(M_locRes[q].data());
    }

    value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return M_locRes[q]( c1,c2 );
    }
    ret_type
    evalq( uint16_type q ) const
    {
        return ret_type(M_locRes[q].data());
    }


    //private:
protected :
    expr_type const& M_expr;
    pc_disp_ptrtype M_pcDisp;
    ctx_disp_ptrtype M_ctxDisp;
    pc_pressure_ptrtype M_pcPressure;
    ctx_pressure_ptrtype M_ctxPressure;

    typename super_type::array_shape_type M_locRes;
    typename super_type::array_tensor2_type M_locGradDisplacement;
    //typename super_type::array_matrix_tensor2_type M_locMatrixGradDisplacement;
    typename super_type::array_scalar_type M_locIdPressure;
};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationMultiplierClassic : public tensorSolidMecPressureFormulationMultiplierBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
{
    typedef tensorSolidMecPressureFormulationMultiplierBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    typedef typename super_type::expr_type expr_type;
    typedef typename expr_type::spec_expr_type SpecificExprType;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;
    
    tensorSolidMecPressureFormulationMultiplierClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( expr, geom, fev, feu)
        {}
    tensorSolidMecPressureFormulationMultiplierClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( expr, geom, fev )
        {}
    tensorSolidMecPressureFormulationMultiplierClassic( expr_type const& expr, Geo_t const& geom )
        :
        super_type( expr, geom )
        {}

    void update( Geo_t const& geom )
    {
        this->setGmc( geom );
        std::fill( this->M_locRes.data(), this->M_locRes.data()+this->M_locRes.num_elements(), super_type::matrix_shape_type::Zero() );
        std::fill( this->M_locGradDisplacement.data(), this->M_locGradDisplacement.data()+this->M_locGradDisplacement.num_elements(), this->M_zeroLocTensor2/*super_type::loc_tensor2_type::Zero()*/ );
        std::fill( this->M_locIdPressure.data(), this->M_locIdPressure.data()+this->M_locIdPressure.num_elements(), this->M_zeroLocScalar/*super_type::loc_scalar_type::Zero()*/ );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            this->M_pcDisp->update( this->gmc()->pc()->nodes() );
        this->M_ctxDisp->update( this->gmc(),  (typename super_type::pc_disp_ptrtype const&) this->M_pcDisp );
        this->M_expr.disp().grad( *this->M_ctxDisp, this->M_locGradDisplacement );
        if ( SpecificExprType::value != ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES )
        {
            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
                this->M_pcPressure->update( this->gmc()->pc()->nodes() );
            this->M_ctxPressure->update( this->gmc(),  (typename super_type::pc_pressure_ptrtype const&) this->M_pcPressure );
            this->M_expr.pressure().id( *this->M_ctxPressure, this->M_locIdPressure );
        }
#if 0
        for ( uint16_type q=0;q<this->gmc()->nPoints();++q )
        {
            Eigen::Map< Eigen::Matrix<typename super_type::value_type,Eigen::Dynamic,Eigen::Dynamic/*,Eigen::ColMajor*/ > >
                gradDisplacementEval( this->M_locGradDisplacement[q].data(), super_type::shape_tensor2::M,super_type::shape_tensor2::N);
            this->M_locMatrixGradDisplacement[q] = gradDisplacementEval;
        }
#endif
        updateImpl( mpl::int_<SpecificExprType::value>() );
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i,j,q, mpl::int_<SpecificExprType::value>() );
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return evalijq( i,j,c1,c2,q, mpl::int_<SpecificExprType::value>() );
    }
    /*value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return evaliq( i,c1,c2,q, mpl::int_<SpecificExprType::value>() );
     }*/

private:


    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL> )
    {
        updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>(), mpl::int_<super_type::gmc_type::nDim>() );
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>, mpl::int_<2> /*Dim*/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // compute : Fv*J*Cv^{-1}*idv(p) = Fv*J*(Fv^T*Fv)^{-1}*idv(p) = J*Fv^{-T}*idv(p) (with J=det(F))
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1);
            const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1);
            typename super_type::matrix_shape_type/*loc_tensor2_type*/ & theLocRes = this->M_locRes[q];
            theLocRes(0,0) = 1+du2vdy;
            theLocRes(0,1) = -du2vdx;
            theLocRes(1,0) = -du1vdy;
            theLocRes(1,1) = 1+du1vdx;
            if ( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL )
            {
                const value_type idPressureEval = this->M_locIdPressure[q](0,0);
                theLocRes *= /*-*/idPressureEval;
            }
        }
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>, mpl::int_<3> /*Dim*/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // compute : Fv*J*Cv^{-1}*idv(p) = Fv*J*(Fv^T*Fv)^{-1}*idv(p) = J*Fv^{-T}*idv(p) (with J=det(F))
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1), du1vdz = gradDisplacementEval(0,2);
            const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1), du2vdz = gradDisplacementEval(1,2);
            const value_type du3vdx = gradDisplacementEval(2,0), du3vdy = gradDisplacementEval(2,1), du3vdz = gradDisplacementEval(2,2);
            typename super_type::matrix_shape_type/*loc_tensor2_type*/ & theLocRes = this->M_locRes[q];
            theLocRes(0,0) = (1+du2vdy)*(1+du3vdz) - du2vdz*du3vdy;
            theLocRes(0,1) = du2vdz*du3vdx - du2vdx*(1+du3vdz);
            theLocRes(0,2) = du2vdx*du3vdy - (1+du2vdy)*du3vdx;
            theLocRes(1,0) = du1vdz*du3vdy - du1vdy*(1+du3vdz);
            theLocRes(1,1) = (1+du1vdx)*(1+du3vdz) - du1vdz*du3vdx;
            theLocRes(1,2) = du1vdy*du3vdx - (1+du1vdx)*du3vdy;
            theLocRes(2,0) = du1vdy*du2vdz - du1vdz*(1+du2vdy);
            theLocRes(2,1) = du1vdz*du2vdx - (1+du1vdx)*du2vdz;
            theLocRes(2,2) = (1+du1vdx)*(1+du2vdy) - du1vdy*du2vdx;
            if ( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL )
            {
                const value_type idPressureEval = this->M_locIdPressure[q](0,0);
                theLocRes *= /*-*/idPressureEval;
            }
        }
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) {}
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> )
    {
        updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>() );
    }
    //---------------------------------------------------------//
    /*value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::EVAL> ) const
    {
        return this->evalq( c1,c2,q );
    }
    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) const
    {
        CHECK( false ) << "not allow"; return 0;
    }
    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> ) const
    {
        CHECK( false ) << "not allow"; return 0;
     }*/
    //---------------------------------------------------------//
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::EVAL> ) const
    {
        return super_type::evalijq( i,j,q); // not allow
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) const
    {
        return evalijq( i,j,q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP>(),mpl::int_<super_type::gmc_type::nDim>() );
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> ) const
    {
        // compute -idt(p)*(Fv*Cv^{-1}) and Fv*Cv^{-1} is computed in update
        const value_type idTrialPressure = this->fecTrial()->id( j, 0, 0, q );
        matrix_shape_type & thelocMat = this->locMatrixShape();
        thelocMat = /*-*/idTrialPressure*this->M_locRes[q];
        return ret_type(thelocMat.data());
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ,mpl::int_<2> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 = gradTrial( 0, 0, 0 ), dF12 = gradTrial( 0, 1, 0 );
        const value_type dF21 = gradTrial( 1, 0, 0 ), dF22 = gradTrial( 1, 1, 0 );
        matrix_shape_type & thelocRes = this->locMatrixShape();
        thelocRes(0,0) =  dF22;
        thelocRes(0,1) = -dF21;
        thelocRes(1,0) = -dF12;
        thelocRes(1,1) =  dF11;
        const value_type idPressureEval = this->M_locIdPressure[q](0,0);
        thelocRes *= /*-*/idPressureEval;
        return ret_type(this->locMatrixShape().data());
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ,mpl::int_<3> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 = gradTrial( 0, 0, 0 ), dF12 = gradTrial( 0, 1, 0 ), dF13 = gradTrial( 0, 2, 0 );
        const value_type dF21 = gradTrial( 1, 0, 0 ), dF22 = gradTrial( 1, 1, 0 ), dF23 = gradTrial( 1, 2, 0 );
        const value_type dF31 = gradTrial( 2, 0, 0 ), dF32 = gradTrial( 2, 1, 0 ), dF33 = gradTrial( 2, 2, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type Fv11 = 1+gradDisplacementEval(0,0), Fv12 =   gradDisplacementEval(0,1), Fv13 =   gradDisplacementEval(0,2);
        const value_type Fv21 =   gradDisplacementEval(1,0), Fv22 = 1+gradDisplacementEval(1,1), Fv23 =   gradDisplacementEval(1,2);
        const value_type Fv31 =   gradDisplacementEval(2,0), Fv32 =   gradDisplacementEval(2,1), Fv33 = 1+gradDisplacementEval(2,2);

        matrix_shape_type & thelocRes = this->locMatrixShape();
        thelocRes(0,0) = dF22+dF33+dF22*Fv33+Fv22*dF33-dF23*Fv32-Fv23*dF32;
        thelocRes(0,1) = dF23*Fv31+Fv23*dF31-dF21-dF21*Fv33-Fv21*dF33;
        thelocRes(0,2) = dF21*Fv32+Fv21*dF32-dF31-dF31*Fv22-Fv31*dF22;
        thelocRes(1,0) = dF13*Fv32+Fv13*dF32-dF12-dF12*Fv33-Fv12*dF33;
        thelocRes(1,1) = dF11+dF33+dF11*Fv33+Fv11*dF33-dF13*Fv31-Fv13*dF31;
        thelocRes(1,2) = dF12*Fv31+Fv12*dF31-dF32-dF32*Fv11-Fv32*dF11;
        thelocRes(2,0) = dF12*Fv23+Fv12*dF23-dF13-dF13*Fv22-Fv13*dF22;
        thelocRes(2,1) = dF13*Fv21+Fv13*dF21-dF23-dF23*Fv11-Fv23*dF11;
        thelocRes(2,2) = dF11+dF22+dF11*Fv22+Fv11*dF22-dF12*Fv21-Fv12*dF21;

        const value_type idPressureEval = this->M_locIdPressure[q](0,0);
        thelocRes *= /*-*/idPressureEval;
        return ret_type(thelocRes.data());
    }

    value_type
    evalijq( uint16_type i,uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::EVAL> ) const
    {
        return this->evalq( c1,c2,q );
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) const
    {
        return evalijq( i,j,c1,c2,q,mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP>(),mpl::int_<super_type::gmc_type::nDim>() );
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> ) const
    {
        //CHECK( false ) << "TODO";return 0;
        // compute -idt(p)*(Fv*Cv^{-1}) and Fv*Cv^{-1} is computed in update
        const value_type idTrialPressure = this->fecTrial()->id( j, 0, 0, q );
        return  /*-*/idTrialPressure*this->M_locRes[q](c1,c2);
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ,mpl::int_<2> /*Dim*/ ) const
    {
        //CHECK( false ) << "TO UPGRADE WITH PRESSURE";
        const value_type idPressureEval = this->M_locIdPressure[q](0,0);
        if (c1==0)
        {
            if (c2==0)
                return idPressureEval*this->fecTrial()->grad( j, 1, 1, q );
            else
                return -idPressureEval*this->fecTrial()->grad( j, 1, 0, q );
        }
        else
        {
            if (c2==0)
                return -idPressureEval*this->fecTrial()->grad( j, 0, 1, q );
            else
                return idPressureEval*this->fecTrial()->grad( j, 0, 0, q );
        }
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ,mpl::int_<3> /*Dim*/ ) const
    {
        //CHECK( false ) << "TO UPGRADE WITH PRESSURE";
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 = gradTrial( 0, 0, 0 ), dF12 = gradTrial( 0, 1, 0 ), dF13 = gradTrial( 0, 2, 0 );
        const value_type dF21 = gradTrial( 1, 0, 0 ), dF22 = gradTrial( 1, 1, 0 ), dF23 = gradTrial( 1, 2, 0 );
        const value_type dF31 = gradTrial( 2, 0, 0 ), dF32 = gradTrial( 2, 1, 0 ), dF33 = gradTrial( 2, 2, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type Fv11 = 1+gradDisplacementEval(0,0), Fv12 = gradDisplacementEval(0,1), Fv13 = gradDisplacementEval(0,2);
        const value_type Fv21 = gradDisplacementEval(1,0), Fv22 = 1+gradDisplacementEval(1,1), Fv23 = gradDisplacementEval(1,2);
        const value_type Fv31 = gradDisplacementEval(2,0), Fv32 = gradDisplacementEval(2,1), Fv33 = 1+gradDisplacementEval(2,2);

        const value_type idPressureEval = this->M_locIdPressure[q](0,0);

        if (c1==0)
        {
            if (c2==0)
                return idPressureEval*(dF22+dF33+dF22*Fv33+Fv22*dF33-dF23*Fv32-Fv23*dF32);
            else if (c2==1)
                return idPressureEval*(dF23*Fv31+Fv23*dF31-dF21-dF21*Fv33-Fv21*dF33);
            else
                return idPressureEval*(dF21*Fv32+Fv21*dF32-dF31-dF31*Fv22-Fv31*dF22);
        }
        else if (c1==1)
        {
            if (c2==0)
                return idPressureEval*(dF13*Fv32+Fv13*dF32-dF12-dF12*Fv33-Fv12*dF33);
            else if (c2==1)
                return idPressureEval*(dF11+dF33+dF11*Fv33+Fv11*dF33-dF13*Fv31-Fv13*dF31);
            else
                return idPressureEval*(dF12*Fv31+Fv12*dF31-dF32-dF32*Fv11-Fv32*dF11);
        }
        else
        {
            if (c2==0)
                return idPressureEval*(dF12*Fv23+Fv12*dF23-dF13-dF13*Fv22-Fv13*dF22);
            else if (c2==1)
                return idPressureEval*(dF13*Fv21+Fv13*dF21-dF23-dF23*Fv11-Fv23*dF11);
            else
                return idPressureEval*(dF11+dF22+dF11*Fv22+Fv11*dF22-dF12*Fv21-Fv12*dF21);
        }
    }

};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff : public tensorSolidMecPressureFormulationMultiplierBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
{
    typedef tensorSolidMecPressureFormulationMultiplierBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    typedef typename super_type::expr_type expr_type;
    typedef typename expr_type::spec_expr_type SpecificExprType;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;
    
    tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( expr, geom, fev, feu)
        {}
    tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( expr, geom, fev )
        {}
    tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff( expr_type const& expr, Geo_t const& geom )
        :
        super_type( expr, geom )
        {}

    void update( Geo_t const& geom )
    {
        this->setGmc( geom );
        std::fill( this->M_locRes.data(), this->M_locRes.data()+this->M_locRes.num_elements(), super_type::matrix_shape_type::Zero() );
        std::fill( this->M_locGradDisplacement.data(), this->M_locGradDisplacement.data()+this->M_locGradDisplacement.num_elements(), this->M_zeroLocTensor2/*super_type::loc_tensor2_type::Zero()*/ );
        std::fill( this->M_locIdPressure.data(), this->M_locIdPressure.data()+this->M_locIdPressure.num_elements(), this->M_zeroLocScalar/*super_type::loc_scalar_type::Zero()*/ );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            this->M_pcDisp->update( this->gmc()->pc()->nodes() );
        this->M_ctxDisp->update( this->gmc(),  (typename super_type::pc_disp_ptrtype const&) this->M_pcDisp );
        this->M_expr.disp().grad( *this->M_ctxDisp, this->M_locGradDisplacement );
        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            this->M_pcPressure->update( this->gmc()->pc()->nodes() );
        this->M_ctxPressure->update( this->gmc(),  (typename super_type::pc_pressure_ptrtype const&) this->M_pcPressure );
        this->M_expr.pressure().id( *this->M_ctxPressure, this->M_locIdPressure );
        updateImpl( mpl::int_<SpecificExprType::value>() );
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i,j,q, mpl::int_<SpecificExprType::value>() );
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        //CHECK( false ) << "TODO";
        //return 0;
        //LOG(WARNING) << "evalijq non optimized";
        return this->evalijq(i,j,q)(c1,c2);
    }
    /*value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return evaliq( i,c1,c2,q, mpl::int_<SpecificExprType::value>() );
     }*/

private:


    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL> )
    {
        updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>(), mpl::int_<super_type::gmc_type::nDim>() );
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>, mpl::int_<2> /*Dim*/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // compute : F*(-p*I)
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            typename super_type::matrix_shape_type/*loc_tensor2_type*/ & theLocRes = this->M_locRes[q];
            theLocRes(0,0) = 1+gradDisplacementEval(0,0);
            theLocRes(0,1) =   gradDisplacementEval(0,1);
            theLocRes(1,0) =   gradDisplacementEval(1,0);
            theLocRes(1,1) = 1+gradDisplacementEval(1,1);
            if ( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL )
            {
                const value_type idPressureEval = this->M_locIdPressure[q](0,0);
                theLocRes *= /*-*/idPressureEval;
            }
        }
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>, mpl::int_<3> /*Dim*/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // compute : F*(-p*I)
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            typename super_type::matrix_shape_type/*loc_tensor2_type*/ & theLocRes = this->M_locRes[q];
            theLocRes(0,0) = 1+gradDisplacementEval(0,0);
            theLocRes(0,1) =   gradDisplacementEval(0,1);
            theLocRes(0,2) =   gradDisplacementEval(0,2);
            theLocRes(1,0) =   gradDisplacementEval(1,0);
            theLocRes(1,1) = 1+gradDisplacementEval(1,1);
            theLocRes(1,2) =   gradDisplacementEval(1,2);
            theLocRes(2,0) =   gradDisplacementEval(2,0);
            theLocRes(2,1) =   gradDisplacementEval(2,1);
            theLocRes(2,2) = 1+gradDisplacementEval(2,2);
            if ( SpecificExprType::value == ExprApplySolidMecPresFormType::EVAL )
            {
                const value_type idPressureEval = this->M_locIdPressure[q](0,0);
                theLocRes *= /*-*/idPressureEval;
            }
        }
    }
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) {}
    void updateImpl( mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> )
    {
        updateImpl( mpl::int_<ExprApplySolidMecPresFormType::EVAL>() );
    }
    //---------------------------------------------------------//
    /*value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::EVAL> ) const
    {
        return this->evalq( c1,c2,q );
    }
    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) const
    {
        CHECK( false ) << "not allow"; return 0;
    }
    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> ) const
    {
        CHECK( false ) << "not allow"; return 0;
     }*/
    //---------------------------------------------------------//
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::EVAL> ) const
    {
        CHECK( false ) << "TODO";
        return ret_type(this->locMatrixShape().data());
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> ) const
    {
        // compute -idv(p)*dF and Fv is computed in update
        auto const& _gradTrial = this->fecTrial()->grad( j, q );
        Eigen::Map< const Eigen::Matrix<typename super_type::value_type,Eigen::Dynamic,Eigen::Dynamic/*,Eigen::ColMajor*/ > >
            gradTrial( _gradTrial.data(), super_type::shape_tensor2::M,super_type::shape_tensor2::N);
        const value_type idPressureEval = this->M_locIdPressure[q](0,0);
        matrix_shape_type & thelocRes = this->locMatrixShape();
        thelocRes = /*-*/idPressureEval*gradTrial;
        return ret_type(this->locMatrixShape().data());
    }
    ret_type
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> ) const
    {
        // compute -idt(p)*Fv and Fv is computed in update
        const value_type idTrialPressure = this->fecTrial()->id( j, 0, 0, q );
        matrix_shape_type & thelocMat = this->locMatrixShape();
        thelocMat = /*-*/idTrialPressure*this->M_locRes[q];
        return ret_type(thelocMat.data());
    }

};









template<typename ElementDispType, typename ElementPressureType, typename ModelPhysicSolidType,typename SpecificExprType>
class SolidMecPressureFormulationMultiplier
{
public:

    typedef SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,SpecificExprType> self_type;

    static const size_type context_disp = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context_pressure = vm::JACOBIAN;
    static const size_type context = context_disp;

    typedef ElementDispType element_disp_type;
    typedef ElementPressureType element_pressure_type;
    typedef SpecificExprType spec_expr_type;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_disp_type::functionspace_type functionspace_disp_type;
    typedef typename functionspace_disp_type::reference_element_type* fe_disp_ptrtype;
    typedef typename functionspace_disp_type::reference_element_type fe_disp_type;

    typedef typename element_pressure_type::functionspace_type functionspace_pressure_type;
    typedef typename functionspace_pressure_type::reference_element_type* fe_pressure_ptrtype;
    typedef typename functionspace_pressure_type::reference_element_type fe_pressure_type;


    typedef typename functionspace_disp_type::geoelement_type geoelement_type;
    typedef typename functionspace_disp_type::gm_type gm_type;
    typedef typename functionspace_disp_type::value_type value_type;
    typedef value_type evaluate_type;
    //static const uint16_type rank = fe_disp_type::rank;
    //static const uint16_type nComponents1 = fe_disp_type::nComponents1;
    //static const uint16_type nComponents2 = fe_disp_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type orderdisplacement = functionspace_disp_type::basis_type::nOrder;
    static const uint16_type orderpressure = functionspace_pressure_type::basis_type::nOrder;
    static const uint16_type nDim = functionspace_disp_type::nDim;
    //------------------------------------------------------------------------------//

    template<typename Func>
    struct HasTestFunction
    {
      static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result_disp = mpl::if_<  mpl::bool_< ( SpecificExprType::value == ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP ) >,
                                                   typename mpl::if_<boost::is_same<Func,fe_disp_type>,
                                                                     mpl::bool_<true>, mpl::bool_<false> >::type,
                                                   mpl::bool_<false> >::type::value;
        static const bool result_pressure = mpl::if_<  mpl::bool_< ( SpecificExprType::value == ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES ) >,
                                                       typename mpl::if_<boost::is_same<Func,fe_pressure_type>,
                                                                         mpl::bool_<true>, mpl::bool_<false> >::type,
                                                       mpl::bool_<false> >::type::value;
        static const bool result = (SpecificExprType::value == ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP)? result_disp : result_pressure;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    SolidMecPressureFormulationMultiplier( element_disp_type const& disp, element_pressure_type const& p, ModelPhysicSolidType const& physicSolidData )
        :
        M_disp( disp ),
        M_pressure( p ),
        M_physicSolidData( physicSolidData )
    {}
    SolidMecPressureFormulationMultiplier( SolidMecPressureFormulationMultiplier const & op ) = default;

    ~SolidMecPressureFormulationMultiplier()
    {}

    //! polynomial order
    uint16_type polynomialOrder() const { return (orderdisplacement-1)*(nDim-1)+orderpressure; }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    element_disp_type const& disp() const { return M_disp; }
    element_pressure_type const& pressure() const { return M_pressure; }
    ModelPhysicSolidType const& physicSolidData() const { return M_physicSolidData; }

    typedef Shape<nDim, Tensor2, false, false> my_shape_type;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        //typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,my_shape_type,value_type > tensorbase_type;
        typedef tensorSolidMecPressureFormulationMultiplierBase<Geo_t,Basis_i_t,Basis_j_t,self_type > tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        typedef typename tensorbase_type::value_type value_type;
        typedef typename tensorbase_type::shape_type shape;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename tensorbase_type::ret_type;
        struct is_zero { static const bool value = false; };

        tensor( self_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev,feu) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev,feu) );
        }
        tensor( self_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev) );
        }
        tensor( self_type const& expr, Geo_t const& geom )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationMultiplierClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom) );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensorbase->update(geom,fev,feu);
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensorbase->update(geom,fev);
        }
        void update( Geo_t const& geom )
        {
            M_tensorbase->update( geom );
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

    private:
        tensorbase_ptrtype M_tensorbase;
    };

private:
    element_disp_type const& M_disp;
    element_pressure_type const& M_pressure;
    ModelPhysicSolidType const& M_physicSolidData;
};

#endif

 
  
    
         
               
                    
               
         
    
  
 

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationConstraintBase : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,typename ExprType::my_shape_type,typename ExprType::value_type >
{
    typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,typename ExprType::my_shape_type,typename ExprType::value_type > super_type;
    typedef ExprType expr_type;
    typedef typename expr_type::value_type value_type;
    typedef typename expr_type::geoelement_type geoelement_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::shape_type shape;
    using ret_type = typename super_type::ret_type;
    // fe disp context
    typedef typename expr_type::fe_disp_type::PreCompute pc_disp_type;
    typedef std::shared_ptr<pc_disp_type> pc_disp_ptrtype;
    typedef typename expr_type::fe_disp_type::template Context<expr_type::context_disp, typename expr_type::fe_disp_type,
                                                               gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_disp_type;
    typedef std::shared_ptr<ctx_disp_type> ctx_disp_ptrtype;

    tensorSolidMecPressureFormulationConstraintBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( geom, fev, feu ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_locRes( expr.disp().gradExtents(*this->gmc()) ),
        M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) )
        {}
    tensorSolidMecPressureFormulationConstraintBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( geom, fev ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_locRes( expr.disp().gradExtents(*this->gmc()) ),
        M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) )
        {}
    tensorSolidMecPressureFormulationConstraintBase( expr_type const& expr, Geo_t const& geom )
        :
        super_type( geom ),
        M_expr( expr ),
        M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
        M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
        M_locRes( expr.disp().gradExtents(*this->gmc()) ),
        M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) )
        {}

    void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) override { this->update(geom); }
    void update( Geo_t const& geom, Basis_i_t const& /*fev*/ ) override { this->update(geom); }
    virtual void update( Geo_t const& geom ) override = 0;

    value_type evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const override
    {
        DCHECK( expr_type::is_applied_as_eval ) << "only for EVAL expression";
        return M_locRes[q]( c1,c2 );
    }
    ret_type
    evaliq( uint16_type i, uint16_type q ) const override
    {
        DCHECK( expr_type::is_applied_as_eval ) << "only for EVAL expression";
        return ret_type(M_locRes[q].data());
    }

    value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const override
    {
        return M_locRes[q]( c1,c2 );
    }
    ret_type
    evalq( uint16_type q ) const override
    {
        return ret_type(M_locRes[q].data());
    }


protected:
    expr_type const& M_expr;
    pc_disp_ptrtype M_pcDisp;
    ctx_disp_ptrtype M_ctxDisp;

    typename super_type::array_shape_type M_locRes;
    typename super_type::array_tensor2_type M_locGradDisplacement;
};
template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationConstraintClassic : public tensorSolidMecPressureFormulationConstraintBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
{
    typedef tensorSolidMecPressureFormulationConstraintBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    typedef typename super_type::expr_type expr_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    using ret_type = typename super_type::ret_type;
    tensorSolidMecPressureFormulationConstraintClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( expr, geom, fev, feu)
        {}
    tensorSolidMecPressureFormulationConstraintClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( expr, geom, fev )
        {}
    tensorSolidMecPressureFormulationConstraintClassic( expr_type const& expr, Geo_t const& geom )
        :
        super_type( expr, geom )
        {}

    void update( Geo_t const& geom ) override
    {
        this->setGmc( geom );
        std::fill( this->M_locRes.data(), this->M_locRes.data()+this->M_locRes.num_elements(), super_type::matrix_shape_type::Zero() );
        std::fill( this->M_locGradDisplacement.data(), this->M_locGradDisplacement.data()+this->M_locGradDisplacement.num_elements(), this->M_zeroLocTensor2/*super_type::loc_tensor2_type::Zero()*/ );
        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            this->M_pcDisp->update( this->gmc()->pc()->nodes() );
        this->M_ctxDisp->update( this->gmc(),  (typename super_type::pc_disp_ptrtype const&) this->M_pcDisp );
        this->M_expr.disp().grad( *this->M_ctxDisp, this->M_locGradDisplacement );

        if constexpr ( expr_type::is_applied_as_eval )
        {
            this->updateImpl();
        }
    }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
    {
        if constexpr ( expr_type::is_applied_as_eval )
            return this->evalq( c1,c2,q );
        else
            return evalijq( i,j,c1,c2,q,mpl::int_<super_type::gmc_type::nDim>() );
    }

private:

    //----------------------------------------------------------------------------------------//

    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            if constexpr ( super_type::gmc_type::nDim == 2 )
            {
                const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1);
                const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1);
                this->M_locRes[q](0,0) = du1vdx + du2vdy + du1vdx*du2vdy - du1vdy*du2vdx;
            }
            else
            {
                const value_type F11 = 1+gradDisplacementEval(0,0), F12 =   gradDisplacementEval(0,1), F13 =   gradDisplacementEval(0,2);
                const value_type F21 =   gradDisplacementEval(1,0), F22 = 1+gradDisplacementEval(1,1), F23 =   gradDisplacementEval(1,2);
                const value_type F31 =   gradDisplacementEval(2,0), F32 =   gradDisplacementEval(2,1), F33 = 1+gradDisplacementEval(2,2);
                const value_type detF = F11*(F22*F33-F23*F32) - F21*(F12*F33-F13*F32) + F31*(F12*F23 - F13*F22);
                this->M_locRes[q](0,0) = detF-1;
            }
        }
    }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type Fat11 = gradTrial( 0, 0, 0 ), Fat12 = gradTrial( 0, 1, 0 );
        const value_type Fat21 = gradTrial( 1, 0, 0 ), Fat22 = gradTrial( 1, 1, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type Fav11 = gradDisplacementEval(0,0), Fav12 = gradDisplacementEval(0,1);
        const value_type Fav21 = gradDisplacementEval(1,0), Fav22 = gradDisplacementEval(1,1);

        const value_type detJm1 = Fat22+Fat11+Fat11*Fav22+Fav11*Fat22 -Fav21*Fat12 - Fat21*Fav12;
        return detJm1;
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<3> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type Fat11 = gradTrial( 0, 0, 0 ), Fat12 = gradTrial( 0, 1, 0 ), Fat13 = gradTrial( 0, 2, 0 );
        const value_type Fat21 = gradTrial( 1, 0, 0 ), Fat22 = gradTrial( 1, 1, 0 ), Fat23 = gradTrial( 1, 2, 0 );
        const value_type Fat31 = gradTrial( 2, 0, 0 ), Fat32 = gradTrial( 2, 1, 0 ), Fat33 = gradTrial( 2, 2, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type Fav11 = gradDisplacementEval(0,0), Fav12 = gradDisplacementEval(0,1), Fav13 = gradDisplacementEval(0,2);
        const value_type Fav21 = gradDisplacementEval(1,0), Fav22 = gradDisplacementEval(1,1), Fav23 = gradDisplacementEval(1,2);
        const value_type Fav31 = gradDisplacementEval(2,0), Fav32 = gradDisplacementEval(2,1), Fav33 = gradDisplacementEval(2,2);

        const value_type detJm1 = Fat11 + Fat22 + Fat33 + (Fat22*Fav33 + Fav22*Fat33) + (Fat11*Fav22 + Fav11*Fat22) + (Fat11*Fav33 + Fav11*Fat33)
            + Fat11*Fav22*Fav33 + Fav11*Fat22*Fav33 + Fav11*Fav22*Fat33
            + Fat12*Fav23*Fav31 + Fav12*Fat23*Fav31 + Fav12*Fav23*Fat31
            + Fat13*Fav21*Fav32 + Fav13*Fat21*Fav32 + Fav13*Fav21*Fat32
            - (Fat23*Fav32 + Fav23*Fat32) - (Fat12*Fav21 + Fav12*Fat21) - (Fat13*Fav31 + Fav13*Fat31)
            - (Fat11*Fav23*Fav32 + Fav11*Fat23*Fav32 + Fav11*Fav23*Fat32)
            - (Fat12*Fav21*Fav33 + Fav12*Fat21*Fav33 + Fav12*Fav21*Fat33)
            - (Fat13*Fav22*Fav31 + Fav13*Fat22*Fav31 + Fav13*Fav22*Fat31) ;
        return detJm1;
    }

};

template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType >
struct tensorSolidMecPressureFormulationConstraintStVenantKirchhoff : public tensorSolidMecPressureFormulationConstraintBase<Geo_t,Basis_i_t,Basis_j_t,ExprType>
{
    typedef tensorSolidMecPressureFormulationConstraintBase<Geo_t,Basis_i_t,Basis_j_t,ExprType> super_type;
    typedef typename super_type::expr_type expr_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::matrix_shape_type matrix_shape_type;

    tensorSolidMecPressureFormulationConstraintStVenantKirchhoff( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        :
        super_type( expr, geom, fev, feu)
        {}
    tensorSolidMecPressureFormulationConstraintStVenantKirchhoff( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        :
        super_type( expr, geom, fev )
        {}
    tensorSolidMecPressureFormulationConstraintStVenantKirchhoff( expr_type const& expr, Geo_t const& geom )
        :
        super_type( expr, geom )
        {}

    void update( Geo_t const& geom ) override
    {
        this->setGmc( geom );
        std::fill( this->M_locRes.data(), this->M_locRes.data()+this->M_locRes.num_elements(), super_type::matrix_shape_type::Zero() );
        std::fill( this->M_locGradDisplacement.data(), this->M_locGradDisplacement.data()+this->M_locGradDisplacement.num_elements(), this->M_zeroLocTensor2/*super_type::loc_tensor2_type::Zero()*/ );
        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            this->M_pcDisp->update( this->gmc()->pc()->nodes() );
        this->M_ctxDisp->update( this->gmc(),  (typename super_type::pc_disp_ptrtype const&) this->M_pcDisp );
        this->M_expr.disp().grad( *this->M_ctxDisp, this->M_locGradDisplacement );

        if constexpr ( expr_type::is_applied_as_eval )
        {
            this->updateImpl();
        }
    }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const override
    {
        //return evalijq( i,j,c1,c2,q, mpl::int_<SpecificExprType::value>() );
        if constexpr ( expr_type::is_applied_as_eval )
            return this->evalq( c1,c2,q );
        else
            return evalijq( i,j,c1,c2,q,mpl::int_<super_type::gmc_type::nDim>() );
    }
    /*value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return evaliq( i,c1,c2,q, mpl::int_<SpecificExprType::value>() );
     }*/

private:

    //----------------------------------------------------------------------------------------//
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // trace E
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            if constexpr ( super_type::gmc_type::nDim == 2 )
            {
                const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1);
                const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1);
                const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2));
                const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2));
                const value_type E11 = du1vdx + subtraceE1;
                const value_type E22 = du2vdy + subtraceE2;
                this->M_locRes[q](0,0) = E11 + E22;
            }
            else
            {
                const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1), du1vdz = gradDisplacementEval(0,2);
                const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1), du2vdz = gradDisplacementEval(1,2);
                const value_type du3vdx = gradDisplacementEval(2,0), du3vdy = gradDisplacementEval(2,1), du3vdz = gradDisplacementEval(2,2);
                const value_type subtraceE1 = 0.5*(std::pow(du1vdx,2)+std::pow(du2vdx,2)+std::pow(du3vdx,2));
                const value_type subtraceE2 = 0.5*(std::pow(du1vdy,2)+std::pow(du2vdy,2)+std::pow(du3vdy,2));
                const value_type subtraceE3 = 0.5*(std::pow(du1vdz,2)+std::pow(du2vdz,2)+std::pow(du3vdz,2));
                const value_type E11 = du1vdx + subtraceE1;
                const value_type E22 = du2vdy + subtraceE2;
                const value_type E33 = du3vdz + subtraceE3;
                this->M_locRes[q](0,0) = E11+E22+E33;
            }
        }
    }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<2> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type du1tdx = gradTrial( 0, 0, 0 ), du1tdy = gradTrial( 0, 1, 0 );
        const value_type du2tdx = gradTrial( 1, 0, 0 ), du2tdy = gradTrial( 1, 1, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1);
        const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1);

        //auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
        const value_type dE11 = du1tdx + /*0.5**/(du1vdx*du1tdx + du2vdx*du2tdx);
        const value_type dE22 = du2tdy + du1vdy*du1tdy + du2vdy*du2tdy;
        const value_type tracedE = dE11 + dE22;
        return tracedE;
    }
    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<3> /*Dim*/ ) const
    {
        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type du1tdx = gradTrial( 0, 0, 0 ), du1tdy = gradTrial( 0, 1, 0 ), du1tdz=gradTrial( 0, 2, 0 );
        const value_type du2tdx = gradTrial( 1, 0, 0 ), du2tdy = gradTrial( 1, 1, 0 ), du2tdz=gradTrial( 1, 2, 0 );
        const value_type du3tdx = gradTrial( 2, 0, 0 ), du3tdy = gradTrial( 2, 1, 0 ), du3tdz=gradTrial( 2, 2, 0 );
        auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
        const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1), du1vdz = gradDisplacementEval(0,2);
        const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1), du2vdz = gradDisplacementEval(1,2);
        const value_type du3vdx = gradDisplacementEval(2,0), du3vdy = gradDisplacementEval(2,1), du3vdz = gradDisplacementEval(2,2);

        //auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
        const value_type dE11 = du1tdx + /*0.5**/(du1vdx*du1tdx + du2vdx*du2tdx + du3vdx*du3tdx);
        const value_type dE22 = du2tdy + du1vdy*du1tdy + du2vdy*du2tdy + du3vdy*du3tdy;
        const value_type dE33 = du3tdz + du1vdz*du1tdz + du2vdz*du2tdz + du3vdz*du3tdz;
        const value_type tracedE = dE11 + dE22 + dE33;
        return tracedE;
    }

};


template<typename ElementDispType, typename ModelPhysicSolidType, ExprApplyType ExprApplied>
class SolidMecPressureFormulationConstraint
{
public:

    typedef SolidMecPressureFormulationConstraint<ElementDispType,ModelPhysicSolidType,ExprApplied> self_type;

    static const size_type context_disp = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_disp;

    typedef ElementDispType element_disp_type;

    static constexpr bool is_applied_as_eval = ExprApplied == ExprApplyType::EVAL;
    static constexpr bool is_applied_as_jacobian = ExprApplied == ExprApplyType::JACOBIAN;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_disp_type::functionspace_type functionspace_disp_type;
    typedef typename functionspace_disp_type::reference_element_type* fe_disp_ptrtype;
    typedef typename functionspace_disp_type::reference_element_type fe_disp_type;

    typedef typename functionspace_disp_type::geoelement_type geoelement_type;
    typedef typename functionspace_disp_type::gm_type gm_type;
    typedef typename functionspace_disp_type::value_type value_type;
    typedef value_type evaluate_type;
    //static const uint16_type rank = fe_disp_type::rank;
    //static const uint16_type nComponents1 = fe_disp_type::nComponents1;
    //static const uint16_type nComponents2 = fe_disp_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type orderdisplacement = functionspace_disp_type::basis_type::nOrder;
    static const uint16_type nDim = functionspace_disp_type::nDim;
    //------------------------------------------------------------------------------//

    template<typename Func>
    struct HasTestFunction
    {
      static const bool result = false;
    };

    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = is_applied_as_jacobian && std::is_same_v<Func,fe_disp_type>;
    };

    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;


    SolidMecPressureFormulationConstraint( element_disp_type const& disp, ModelPhysicSolidType const& physicSolidData )
        :
        M_disp( disp ),
        M_physicSolidData( physicSolidData )
    {}

    SolidMecPressureFormulationConstraint( SolidMecPressureFormulationConstraint const & op ) = default;

    ~SolidMecPressureFormulationConstraint()
    {}

        //! polynomial order
    uint16_type polynomialOrder() const { return (orderdisplacement-1)*nDim; }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    element_disp_type const& disp() const { return M_disp; }
    ModelPhysicSolidType const& physicSolidData() const { return M_physicSolidData; }

    typedef Shape<nDim, Scalar, false, false> my_shape_type;

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef tensorSolidMecPressureFormulationConstraintBase<Geo_t,Basis_i_t,Basis_j_t,self_type > tensorbase_type;
        typedef std::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        typedef typename tensorbase_type::value_type value_type;
        typedef typename tensorbase_type::shape_type shape;
        typedef typename tensorbase_type::matrix_shape_type matrix_shape_type;

        struct is_zero { static const bool value = false; };

        tensor( self_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev,feu) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev,feu) );
        }
        tensor( self_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom,fev) );
        }
        tensor( self_type const& expr, Geo_t const& geom )
        {
            if ( expr.physicSolidData().materialModel() == "StVenantKirchhoff" )
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom) );
            else
                M_tensorbase.reset( new tensorSolidMecPressureFormulationConstraintClassic<Geo_t, Basis_i_t, Basis_j_t,self_type>(expr,geom) );
        }
        void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            M_tensorbase->update(geom,fev,feu);
        }
        void update( Geo_t const& geom, Basis_i_t const& fev )
        {
            M_tensorbase->update(geom,fev);
        }
        void update( Geo_t const& geom )
        {
            M_tensorbase->update( geom );
        }
#if 0
        typename tensorbase_type::matrix_shape_type const&
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i,j,q );
        }
#endif
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
    element_disp_type const& M_disp;
    ModelPhysicSolidType const& M_physicSolidData;
};














/**
 * keywords
 */
#if 0
template<typename ElementDispType,typename ElementPressureType,typename ModelPhysicSolidType>
inline
Expr< SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,mpl::int_<ExprApplySolidMecPresFormType::EVAL> > >
solidMecPressureFormulationMultiplier( ElementDispType const& v, ElementPressureType const& p,
                                       ModelPhysicSolidType const& physicSolidData )
{
    typedef SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,mpl::int_<ExprApplySolidMecPresFormType::EVAL> > myexpr_type;
    return Expr< myexpr_type >( myexpr_type( v,p,physicSolidData ) );
}

template<typename ElementDispType,typename ElementPressureType,typename ModelPhysicSolidType>
inline
Expr< SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,
                                            mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> > >
solidMecPressureFormulationMultiplierJacobianTrialDisp( ElementDispType const& v,ElementPressureType const& p,
                                                        ModelPhysicSolidType const& physicSolidData )
{
    typedef SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,
                                                  mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_DISP> > myexpr_type;
    return Expr< myexpr_type >( myexpr_type( v,p,physicSolidData ) );
}

template<typename ElementDispType,typename ElementPressureType,typename ModelPhysicSolidType>
inline
Expr< SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,
                                            mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> > >
solidMecPressureFormulationMultiplierJacobianTrialPressure( ElementDispType const& v,ElementPressureType const& p,
                                                            ModelPhysicSolidType const& physicSolidData )
{
    typedef SolidMecPressureFormulationMultiplier<ElementDispType,ElementPressureType,ModelPhysicSolidType,
                                                  mpl::int_<ExprApplySolidMecPresFormType::JACOBIAN_TRIAL_PRES> > myexpr_type;
    return Expr< myexpr_type >(  myexpr_type( v,p,physicSolidData ) );
}
#endif


template<typename ElementDispType,typename ModelPhysicSolidType>
inline
auto
solidMecPressureFormulationConstraint( ElementDispType const& v,
                                       ModelPhysicSolidType const& physicSolidData )
{
    typedef SolidMecPressureFormulationConstraint<unwrap_ptr_t<ElementDispType>,ModelPhysicSolidType,ExprApplyType::EVAL > myexpr_type;
    return Expr< myexpr_type >( myexpr_type( unwrap_ptr(v),physicSolidData ) );
}

template<typename ElementDispType,typename ModelPhysicSolidType>
inline
auto
solidMecPressureFormulationConstraintJacobian( ElementDispType const& v,
                                               ModelPhysicSolidType const& physicSolidData )
{
    typedef SolidMecPressureFormulationConstraint<unwrap_ptr_t<ElementDispType>,ModelPhysicSolidType,ExprApplyType::JACOBIAN > myexpr_type;
    return Expr< myexpr_type >(  myexpr_type( unwrap_ptr(v),physicSolidData ) );
}

} // namespace FeelModels
} // namespace Feel
#endif /* __SOLIDMECINCOMPRESSIBILITY_H */
