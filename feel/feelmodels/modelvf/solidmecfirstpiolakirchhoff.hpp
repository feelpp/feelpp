/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 2014-08-19

  Copyright (C) 2014 Universite Joseph Fourier (Grenoble I)

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
   \file solidmecfirstpiolakirchhoff.hpp
   \author Vincent Chabannes
   \date 2014-08-19
 */
#ifndef __FEELPP_MODELS_VF_SOLIDMEC_FIRSTPIOLAKIRCHHOFF_H
#define __FEELPP_MODELS_VF_SOLIDMEC_FIRSTPIOLAKIRCHHOFF_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>
#include <feel/feelmodels/solid/mechanicalpropertiesdescription.hpp>

namespace Feel
{
namespace vf
{
namespace FeelModels
{

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic;
template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985;

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffCompressibleVolumicPartBase : public tensorBase<Geo_t, Basis_i_t, Basis_j_t, typename ExprType::shape_type, typename ExprType::value_type>
{
    typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t, typename ExprType::shape_type, typename ExprType::value_type> super_type;
    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> self_type;
    typedef ExprType expr_type;
    typedef typename expr_type::geoelement_type geoelement_type;
    typedef typename expr_type::value_type value_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::gmc_type gmc_type;

    typedef typename expr_type::fe_mechprop_scalar_type fe_mechprop_scalar_type;
    static const size_type context_mechprop_scalar = expr_type::context_mechprop_scalar;
    typedef typename fe_mechprop_scalar_type::template Context<context_mechprop_scalar, fe_mechprop_scalar_type, gm_type, geoelement_type, gmc_type::context> ctx_mechprop_scalar_type;

    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic<Geo_t, Basis_i_t, Basis_j_t, ExprType> tensor_volumic_part_classic_type;
    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985<Geo_t, Basis_i_t, Basis_j_t, ExprType> tensor_volumic_part_simo1985_type;

    tensorFirstPiolaKirchhoffCompressibleVolumicPartBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( geom, fev, feu ),
          M_expr( expr ),
          M_locEvalBulkModulus( M_expr.mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locRes( expr.displacement().gradExtents( *this->gmc() ) )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( geom, fev ),
          M_expr( expr ),
          M_locEvalBulkModulus( M_expr.mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locRes( M_expr.displacement().gradExtents( *this->gmc() ) )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartBase( expr_type const& expr, Geo_t const& geom )
        : super_type( geom ),
          M_expr( expr ),
          M_locEvalBulkModulus( M_expr.mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locRes( M_expr.displacement().gradExtents( *this->gmc() ) )
    {
    }

    void update( Geo_t const& geom ) { CHECK( false ) << "TODO"; }
    void update( Geo_t const& geom, uint16_type face ) { CHECK( false ) << "TODO"; }
    virtual void updateImpl( ctx_mechprop_scalar_type const& ctxMechProp, typename super_type::array_tensor2_type const& locGradEval )
    {
        std::fill( M_locEvalBulkModulus.data(), M_locEvalBulkModulus.data() + M_locEvalBulkModulus.num_elements(), super_type::loc_scalar_type::Zero() );
        M_expr.mechanicalPropertiesDesc().fieldBulkModulus().id( ctxMechProp, M_locEvalBulkModulus );
    }

    using super_type::evalijq; // fix clang warning
    virtual
        typename super_type::matrix_shape_type const&
        evalijq( uint16_type i, uint16_type j, uint16_type q,
                 const value_type detFv, const value_type ddetF,
                 typename super_type::loc_tensor2_type const& dFmt ) const = 0;

    expr_type const& expr() const { return M_expr; }
    typename super_type::array_shape_type const& locRes() const { return M_locRes; }
    typename super_type::matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }

    static boost::shared_ptr<self_type>
    New( std::string const& name, expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
    {
        boost::shared_ptr<self_type> res;
        if ( name == "classic" )
            res.reset( new tensor_volumic_part_classic_type( expr, geom, fev, feu ) );
        else if ( name == "simo1985" )
            res.reset( new tensor_volumic_part_simo1985_type( expr, geom, fev, feu ) );
        return res;
    }
    static boost::shared_ptr<self_type>
    New( std::string const& name, expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
    {
        boost::shared_ptr<self_type> res;
        if ( name == "classic" )
            res.reset( new tensor_volumic_part_classic_type( expr, geom, fev ) );
        else if ( name == "simo1985" )
            res.reset( new tensor_volumic_part_simo1985_type( expr, geom, fev ) );
        return res;
    }
    static boost::shared_ptr<self_type>
    New( std::string const& name, expr_type const& expr, Geo_t const& geom )
    {
        boost::shared_ptr<self_type> res;
        if ( name == "classic" )
            res.reset( new tensor_volumic_part_classic_type( expr, geom ) );
        else if ( name == "simo1985" )
            res.reset( new tensor_volumic_part_simo1985_type( expr, geom ) );
        return res;
    }

  private:
    expr_type const& M_expr;

  protected:
    typename super_type::array_scalar_type M_locEvalBulkModulus;
    typename super_type::array_shape_type M_locRes;
};

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic : public tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
    typedef ExprType expr_type;
    typedef typename super_type::value_type value_type;

    tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }

    void updateImpl( typename super_type::ctx_mechprop_scalar_type const& ctxMechProp, typename super_type::array_tensor2_type const& locGradEval )
    {
        super_type::updateImpl( ctxMechProp, locGradEval );
        updateImpl( locGradEval );
    }
    void updateImpl( typename super_type::array_tensor2_type const& locGradEval )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // F*\kappa*(J-1)*J*C^{-1} = \kappa*(J-1)*J*F^{-T}
            const value_type bulkModulus = this->M_locEvalBulkModulus[q]( 0, 0 );
            auto const& gradDisplacementEval = locGradEval[q];
            auto Id = super_type::loc_tensor2_type::Identity();
            auto F = Id + gradDisplacementEval;
            const value_type detFv = F.determinant();

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                auto InvFvmt = F.inverse().transpose();
                this->M_locRes[q] = bulkModulus * detFv * ( detFv - 1 ) * InvFvmt;
            }
            else
            {
                M_locEvalFmt[q] = F.inverse().transpose();
            }
        }
    }
    using super_type::evalijq; // fix clang warning
    typename super_type::matrix_shape_type const&
    evalijq( uint16_type i, uint16_type j, uint16_type q,
             const value_type detFv, const value_type ddetF,
             typename super_type::loc_tensor2_type const& dFmt ) const
    {
        const value_type bulkModulus = this->M_locEvalBulkModulus[q]( 0, 0 );
        this->locMatrixShape() = bulkModulus * ( detFv * ( detFv - 1 ) * dFmt + ddetF * ( 2 * detFv - 1 ) * M_locEvalFmt[q] );
        return this->locMatrixShape();
    }

  private:
    typename super_type::array_tensor2_type M_locEvalFmt;
};

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985 : public tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;
    typedef ExprType expr_type;
    typedef typename super_type::value_type value_type;

    tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalFmt( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }

    void updateImpl( typename super_type::ctx_mechprop_scalar_type const& ctxMechProp, typename super_type::array_tensor2_type const& locGradEval )
    {
        super_type::updateImpl( ctxMechProp, locGradEval );
        updateImpl( locGradEval );
    }
    void updateImpl( typename super_type::array_tensor2_type const& locGradEval )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            // F*\kappa*ln(J)*C^{-1} = \kappa*ln(J)*F^{-T}
            const value_type bulkModulus = this->M_locEvalBulkModulus[q]( 0, 0 );
            auto const& gradDisplacementEval = locGradEval[q];
            auto Id = super_type::loc_tensor2_type::Identity();
            auto F = Id + gradDisplacementEval;
            const value_type detFv = F.determinant();

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                auto InvFvmt = F.inverse().transpose();
                this->M_locRes[q] = bulkModulus * math::log( detFv ) * InvFvmt;
            }
            else
            {
                M_locEvalFmt[q] = F.inverse().transpose();
                M_locEvalPrecomputeLogDetF[q]( 0, 0 ) = bulkModulus * math::log( detFv );
            }
        }
    }
    using super_type::evalijq; // fix clang warning
    typename super_type::matrix_shape_type const&
    evalijq( uint16_type i, uint16_type j, uint16_type q,
             const value_type detFv, const value_type ddetF,
             typename super_type::loc_tensor2_type const& dFmt ) const
    {
        typename super_type::matrix_shape_type& locMat = this->locMatrixShape();
        const value_type bulkModulus = this->M_locEvalBulkModulus[q]( 0, 0 );
        const value_type precomputeLog = M_locEvalPrecomputeLogDetF[q]( 0, 0 );
        locMat = precomputeLog * dFmt;
        const value_type factorOther2 = bulkModulus * ddetF / detFv;
        typename super_type::loc_tensor2_type const& Fmt = M_locEvalFmt[q];
        locMat += factorOther2 * Fmt;
        return this->locMatrixShape();
    }

  private:
    typename super_type::array_scalar_type M_locEvalPrecomputeLogDetF;
    typename super_type::array_tensor2_type M_locEvalFmt;
};

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorFirstPiolaKirchhoffBase : public tensorBase<Geo_t, Basis_i_t, Basis_j_t, typename ExprType::shape_type, typename ExprType::value_type>
{
    typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t, typename ExprType::shape_type, typename ExprType::value_type> super_type;

    typedef ExprType expr_type;
    typedef typename expr_type::fe_displacement_type fe_displacement_type;
    typedef typename expr_type::geoelement_type geoelement_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::gmc_type gmc_type;
    // fe displacement context
    typedef typename fe_displacement_type::PreCompute pc_displacement_type;
    typedef boost::shared_ptr<pc_displacement_type> pc_displacement_ptrtype;
    typedef typename fe_displacement_type::template Context<expr_type::context_displacement, fe_displacement_type, gm_type, geoelement_type, gmc_type::context> ctx_displacement_type;
    typedef boost::shared_ptr<ctx_displacement_type> ctx_displacement_ptrtype;
    typedef typename expr_type::value_type value_type;

    typedef typename super_type::shape_type shape;
    typedef typename super_type::matrix_shape_type matrix_shape_type;
    typedef typename super_type::array_shape_type array_shape_type;

    tensorFirstPiolaKirchhoffBase( expr_type const& expr,
                                   Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( geom, fev, feu ),
          M_expr( expr ),
          M_pcDisplacement( new pc_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxDisplacement( new ctx_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc(), (pc_displacement_ptrtype const&)M_pcDisplacement ) ),
          M_locRes( expr.displacement().gradExtents( *this->gmc() ) ),
          M_locGradDisplacement( expr.displacement().gradExtents( *this->gmc() ) )
    {
    }
    tensorFirstPiolaKirchhoffBase( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( geom, fev ),
          M_expr( expr ),
          M_pcDisplacement( new pc_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxDisplacement( new ctx_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc(), (pc_displacement_ptrtype const&)M_pcDisplacement ) ),
          M_locRes( expr.displacement().gradExtents( *this->gmc() ) ),
          M_locGradDisplacement( expr.displacement().gradExtents( *this->gmc() ) )
    {
    }

    tensorFirstPiolaKirchhoffBase( expr_type const& expr, Geo_t const& geom )
        : super_type( geom ),
          M_expr( expr ),
          M_pcDisplacement( new pc_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxDisplacement( new ctx_displacement_type( expr.displacement().functionSpace()->fe(), this->gmc(), (pc_displacement_ptrtype const&)M_pcDisplacement ) ),
          M_locRes( expr.displacement().gradExtents( *this->gmc() ) ),
          M_locGradDisplacement( expr.displacement().gradExtents( *this->gmc() ) )
    {
    }

    expr_type const& expr() const { return M_expr; }

    array_shape_type const& locGradDisplacement() const { return M_locGradDisplacement; }
    matrix_shape_type const& locGradDisplacement( uint16_type q ) const { return M_locGradDisplacement[q]; }

    array_shape_type& locRes() { return M_locRes; }
    array_shape_type const& locRes() const { return M_locRes; }
    matrix_shape_type& locRes( uint16_type q ) { return M_locRes[q]; }
    matrix_shape_type const& locRes( uint16_type q ) const { return M_locRes[q]; }

    void update( Geo_t const& geom )
    {
        M_ctxDisplacement->update( this->gmc(), (pc_displacement_ptrtype const&)M_pcDisplacement );
        std::fill( M_locGradDisplacement.data(), M_locGradDisplacement.data() + M_locGradDisplacement.num_elements(), matrix_shape_type::Zero() );
        this->expr().displacement().grad( *M_ctxDisplacement, M_locGradDisplacement );
    }
    void update( Geo_t const& geom, uint16_type face )
    {
        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcDisplacement->update( this->gmc()->pc()->nodes() );
        M_ctxDisplacement->update( this->gmc(), (pc_displacement_ptrtype const&)M_pcDisplacement );
        std::fill( M_locGradDisplacement.data(), M_locGradDisplacement.data() + M_locGradDisplacement.num_elements(), matrix_shape_type::Zero() );
        this->expr().displacement().grad( *M_ctxDisplacement, M_locGradDisplacement );
    }

    matrix_shape_type const&
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        CHECK( false ) << "not allow\n";
        return this->M_locMatrixShape;
    }

    value_type
    evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        CHECK( false ) << "not allow\n";
        return 0;
    }

    value_type
    evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return evalq( c1, c2, q );
    }

    value_type
    evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
    {
        return M_locRes[q]( c1, c2 );
    }

  private:
    expr_type const& M_expr;

    pc_displacement_ptrtype M_pcDisplacement;
    ctx_displacement_ptrtype M_ctxDisplacement;

    array_shape_type M_locRes;
    array_shape_type M_locGradDisplacement;
};

/**
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     * StVenantKirchhoff :
     *   S = \lambda trace(E)*Id + 2*\mu*E
     *   S = 2*\mu*E ( + pressure part : p*Id )
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     */
template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType, bool useDispPresForm>
struct tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff : public tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;

  public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::value_type value_type;

    // fe material scalar context
    typedef typename expr_type::fe_mechprop_scalar_type fe_mechprop_scalar_type;
    static const size_type context_mechprop_scalar = expr_type::context_mechprop_scalar;
    typedef typename fe_mechprop_scalar_type::PreCompute pc_mechprop_scalar_type;
    typedef boost::shared_ptr<pc_mechprop_scalar_type> pc_mechprop_scalar_ptrtype;
    typedef typename fe_mechprop_scalar_type::template Context<context_mechprop_scalar, fe_mechprop_scalar_type, gm_type, geoelement_type, gmc_type::context> ctx_mechprop_scalar_type;
    typedef boost::shared_ptr<ctx_mechprop_scalar_type> ctx_mechprop_scalar_ptrtype;

    tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff( expr_type const& expr,
                                                        Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame1( this->expr().mechanicalPropertiesDesc().fieldCoeffLame1().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) )
    {
    }
    tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame1( this->expr().mechanicalPropertiesDesc().fieldCoeffLame1().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) )
    {
    }
    tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame1( this->expr().mechanicalPropertiesDesc().fieldCoeffLame1().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) )
    {
    }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        updateImpl( geom );
    }
    void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcMechPropField->update( this->gmc()->pc()->nodes() );
        updateImpl( geom );
    }
    void updateImpl( Geo_t const& geom )
    {
        M_ctxMechPropField->update( this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField );
        if ( !useDispPresForm )
        {
            std::fill( M_locEvalFieldCoefflame1.data(), M_locEvalFieldCoefflame1.data() + M_locEvalFieldCoefflame1.num_elements(), /*typename*/ super_type::loc_scalar_type::Zero() );
            this->expr().mechanicalPropertiesDesc().fieldCoeffLame1().id( *M_ctxMechPropField, M_locEvalFieldCoefflame1 );
        }
        std::fill( M_locEvalFieldCoefflame2.data(), M_locEvalFieldCoefflame2.data() + M_locEvalFieldCoefflame2.num_elements(), /*typename*/ super_type::loc_scalar_type::Zero() );
        this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().id( *M_ctxMechPropField, M_locEvalFieldCoefflame2 );

        std::fill( this->locRes().data(), this->locRes().data() + this->locRes().num_elements(), super_type::matrix_shape_type /*loc_res_type*/ ::Zero() );
        updateImpl( mpl::int_<expr_type::nRealDim>() );
    }
    void updateImpl( mpl::int_<2> /**/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type du1vdx = gradDisplacementEval( 0, 0 ), du1vdy = gradDisplacementEval( 0, 1 );
            const value_type du2vdx = gradDisplacementEval( 1, 0 ), du2vdy = gradDisplacementEval( 1, 1 );

            const value_type subtraceE1 = 0.5 * ( std::pow( du1vdx, 2 ) + std::pow( du2vdx, 2 ) );
            const value_type subtraceE2 = 0.5 * ( std::pow( du1vdy, 2 ) + std::pow( du2vdy, 2 ) );

            const value_type E11 = du1vdx + subtraceE1;
            const value_type E12 = 0.5 * ( du1vdy + du2vdx + du1vdx * du1vdy + du2vdx * du2vdy );
            const value_type E21 = 0.5 * ( du2vdx + du1vdy + du1vdy * du1vdx + du2vdy * du2vdx );
            const value_type E22 = du2vdy + subtraceE2;

            typename super_type::loc_tensor2_type& theLocRes = this->locRes( q );

            if ( useDispPresForm )
            {
                const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
                const value_type S11 = 2 * coefflame2 * E11;
                const value_type S12 = 2 * coefflame2 * E12;
                const value_type S21 = 2 * coefflame2 * E21;
                const value_type S22 = 2 * coefflame2 * E22;
                if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
                {
                    // F*S
                    theLocRes( 0, 0 ) = ( 1 + du1vdx ) * S11 + du1vdy * S21;
                    theLocRes( 0, 1 ) = ( 1 + du1vdx ) * S12 + du1vdy * S22;
                    theLocRes( 1, 0 ) = du2vdx * S11 + ( 1 + du2vdy ) * S21;
                    theLocRes( 1, 1 ) = du2vdx * S12 + ( 1 + du2vdy ) * S22;
                }
                else
                {
                    theLocRes( 0, 0 ) = S11;
                    theLocRes( 0, 1 ) = S12;
                    theLocRes( 1, 0 ) = S21;
                    theLocRes( 1, 1 ) = S22;
                }
            }
            else
            {
                const value_type coefflame1 = M_locEvalFieldCoefflame1[q]( 0, 0 );
                const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
                const value_type traceE = E11 + E22;
                const value_type S11 = coefflame1 * traceE + 2 * coefflame2 * E11;
                const value_type S12 = 2 * coefflame2 * E12;
                const value_type S21 = 2 * coefflame2 * E21;
                const value_type S22 = coefflame1 * traceE + 2 * coefflame2 * E22;

                if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
                {
                    // F*S
                    theLocRes( 0, 0 ) = ( 1 + du1vdx ) * S11 + du1vdy * S21;
                    theLocRes( 0, 1 ) = ( 1 + du1vdx ) * S12 + du1vdy * S22;
                    theLocRes( 1, 0 ) = du2vdx * S11 + ( 1 + du2vdy ) * S21;
                    theLocRes( 1, 1 ) = du2vdx * S12 + ( 1 + du2vdy ) * S22;
                }
                else
                {
                    theLocRes( 0, 0 ) = S11;
                    theLocRes( 0, 1 ) = S12;
                    theLocRes( 1, 0 ) = S21;
                    theLocRes( 1, 1 ) = S22;
                }
            }
        }
    }

    void updateImpl( mpl::int_<3> /**/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type du1vdx = gradDisplacementEval( 0, 0 ), du1vdy = gradDisplacementEval( 0, 1 ), du1vdz = gradDisplacementEval( 0, 2 );
            const value_type du2vdx = gradDisplacementEval( 1, 0 ), du2vdy = gradDisplacementEval( 1, 1 ), du2vdz = gradDisplacementEval( 1, 2 );
            const value_type du3vdx = gradDisplacementEval( 2, 0 ), du3vdy = gradDisplacementEval( 2, 1 ), du3vdz = gradDisplacementEval( 2, 2 );

            const value_type subtraceE1 = 0.5 * ( std::pow( du1vdx, 2 ) + std::pow( du2vdx, 2 ) + std::pow( du3vdx, 2 ) );
            const value_type subtraceE2 = 0.5 * ( std::pow( du1vdy, 2 ) + std::pow( du2vdy, 2 ) + std::pow( du3vdy, 2 ) );
            const value_type subtraceE3 = 0.5 * ( std::pow( du1vdz, 2 ) + std::pow( du2vdz, 2 ) + std::pow( du3vdz, 2 ) );

            const value_type E11 = du1vdx + subtraceE1;
            const value_type E12 = 0.5 * ( du1vdy + du2vdx + du1vdx * du1vdy + du2vdx * du2vdy + du3vdx * du3vdy );
            const value_type E13 = 0.5 * ( du1vdz + du3vdx + du1vdx * du1vdz + du2vdx * du2vdz + du3vdx * du3vdz );
            const value_type E21 = 0.5 * ( du2vdx + du1vdy + du1vdy * du1vdx + du2vdy * du2vdx + du3vdy * du3vdx );
            const value_type E22 = du2vdy + subtraceE2;
            const value_type E23 = 0.5 * ( du2vdz + du3vdy + du1vdy * du1vdz + du2vdy * du2vdz + du3vdy * du3vdz );
            const value_type E31 = 0.5 * ( du3vdx + du1vdz + du1vdz * du1vdx + du2vdz * du2vdx + du3vdz * du3vdx );
            const value_type E32 = 0.5 * ( du3vdy + du2vdz + du1vdz * du1vdy + du2vdz * du2vdy + du3vdz * du3vdy );
            const value_type E33 = du3vdz + subtraceE3;

            typename super_type::loc_tensor2_type& theLocRes = this->locRes( q );
            if ( useDispPresForm )
            {
                const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
                const value_type S11 = 2 * coefflame2 * E11;
                const value_type S12 = 2 * coefflame2 * E12;
                const value_type S13 = 2 * coefflame2 * E13;
                const value_type S21 = 2 * coefflame2 * E21;
                const value_type S22 = 2 * coefflame2 * E22;
                const value_type S23 = 2 * coefflame2 * E23;
                const value_type S31 = 2 * coefflame2 * E31;
                const value_type S32 = 2 * coefflame2 * E32;
                const value_type S33 = 2 * coefflame2 * E33;

                if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
                {
                    // F*S
                    theLocRes( 0, 0 ) = ( 1 + du1vdx ) * S11 + du1vdy * S21 + du1vdz * S31;
                    theLocRes( 0, 1 ) = ( 1 + du1vdx ) * S12 + du1vdy * S22 + du1vdz * S32;
                    theLocRes( 0, 2 ) = ( 1 + du1vdx ) * S13 + du1vdy * S23 + du1vdz * S33;
                    theLocRes( 1, 0 ) = du2vdx * S11 + ( 1 + du2vdy ) * S21 + du2vdz * S31;
                    theLocRes( 1, 1 ) = du2vdx * S12 + ( 1 + du2vdy ) * S22 + du2vdz * S32;
                    theLocRes( 1, 2 ) = du2vdx * S13 + ( 1 + du2vdy ) * S23 + du2vdz * S33;
                    theLocRes( 2, 0 ) = du3vdx * S11 + du3vdy * S21 + ( 1 + du3vdz ) * S31;
                    theLocRes( 2, 1 ) = du3vdx * S12 + du3vdy * S22 + ( 1 + du3vdz ) * S32;
                    theLocRes( 2, 2 ) = du3vdx * S13 + du3vdy * S23 + ( 1 + du3vdz ) * S33;
                }
                else
                {
                    // S
                    theLocRes( 0, 0 ) = S11;
                    theLocRes( 0, 1 ) = S12;
                    theLocRes( 0, 2 ) = S13;
                    theLocRes( 1, 0 ) = S21;
                    theLocRes( 1, 1 ) = S22;
                    theLocRes( 1, 2 ) = S23;
                    theLocRes( 2, 0 ) = S31;
                    theLocRes( 2, 1 ) = S32;
                    theLocRes( 2, 2 ) = S33;
                }
            }
            else
            {
                const value_type coefflame1 = M_locEvalFieldCoefflame1[q]( 0, 0 );
                const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
                const value_type traceE = E11 + E22 + E33;
                const value_type S11 = coefflame1 * traceE + 2 * coefflame2 * E11;
                const value_type S12 = 2 * coefflame2 * E12;
                const value_type S13 = 2 * coefflame2 * E13;
                const value_type S21 = 2 * coefflame2 * E21;
                const value_type S22 = coefflame1 * traceE + 2 * coefflame2 * E22;
                const value_type S23 = 2 * coefflame2 * E23;
                const value_type S31 = 2 * coefflame2 * E31;
                const value_type S32 = 2 * coefflame2 * E32;
                const value_type S33 = coefflame1 * traceE + 2 * coefflame2 * E33;

                if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
                {
                    // F*S
                    theLocRes( 0, 0 ) = ( 1 + du1vdx ) * S11 + du1vdy * S21 + du1vdz * S31;
                    theLocRes( 0, 1 ) = ( 1 + du1vdx ) * S12 + du1vdy * S22 + du1vdz * S32;
                    theLocRes( 0, 2 ) = ( 1 + du1vdx ) * S13 + du1vdy * S23 + du1vdz * S33;
                    theLocRes( 1, 0 ) = du2vdx * S11 + ( 1 + du2vdy ) * S21 + du2vdz * S31;
                    theLocRes( 1, 1 ) = du2vdx * S12 + ( 1 + du2vdy ) * S22 + du2vdz * S32;
                    theLocRes( 1, 2 ) = du2vdx * S13 + ( 1 + du2vdy ) * S23 + du2vdz * S33;
                    theLocRes( 2, 0 ) = du3vdx * S11 + du3vdy * S21 + ( 1 + du3vdz ) * S31;
                    theLocRes( 2, 1 ) = du3vdx * S12 + du3vdy * S22 + ( 1 + du3vdz ) * S32;
                    theLocRes( 2, 2 ) = du3vdx * S13 + du3vdy * S23 + ( 1 + du3vdz ) * S33;
                }
                else
                {
                    // S
                    theLocRes( 0, 0 ) = S11;
                    theLocRes( 0, 1 ) = S12;
                    theLocRes( 0, 2 ) = S13;
                    theLocRes( 1, 0 ) = S21;
                    theLocRes( 1, 1 ) = S22;
                    theLocRes( 1, 2 ) = S23;
                    theLocRes( 2, 0 ) = S31;
                    theLocRes( 2, 1 ) = S32;
                    theLocRes( 2, 2 ) = S33;
                }
            }
        }
    }
    using super_type::evalijq; // fix clang warning

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i, j, q, mpl::int_<expr_type::specific_expr_type::value>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::EVAL> /**/ ) const
    {
        CHECK( false ) << "not allow";
        return this->locMatrixShape();
    }

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/ ) const
    {
        return evalijq( i, j, q, mpl::int_<ExprApplyType::JACOBIAN>(), mpl::int_<expr_type::nRealDim>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<2> /**/ ) const
    {
        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type du1vdx = gradDisplacementEval( 0, 0 ), du1vdy = gradDisplacementEval( 0, 1 );
        const value_type du2vdx = gradDisplacementEval( 1, 0 ), du2vdy = gradDisplacementEval( 1, 1 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type du1tdx = gradTrial( 0, 0 ), du1tdy = gradTrial( 0, 1 );
        const value_type du2tdx = gradTrial( 1, 0 ), du2tdy = gradTrial( 1, 1 );

        // dFSv = dF*val(Sv)
        auto const& localEval = this->locRes( q );
        const value_type dFSv11 = du1tdx * localEval( 0, 0 ) + du1tdy * localEval( 1, 0 );
        const value_type dFSv12 = du1tdx * localEval( 0, 1 ) + du1tdy * localEval( 1, 1 );
        const value_type dFSv21 = du2tdx * localEval( 0, 0 ) + du2tdy * localEval( 1, 0 );
        const value_type dFSv22 = du2tdx * localEval( 0, 1 ) + du2tdy * localEval( 1, 1 );

        //auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
        const value_type dE11 = du1tdx + /*0.5**/ ( du1vdx * du1tdx + du2vdx * du2tdx );
        const value_type dE12 = 0.5 * ( du1tdy + du2tdx + du1tdx * du1vdy + du2tdx * du2vdy + du1vdx * du1tdy + du2vdx * du2tdy );
        const value_type dE21 = 0.5 * ( du2tdx + du1tdy + du1tdy * du1vdx + du2tdy * du2vdx + du1vdy * du1tdx + du2vdy * du2tdx );
        const value_type dE22 = du2tdy + du1vdy * du1tdy + du2vdy * du2tdy;

        typename super_type::matrix_shape_type& thelocMat = this->locMatrixShape();
        if ( useDispPresForm )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
            // auto dS = 2*idv(CoeffLame2)*dE;
            const value_type dS11 = 2 * coefflame2 * dE11;
            const value_type dS12 = 2 * coefflame2 * dE12;
            const value_type dS21 = 2 * coefflame2 * dE21;
            const value_type dS22 = 2 * coefflame2 * dE22;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = ( 1. + du1vdx ) * dS11 + du1vdy * dS21;
            const value_type FvdS12 = ( 1. + du1vdx ) * dS12 + du1vdy * dS22;
            const value_type FvdS21 = du2vdx * dS11 + ( 1. + du2vdy ) * dS21;
            const value_type FvdS22 = du2vdx * dS12 + ( 1. + du2vdy ) * dS22;

            // dF*val(Sv) + val(Fv)*dS
            thelocMat( 0, 0 ) = dFSv11 + FvdS11;
            thelocMat( 0, 1 ) = dFSv12 + FvdS12;
            thelocMat( 1, 0 ) = dFSv21 + FvdS21;
            thelocMat( 1, 1 ) = dFSv22 + FvdS22;
        }
        else
        {
            const value_type coefflame1 = M_locEvalFieldCoefflame1[q]( 0, 0 );
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
            const value_type tracedE = dE11 + dE22;

            // auto dS = idv(CoeffLame1)*trace(dE)*Id + 2*idv(CoeffLame2)*dE;
            const value_type dS11 = coefflame1 * tracedE + 2 * coefflame2 * dE11;
            const value_type dS12 = 2 * coefflame2 * dE12;
            const value_type dS21 = 2 * coefflame2 * dE21;
            const value_type dS22 = coefflame1 * tracedE + 2 * coefflame2 * dE22;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = ( 1. + du1vdx ) * dS11 + du1vdy * dS21;
            const value_type FvdS12 = ( 1. + du1vdx ) * dS12 + du1vdy * dS22;
            const value_type FvdS21 = du2vdx * dS11 + ( 1. + du2vdy ) * dS21;
            const value_type FvdS22 = du2vdx * dS12 + ( 1. + du2vdy ) * dS22;

            // dF*val(Sv) + val(Fv)*dS
            thelocMat( 0, 0 ) = dFSv11 + FvdS11;
            thelocMat( 0, 1 ) = dFSv12 + FvdS12;
            thelocMat( 1, 0 ) = dFSv21 + FvdS21;
            thelocMat( 1, 1 ) = dFSv22 + FvdS22;
        }
        return this->locMatrixShape();
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<3> /**/ ) const
    {
        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type du1vdx = gradDisplacementEval( 0, 0 ), du1vdy = gradDisplacementEval( 0, 1 ), du1vdz = gradDisplacementEval( 0, 2 );
        const value_type du2vdx = gradDisplacementEval( 1, 0 ), du2vdy = gradDisplacementEval( 1, 1 ), du2vdz = gradDisplacementEval( 1, 2 );
        const value_type du3vdx = gradDisplacementEval( 2, 0 ), du3vdy = gradDisplacementEval( 2, 1 ), du3vdz = gradDisplacementEval( 2, 2 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type du1tdx = gradTrial( 0, 0 ), du1tdy = gradTrial( 0, 1 ), du1tdz = gradTrial( 0, 2 );
        const value_type du2tdx = gradTrial( 1, 0 ), du2tdy = gradTrial( 1, 1 ), du2tdz = gradTrial( 1, 2 );
        const value_type du3tdx = gradTrial( 2, 0 ), du3tdy = gradTrial( 2, 1 ), du3tdz = gradTrial( 2, 2 );

        // dFSv = dF*val(Sv)
        auto const& localEval = this->locRes( q );
        const value_type dFSv11 = du1tdx * localEval( 0, 0 ) + du1tdy * localEval( 1, 0 ) + du1tdz * localEval( 2, 0 );
        const value_type dFSv12 = du1tdx * localEval( 0, 1 ) + du1tdy * localEval( 1, 1 ) + du1tdz * localEval( 2, 1 );
        const value_type dFSv13 = du1tdx * localEval( 0, 2 ) + du1tdy * localEval( 1, 2 ) + du1tdz * localEval( 2, 2 );

        const value_type dFSv21 = du2tdx * localEval( 0, 0 ) + du2tdy * localEval( 1, 0 ) + du2tdz * localEval( 2, 0 );
        const value_type dFSv22 = du2tdx * localEval( 0, 1 ) + du2tdy * localEval( 1, 1 ) + du2tdz * localEval( 2, 1 );
        const value_type dFSv23 = du2tdx * localEval( 0, 2 ) + du2tdy * localEval( 1, 2 ) + du2tdz * localEval( 2, 2 );

        const value_type dFSv31 = du3tdx * localEval( 0, 0 ) + du3tdy * localEval( 1, 0 ) + du3tdz * localEval( 2, 0 );
        const value_type dFSv32 = du3tdx * localEval( 0, 1 ) + du3tdy * localEval( 1, 1 ) + du3tdz * localEval( 2, 1 );
        const value_type dFSv33 = du3tdx * localEval( 0, 2 ) + du3tdy * localEval( 1, 2 ) + du3tdz * localEval( 2, 2 );

        //auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
        const value_type dE11 = du1tdx + /*0.5**/ ( du1vdx * du1tdx + du2vdx * du2tdx + du3vdx * du3tdx );
        const value_type dE12 = 0.5 * ( du1tdy + du2tdx + du1tdx * du1vdy + du2tdx * du2vdy + du3tdx * du3vdy + du1vdx * du1tdy + du2vdx * du2tdy + du3vdx * du3tdy );
        const value_type dE13 = 0.5 * ( du1tdz + du3tdx + du1tdx * du1vdz + du2tdx * du2vdz + du3tdx * du3vdz + du1vdx * du1tdz + du2vdx * du2tdz + du3vdx * du3tdz );
        const value_type dE21 = 0.5 * ( du2tdx + du1tdy + du1tdy * du1vdx + du2tdy * du2vdx + du3tdy * du3vdx + du1vdy * du1tdx + du2vdy * du2tdx + du3vdy * du3tdx );
        const value_type dE22 = du2tdy + du1vdy * du1tdy + du2vdy * du2tdy + du3vdy * du3tdy;
        const value_type dE23 = 0.5 * ( du2tdz + du3tdy + du1tdy * du1vdz + du2tdy * du2vdz + du3tdy * du3vdz + du1vdy * du1tdz + du2vdy * du2tdz + du3vdy * du3tdz );
        const value_type dE31 = 0.5 * ( du3tdx + du1tdz + du1tdz * du1vdx + du2tdz * du2vdx + du3tdz * du3vdx + du1vdz * du1tdx + du2vdz * du2tdx + du3vdz * du3tdx );
        const value_type dE32 = 0.5 * ( du3tdy + du2tdz + du1tdz * du1vdy + du2tdz * du2vdy + du3tdz * du3vdy + du1vdz * du1tdy + du2vdz * du2tdy + du3vdz * du3tdy );
        const value_type dE33 = du3tdz + du1vdz * du1tdz + du2vdz * du2tdz + du3vdz * du3tdz;

        typename super_type::matrix_shape_type& thelocMat = this->locMatrixShape();
        if ( useDispPresForm )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

            // auto dS = 2*idv(CoeffLame2)*dE;
            const value_type dS11 = 2 * coefflame2 * dE11;
            const value_type dS12 = 2 * coefflame2 * dE12;
            const value_type dS13 = 2 * coefflame2 * dE13;
            const value_type dS21 = 2 * coefflame2 * dE21;
            const value_type dS22 = 2 * coefflame2 * dE22;
            const value_type dS23 = 2 * coefflame2 * dE23;
            const value_type dS31 = 2 * coefflame2 * dE31;
            const value_type dS32 = 2 * coefflame2 * dE32;
            const value_type dS33 = 2 * coefflame2 * dE33;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = ( 1. + du1vdx ) * dS11 + du1vdy * dS21 + du1vdz * dS31;
            const value_type FvdS12 = ( 1. + du1vdx ) * dS12 + du1vdy * dS22 + du1vdz * dS32;
            const value_type FvdS13 = ( 1. + du1vdx ) * dS13 + du1vdy * dS23 + du1vdz * dS33;
            const value_type FvdS21 = du2vdx * dS11 + ( 1. + du2vdy ) * dS21 + du2vdz * dS31;
            const value_type FvdS22 = du2vdx * dS12 + ( 1. + du2vdy ) * dS22 + du2vdz * dS32;
            const value_type FvdS23 = du2vdx * dS13 + ( 1. + du2vdy ) * dS23 + du2vdz * dS33;
            const value_type FvdS31 = du3vdx * dS11 + du3vdy * dS21 + ( 1. + du3vdz ) * dS31;
            const value_type FvdS32 = du3vdx * dS12 + du3vdy * dS22 + ( 1. + du3vdz ) * dS32;
            const value_type FvdS33 = du3vdx * dS13 + du3vdy * dS23 + ( 1. + du3vdz ) * dS33;

            // dF*val(Sv) + val(Fv)*dS
            thelocMat( 0, 0 ) = dFSv11 + FvdS11;
            thelocMat( 0, 1 ) = dFSv12 + FvdS12;
            thelocMat( 0, 2 ) = dFSv13 + FvdS13;
            thelocMat( 1, 0 ) = dFSv21 + FvdS21;
            thelocMat( 1, 1 ) = dFSv22 + FvdS22;
            thelocMat( 1, 2 ) = dFSv23 + FvdS23;
            thelocMat( 2, 0 ) = dFSv31 + FvdS31;
            thelocMat( 2, 1 ) = dFSv32 + FvdS32;
            thelocMat( 2, 2 ) = dFSv32 + FvdS33;
        }
        else
        {
            const value_type coefflame1 = M_locEvalFieldCoefflame1[q]( 0, 0 );
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
            const value_type tracedE = dE11 + dE22 + dE33;

            // auto dS = idv(CoeffLame1)*trace(dE)*Id + 2*idv(CoeffLame2)*dE;
            const value_type dS11 = coefflame1 * tracedE + 2 * coefflame2 * dE11;
            const value_type dS12 = 2 * coefflame2 * dE12;
            const value_type dS13 = 2 * coefflame2 * dE13;
            const value_type dS21 = 2 * coefflame2 * dE21;
            const value_type dS22 = coefflame1 * tracedE + 2 * coefflame2 * dE22;
            const value_type dS23 = 2 * coefflame2 * dE23;
            const value_type dS31 = 2 * coefflame2 * dE31;
            const value_type dS32 = 2 * coefflame2 * dE32;
            const value_type dS33 = coefflame1 * tracedE + 2 * coefflame2 * dE33;

            // FvdS = val(Fv)*dS
            const value_type FvdS11 = ( 1. + du1vdx ) * dS11 + du1vdy * dS21 + du1vdz * dS31;
            const value_type FvdS12 = ( 1. + du1vdx ) * dS12 + du1vdy * dS22 + du1vdz * dS32;
            const value_type FvdS13 = ( 1. + du1vdx ) * dS13 + du1vdy * dS23 + du1vdz * dS33;
            const value_type FvdS21 = du2vdx * dS11 + ( 1. + du2vdy ) * dS21 + du2vdz * dS31;
            const value_type FvdS22 = du2vdx * dS12 + ( 1. + du2vdy ) * dS22 + du2vdz * dS32;
            const value_type FvdS23 = du2vdx * dS13 + ( 1. + du2vdy ) * dS23 + du2vdz * dS33;
            const value_type FvdS31 = du3vdx * dS11 + du3vdy * dS21 + ( 1. + du3vdz ) * dS31;
            const value_type FvdS32 = du3vdx * dS12 + du3vdy * dS22 + ( 1. + du3vdz ) * dS32;
            const value_type FvdS33 = du3vdx * dS13 + du3vdy * dS23 + ( 1. + du3vdz ) * dS33;

            // dF*val(Sv) + val(Fv)*dS
            thelocMat( 0, 0 ) = dFSv11 + FvdS11;
            thelocMat( 0, 1 ) = dFSv12 + FvdS12;
            thelocMat( 0, 2 ) = dFSv13 + FvdS13;
            thelocMat( 1, 0 ) = dFSv21 + FvdS21;
            thelocMat( 1, 1 ) = dFSv22 + FvdS22;
            thelocMat( 1, 2 ) = dFSv23 + FvdS23;
            thelocMat( 2, 0 ) = dFSv31 + FvdS31;
            thelocMat( 2, 1 ) = dFSv32 + FvdS32;
            thelocMat( 2, 2 ) = dFSv32 + FvdS33;
        }
        return this->locMatrixShape();
    }

  private:
    pc_mechprop_scalar_ptrtype M_pcMechPropField;
    ctx_mechprop_scalar_ptrtype M_ctxMechPropField;

    typename super_type::array_scalar_type M_locEvalFieldCoefflame1, M_locEvalFieldCoefflame2;

}; // tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff

/**
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     * NeoHookeanCompressible :
     *   S = \mu*J^{-2/3}*(I - (1/3)*trace(C)*C^{-1}) (+ pressure part : pJC^{-1})
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     */
template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType, bool useDispPresForm>
struct tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible : public tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;

  public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::value_type value_type;

    // fe mech prop context
    typedef typename expr_type::fe_mechprop_scalar_type fe_mechprop_scalar_type;
    static const size_type context_mechprop_scalar = expr_type::context_mechprop_scalar;
    typedef typename fe_mechprop_scalar_type::PreCompute pc_mechprop_scalar_type;
    typedef boost::shared_ptr<pc_mechprop_scalar_type> pc_mechprop_scalar_ptrtype;
    typedef typename fe_mechprop_scalar_type::template Context<context_mechprop_scalar, fe_mechprop_scalar_type, gm_type, geoelement_type, gmc_type::context> ctx_mechprop_scalar_type;
    typedef boost::shared_ptr<ctx_mechprop_scalar_type> ctx_mechprop_scalar_ptrtype;

    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> tensor_volumic_part_type;

    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(),
                                                            this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeDetF1( boost::extents[this->gmc()->xRefs().size2()] ),
          //M_locEvalPrecomputeDetF2( boost::extents[ this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF3( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeInvDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTraceC( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTrialDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom, fev, feu );
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(),
                                                            this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeDetF1( boost::extents[this->gmc()->xRefs().size2()] ),
          //M_locEvalPrecomputeDetF2( boost::extents[ this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF3( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeInvDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTraceC( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTrialDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom, fev );
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(),
                                                            this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeDetF1( boost::extents[this->gmc()->xRefs().size2()] ),
          //M_locEvalPrecomputeDetF2( boost::extents[ this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF3( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeInvDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTraceC( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeTrialDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom );
    }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        updateImpl( geom );
    }
    void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcMechPropField->update( this->gmc()->pc()->nodes() );
        updateImpl( geom );
    }
#if 1
    void updateImpl( Geo_t const& geom )
    {
        M_ctxMechPropField->update( this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField );
        std::fill( M_locEvalFieldCoefflame2.data(), M_locEvalFieldCoefflame2.data() + M_locEvalFieldCoefflame2.num_elements(), super_type::loc_scalar_type::Zero() );
        this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().id( *M_ctxMechPropField, M_locEvalFieldCoefflame2 );

        if ( !useDispPresForm )
            M_tensorVolumicPart->updateImpl( *M_ctxMechPropField, this->locGradDisplacement() );

        std::fill( this->locRes().data(), this->locRes().data() + this->locRes().num_elements(), super_type::matrix_shape_type /*loc_res_type*/ ::Zero() );
        //updateImpl( /*mpl::int_<expr_type::nRealDim>()*/ );
        updateImpl( mpl::int_<expr_type::nRealDim>() );
    }
#else
    // generic update but seems more slow
    void updateImpl( Geo_t const& geom )
    {
        M_ctxMechPropField->update( this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField );
        std::fill( M_locEvalFieldCoefflame2.data(), M_locEvalFieldCoefflame2.data() + M_locEvalFieldCoefflame2.num_elements(), super_type::loc_scalar_type::Zero() );
        this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().id( *M_ctxMechPropField, M_locEvalFieldCoefflame2 );

        if ( !useDispPresForm )
            M_tensorVolumicPart->updateImpl( *M_ctxMechPropField, this->locGradDisplacement() );

        std::fill( this->locRes().data(), this->locRes().data() + this->locRes().num_elements(), super_type::matrix_shape_type /*loc_res_type*/ ::Zero() );

        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
            auto Id = super_type::loc_tensor2_type::Identity();
            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            auto F = Id + gradDisplacementEval;
            const value_type detFv = F.determinant();
            //const value_type traceCv = std::pow(Fv11,2) + std::pow(Fv21,2) + std::pow(Fv12,2) + std::pow(Fv22,2);
            const value_type traceCv = F.diagonal().squaredNorm();
            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                auto InvFvmt = F.inverse().transpose();
                //const value_type A10 = coefflame2/2;
                this->locRes( q ) = coefflame2 * std::pow( detFv, -2. / expr_type::nRealDim ) * ( F - ( 1. / expr_type::nRealDim ) * traceCv * InvFvmt );
                if ( !useDispPresForm )
                    this->locRes( q ) += M_tensorVolumicPart->locRes( q );
            }
            else
            {
                const value_type factorWithDetFv1 = 1. / math::pow( detFv, 2 );
                const value_type factorWithDetFv2 = -( 2. / expr_type::nRealDim ) * coefflame2 * math::pow( detFv, ( ( expr_type::nRealDim == 2 ) ? -4. : -5. ) / expr_type::nRealDim );
                //const value_type factorWithDetFv2 = -(2./3.)*coefflame2*math::pow(detFv,-5./3.);
                const value_type factorWithDetFv3 = coefflame2 * math::pow( detFv, -2. / expr_type::nRealDim );
                M_locEvalPrecomputeInvDetF[q]( 0, 0 ) = 1. / detFv;
                M_locEvalPrecomputeDetF1[q]( 0, 0 ) = factorWithDetFv1;
                //M_locEvalPrecomputeDetF2[q](0,0) = factorWithDetFv2;
                M_locEvalPrecomputeDetF3[q]( 0, 0 ) = factorWithDetFv3;
                M_locEvalPrecomputeTraceC[q]( 0, 0 ) = traceCv;
                //TODO here
                const value_type factorTrialDetF = ( 1. / expr_type::nRealDim ) * traceCv;
                typename super_type::loc_tensor2_type& matLocTrialDetF = M_locEvalPrecomputeTrialDetF[q];
                matLocTrialDetF = factorWithDetFv2 * ( F - factorTrialDetF * F.inverse().transpose() );
            }
        }
    }
#endif
    void updateImpl( mpl::int_<2> /**/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 );
            const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 );
            const value_type detFv = Fv11 * Fv22 - Fv12 * Fv21;
            const value_type traceCv = std::pow( Fv11, 2 ) + std::pow( Fv21, 2 ) + std::pow( Fv12, 2 ) + std::pow( Fv22, 2 );

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                const value_type scaleUsedWithInvFv = 1. / detFv;
                const value_type factorWithDetFv = coefflame2 * std::pow( detFv, -2. / 2. );
                const value_type factorWithTraceCv = ( 1. / 2. ) * traceCv;
                const value_type factorOther = factorWithTraceCv * scaleUsedWithInvFv;
                //auto InvFv = (cst(1.)/detFv)*mat<2,2>( Fv22,-Fv12,-Fv21,Fv11);
                //auto FSv_neohookean = idv(CoeffLame2)*pow(detFv,-2./2.)*( Fv - (1./2.)*traceCv*trans(InvFv) );
                typename super_type::matrix_shape_type& matLoc = this->locRes( q );
                matLoc( 0, 0 ) = factorWithDetFv * ( Fv11 - factorOther * Fv22 );
                matLoc( 1, 0 ) = factorWithDetFv * ( Fv21 + factorOther * Fv12 );
                matLoc( 0, 1 ) = factorWithDetFv * ( Fv12 + factorOther * Fv21 );
                matLoc( 1, 1 ) = factorWithDetFv * ( Fv22 - factorOther * Fv11 );
                if ( !useDispPresForm )
                    matLoc += M_tensorVolumicPart->locRes( q );
            }
            else
            {
                const value_type factorWithDetFv1 = 1. / math::pow( detFv, 2 );
                const value_type factorWithDetFv2 = -( 2. / 2. ) * coefflame2 * math::pow( detFv, -4. / 2. );
                const value_type factorWithDetFv3 = coefflame2 * math::pow( detFv, -2. / 2. );
                M_locEvalPrecomputeInvDetF[q]( 0, 0 ) = 1. / detFv;
                M_locEvalPrecomputeDetF1[q]( 0, 0 ) = factorWithDetFv1;
                //M_locEvalPrecomputeDetF2[q](0,0) = factorWithDetFv2;
                M_locEvalPrecomputeDetF3[q]( 0, 0 ) = factorWithDetFv3;
                M_locEvalPrecomputeTraceC[q]( 0, 0 ) = traceCv;
#if 1 // NEW
                const value_type factorTrialDetF = ( 1. / 2. ) * traceCv * ( 1. / detFv );
                typename super_type::loc_tensor2_type& matLocTrialDetF = M_locEvalPrecomputeTrialDetF[q];
                matLocTrialDetF( 0, 0 ) = factorWithDetFv2 * ( Fv11 - factorTrialDetF * Fv22 );
                matLocTrialDetF( 0, 1 ) = factorWithDetFv2 * ( Fv12 + factorTrialDetF * Fv21 );
                matLocTrialDetF( 1, 0 ) = factorWithDetFv2 * ( Fv21 + factorTrialDetF * Fv12 );
                matLocTrialDetF( 1, 1 ) = factorWithDetFv2 * ( Fv22 - factorTrialDetF * Fv11 );
#endif
            }
        }
    }

    void updateImpl( mpl::int_<3> /**/ )
    {

        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 ), Fv13 = gradDisplacementEval( 0, 2 );
            const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 ), Fv23 = gradDisplacementEval( 1, 2 );
            const value_type Fv31 = gradDisplacementEval( 2, 0 ), Fv32 = gradDisplacementEval( 2, 1 ), Fv33 = 1. + gradDisplacementEval( 2, 2 );

            const value_type detFv = Fv11 * ( Fv22 * Fv33 - Fv23 * Fv32 ) - Fv21 * ( Fv12 * Fv33 - Fv13 * Fv32 ) + Fv31 * ( Fv12 * Fv23 - Fv13 * Fv22 );
#if 0
                //const value_type detFv = Fv11*(Fv22*Fv33-Fv23*Fv32) - Fv12*(Fv21*Fv33-Fv31*Fv23) + Fv13*(Fv21*Fv32 - Fv31*Fv22);
                if ( detFv < 0 )
                {
                    std::cout << "Negative : " << detFv << "\n"
                              << Fv11 << " " << Fv12 << " " << Fv13 << "\n"
                              << Fv21 << " " << Fv22 << " " << Fv23 << "\n"
                              << Fv31 << " " << Fv32 << " " << Fv33 << "\n";
                }
#endif
            const value_type traceCv = math::pow( Fv11, 2 ) + math::pow( Fv21, 2 ) + math::pow( Fv31, 2 ) + math::pow( Fv12, 2 ) +
                                       math::pow( Fv22, 2 ) + math::pow( Fv32, 2 ) + math::pow( Fv13, 2 ) + math::pow( Fv23, 2 ) + math::pow( Fv33, 2 );

            /*auto InvFv = (cst(1.)/detFv)*mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
                 Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
                 Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21 );*/
            // InvFv ( without 1/detF )
            const value_type InvFv11 = Fv22 * Fv33 - Fv23 * Fv32;
            const value_type InvFv12 = Fv13 * Fv32 - Fv12 * Fv33;
            const value_type InvFv13 = Fv12 * Fv23 - Fv13 * Fv22;
            const value_type InvFv21 = Fv23 * Fv31 - Fv21 * Fv33;
            const value_type InvFv22 = Fv11 * Fv33 - Fv13 * Fv31;
            const value_type InvFv23 = Fv13 * Fv21 - Fv11 * Fv23;
            const value_type InvFv31 = Fv21 * Fv32 - Fv22 * Fv31;
            const value_type InvFv32 = Fv12 * Fv31 - Fv11 * Fv32;
            const value_type InvFv33 = Fv11 * Fv22 - Fv12 * Fv21;

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {

                //auto FSv_neohookean = idv(CoeffLame2)*pow(detFv,-2./3.)*( Fv - (1./3.)*traceCv*trans(InvFv) );
                const value_type scaleUsedWithInvFv = 1. / detFv;
                const value_type factorWithDetFv = coefflame2 * math::pow( detFv, -2. / 3. );
                const value_type factorWithTraceCv = ( 1. / 3. ) * traceCv;
                const value_type factorOther = factorWithTraceCv * scaleUsedWithInvFv;
                typename super_type::matrix_shape_type& matLoc = this->locRes( q );
                matLoc( 0, 0 ) = factorWithDetFv * ( Fv11 - factorOther * InvFv11 );
                matLoc( 1, 0 ) = factorWithDetFv * ( Fv21 - factorOther * InvFv12 );
                matLoc( 2, 0 ) = factorWithDetFv * ( Fv31 - factorOther * InvFv13 );
                matLoc( 0, 1 ) = factorWithDetFv * ( Fv12 - factorOther * InvFv21 );
                matLoc( 1, 1 ) = factorWithDetFv * ( Fv22 - factorOther * InvFv22 );
                matLoc( 2, 1 ) = factorWithDetFv * ( Fv32 - factorOther * InvFv23 );
                matLoc( 0, 2 ) = factorWithDetFv * ( Fv13 - factorOther * InvFv31 );
                matLoc( 1, 2 ) = factorWithDetFv * ( Fv23 - factorOther * InvFv32 );
                matLoc( 2, 2 ) = factorWithDetFv * ( Fv33 - factorOther * InvFv33 );
                if ( !useDispPresForm )
                    matLoc += M_tensorVolumicPart->locRes( q );
            }
            else
            {
                const value_type factorWithDetFv1 = 1. / math::pow( detFv, 2 );
                const value_type factorWithDetFv2 = -( 2. / 3. ) * coefflame2 * math::pow( detFv, -5. / 3. );
                const value_type factorWithDetFv3 = coefflame2 * math::pow( detFv, -2. / 3. );
                M_locEvalPrecomputeInvDetF[q]( 0, 0 ) = 1. / detFv;
                M_locEvalPrecomputeDetF1[q]( 0, 0 ) = factorWithDetFv1;
                //M_locEvalPrecomputeDetF2[q](0,0) = factorWithDetFv2;
                M_locEvalPrecomputeDetF3[q]( 0, 0 ) = factorWithDetFv3;
                M_locEvalPrecomputeTraceC[q]( 0, 0 ) = traceCv;
#if 1 // NEW
                const value_type factorTrialDetF = ( 1. / 3. ) * traceCv * ( 1. / detFv );
                typename super_type::loc_tensor2_type& matLocTrialDetF = M_locEvalPrecomputeTrialDetF[q];
                matLocTrialDetF( 0, 0 ) = factorWithDetFv2 * ( Fv11 - factorTrialDetF * InvFv11 );
                matLocTrialDetF( 0, 1 ) = factorWithDetFv2 * ( Fv12 - factorTrialDetF * InvFv21 );
                matLocTrialDetF( 0, 2 ) = factorWithDetFv2 * ( Fv13 - factorTrialDetF * InvFv31 );
                matLocTrialDetF( 1, 0 ) = factorWithDetFv2 * ( Fv21 - factorTrialDetF * InvFv12 );
                matLocTrialDetF( 1, 1 ) = factorWithDetFv2 * ( Fv22 - factorTrialDetF * InvFv22 );
                matLocTrialDetF( 1, 2 ) = factorWithDetFv2 * ( Fv23 - factorTrialDetF * InvFv32 );
                matLocTrialDetF( 2, 0 ) = factorWithDetFv2 * ( Fv31 - factorTrialDetF * InvFv13 );
                matLocTrialDetF( 2, 1 ) = factorWithDetFv2 * ( Fv32 - factorTrialDetF * InvFv23 );
                matLocTrialDetF( 2, 2 ) = factorWithDetFv2 * ( Fv33 - factorTrialDetF * InvFv33 );
#endif
            }
        }
    }

    using super_type::evalijq; // fix clang warning

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i, j, q, mpl::int_<expr_type::specific_expr_type::value>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::EVAL> /**/ ) const
    {
        CHECK( false ) << "not allow";
        return this->locMatrixShape();
    }

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/ ) const
    {
        return evalijq( i, j, q, mpl::int_<ExprApplyType::JACOBIAN>(), mpl::int_<expr_type::nRealDim>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<2> /**/ ) const
    {
        //const value_type coefflame2 = M_locEvalFieldCoefflame2[q](0,0);

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 );

        const value_type ddetF = dF11 * Fv22 + Fv11 * dF22 - Fv21 * dF12 - dF21 * Fv12;
        const value_type dtraceC = 2 * Fv11 * dF11 + 2 * Fv21 * dF21 + 2 * Fv12 * dF12 + 2 * Fv22 * dF22;

        const value_type traceCv = M_locEvalPrecomputeTraceC[q]( 0, 0 );
        // 1./detFv
        const value_type precomputeInvDetF = M_locEvalPrecomputeInvDetF[q]( 0, 0 );
        // 1./math::pow(detFv,2);
        const value_type precomputeDetF1 = M_locEvalPrecomputeDetF1[q]( 0, 0 );

        // F^{-T} trial
        const value_type dFmt11 = -ddetF * precomputeDetF1 * Fv22 + precomputeInvDetF * dF22;
        const value_type dFmt12 = ddetF * precomputeDetF1 * Fv21 - precomputeInvDetF * dF21;
        const value_type dFmt21 = ddetF * precomputeDetF1 * Fv12 - precomputeInvDetF * dF12;
        const value_type dFmt22 = -ddetF * precomputeDetF1 * Fv11 + precomputeInvDetF * dF11;
#if 0 // NEW
            const value_type precomputeDetF2 = M_locEvalPrecomputeDetF2[q](0,0);
            const value_type dFS_neohookean_a11 = precomputeDetF2*ddetF*( Fv11 - (1./2.)*traceCv*precomputeInvDetF*Fv22 );
            const value_type dFS_neohookean_a12 = precomputeDetF2*ddetF*( Fv12 + (1./2.)*traceCv*precomputeInvDetF*Fv21 );
            const value_type dFS_neohookean_a21 = precomputeDetF2*ddetF*( Fv21 + (1./2.)*traceCv*precomputeInvDetF*Fv12 );
            const value_type dFS_neohookean_a22 = precomputeDetF2*ddetF*( Fv22 - (1./2.)*traceCv*precomputeInvDetF*Fv11 );

            const value_type precomputeDetF3 = M_locEvalPrecomputeDetF3[q](0,0);
            const value_type dFS_neohookean_b11 = precomputeDetF3*( dF11 - (1./2.)*(dtraceC*precomputeInvDetF*Fv22+traceCv*dFmt11) );
            const value_type dFS_neohookean_b12 = precomputeDetF3*( dF12 - (1./2.)*(-dtraceC*precomputeInvDetF*Fv21+traceCv*dFmt12) );
            const value_type dFS_neohookean_b21 = precomputeDetF3*( dF21 - (1./2.)*(-dtraceC*precomputeInvDetF*Fv12+traceCv*dFmt21) );
            const value_type dFS_neohookean_b22 = precomputeDetF3*( dF22 - (1./2.)*(dtraceC*precomputeInvDetF*Fv11+traceCv*dFmt22) );

            typename super_type::matrix_shape_type & matLoc = this->locMatrixShape();
            matLoc(0,0) = dFS_neohookean_a11+dFS_neohookean_b11;
            matLoc(1,0) = dFS_neohookean_a21+dFS_neohookean_b21;
            matLoc(0,1) = dFS_neohookean_a12+dFS_neohookean_b12;
            matLoc(1,1) = dFS_neohookean_a22+dFS_neohookean_b22;
#else
        typename super_type::matrix_shape_type& matLoc = this->locMatrixShape();
        matLoc = ddetF * M_locEvalPrecomputeTrialDetF[q];
        const value_type precomputeDetF3 = M_locEvalPrecomputeDetF3[q]( 0, 0 );
        matLoc( 0, 0 ) += precomputeDetF3 * ( dF11 - ( 1. / 2. ) * ( dtraceC * precomputeInvDetF * Fv22 + traceCv * dFmt11 ) );
        matLoc( 0, 1 ) += precomputeDetF3 * ( dF12 - ( 1. / 2. ) * ( -dtraceC * precomputeInvDetF * Fv21 + traceCv * dFmt12 ) );
        matLoc( 1, 0 ) += precomputeDetF3 * ( dF21 - ( 1. / 2. ) * ( -dtraceC * precomputeInvDetF * Fv12 + traceCv * dFmt21 ) );
        matLoc( 1, 1 ) += precomputeDetF3 * ( dF22 - ( 1. / 2. ) * ( dtraceC * precomputeInvDetF * Fv11 + traceCv * dFmt22 ) );
#endif
        if ( !useDispPresForm )
        {
            typename super_type::loc_tensor2_type dFmt;
            dFmt( 0, 0 ) = dFmt11;
            dFmt( 0, 1 ) = dFmt12;
            dFmt( 1, 0 ) = dFmt21;
            dFmt( 1, 1 ) = dFmt22;
            const value_type detFv = 1. / precomputeInvDetF;
            matLoc += M_tensorVolumicPart->evalijq( i, j, q, detFv, ddetF, dFmt );
        }
        return this->locMatrixShape();
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<3> /**/ ) const
    {
        //const value_type coefflame2 = M_locEvalFieldCoefflame2[q](0,0);

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 ), Fv13 = gradDisplacementEval( 0, 2 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 ), Fv23 = gradDisplacementEval( 1, 2 );
        const value_type Fv31 = gradDisplacementEval( 2, 0 ), Fv32 = gradDisplacementEval( 2, 1 ), Fv33 = 1. + gradDisplacementEval( 2, 2 );

        const value_type Fav11 = gradDisplacementEval( 0, 0 ), Fav12 = Fv12, Fav13 = Fv13;
        const value_type Fav21 = Fv21, Fav22 = gradDisplacementEval( 1, 1 ), Fav23 = Fv23;
        const value_type Fav31 = Fv31, Fav32 = Fv32, Fav33 = gradDisplacementEval( 2, 2 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 ), dF13 /*du1tdz*/ = gradTrial( 0, 2 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 ), dF23 /*du2tdz*/ = gradTrial( 1, 2 );
        const value_type dF31 /*du3tdx*/ = gradTrial( 2, 0 ), dF32 /*du3tdy*/ = gradTrial( 2, 1 ), dF33 /*du3tdz*/ = gradTrial( 2, 2 );

        /*auto InvFvMat = mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
              Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
              Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21 );*/
        // InvFv ( without 1/detF )
        const value_type InvFv11 = Fv22 * Fv33 - Fv23 * Fv32;
        const value_type InvFv12 = Fv13 * Fv32 - Fv12 * Fv33;
        const value_type InvFv13 = Fv12 * Fv23 - Fv13 * Fv22;
        const value_type InvFv21 = Fv23 * Fv31 - Fv21 * Fv33;
        const value_type InvFv22 = Fv11 * Fv33 - Fv13 * Fv31;
        const value_type InvFv23 = Fv13 * Fv21 - Fv11 * Fv23;
        const value_type InvFv31 = Fv21 * Fv32 - Fv22 * Fv31;
        const value_type InvFv32 = Fv12 * Fv31 - Fv11 * Fv32;
        const value_type InvFv33 = Fv11 * Fv22 - Fv12 * Fv21;

        const value_type FmtNL11 = dF22 + dF33 + dF22 * Fav33 + Fav22 * dF33 - dF23 * Fav32 - Fav23 * dF32;
        const value_type FmtNL12 = dF23 * Fav31 + Fav23 * dF31 - dF21 - dF21 * Fav33 - Fav21 * dF33;
        const value_type FmtNL13 = dF21 * Fav32 + Fav21 * dF32 - dF31 - dF31 * Fav22 - Fav31 * dF22;
        const value_type FmtNL21 = dF13 * Fav32 + Fav13 * dF32 - dF12 - dF12 * Fav33 - Fav12 * dF33;
        const value_type FmtNL22 = dF11 + dF33 + dF11 * Fav33 + Fav11 * dF33 - dF13 * Fav31 - Fav13 * dF31;
        const value_type FmtNL23 = dF12 * Fav31 + Fav12 * dF31 - dF32 - dF32 * Fav11 - Fav32 * dF11;
        const value_type FmtNL31 = dF12 * Fav23 + Fav12 * dF23 - dF13 - dF13 * Fav22 - Fav13 * dF22;
        const value_type FmtNL32 = dF13 * Fav21 + Fav13 * dF21 - dF23 - dF23 * Fav11 - Fav23 * dF11;
        const value_type FmtNL33 = dF11 + dF22 + dF11 * Fav22 + Fav11 * dF22 - dF12 * Fav21 - Fav12 * dF21;

        const value_type ddetF = dF11 + dF22 + dF33 + ( dF22 * Fav33 + Fav22 * dF33 ) + ( dF11 * Fav22 + Fav11 * dF22 ) + ( dF11 * Fav33 + Fav11 * dF33 ) + dF11 * Fav22 * Fav33 + Fav11 * dF22 * Fav33 + Fav11 * Fav22 * dF33 + dF12 * Fav23 * Fav31 + Fav12 * dF23 * Fav31 + Fav12 * Fav23 * dF31 + dF13 * Fav21 * Fav32 + Fav13 * dF21 * Fav32 + Fav13 * Fav21 * dF32 - ( dF23 * Fav32 + Fav23 * dF32 ) - ( dF12 * Fav21 + Fav12 * dF21 ) - ( dF13 * Fav31 + Fav13 * dF31 ) - ( dF11 * Fav23 * Fav32 + Fav11 * dF23 * Fav32 + Fav11 * Fav23 * dF32 ) - ( dF12 * Fav21 * Fav33 + Fav12 * dF21 * Fav33 + Fav12 * Fav21 * dF33 ) - ( dF13 * Fav22 * Fav31 + Fav13 * dF22 * Fav31 + Fav13 * Fav22 * dF31 );

        const value_type dtraceC = 2 * Fv11 * dF11 + 2 * Fv21 * dF21 + 2 * Fv31 * dF31 + 2 * Fv12 * dF12 + 2 * Fv22 * dF22 + 2 * Fv32 * dF32 + 2 * Fv13 * dF13 + 2 * Fv23 * dF23 + 2 * Fv33 * dF33;

        /*auto InvFtMat =  mat<3,3>(FmtNL11,FmtNL12,FmtNL13,
                                  FmtNL21,FmtNL22,FmtNL23,
                                  FmtNL31,FmtNL32,FmtNL33 );*/
        //auto dFmt = -(ddetF/pow(detFv,2))*trans(InvFvMat)+(cst(1.)/detFv)*/*trans*/(InvFtMat) ;

        //auto dFS_neohookean_a = -(2./3.)*idv(CoeffLame2)*ddetF*pow(detFv,-5./3.) *( Fv - (1./3.)*traceCv*trans(InvFv) );
        //auto dFS_neohookean_b = idv(CoeffLame2)*pow(detFv,-2./3.)*( dF - (1./3.)*(dtraceC*trans(InvFv)+traceCv*dFmt) );
        //auto dFS_neohookean = dFS_neohookean_a+dFS_neohookean_b;

        const value_type traceCv = M_locEvalPrecomputeTraceC[q]( 0, 0 );
        const value_type precomputeInvDetF = M_locEvalPrecomputeInvDetF[q]( 0, 0 ); // 1./detFv;
        const value_type precomputeDetF1 = M_locEvalPrecomputeDetF1[q]( 0, 0 );     // 1./math::pow(detFv,2);
        const value_type dFmt11 = -ddetF * precomputeDetF1 * InvFv11 + precomputeInvDetF * FmtNL11;
        const value_type dFmt12 = -ddetF * precomputeDetF1 * InvFv21 + precomputeInvDetF * FmtNL12;
        const value_type dFmt13 = -ddetF * precomputeDetF1 * InvFv31 + precomputeInvDetF * FmtNL13;
        const value_type dFmt21 = -ddetF * precomputeDetF1 * InvFv12 + precomputeInvDetF * FmtNL21;
        const value_type dFmt22 = -ddetF * precomputeDetF1 * InvFv22 + precomputeInvDetF * FmtNL22;
        const value_type dFmt23 = -ddetF * precomputeDetF1 * InvFv32 + precomputeInvDetF * FmtNL23;
        const value_type dFmt31 = -ddetF * precomputeDetF1 * InvFv13 + precomputeInvDetF * FmtNL31;
        const value_type dFmt32 = -ddetF * precomputeDetF1 * InvFv23 + precomputeInvDetF * FmtNL32;
        const value_type dFmt33 = -ddetF * precomputeDetF1 * InvFv33 + precomputeInvDetF * FmtNL33;
#if 0 // NEW
            const value_type precomputeDetF2 = M_locEvalPrecomputeDetF2[q](0,0);
            const value_type dFS_neohookean_a11 = precomputeDetF2*ddetF*( Fv11 - (1./3.)*traceCv*precomputeInvDetF*InvFv11 );
            const value_type dFS_neohookean_a12 = precomputeDetF2*ddetF*( Fv12 - (1./3.)*traceCv*precomputeInvDetF*InvFv21 );
            const value_type dFS_neohookean_a13 = precomputeDetF2*ddetF*( Fv13 - (1./3.)*traceCv*precomputeInvDetF*InvFv31 );
            const value_type dFS_neohookean_a21 = precomputeDetF2*ddetF*( Fv21 - (1./3.)*traceCv*precomputeInvDetF*InvFv12 );
            const value_type dFS_neohookean_a22 = precomputeDetF2*ddetF*( Fv22 - (1./3.)*traceCv*precomputeInvDetF*InvFv22 );
            const value_type dFS_neohookean_a23 = precomputeDetF2*ddetF*( Fv23 - (1./3.)*traceCv*precomputeInvDetF*InvFv32 );
            const value_type dFS_neohookean_a31 = precomputeDetF2*ddetF*( Fv31 - (1./3.)*traceCv*precomputeInvDetF*InvFv13 );
            const value_type dFS_neohookean_a32 = precomputeDetF2*ddetF*( Fv32 - (1./3.)*traceCv*precomputeInvDetF*InvFv23 );
            const value_type dFS_neohookean_a33 = precomputeDetF2*ddetF*( Fv33 - (1./3.)*traceCv*precomputeInvDetF*InvFv33 );

            const value_type precomputeDetF3 = M_locEvalPrecomputeDetF3[q](0,0);
            const value_type dFS_neohookean_b11 = precomputeDetF3*( dF11 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv11+traceCv*dFmt11) );
            const value_type dFS_neohookean_b12 = precomputeDetF3*( dF12 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv21+traceCv*dFmt12) );
            const value_type dFS_neohookean_b13 = precomputeDetF3*( dF13 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv31+traceCv*dFmt13) );
            const value_type dFS_neohookean_b21 = precomputeDetF3*( dF21 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv12+traceCv*dFmt21) );
            const value_type dFS_neohookean_b22 = precomputeDetF3*( dF22 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv22+traceCv*dFmt22) );
            const value_type dFS_neohookean_b23 = precomputeDetF3*( dF23 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv32+traceCv*dFmt23) );
            const value_type dFS_neohookean_b31 = precomputeDetF3*( dF31 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv13+traceCv*dFmt31) );
            const value_type dFS_neohookean_b32 = precomputeDetF3*( dF32 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv23+traceCv*dFmt32) );
            const value_type dFS_neohookean_b33 = precomputeDetF3*( dF33 - (1./3.)*(dtraceC*precomputeInvDetF*InvFv33+traceCv*dFmt33) );

            typename super_type::matrix_shape_type & matLoc = this->locMatrixShape();
            matLoc(0,0) = dFS_neohookean_a11 + dFS_neohookean_b11;
            matLoc(1,0) = dFS_neohookean_a21 + dFS_neohookean_b21;
            matLoc(2,0) = dFS_neohookean_a31 + dFS_neohookean_b31;
            matLoc(0,1) = dFS_neohookean_a12 + dFS_neohookean_b12;
            matLoc(1,1) = dFS_neohookean_a22 + dFS_neohookean_b22;
            matLoc(2,1) = dFS_neohookean_a32 + dFS_neohookean_b32;
            matLoc(0,2) = dFS_neohookean_a13 + dFS_neohookean_b13;
            matLoc(1,2) = dFS_neohookean_a23 + dFS_neohookean_b23;
            matLoc(2,2) = dFS_neohookean_a33 + dFS_neohookean_b33;
#else
        typename super_type::matrix_shape_type& matLoc = this->locMatrixShape();
        matLoc = ddetF * M_locEvalPrecomputeTrialDetF[q];
        const value_type precomputeDetF3 = M_locEvalPrecomputeDetF3[q]( 0, 0 );
        matLoc( 0, 0 ) += precomputeDetF3 * ( dF11 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv11 + traceCv * dFmt11 ) );
        matLoc( 0, 1 ) += precomputeDetF3 * ( dF12 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv21 + traceCv * dFmt12 ) );
        matLoc( 0, 2 ) += precomputeDetF3 * ( dF13 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv31 + traceCv * dFmt13 ) );
        matLoc( 1, 0 ) += precomputeDetF3 * ( dF21 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv12 + traceCv * dFmt21 ) );
        matLoc( 1, 1 ) += precomputeDetF3 * ( dF22 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv22 + traceCv * dFmt22 ) );
        matLoc( 1, 2 ) += precomputeDetF3 * ( dF23 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv32 + traceCv * dFmt23 ) );
        matLoc( 2, 0 ) += precomputeDetF3 * ( dF31 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv13 + traceCv * dFmt31 ) );
        matLoc( 2, 1 ) += precomputeDetF3 * ( dF32 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv23 + traceCv * dFmt32 ) );
        matLoc( 2, 2 ) += precomputeDetF3 * ( dF33 - ( 1. / 3. ) * ( dtraceC * precomputeInvDetF * InvFv33 + traceCv * dFmt33 ) );
#endif
        if ( !useDispPresForm )
        {
            typename super_type::loc_tensor2_type dFmt;
            dFmt( 0, 0 ) = dFmt11;
            dFmt( 0, 1 ) = dFmt12;
            dFmt( 0, 2 ) = dFmt13;
            dFmt( 1, 0 ) = dFmt21;
            dFmt( 1, 1 ) = dFmt22;
            dFmt( 1, 2 ) = dFmt23;
            dFmt( 2, 0 ) = dFmt31;
            dFmt( 2, 1 ) = dFmt32;
            dFmt( 2, 2 ) = dFmt33;
            const value_type detFv = 1. / precomputeInvDetF;
            matLoc += M_tensorVolumicPart->evalijq( i, j, q, detFv, ddetF, dFmt );
        }

        return this->locMatrixShape();
    }

  private:
    pc_mechprop_scalar_ptrtype M_pcMechPropField;
    ctx_mechprop_scalar_ptrtype M_ctxMechPropField;

    typename super_type::array_scalar_type M_locEvalFieldCoefflame2;
    typename super_type::array_scalar_type M_locEvalPrecomputeDetF1 /*,M_locEvalPrecomputeDetF2*/, M_locEvalPrecomputeDetF3, M_locEvalPrecomputeInvDetF, M_locEvalPrecomputeTraceC;
    typename super_type::array_tensor2_type M_locEvalPrecomputeTrialDetF;

    boost::shared_ptr<tensor_volumic_part_type> M_tensorVolumicPart;
};

template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType, bool useDispPresForm>
struct tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory : public tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;

  public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::value_type value_type;

    // fe mech prop context
    typedef typename expr_type::fe_mechprop_scalar_type fe_mechprop_scalar_type;
    static const size_type context_mechprop_scalar = expr_type::context_mechprop_scalar;
    typedef typename fe_mechprop_scalar_type::PreCompute pc_mechprop_scalar_type;
    typedef boost::shared_ptr<pc_mechprop_scalar_type> pc_mechprop_scalar_ptrtype;
    typedef typename fe_mechprop_scalar_type::template Context<context_mechprop_scalar, fe_mechprop_scalar_type, gm_type, geoelement_type, gmc_type::context> ctx_mechprop_scalar_type;
    typedef boost::shared_ptr<ctx_mechprop_scalar_type> ctx_mechprop_scalar_ptrtype;

    //
    typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> tensor_volumic_part_type;
    //typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartClassic<Geo_t,Basis_i_t,Basis_j_t,ExprType> tensor_volumic_part_classic_type;
    //typedef tensorFirstPiolaKirchhoffCompressibleVolumicPartSimo1985<Geo_t,Basis_i_t,Basis_j_t,ExprType> tensor_volumic_part_simo1985_type;
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom, fev, feu );
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom, fev );
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
        if ( !useDispPresForm )
            M_tensorVolumicPart = tensor_volumic_part_type::New( this->expr().mechanicalPropertiesDesc().decouplingEnergyVolumicLaw(), expr, geom );
    }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        updateImpl( geom );
    }
    void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcMechPropField->update( this->gmc()->pc()->nodes() );
        updateImpl( geom );
    }
    void updateImpl( Geo_t const& geom )
    {
        M_ctxMechPropField->update( this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField );
        std::fill( M_locEvalFieldCoefflame2.data(), M_locEvalFieldCoefflame2.data() + M_locEvalFieldCoefflame2.num_elements(), super_type::loc_scalar_type::Zero() );
        this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().id( *M_ctxMechPropField, M_locEvalFieldCoefflame2 );
        //std::fill( this->locRes().data(), this->locRes().data()+this->locRes().num_elements(), super_type::loc_res_type::Zero() );

        if ( !useDispPresForm )
            M_tensorVolumicPart->updateImpl( *M_ctxMechPropField, this->locGradDisplacement() );

        updateImpl();
    }
    void updateImpl()
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );
            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            auto Id = super_type::loc_tensor2_type::Identity();
            auto F = Id + gradDisplacementEval;
            const value_type detFv = F.determinant();
            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                auto InvFvmt = F.inverse().transpose();
                this->locRes( q ) = coefflame2 * ( F - InvFvmt );
                if ( !useDispPresForm )
                    this->locRes( q ) += M_tensorVolumicPart->locRes( q );
            }
            else
            {
                M_locEvalPrecomputePowDetF[q]( 0, 0 ) = math::pow( detFv, 2 );
                M_locEvalPrecomputeDetF[q]( 0, 0 ) = detFv;
            }
        }
    }

    using super_type::evalijq; // fix clang warning

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i, j, q, mpl::int_<expr_type::specific_expr_type::value>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::EVAL> /**/ ) const
    {
        CHECK( false ) << "not allow";
        return this->locMatrixShape();
    }

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/ ) const
    {
        return evalijq( i, j, q, mpl::int_<ExprApplyType::JACOBIAN>(), mpl::int_<expr_type::nRealDim>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<2> /**/ ) const
    {
        const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 );

        const value_type ddetF = dF11 * Fv22 + Fv11 * dF22 - Fv21 * dF12 - dF21 * Fv12;
        const value_type detFv = M_locEvalPrecomputeDetF[q]( 0, 0 );

        const value_type scaleUsedWithInvFv = 1. / detFv;
        const value_type precomputePowDetF = M_locEvalPrecomputePowDetF[q]( 0, 0 );

        const value_type factorOther = ddetF / precomputePowDetF;

        typename super_type::loc_tensor2_type dFmt;
        dFmt( 0, 0 ) = -factorOther * Fv22 + scaleUsedWithInvFv * dF22;
        dFmt( 0, 1 ) = factorOther * Fv21 - scaleUsedWithInvFv * dF21;
        dFmt( 1, 0 ) = factorOther * Fv12 - scaleUsedWithInvFv * dF12;
        dFmt( 1, 1 ) = -factorOther * Fv11 + scaleUsedWithInvFv * dF11;

        //typename super_type::matrix_shape_type& locMat = this->locMatrixShape();
        this->locMatrixShape() = coefflame2 * ( gradTrial - dFmt );
        if ( !useDispPresForm )
            this->locMatrixShape() += M_tensorVolumicPart->evalijq( i, j, q, detFv, ddetF, dFmt );

        return this->locMatrixShape();
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<3> /**/ ) const
    {
        const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 ), Fv13 = gradDisplacementEval( 0, 2 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 ), Fv23 = gradDisplacementEval( 1, 2 );
        const value_type Fv31 = gradDisplacementEval( 2, 0 ), Fv32 = gradDisplacementEval( 2, 1 ), Fv33 = 1. + gradDisplacementEval( 2, 2 );

        const value_type Fav11 = gradDisplacementEval( 0, 0 ), Fav12 = Fv12, Fav13 = Fv13;
        const value_type Fav21 = Fv21, Fav22 = gradDisplacementEval( 1, 1 ), Fav23 = Fv23;
        const value_type Fav31 = Fv31, Fav32 = Fv32, Fav33 = gradDisplacementEval( 2, 2 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 ), dF13 /*du1tdz*/ = gradTrial( 0, 2 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 ), dF23 /*du2tdz*/ = gradTrial( 1, 2 );
        const value_type dF31 /*du3tdx*/ = gradTrial( 2, 0 ), dF32 /*du3tdy*/ = gradTrial( 2, 1 ), dF33 /*du3tdz*/ = gradTrial( 2, 2 );

        const value_type ddetF = dF11 + dF22 + dF33 + ( dF22 * Fav33 + Fav22 * dF33 ) + ( dF11 * Fav22 + Fav11 * dF22 ) + ( dF11 * Fav33 + Fav11 * dF33 ) + dF11 * Fav22 * Fav33 + Fav11 * dF22 * Fav33 + Fav11 * Fav22 * dF33 + dF12 * Fav23 * Fav31 + Fav12 * dF23 * Fav31 + Fav12 * Fav23 * dF31 + dF13 * Fav21 * Fav32 + Fav13 * dF21 * Fav32 + Fav13 * Fav21 * dF32 - ( dF23 * Fav32 + Fav23 * dF32 ) - ( dF12 * Fav21 + Fav12 * dF21 ) - ( dF13 * Fav31 + Fav13 * dF31 ) - ( dF11 * Fav23 * Fav32 + Fav11 * dF23 * Fav32 + Fav11 * Fav23 * dF32 ) - ( dF12 * Fav21 * Fav33 + Fav12 * dF21 * Fav33 + Fav12 * Fav21 * dF33 ) - ( dF13 * Fav22 * Fav31 + Fav13 * dF22 * Fav31 + Fav13 * Fav22 * dF31 );

        /*auto InvFvMat = mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
              Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
              Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21 );*/
        // InvFv ( without 1/detF )
        const value_type InvFv11 = Fv22 * Fv33 - Fv23 * Fv32;
        const value_type InvFv12 = Fv13 * Fv32 - Fv12 * Fv33;
        const value_type InvFv13 = Fv12 * Fv23 - Fv13 * Fv22;
        const value_type InvFv21 = Fv23 * Fv31 - Fv21 * Fv33;
        const value_type InvFv22 = Fv11 * Fv33 - Fv13 * Fv31;
        const value_type InvFv23 = Fv13 * Fv21 - Fv11 * Fv23;
        const value_type InvFv31 = Fv21 * Fv32 - Fv22 * Fv31;
        const value_type InvFv32 = Fv12 * Fv31 - Fv11 * Fv32;
        const value_type InvFv33 = Fv11 * Fv22 - Fv12 * Fv21;

        const value_type FmtNL11 = dF22 + dF33 + dF22 * Fav33 + Fav22 * dF33 - dF23 * Fav32 - Fav23 * dF32;
        const value_type FmtNL12 = dF23 * Fav31 + Fav23 * dF31 - dF21 - dF21 * Fav33 - Fav21 * dF33;
        const value_type FmtNL13 = dF21 * Fav32 + Fav21 * dF32 - dF31 - dF31 * Fav22 - Fav31 * dF22;
        const value_type FmtNL21 = dF13 * Fav32 + Fav13 * dF32 - dF12 - dF12 * Fav33 - Fav12 * dF33;
        const value_type FmtNL22 = dF11 + dF33 + dF11 * Fav33 + Fav11 * dF33 - dF13 * Fav31 - Fav13 * dF31;
        const value_type FmtNL23 = dF12 * Fav31 + Fav12 * dF31 - dF32 - dF32 * Fav11 - Fav32 * dF11;
        const value_type FmtNL31 = dF12 * Fav23 + Fav12 * dF23 - dF13 - dF13 * Fav22 - Fav13 * dF22;
        const value_type FmtNL32 = dF13 * Fav21 + Fav13 * dF21 - dF23 - dF23 * Fav11 - Fav23 * dF11;
        const value_type FmtNL33 = dF11 + dF22 + dF11 * Fav22 + Fav11 * dF22 - dF12 * Fav21 - Fav12 * dF21;

        /*auto InvFtMat =  mat<3,3>(FmtNL11,FmtNL12,FmtNL13,
                                  FmtNL21,FmtNL22,FmtNL23,
                                  FmtNL31,FmtNL32,FmtNL33 );*/
        //auto dFmt = -(ddetF/pow(detFv,2))*trans(InvFvMat)+(cst(1.)/detFv)*/*trans*/(InvFtMat) ;

        const value_type detFv = M_locEvalPrecomputeDetF[q]( 0, 0 );
        const value_type scaleUsedWithInvFv = 1. / detFv;
        const value_type precomputePowDetF = M_locEvalPrecomputePowDetF[q]( 0, 0 ); //pow(detFv,2)  //idv(CoeffLame1)/pow(detFv,2)

        const value_type factorOther = ddetF / precomputePowDetF; //math::pow(detFv,2);

        typename super_type::loc_tensor2_type dFmt;
        dFmt( 0, 0 ) = -factorOther * InvFv11 + scaleUsedWithInvFv * FmtNL11;
        dFmt( 0, 1 ) = -factorOther * InvFv21 + scaleUsedWithInvFv * FmtNL12;
        dFmt( 0, 2 ) = -factorOther * InvFv31 + scaleUsedWithInvFv * FmtNL13;
        dFmt( 1, 0 ) = -factorOther * InvFv12 + scaleUsedWithInvFv * FmtNL21;
        dFmt( 1, 1 ) = -factorOther * InvFv22 + scaleUsedWithInvFv * FmtNL22;
        dFmt( 1, 2 ) = -factorOther * InvFv32 + scaleUsedWithInvFv * FmtNL23;
        dFmt( 2, 0 ) = -factorOther * InvFv13 + scaleUsedWithInvFv * FmtNL31;
        dFmt( 2, 1 ) = -factorOther * InvFv23 + scaleUsedWithInvFv * FmtNL32;
        dFmt( 2, 2 ) = -factorOther * InvFv33 + scaleUsedWithInvFv * FmtNL33;

        this->locMatrixShape() = coefflame2 * ( gradTrial - dFmt );
        if ( !useDispPresForm )
            this->locMatrixShape() += M_tensorVolumicPart->evalijq( i, j, q, detFv, ddetF, dFmt );

        return this->locMatrixShape();
    }

  private:
    pc_mechprop_scalar_ptrtype M_pcMechPropField;
    ctx_mechprop_scalar_ptrtype M_ctxMechPropField;
    typename super_type::array_scalar_type M_locEvalFieldCoefflame2;
    typename super_type::array_scalar_type M_locEvalPrecomputePowDetF, M_locEvalPrecomputeDetF;
    boost::shared_ptr<tensor_volumic_part_type> M_tensorVolumicPart;
};

/**
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     * NeoHookeanCompressible :
     *   S = \mu*(Id-C^{-1}) + \kappa*ln(J)*C^{-1}
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     */
template <typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ExprType>
struct tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985 : public tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType>
{
    typedef tensorFirstPiolaKirchhoffBase<Geo_t, Basis_i_t, Basis_j_t, ExprType> super_type;

  public:
    typedef ExprType expr_type;
    typedef typename super_type::geoelement_type geoelement_type;
    typedef typename super_type::gmc_type gmc_type;
    typedef typename super_type::gm_type gm_type;
    typedef typename super_type::value_type value_type;

    // fe mech prop context
    typedef typename expr_type::fe_mechprop_scalar_type fe_mechprop_scalar_type;
    static const size_type context_mechprop_scalar = expr_type::context_mechprop_scalar;
    typedef typename fe_mechprop_scalar_type::PreCompute pc_mechprop_scalar_type;
    typedef boost::shared_ptr<pc_mechprop_scalar_type> pc_mechprop_scalar_ptrtype;
    typedef typename fe_mechprop_scalar_type::template Context<context_mechprop_scalar, fe_mechprop_scalar_type, gm_type, geoelement_type, gmc_type::context> ctx_mechprop_scalar_type;
    typedef boost::shared_ptr<ctx_mechprop_scalar_type> ctx_mechprop_scalar_ptrtype;

    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        : super_type( expr, geom, fev, feu ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalBulkModulus( this->expr().mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985( expr_type const& expr, Geo_t const& geom, Basis_i_t const& fev )
        : super_type( expr, geom, fev ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalBulkModulus( this->expr().mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }
    tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985( expr_type const& expr, Geo_t const& geom )
        : super_type( expr, geom ),
          M_pcMechPropField( new pc_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc()->xRefs() ) ),
          M_ctxMechPropField( new ctx_mechprop_scalar_type( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().functionSpace()->fe(), this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField ) ),
          M_locEvalBulkModulus( this->expr().mechanicalPropertiesDesc().fieldBulkModulus().idExtents( *this->gmc() ) ),
          M_locEvalFieldCoefflame2( this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().idExtents( *this->gmc() ) ),
          M_locEvalPrecomputeLogDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputePowDetF( boost::extents[this->gmc()->xRefs().size2()] ),
          M_locEvalPrecomputeDetF( boost::extents[this->gmc()->xRefs().size2()] )
    {
    }

    void update( Geo_t const& geom )
    {
        super_type::update( geom );
        updateImpl( geom );
    }
    void update( Geo_t const& geom, uint16_type face )
    {
        super_type::update( geom, face );

        if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
            M_pcMechPropField->update( this->gmc()->pc()->nodes() );
        updateImpl( geom );
    }
    void updateImpl( Geo_t const& geom )
    {
        M_ctxMechPropField->update( this->gmc(), (pc_mechprop_scalar_ptrtype const&)M_pcMechPropField );
        std::fill( M_locEvalBulkModulus.data(), M_locEvalBulkModulus.data() + M_locEvalBulkModulus.num_elements(), /*typename*/ super_type::loc_scalar_type::Zero() );
        std::fill( M_locEvalFieldCoefflame2.data(), M_locEvalFieldCoefflame2.data() + M_locEvalFieldCoefflame2.num_elements(), /*typename*/ super_type::loc_scalar_type::Zero() );
        this->expr().mechanicalPropertiesDesc().fieldBulkModulus().id( *M_ctxMechPropField, M_locEvalBulkModulus );
        this->expr().mechanicalPropertiesDesc().fieldCoeffLame2().id( *M_ctxMechPropField, M_locEvalFieldCoefflame2 );

        std::fill( this->locRes().data(), this->locRes().data() + this->locRes().num_elements(), super_type::matrix_shape_type /*loc_res_type*/ ::Zero() );
        updateImpl( mpl::int_<expr_type::nRealDim>() );
    }
    void updateImpl( mpl::int_<2> /**/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type bulkModulus = M_locEvalBulkModulus[q]( 0, 0 );
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 );
            const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 );
            const value_type detFv = Fv11 * Fv22 - Fv12 * Fv21;
            //const value_type traceCv = std::pow(Fv11,2) + std::pow(Fv21,2) + std::pow(Fv12,2) + std::pow(Fv22,2);

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                //auto InvFv = (cst(1.)/detFv)*mat<2,2>( Fv22,-Fv12,-Fv21,Fv11);
                const value_type scaleUsedWithInvFv = 1. / detFv;
                const value_type InvFv11 = scaleUsedWithInvFv * Fv22, InvFv12 = -scaleUsedWithInvFv * Fv12;
                const value_type InvFv21 = -scaleUsedWithInvFv * Fv21, InvFv22 = scaleUsedWithInvFv * Fv11;

                /*const value_type factorWithDetFv = coefflame2*std::pow(detFv,-2./3.);
                      const value_type factorWithTraceCv = (1./3.)*traceCv;
                      const value_type factorOther = factorWithDetFv*factorWithTraceCv*scaleUsedWithInvFv;*/

                //auto FSv_neohookean = idv(CoeffLame2)*(Fv-trans(InvFv)) + idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv))*trans(InvFv);
                const value_type factorWithLogDetFv = bulkModulus * math::log( detFv );
                this->locRes( q )( 0, 0 ) = coefflame2 * ( Fv11 - InvFv11 ) + factorWithLogDetFv * InvFv11;
                this->locRes( q )( 0, 1 ) = coefflame2 * ( Fv12 - InvFv21 ) + factorWithLogDetFv * InvFv21;
                this->locRes( q )( 1, 0 ) = coefflame2 * ( Fv21 - InvFv12 ) + factorWithLogDetFv * InvFv12;
                this->locRes( q )( 1, 1 ) = coefflame2 * ( Fv22 - InvFv22 ) + factorWithLogDetFv * InvFv22;
            }
            else
            {
                //const value_type precomputeLog = M_locEvalPrecomputeLog[q](0,0); // idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv)) -idv(CoeffLame2)
                //const value_type precomputePowDetF = M_locEvalPrecomputePowDetF[q](0,0);//idv(CoeffLame1)/pow(detFv,2)
                M_locEvalPrecomputeLogDetF[q]( 0, 0 ) = bulkModulus * math::log( detFv ) - coefflame2;
                M_locEvalPrecomputePowDetF[q]( 0, 0 ) = /*coefflame1/*/ math::pow( detFv, 2 );
                M_locEvalPrecomputeDetF[q]( 0, 0 ) = detFv;
            }
        }
    }

    void updateImpl( mpl::int_<3> /**/ )
    {
        for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
        {
            const value_type bulkModulus = M_locEvalBulkModulus[q]( 0, 0 );
            const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

            auto const& gradDisplacementEval = this->locGradDisplacement( q );
            const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 ), Fv13 = gradDisplacementEval( 0, 2 );
            const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 ), Fv23 = gradDisplacementEval( 1, 2 );
            const value_type Fv31 = gradDisplacementEval( 2, 0 ), Fv32 = gradDisplacementEval( 2, 1 ), Fv33 = 1. + gradDisplacementEval( 2, 2 );

            const value_type detFv = Fv11 * ( Fv22 * Fv33 - Fv23 * Fv32 ) - Fv21 * ( Fv12 * Fv33 - Fv13 * Fv32 ) + Fv31 * ( Fv12 * Fv23 - Fv13 * Fv22 );
            //const value_type traceCv = math::pow(Fv11,2) + math::pow(Fv21,2) + math::pow(Fv31,2) + math::pow(Fv12,2) +
            //    math::pow(Fv22,2) + math::pow(Fv32,2) + math::pow(Fv13,2) + math::pow(Fv23,2) + math::pow(Fv33,2);

            if ( expr_type::specific_expr_type::value == ExprApplyType::EVAL )
            {
                // InvFv
                /*auto InvFv = (cst(1.)/detFv)*mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
                      Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
                      Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21 );*/
                const value_type scaleUsedWithInvFv = 1. / detFv;
                const value_type InvFv11 = scaleUsedWithInvFv * ( Fv22 * Fv33 - Fv23 * Fv32 );
                const value_type InvFv12 = scaleUsedWithInvFv * ( Fv13 * Fv32 - Fv12 * Fv33 );
                const value_type InvFv13 = scaleUsedWithInvFv * ( Fv12 * Fv23 - Fv13 * Fv22 );
                const value_type InvFv21 = scaleUsedWithInvFv * ( Fv23 * Fv31 - Fv21 * Fv33 );
                const value_type InvFv22 = scaleUsedWithInvFv * ( Fv11 * Fv33 - Fv13 * Fv31 );
                const value_type InvFv23 = scaleUsedWithInvFv * ( Fv13 * Fv21 - Fv11 * Fv23 );
                const value_type InvFv31 = scaleUsedWithInvFv * ( Fv21 * Fv32 - Fv22 * Fv31 );
                const value_type InvFv32 = scaleUsedWithInvFv * ( Fv12 * Fv31 - Fv11 * Fv32 );
                const value_type InvFv33 = scaleUsedWithInvFv * ( Fv11 * Fv22 - Fv12 * Fv21 );

                //auto FSv_neohookean = idv(CoeffLame2)*(Fv-trans(InvFv)) + idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv))*trans(InvFv);
                const value_type factorWithLogDetFv = bulkModulus * math::log( detFv );
                this->locRes( q )( 0, 0 ) = coefflame2 * ( Fv11 - InvFv11 ) + factorWithLogDetFv * InvFv11;
                this->locRes( q )( 0, 1 ) = coefflame2 * ( Fv12 - InvFv21 ) + factorWithLogDetFv * InvFv21;
                this->locRes( q )( 0, 2 ) = coefflame2 * ( Fv13 - InvFv31 ) + factorWithLogDetFv * InvFv31;
                this->locRes( q )( 1, 0 ) = coefflame2 * ( Fv21 - InvFv12 ) + factorWithLogDetFv * InvFv12;
                this->locRes( q )( 1, 1 ) = coefflame2 * ( Fv22 - InvFv22 ) + factorWithLogDetFv * InvFv22;
                this->locRes( q )( 1, 2 ) = coefflame2 * ( Fv23 - InvFv32 ) + factorWithLogDetFv * InvFv32;
                this->locRes( q )( 2, 0 ) = coefflame2 * ( Fv31 - InvFv13 ) + factorWithLogDetFv * InvFv13;
                this->locRes( q )( 2, 1 ) = coefflame2 * ( Fv32 - InvFv23 ) + factorWithLogDetFv * InvFv23;
                this->locRes( q )( 2, 2 ) = coefflame2 * ( Fv33 - InvFv33 ) + factorWithLogDetFv * InvFv33;
            }
            else
            {
                M_locEvalPrecomputeLogDetF[q]( 0, 0 ) = bulkModulus * math::log( detFv ) - coefflame2;
                M_locEvalPrecomputePowDetF[q]( 0, 0 ) = /*coefflame1/*/ math::pow( detFv, 2 );
                M_locEvalPrecomputeDetF[q]( 0, 0 ) = detFv;
            }
        }
    }
    using super_type::evalijq; // fix clang warning

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q ) const
    {
        return evalijq( i, j, q, mpl::int_<expr_type::specific_expr_type::value>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::EVAL> /**/ ) const
    {
        CHECK( false ) << "not allow";
        return this->locMatrixShape();
    }

    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/ ) const
    {
        return evalijq( i, j, q, mpl::int_<ExprApplyType::JACOBIAN>(), mpl::int_<expr_type::nRealDim>() );
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<2> /**/ ) const
    {
        const value_type bulkModulus = M_locEvalBulkModulus[q]( 0, 0 );
        const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 );

        const value_type ddetF = dF11 * Fv22 + Fv11 * dF22 - Fv21 * dF12 - dF21 * Fv12;
        //const value_type ddetF = dF22+dF11+dF11*Fv22+Fv11*dF22 -Fv21*dF12 - dF21*Fv12;
        //const value_type dtraceC = 2*Fv11*dF11 + 2*Fv21*dF21 + 2*Fv12*dF12 + 2*Fv22*dF22;

        //auto InvFtMat = mat<2,2>( dF22,-dF21,-dF12,dF11);//already transpose
        //auto InvFvMat = mat<2,2>( Fv22,-Fv12,-Fv21,Fv11);
        //auto dFmt = -(ddetF/pow(detFv,2))*trans(InvFvMat)+(cst(1.)/detFv)*/*trans*/(InvFtMat);
        const value_type detFv = M_locEvalPrecomputeDetF[q]( 0, 0 );
        const value_type scaleUsedWithInvFv = 1. / detFv;
        const value_type precomputePowDetF = M_locEvalPrecomputePowDetF[q]( 0, 0 ); //pow(detFv,2)  //idv(CoeffLame1)/pow(detFv,2)

        const value_type factorOther = ddetF / precomputePowDetF; //math::pow(detFv,2);
        const value_type dFmt11 = -factorOther * Fv22 + scaleUsedWithInvFv * dF22;
        const value_type dFmt12 = factorOther * Fv21 - scaleUsedWithInvFv * dF21;
        const value_type dFmt21 = factorOther * Fv12 - scaleUsedWithInvFv * dF12;
        const value_type dFmt22 = -factorOther * Fv11 + scaleUsedWithInvFv * dF11;

        //idv(CoeffLame2)*dF
        const value_type dFS_neohookean_a11 = coefflame2 * dF11;
        const value_type dFS_neohookean_a12 = coefflame2 * dF12;
        const value_type dFS_neohookean_a21 = coefflame2 * dF21;
        const value_type dFS_neohookean_a22 = coefflame2 * dF22;

        // _expr= val( idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv)) -idv(CoeffLame2) )*trace(dFmt*trans(grad(v)) ),
        const value_type precomputeLog = M_locEvalPrecomputeLogDetF[q]( 0, 0 ); // idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv)) -idv(CoeffLame2)
        const value_type dFS_neohookean_b11 = precomputeLog * dFmt11;
        const value_type dFS_neohookean_b12 = precomputeLog * dFmt12;
        const value_type dFS_neohookean_b21 = precomputeLog * dFmt21;
        const value_type dFS_neohookean_b22 = precomputeLog * dFmt22;

        //_expr= ddetF*trace( val(idv(CoeffLame1)*trans(InvFvMat)/pow(detFv,2) )*trans(grad(v)) ),
        const value_type factorOther2 = ddetF * bulkModulus / precomputePowDetF; //        ddetF/precomputePowDetF;//math::pow(detFv,2);
        const value_type dFS_neohookean_c11 = factorOther2 * Fv22;
        const value_type dFS_neohookean_c12 = -factorOther2 * Fv21;
        const value_type dFS_neohookean_c21 = -factorOther2 * Fv12;
        const value_type dFS_neohookean_c22 = factorOther2 * Fv11;

        /*const value_type traceCv = M_locEvalPrecomputeTraceC[q](0,0);
            const value_type precomputeInvDetF = M_locEvalPrecomputeInvDetF[q](0,0);
            //const value_type precomputeInvDetF = 1./detFv;
            const value_type precomputeDetF1 = M_locEvalPrecomputeDetF1[q](0,0);*/

        this->locMatrixShape()( 0, 0 ) = dFS_neohookean_a11 + dFS_neohookean_b11 + dFS_neohookean_c11;
        this->locMatrixShape()( 1, 0 ) = dFS_neohookean_a21 + dFS_neohookean_b21 + dFS_neohookean_c21;
        this->locMatrixShape()( 0, 1 ) = dFS_neohookean_a12 + dFS_neohookean_b12 + dFS_neohookean_c12;
        this->locMatrixShape()( 1, 1 ) = dFS_neohookean_a22 + dFS_neohookean_b22 + dFS_neohookean_c22;

        return this->locMatrixShape();
    }
    Eigen::Matrix<value_type, super_type::shape::M, super_type::shape::N> const&
    evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> /**/, mpl::int_<3> /**/ ) const
    {
        const value_type bulkModulus = M_locEvalBulkModulus[q]( 0, 0 );
        const value_type coefflame2 = M_locEvalFieldCoefflame2[q]( 0, 0 );

        auto const& gradDisplacementEval = this->locGradDisplacement( q );
        const value_type Fv11 = 1. + gradDisplacementEval( 0, 0 ), Fv12 = gradDisplacementEval( 0, 1 ), Fv13 = gradDisplacementEval( 0, 2 );
        const value_type Fv21 = gradDisplacementEval( 1, 0 ), Fv22 = 1. + gradDisplacementEval( 1, 1 ), Fv23 = gradDisplacementEval( 1, 2 );
        const value_type Fv31 = gradDisplacementEval( 2, 0 ), Fv32 = gradDisplacementEval( 2, 1 ), Fv33 = 1. + gradDisplacementEval( 2, 2 );

        const value_type Fav11 = gradDisplacementEval( 0, 0 ), Fav12 = Fv12, Fav13 = Fv13;
        const value_type Fav21 = Fv21, Fav22 = gradDisplacementEval( 1, 1 ), Fav23 = Fv23;
        const value_type Fav31 = Fv31, Fav32 = Fv32, Fav33 = gradDisplacementEval( 2, 2 );

        auto const& gradTrial = this->fecTrial()->grad( j, q );
        const value_type dF11 /*du1tdx*/ = gradTrial( 0, 0 ), dF12 /*du1tdy*/ = gradTrial( 0, 1 ), dF13 /*du1tdz*/ = gradTrial( 0, 2 );
        const value_type dF21 /*du2tdx*/ = gradTrial( 1, 0 ), dF22 /*du2tdy*/ = gradTrial( 1, 1 ), dF23 /*du2tdz*/ = gradTrial( 1, 2 );
        const value_type dF31 /*du3tdx*/ = gradTrial( 2, 0 ), dF32 /*du3tdy*/ = gradTrial( 2, 1 ), dF33 /*du3tdz*/ = gradTrial( 2, 2 );

        const value_type ddetF = dF11 + dF22 + dF33 + ( dF22 * Fav33 + Fav22 * dF33 ) + ( dF11 * Fav22 + Fav11 * dF22 ) + ( dF11 * Fav33 + Fav11 * dF33 ) + dF11 * Fav22 * Fav33 + Fav11 * dF22 * Fav33 + Fav11 * Fav22 * dF33 + dF12 * Fav23 * Fav31 + Fav12 * dF23 * Fav31 + Fav12 * Fav23 * dF31 + dF13 * Fav21 * Fav32 + Fav13 * dF21 * Fav32 + Fav13 * Fav21 * dF32 - ( dF23 * Fav32 + Fav23 * dF32 ) - ( dF12 * Fav21 + Fav12 * dF21 ) - ( dF13 * Fav31 + Fav13 * dF31 ) - ( dF11 * Fav23 * Fav32 + Fav11 * dF23 * Fav32 + Fav11 * Fav23 * dF32 ) - ( dF12 * Fav21 * Fav33 + Fav12 * dF21 * Fav33 + Fav12 * Fav21 * dF33 ) - ( dF13 * Fav22 * Fav31 + Fav13 * dF22 * Fav31 + Fav13 * Fav22 * dF31 );

        /*auto InvFvMat = mat<3,3>( Fv22*Fv33-Fv23*Fv32 , Fv13*Fv32-Fv12*Fv33 , Fv12*Fv23-Fv13*Fv22,
              Fv23*Fv31-Fv21*Fv33 , Fv11*Fv33-Fv13*Fv31 , Fv13*Fv21-Fv11*Fv23,
              Fv21*Fv32-Fv22*Fv31 , Fv12*Fv31-Fv11*Fv32 , Fv11*Fv22-Fv12*Fv21 );*/
        // InvFv ( without 1/detF )
        const value_type InvFv11 = Fv22 * Fv33 - Fv23 * Fv32;
        const value_type InvFv12 = Fv13 * Fv32 - Fv12 * Fv33;
        const value_type InvFv13 = Fv12 * Fv23 - Fv13 * Fv22;
        const value_type InvFv21 = Fv23 * Fv31 - Fv21 * Fv33;
        const value_type InvFv22 = Fv11 * Fv33 - Fv13 * Fv31;
        const value_type InvFv23 = Fv13 * Fv21 - Fv11 * Fv23;
        const value_type InvFv31 = Fv21 * Fv32 - Fv22 * Fv31;
        const value_type InvFv32 = Fv12 * Fv31 - Fv11 * Fv32;
        const value_type InvFv33 = Fv11 * Fv22 - Fv12 * Fv21;

        const value_type FmtNL11 = dF22 + dF33 + dF22 * Fav33 + Fav22 * dF33 - dF23 * Fav32 - Fav23 * dF32;
        const value_type FmtNL12 = dF23 * Fav31 + Fav23 * dF31 - dF21 - dF21 * Fav33 - Fav21 * dF33;
        const value_type FmtNL13 = dF21 * Fav32 + Fav21 * dF32 - dF31 - dF31 * Fav22 - Fav31 * dF22;
        const value_type FmtNL21 = dF13 * Fav32 + Fav13 * dF32 - dF12 - dF12 * Fav33 - Fav12 * dF33;
        const value_type FmtNL22 = dF11 + dF33 + dF11 * Fav33 + Fav11 * dF33 - dF13 * Fav31 - Fav13 * dF31;
        const value_type FmtNL23 = dF12 * Fav31 + Fav12 * dF31 - dF32 - dF32 * Fav11 - Fav32 * dF11;
        const value_type FmtNL31 = dF12 * Fav23 + Fav12 * dF23 - dF13 - dF13 * Fav22 - Fav13 * dF22;
        const value_type FmtNL32 = dF13 * Fav21 + Fav13 * dF21 - dF23 - dF23 * Fav11 - Fav23 * dF11;
        const value_type FmtNL33 = dF11 + dF22 + dF11 * Fav22 + Fav11 * dF22 - dF12 * Fav21 - Fav12 * dF21;

        /*auto InvFtMat =  mat<3,3>(FmtNL11,FmtNL12,FmtNL13,
                                  FmtNL21,FmtNL22,FmtNL23,
                                  FmtNL31,FmtNL32,FmtNL33 );*/
        //auto dFmt = -(ddetF/pow(detFv,2))*trans(InvFvMat)+(cst(1.)/detFv)*/*trans*/(InvFtMat) ;

        const value_type detFv = M_locEvalPrecomputeDetF[q]( 0, 0 );
        const value_type scaleUsedWithInvFv = 1. / detFv;
        const value_type precomputePowDetF = M_locEvalPrecomputePowDetF[q]( 0, 0 ); //pow(detFv,2)  //idv(CoeffLame1)/pow(detFv,2)

        const value_type factorOther = ddetF / precomputePowDetF; //math::pow(detFv,2);
        const value_type dFmt11 = -factorOther * InvFv11 + scaleUsedWithInvFv * FmtNL11;
        const value_type dFmt12 = -factorOther * InvFv21 + scaleUsedWithInvFv * FmtNL12;
        const value_type dFmt13 = -factorOther * InvFv31 + scaleUsedWithInvFv * FmtNL13;
        const value_type dFmt21 = -factorOther * InvFv12 + scaleUsedWithInvFv * FmtNL21;
        const value_type dFmt22 = -factorOther * InvFv22 + scaleUsedWithInvFv * FmtNL22;
        const value_type dFmt23 = -factorOther * InvFv32 + scaleUsedWithInvFv * FmtNL23;
        const value_type dFmt31 = -factorOther * InvFv13 + scaleUsedWithInvFv * FmtNL31;
        const value_type dFmt32 = -factorOther * InvFv23 + scaleUsedWithInvFv * FmtNL32;
        const value_type dFmt33 = -factorOther * InvFv33 + scaleUsedWithInvFv * FmtNL33;

        //--------------------------------//
        //idv(CoeffLame2)*dF
        const value_type dFS_neohookean_a11 = coefflame2 * dF11;
        const value_type dFS_neohookean_a12 = coefflame2 * dF12;
        const value_type dFS_neohookean_a13 = coefflame2 * dF13;
        const value_type dFS_neohookean_a21 = coefflame2 * dF21;
        const value_type dFS_neohookean_a22 = coefflame2 * dF22;
        const value_type dFS_neohookean_a23 = coefflame2 * dF23;
        const value_type dFS_neohookean_a31 = coefflame2 * dF31;
        const value_type dFS_neohookean_a32 = coefflame2 * dF32;
        const value_type dFS_neohookean_a33 = coefflame2 * dF33;

        // _expr= val( idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv)) -idv(CoeffLame2) )*trace(dFmt*trans(grad(v)) ),
        const value_type precomputeLog = M_locEvalPrecomputeLogDetF[q]( 0, 0 ); // idv(CoeffLame1)*vf::log(/*vf::abs*/(detFv)) -idv(CoeffLame2)
        const value_type dFS_neohookean_b11 = precomputeLog * dFmt11;
        const value_type dFS_neohookean_b12 = precomputeLog * dFmt12;
        const value_type dFS_neohookean_b13 = precomputeLog * dFmt13;
        const value_type dFS_neohookean_b21 = precomputeLog * dFmt21;
        const value_type dFS_neohookean_b22 = precomputeLog * dFmt22;
        const value_type dFS_neohookean_b23 = precomputeLog * dFmt23;
        const value_type dFS_neohookean_b31 = precomputeLog * dFmt31;
        const value_type dFS_neohookean_b32 = precomputeLog * dFmt32;
        const value_type dFS_neohookean_b33 = precomputeLog * dFmt33;

        //_expr= ddetF*trace( val(idv(CoeffLame1)*trans(InvFvMat)/pow(detFv,2) )*trans(grad(v)) ),
        const value_type factorOther2 = ddetF * bulkModulus / precomputePowDetF; //        ddetF/precomputePowDetF;//math::pow(detFv,2);
        const value_type dFS_neohookean_c11 = factorOther2 * InvFv11;
        const value_type dFS_neohookean_c12 = factorOther2 * InvFv21;
        const value_type dFS_neohookean_c13 = factorOther2 * InvFv31;
        const value_type dFS_neohookean_c21 = factorOther2 * InvFv12;
        const value_type dFS_neohookean_c22 = factorOther2 * InvFv22;
        const value_type dFS_neohookean_c23 = factorOther2 * InvFv32;
        const value_type dFS_neohookean_c31 = factorOther2 * InvFv13;
        const value_type dFS_neohookean_c32 = factorOther2 * InvFv23;
        const value_type dFS_neohookean_c33 = factorOther2 * InvFv33;

        this->locMatrixShape()( 0, 0 ) = dFS_neohookean_a11 + dFS_neohookean_b11 + dFS_neohookean_c11;
        this->locMatrixShape()( 1, 0 ) = dFS_neohookean_a21 + dFS_neohookean_b21 + dFS_neohookean_c21;
        this->locMatrixShape()( 2, 0 ) = dFS_neohookean_a31 + dFS_neohookean_b31 + dFS_neohookean_c31;
        this->locMatrixShape()( 0, 1 ) = dFS_neohookean_a12 + dFS_neohookean_b12 + dFS_neohookean_c12;
        this->locMatrixShape()( 1, 1 ) = dFS_neohookean_a22 + dFS_neohookean_b22 + dFS_neohookean_c22;
        this->locMatrixShape()( 2, 1 ) = dFS_neohookean_a32 + dFS_neohookean_b32 + dFS_neohookean_c32;
        this->locMatrixShape()( 0, 2 ) = dFS_neohookean_a13 + dFS_neohookean_b13 + dFS_neohookean_c13;
        this->locMatrixShape()( 1, 2 ) = dFS_neohookean_a23 + dFS_neohookean_b23 + dFS_neohookean_c23;
        this->locMatrixShape()( 2, 2 ) = dFS_neohookean_a33 + dFS_neohookean_b33 + dFS_neohookean_c33;

        return this->locMatrixShape();
    }

  private:
    pc_mechprop_scalar_ptrtype M_pcMechPropField;
    ctx_mechprop_scalar_ptrtype M_ctxMechPropField;

    typename super_type::array_scalar_type M_locEvalBulkModulus, M_locEvalFieldCoefflame2;
    typename super_type::array_scalar_type M_locEvalPrecomputeLogDetF, M_locEvalPrecomputePowDetF, M_locEvalPrecomputeDetF;
};

/**
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     * MooneyRivlinCompressible :
     *   S = TODO
     * -------------------------------------------------------------------------------------------------  *
     * -------------------------------------------------------------------------------------------------  *
     */

/**
 * \class SolidMecStressTensorImpl
 * \brief det of a matrix
 *
 * @author Vincent Chabannes
 * @see
 */
template <typename ElementDisplacementType, typename ElementMechPropType, typename SpecificExprType, int QuadOrder>
class SolidMecStressTensorImpl
{
  public:
    typedef SolidMecStressTensorImpl<ElementDisplacementType, ElementMechPropType, SpecificExprType, QuadOrder> this_type;

    static const size_type context_displacement = vm::JACOBIAN | vm::KB | vm::GRAD;
    static const size_type context_mechprop_scalar = vm::JACOBIAN;
    static const size_type context = context_displacement;

    typedef ElementDisplacementType element_displacement_type;
    typedef ElementMechPropType element_mechprop_type;
    typedef SpecificExprType specific_expr_type;
    typedef Feel::FeelModels::MechanicalPropertiesDescription<element_mechprop_type> mechanical_properties_desc_type;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_displacement_type::functionspace_type functionspace_displacement_type;
    typedef typename functionspace_displacement_type::reference_element_type* fe_displacement_ptrtype;
    typedef typename functionspace_displacement_type::reference_element_type fe_displacement_type;
    //------------------------------------------------------------------------------//
    // mech prop functionspace
    typedef typename element_mechprop_type::functionspace_type functionspace_mechprop_type;
    typedef typename functionspace_mechprop_type::reference_element_type fe_mechprop_scalar_type;
    //------------------------------------------------------------------------------//
    // expression desc
    typedef typename functionspace_displacement_type::geoelement_type geoelement_type;
    //typedef typename functionspace_displacement_type::gm_type gm_type;
    typedef typename functionspace_displacement_type::value_type value_type;
    typedef value_type evaluate_type;
    static const uint16_type nDim = fe_displacement_type::nDim;
    static const uint16_type nRealDim = fe_displacement_type::nRealDim;
    static const uint16_type rank = fe_displacement_type::rank;
    static const uint16_type nComponents1 = fe_displacement_type::nComponents1;
    static const uint16_type nComponents2 = fe_displacement_type::nComponents2;
    static const bool is_terminal = true;

    static const uint16_type orderdisplacement = functionspace_displacement_type::basis_type::nOrder;

    //QuadOrder
    static const uint16_type imorderAuto = mpl::if_<boost::is_same<SpecificExprType, mpl::int_<ExprApplyType::EVAL>>,
                                                    mpl::int_<2 * orderdisplacement>,
                                                    mpl::int_<2 * orderdisplacement>>::type::value;
    static const uint16_type imorder = ( QuadOrder >= 0 ) ? QuadOrder : imorderAuto;
    static const bool imIsPoly = true;

    typedef Shape<nDim, Tensor2, false, false> shape_type;
    typedef Eigen::Matrix<value_type, shape_type::M, shape_type::N> matrix_shape_type;

    template <typename Func>
    struct HasTestFunction
    {
        static const bool result = false;
    };

    template <typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<boost::is_same<SpecificExprType, mpl::int_<ExprApplyType::JACOBIAN>>,
                                            typename mpl::if_<boost::is_same<Func, fe_displacement_type>,
                                                              mpl::bool_<true>,
                                                              mpl::bool_<false>>::type,
                                            mpl::bool_<false>>::type::value;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    SolidMecStressTensorImpl( element_displacement_type const& u,
                              Feel::FeelModels::MechanicalPropertiesDescription<element_mechprop_type> const& mechanicalPropertiesDesc )
        : M_displacement( boost::cref( u ) ),
          M_mechanicalPropertiesDesc( mechanicalPropertiesDesc )
    {
    }
    SolidMecStressTensorImpl( SolidMecStressTensorImpl const& op )
        : M_displacement( op.M_displacement ),
          M_mechanicalPropertiesDesc( op.M_mechanicalPropertiesDesc )
    {
    }
    ~SolidMecStressTensorImpl() {}

    element_displacement_type const& displacement() const { return M_displacement; }
    mechanical_properties_desc_type const& mechanicalPropertiesDesc() const { return M_mechanicalPropertiesDesc; }

    template <typename Geo_t, typename Basis_i_t, typename Basis_j_t>
    struct tensor
    {
        typedef typename element_displacement_type::value_type value_type;
        typedef shape_type shape;
        struct is_zero
        {
            static const bool value = false;
        };

        typedef tensorBase<Geo_t, Basis_i_t, Basis_j_t, shape, value_type> tensorbase_type;
        typedef boost::shared_ptr<tensorbase_type> tensorbase_ptrtype;

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
        {
            if ( expr.mechanicalPropertiesDesc().materialLaw() == "StVenantKirchhoff" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev, feu ) );
                else
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev, feu ) );
            }
            else if ( expr.mechanicalPropertiesDesc().materialLaw() == "NeoHookean" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev, feu ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev, feu ) );
                    else
                        CHECK( false ) << "invalid";
                }
                else
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev, feu ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev, feu ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory-simo1985" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985<Geo_t, Basis_i_t, Basis_j_t, this_type>( expr, geom, fev, feu ) );
                    else
                        CHECK( false ) << "invalid";
                }
            }
            else
                CHECK( false ) << "invalid materialLaw : " << expr.mechanicalPropertiesDesc().materialLaw() << "\n";
        }

        tensor( this_type const& expr,
                Geo_t const& geom, Basis_i_t const& fev )
        {
            if ( expr.mechanicalPropertiesDesc().materialLaw() == "StVenantKirchhoff" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev ) );
                else
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev ) );
            }
            else if ( expr.mechanicalPropertiesDesc().materialLaw() == "NeoHookean" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom, fev ) );
                    else
                        CHECK( false ) << "invalid";
                }
                else
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom, fev ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory-simo1985" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985<Geo_t, Basis_i_t, Basis_j_t, this_type>( expr, geom, fev ) );
                    else
                        CHECK( false ) << "invalid";
                }
            }
            else
                CHECK( false ) << "invalid materialLaw : " << expr.mechanicalPropertiesDesc().materialLaw() << "\n";
        }
        tensor( this_type const& expr, Geo_t const& geom )
        {
            if ( expr.mechanicalPropertiesDesc().materialLaw() == "StVenantKirchhoff" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom ) );
                else
                    M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffStVenantKirchhoff<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom ) );
            }
            else if ( expr.mechanicalPropertiesDesc().materialLaw() == "NeoHookean" )
            {
                if ( expr.mechanicalPropertiesDesc().useDisplacementPressureFormulation() )
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, true>( expr, geom ) );
                    else
                        CHECK( false ) << "invalid";
                }
                else
                {
                    if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "default" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressible<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheory<Geo_t, Basis_i_t, Basis_j_t, this_type, false>( expr, geom ) );
                    else if ( expr.mechanicalPropertiesDesc().compressibleNeoHookeanVariantName() == "molecular-theory-simo1985" )
                        M_tensorbase.reset( new tensorSolidMecFirstPiolaKirchhoffNeoHookeanCompressibleMolecularTheoryAndSimo1985<Geo_t, Basis_i_t, Basis_j_t, this_type>( expr, geom ) );
                    else
                        CHECK( false ) << "invalid";
                }
            }
            else
                CHECK( false ) << "invalid materialLaw : " << expr.mechanicalPropertiesDesc().materialLaw() << "\n";
        }
        tensor( tensor const& t )
            : M_tensorbase( t.M_tensorbase )
        {
        }

        template <typename IM>
        void init( IM const& im )
        {
            //M_tensor_expr.init( im );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ )
        {
            update( geom );
        }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ )
        {
            update( geom );
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
            return M_tensorbase->evalijq( i, j, q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalijq( i, j, c1, c2, q );
        }

        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evaliq( i, c1, c2, q );
        }

        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_tensorbase->evalq( c1, c2, q );
        }

      private:
        tensorbase_ptrtype M_tensorbase;
    };

  private:
    boost::reference_wrapper<const element_displacement_type> M_displacement;
    mechanical_properties_desc_type const& M_mechanicalPropertiesDesc;
};

/**
 * \brief FirstPiolaKirchhoffTensor
 */

template <int QuadOrder = -1, class ElementDisplacementType, class ElementMechPropType>
inline Expr<SolidMecStressTensorImpl<ElementDisplacementType, ElementMechPropType, mpl::int_<ExprApplyType::EVAL>, QuadOrder>>
solidMecFirstPiolaKirchhoffTensor( ElementDisplacementType const& u,
                                   Feel::FeelModels::MechanicalPropertiesDescription<ElementMechPropType> const& mechanicalPropertiesDesc )
{
    typedef SolidMecStressTensorImpl<ElementDisplacementType, ElementMechPropType, mpl::int_<ExprApplyType::EVAL>, QuadOrder> smstresstensor_t;
    return Expr<smstresstensor_t>( smstresstensor_t( u, mechanicalPropertiesDesc ) );
}

template <int QuadOrder = -1, class ElementDisplacementType, class ElementMechPropType>
inline Expr<SolidMecStressTensorImpl<ElementDisplacementType, ElementMechPropType, mpl::int_<ExprApplyType::JACOBIAN>, QuadOrder>>
solidMecFirstPiolaKirchhoffTensorJacobian( ElementDisplacementType const& u,
                                           Feel::FeelModels::MechanicalPropertiesDescription<ElementMechPropType> const& mechanicalPropertiesDesc )
{
    typedef SolidMecStressTensorImpl<ElementDisplacementType, ElementMechPropType, mpl::int_<ExprApplyType::JACOBIAN>, QuadOrder> smstresstensor_t;
    return Expr<smstresstensor_t>( smstresstensor_t( u, mechanicalPropertiesDesc ) );
}

} // namespace FeelModels
} // namespace vf
} // namespace Feel
#endif /* __SOLIDMEC_FIRSTPIOLAKIRCHHOFF_H */
