/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#ifndef __FEELPP_MODELS_VF_SOLIDMECGEOMAPEULERIAN_H
#define __FEELPP_MODELS_VF_SOLIDMECGEOMAPEULERIAN_H 1

#include <feel/feelmodels/modelvf/exprtensorbase.hpp>

namespace Feel
{
namespace FeelModels
{

template<typename ElementDispType,typename SpecificExprType>
class SolidMecGeomapEulerian
{
public:

    typedef SolidMecGeomapEulerian<ElementDispType,SpecificExprType> self_type;

    static const size_type context_disp = vm::JACOBIAN|vm::KB|vm::GRAD;
    static const size_type context = context_disp;

    typedef ElementDispType element_disp_type;
    typedef SpecificExprType spec_expr_type;
    //------------------------------------------------------------------------------//
    // displacement functionspace
    typedef typename element_disp_type::functionspace_type functionspace_disp_type;
    typedef typename functionspace_disp_type::reference_element_type* fe_disp_ptrtype;
    typedef typename functionspace_disp_type::reference_element_type fe_disp_type;

    typedef typename functionspace_disp_type::geoelement_type geoelement_type;
    //typedef typename functionspace_disp_type::gm_type gm_type;
    typedef typename functionspace_disp_type::value_type value_type;
    typedef value_type evaluate_type;
    static const bool is_terminal = true;

    static const uint16_type orderdisplacement = functionspace_disp_type::basis_type::nOrder;
    static const uint16_type nDim = functionspace_disp_type::nDim;

    template<typename Func>
    struct HasTestFunction
    {
      static const bool result = false;
    };
    template<typename Func>
    struct HasTrialFunction
    {
        static const bool result = mpl::if_<  mpl::bool_< ( SpecificExprType::value == ExprApplyType::JACOBIAN ) >,
                                              typename mpl::if_<boost::is_same<Func,fe_disp_type>,
                                                                     mpl::bool_<true>, mpl::bool_<false> >::type,
                                                   mpl::bool_<false> >::type::value;
    };
    using test_basis = std::nullptr_t;
    using trial_basis = std::nullptr_t;

    SolidMecGeomapEulerian( element_disp_type const& disp )
        :
        M_disp( disp )
    {}
    SolidMecGeomapEulerian( SolidMecGeomapEulerian const & op ) = default;

    ~SolidMecGeomapEulerian()
    {}

    //! polynomial order
    uint16_type polynomialOrder() const { return (orderdisplacement-1)*(nDim-1); }

    //! expression is polynomial?
    bool isPolynomial() const { return true; }

    element_disp_type const& disp() const { return M_disp; }

    typedef Shape<nDim, Tensor2, false, false> my_shape_type;


    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t >
    struct tensor : public tensorBase<Geo_t,Basis_i_t,Basis_j_t,my_shape_type,value_type >
    {
        typedef tensorBase<Geo_t,Basis_i_t,Basis_j_t,my_shape_type,value_type > super_type;
        //typedef ExprType expr_type;
        typedef self_type expr_type;
        //typedef typename expr_type::spec_expr_type SpecificExprType;
        //typedef typename expr_type::value_type value_type;
        //typedef typename expr_type::gm_type gm_type;
        typedef typename super_type::gm_type gm_type;
        //typedef typename expr_type::geoelement_type geoelement_type;
        typedef typename super_type::matrix_shape_type matrix_shape_type;
        using ret_type = typename super_type::ret_type;
        typedef typename super_type::gmc_type gmc_type;
        typedef typename super_type::shape_type shape;
        struct is_zero { static const bool value = false; };
        // fe disp context
        typedef typename fe_disp_type::PreCompute pc_disp_type;
        typedef std::shared_ptr<pc_disp_type> pc_disp_ptrtype;
        typedef typename fe_disp_type::template Context<context_disp, fe_disp_type,gm_type,geoelement_type,0, gmc_type::subEntityCoDim> ctx_disp_type;
        typedef std::shared_ptr<ctx_disp_type> ctx_disp_ptrtype;

        tensor( expr_type const& expr,Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            super_type( geom, fev, feu ),
            M_expr( expr ),
            M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(), this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
            M_locRes( expr.disp().gradExtents( *this->gmc()) ),
            M_locGradDisplacement( expr.disp().gradExtents( *this->gmc()) )
            {}

        tensor( expr_type const& expr,Geo_t const& geom, Basis_i_t const& fev )
            :
            super_type( geom, fev ),
            M_expr( expr ),
            M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
            M_locRes( expr.disp().gradExtents(*this->gmc()) ),
            M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) )
            {}

        tensor( expr_type const& expr, Geo_t const& geom )
            :
            super_type( geom ),
            M_expr( expr ),
            M_pcDisp( new pc_disp_type( expr.disp().functionSpace()->fe(), this->gmc()->xRefs() ) ),
            M_ctxDisp( new ctx_disp_type( expr.disp().functionSpace()->fe(),this->gmc(),(pc_disp_ptrtype const&)M_pcDisp ) ),
            M_locRes( expr.disp().gradExtents(*this->gmc()) ),
            M_locGradDisplacement( expr.disp().gradExtents(*this->gmc()) )
            {}

        template<typename IM>
        void init( IM const& im ) {}

        void update( Geo_t const& geom, Basis_i_t const& /*fev*/, Basis_j_t const& /*feu*/ ) { update(geom); }
        void update( Geo_t const& geom, Basis_i_t const& /*fev*/ ) { update(geom); }
        void update( Geo_t const& geom, uint16_type face ) { CHECK( false ) << "TODO"; }

        using super_type::evalijq; // fix clang warning

        value_type evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return M_locRes[q]( c1,c2 );
        }
        ret_type
        evalq( uint16_type q ) const
        {
            return ret_type(M_locRes[q].data());
        }

        void update( Geo_t const& geom )
        {
            this->setGmc( geom );
            std::fill( this->M_locRes.data(), this->M_locRes.data()+this->M_locRes.num_elements(), super_type::matrix_shape_type::Zero() );
            std::fill( this->M_locGradDisplacement.data(), this->M_locGradDisplacement.data()+this->M_locGradDisplacement.num_elements(), this->M_zeroLocTensor2/*super_type::loc_tensor2_type::Zero()*/ );

            if ( this->gmc()->faceId() != invalid_uint16_type_value ) /*face case*/
                this->M_pcDisp->update( this->gmc()->pc()->nodes() );
            this->M_ctxDisp->update( this->gmc(),  (pc_disp_ptrtype const&) this->M_pcDisp );
            this->M_expr.disp().grad( *this->M_ctxDisp, this->M_locGradDisplacement );
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
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            return evaliq( i,c1,c2,q, mpl::int_<SpecificExprType::value>() );
        }
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
            {
            DCHECK( SpecificExprType::value == ExprApplyType::EVAL ) << "only for EVAL expression";
            return ret_type(M_locRes[q].data());
        }

    private:


        void updateImpl( mpl::int_<ExprApplyType::EVAL> )
        {
            updateImpl( mpl::int_<ExprApplyType::EVAL>(), mpl::int_<super_type::gmc_type::nDim>() );
        }
        void updateImpl( mpl::int_<ExprApplyType::EVAL>, mpl::int_<2> /*Dim*/ )
        {
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                // compute : Fv*J*Cv^{-1} = Fv*J*(Fv^T*Fv)^{-1} = J*Fv^{-T} (with J=det(F))
                auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
                const value_type du1vdx = gradDisplacementEval(0,0), du1vdy = gradDisplacementEval(0,1);
                const value_type du2vdx = gradDisplacementEval(1,0), du2vdy = gradDisplacementEval(1,1);
                typename super_type::matrix_shape_type/*loc_tensor2_type*/ & theLocRes = this->M_locRes[q];
                theLocRes(0,0) = 1+du2vdy;
                theLocRes(0,1) = -du2vdx;
                theLocRes(1,0) = -du1vdy;
                theLocRes(1,1) = 1+du1vdx;
            }
        }
        void updateImpl( mpl::int_<ExprApplyType::EVAL>, mpl::int_<3> /*Dim*/ )
        {
            for ( uint16_type q = 0; q < this->gmc()->nPoints(); ++q )
            {
                // compute : Fv*J*Cv^{-1} = Fv*J*(Fv^T*Fv)^{-1} = J*Fv^{-T} (with J=det(F))
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
            }
        }
        void updateImpl( mpl::int_<ExprApplyType::JACOBIAN> ) {}
        //---------------------------------------------------------//
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::EVAL> ) const
        {
            return this->evalq( c1,c2,q );
        }
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ) const
        {
            CHECK( false ) << "not allow"; return 0;
        }
        //---------------------------------------------------------//
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::EVAL> ) const
        {
            return super_type::evalijq( i,j,q); // not allow
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ) const
        {
            return evalijq( i,j,q, mpl::int_<ExprApplyType::JACOBIAN>(),mpl::int_<super_type::gmc_type::nDim>() );
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ,mpl::int_<2> /*Dim*/ ) const
        {
            CHECK( false ) << "TODO mat";
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type dF11 = gradTrial( 0, 0, 0 ), dF12 = gradTrial( 0, 1, 0 );
            const value_type dF21 = gradTrial( 1, 0, 0 ), dF22 = gradTrial( 1, 1, 0 );
            matrix_shape_type & thelocRes = this->locMatrixShape();
            thelocRes(0,0) =  dF22;
            thelocRes(0,1) = -dF21;
            thelocRes(1,0) = -dF12;
            thelocRes(1,1) =  dF11;
            return ret_type(this->locMatrixShape().data());
        }
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ,mpl::int_<3> /*Dim*/ ) const
        {
            CHECK( false ) << "TODO mat";
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
            return ret_type(thelocRes.data());
        }

        value_type
        evalijq( uint16_type i,uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::EVAL> ) const
        {
            return this->evalq( c1,c2,q );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ) const
        {
            return evalijq( i,j,c1,c2,q,mpl::int_<ExprApplyType::JACOBIAN>(),mpl::int_<super_type::gmc_type::nDim>() );
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ,mpl::int_<2> /*Dim*/ ) const
        {
            if (c1==0)
            {
                if (c2==0)
                    return this->fecTrial()->grad( j, 1, 1, q );
                else
                    return -this->fecTrial()->grad( j, 1, 0, q );
            }
            else
            {
                if (c2==0)
                    return -this->fecTrial()->grad( j, 0, 1, q );
                else
                    return this->fecTrial()->grad( j, 0, 0, q );
            }
        }
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q, mpl::int_<ExprApplyType::JACOBIAN> ,mpl::int_<3> /*Dim*/ ) const
        {
            auto const& gradTrial = this->fecTrial()->grad( j, q );
            const value_type dF11 = gradTrial( 0, 0, 0 ), dF12 = gradTrial( 0, 1, 0 ), dF13 = gradTrial( 0, 2, 0 );
            const value_type dF21 = gradTrial( 1, 0, 0 ), dF22 = gradTrial( 1, 1, 0 ), dF23 = gradTrial( 1, 2, 0 );
            const value_type dF31 = gradTrial( 2, 0, 0 ), dF32 = gradTrial( 2, 1, 0 ), dF33 = gradTrial( 2, 2, 0 );
            auto const& gradDisplacementEval = this->M_locGradDisplacement[q];
            const value_type Fv11 = 1+gradDisplacementEval(0,0), Fv12 = gradDisplacementEval(0,1), Fv13 = gradDisplacementEval(0,2);
            const value_type Fv21 = gradDisplacementEval(1,0), Fv22 = 1+gradDisplacementEval(1,1), Fv23 = gradDisplacementEval(1,2);
            const value_type Fv31 = gradDisplacementEval(2,0), Fv32 = gradDisplacementEval(2,1), Fv33 = 1+gradDisplacementEval(2,2);

            if (c1==0)
            {
                if (c2==0)
                    return dF22+dF33+dF22*Fv33+Fv22*dF33-dF23*Fv32-Fv23*dF32;
                else if (c2==1)
                    return dF23*Fv31+Fv23*dF31-dF21-dF21*Fv33-Fv21*dF33;
                else
                    return dF21*Fv32+Fv21*dF32-dF31-dF31*Fv22-Fv31*dF22;
            }
            else if (c1==1)
            {
                if (c2==0)
                    return dF13*Fv32+Fv13*dF32-dF12-dF12*Fv33-Fv12*dF33;
                else if (c2==1)
                    return dF11+dF33+dF11*Fv33+Fv11*dF33-dF13*Fv31-Fv13*dF31;
                else
                    return dF12*Fv31+Fv12*dF31-dF32-dF32*Fv11-Fv32*dF11;
            }
            else
            {
                if (c2==0)
                    return dF12*Fv23+Fv12*dF23-dF13-dF13*Fv22-Fv13*dF22;
                else if (c2==1)
                    return dF13*Fv21+Fv13*dF21-dF23-dF23*Fv11-Fv23*dF11;
                else
                    return dF11+dF22+dF11*Fv22+Fv11*dF22-dF12*Fv21-Fv12*dF21;
            }
        }

    private:
        expr_type const& M_expr;
        pc_disp_ptrtype M_pcDisp;
        ctx_disp_ptrtype M_ctxDisp;

        typename super_type::array_shape_type M_locRes;
        typename super_type::array_tensor2_type M_locGradDisplacement;

    };
private:
    element_disp_type const& M_disp;

}; // class SolidMecGeomapEulerian


template<class ElementDispType>
inline
auto
solidMecGeomapEulerian( ElementDispType const& v )
{
    typedef SolidMecGeomapEulerian<unwrap_ptr_t<ElementDispType>,mpl::int_<ExprApplyType::EVAL> > myexpr_type;
    return Expr< myexpr_type >(  myexpr_type( unwrap_ptr(v) ) );
}
template<class ElementDispType>
inline
auto
solidMecGeomapEulerianJacobian( ElementDispType const& v )
{
    typedef SolidMecGeomapEulerian<unwrap_ptr_t<ElementDispType>,mpl::int_<ExprApplyType::JACOBIAN> > myexpr_type;
    return Expr< myexpr_type >(  myexpr_type( unwrap_ptr(v) ) );
}


} // namespace FeelModels
} // namespace Feel
#endif /* __SOLIDMECGEOMAPEULERIAN_H */
