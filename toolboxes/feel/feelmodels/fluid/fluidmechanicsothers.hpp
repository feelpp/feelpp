#ifndef FEELPP_TOOLBOXES_FLUIDMECHANICS_OTHERS_HPP
#define FEELPP_TOOLBOXES_FLUIDMECHANICS_OTHERS_HPP 1

#include <feel/feelpde/operatorpcd.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateInHousePreconditioner( DataUpdateBase & data, ModelContextType const& mctx ) const
{
    if ( M_preconditionerAttachPMM )
        this->updateInHousePreconditionerPMM( data, mctx );
    if ( M_preconditionerAttachPCD )
        this->updateInHousePreconditionerPCD( data, mctx );
}

template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateInHousePreconditionerPMM( DataUpdateBase & data, ModelContextType const& mctx ) const
{
    if ( !this->algebraicFactory() )
        return; // TODO : should be use with multiphysics toolboxes as heat-fluid

    bool hasAlreadyBuiltPMM = this->algebraicFactory()->hasAuxiliarySparseMatrix( "pmm" );
    if ( hasAlreadyBuiltPMM && !M_pmmNeedUpdate )
        return;
    sparse_matrix_ptrtype pmmMat;
    if ( hasAlreadyBuiltPMM )
        pmmMat = this->algebraicFactory()->auxiliarySparseMatrix( "pmm" );
    else
    {
        pmmMat = this->algebraicBackend()->newMatrix(_trial=this->functionSpacePressure(), _test=this->functionSpacePressure());
        this->algebraicFactory()->attachAuxiliarySparseMatrix( "pmm", pmmMat );
    }
    CHECK( pmmMat ) << "pmmMat is not initialized";

    pmmMat->zero();
    auto massbf = form2( _trial=this->functionSpacePressure(), _test=this->functionSpacePressure(),_matrix=pmmMat);
    //auto u = this->functionSpaceVelocity()->element( vecSol, this->rowStartInVector() );
    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& p = this->fieldPressure();
    auto const& se = mctx.symbolsExpr();

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto coeff = cst(1.)/fluidMecViscosity(gradv(u),*physicFluidData,matProps,se);
            massbf += integrate( _range=range, _expr=coeff*inner( idt(p),id(p) ) );
        }
    }
    pmmMat->close();
    M_pmmNeedUpdate = false;
}


template< typename ConvexType, typename BasisVelocityType, typename BasisPressureType>
template <typename ModelContextType>
void
FluidMechanics<ConvexType,BasisVelocityType,BasisPressureType>::updateInHousePreconditionerPCD( DataUpdateBase & data, ModelContextType const& mctx ) const
{
    this->log("FluidMechanics","updateInHousePreconditionerPCD", "start" );

    CHECK( this->hasOperatorPCD() ) << "operator PCD does not init";

    typedef Feel::Alternatives::OperatorPCD<space_velocity_type,space_pressure_type> op_pcd_type;
    std::shared_ptr<op_pcd_type> myOpPCD =
        std::dynamic_pointer_cast<op_pcd_type>( this->operatorPCD() );

    auto const& u = mctx.field( FieldTag::velocity(this), "velocity" );
    auto const& se = mctx.symbolsExpr();

    myOpPCD->updateStart();

    //CHECK( this->physicsFromCurrentType().size() == 1 ) << "TODO";
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicFluidData = std::static_pointer_cast<ModelPhysicFluid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& therange = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& matProps = this->materialsProperties()->materialProperties( matName );
            auto const& rhoExpr = expr( this->materialsProperties()->density( matName ).template expr<1,1>(), se );
            auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(u),*physicFluidData,matProps,se);

            if ( physicFluidData->equation() == "Stokes" || physicFluidData->equation() == "StokesTransient" )
            {
                if (this->hasMeshMotion() )
                {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, -rhoExpr*idv(this->meshMotionTool()->velocity()), true );
#endif
                }
                else
                {
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, vf::zero<nDim,1>(), false );
                }
            }
            else if ( physicFluidData->equation() == "Navier-Stokes" && this->useSemiImplicitTimeScheme() )
            {
                auto betaU = *M_fieldVelocityExtrapolated;//this->timeStepBDF()->poly();
                if (this->hasMeshMotion() )
                {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*( idv(betaU)-idv(this->meshMotionTool()->velocity()) ), true );
#endif
                }
                else
                {
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*idv(betaU), true );
                }
            }
            else if ( physicFluidData->equation() == "Navier-Stokes" )
            {
                if (this->hasMeshMotion() )
                {
#if defined( FEELPP_MODELS_HAS_MESHALE )
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*( idv(u)-idv(this->meshMotionTool()->velocity()) ), true );
#endif
                }
                else
                {
                    myOpPCD->updateFpDiffusionConvection( therange, muExpr, rhoExpr*idv(u), true );
                }
            }

            if ( !this->isStationaryModel() )
            {
                myOpPCD->updateFpMass( therange, rhoExpr*this->timeStepBDF()->polyDerivCoefficient(0) );
            }
            if ( data.hasInfo( "use-pseudo-transient-continuation" ) )
            {
#if 0
                //Feel::cout << "updsate PCD : use-pseudo-transient-continuation\n";
                //Warning : it's a copy past, should be improve : TODO!
                double pseudoTimeStepDelta = data.doubleInfo("pseudo-transient-continuation.delta");
                auto norm2_uu = this->materialProperties()->fieldRho().functionSpace()->element(); // TODO : improve this (maybe create an expression instead)
                //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(idv(u))/h());
                auto fieldNormu = u.functionSpace()->compSpace()->element( norm2(idv(u)) );
                auto maxu = fieldNormu.max( this->materialProperties()->fieldRho().functionSpace() );
                //auto maxux = u[ComponentType::X].max( this->materialProperties()->fieldRho().functionSpace() );
                //auto maxuy = u[ComponentType::Y].max( this->materialProperties()->fieldRho().functionSpace() );
                //norm2_uu.on(_range=M_rangeMeshElements,_expr=norm2(vec(idv(maxux),idv(maxux)))/h());
                norm2_uu.on(_range=M_rangeMeshElements,_expr=idv(maxu)/h());

                myOpPCD->updateFpMass( therange, (1./pseudoTimeStepDelta)*idv(norm2_uu) );
#else
                CHECK(false) << "TODO : require P0d space, maybe try to find another alternative";
#endif
            }
        }
    } // foreach physic

    if ( !dynamic_cast<DataUpdateJacobian*>(&data) || !boption(_name="pcd.apply-homogeneous-dirichlet-in-newton",_prefix=this->prefix()) )
    {
#if 0 // VINCENT
        // auto const& fieldRho = this->materialProperties()->fieldRho();
        // auto rhoExpr = idv( fieldRho );
        //auto se = this->symbolsExpr();
        auto rhoExpr = this->materialsProperties()->template materialPropertyExpr<1,1>( "density", se );
        for( auto const& d : M_bcDirichlet )
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, name(d), expression(d,se) );
        for( auto const& d : M_bcMovingBoundaryImposed )
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, name(d), idv(M_meshALE->velocity()) );
        for ( auto const& inletbc : M_fluidInletDesc )
        {
            std::string const& marker = std::get<0>( inletbc );
            auto const& inletVel = std::get<0>( M_fluidInletVelocityInterpolated.find(marker)->second );
            myOpPCD->updateFpBoundaryConditionWithDirichlet( rhoExpr, marker, -idv(inletVel)*N() );
        }
#endif
    }
    else
        Feel::cout << "PCD NOT UP BC \n";

    // updated from the outside
    for ( auto const& f : M_addUpdateInHousePreconditionerPCD )
        f.second.second( *myOpPCD, data );

    myOpPCD->updateFinish();

    this->log("FluidMechanics","updateInHousePreconditionerPCD", "finish" );
}



} // namespace FeelModels
} // namespace Feel

#endif
