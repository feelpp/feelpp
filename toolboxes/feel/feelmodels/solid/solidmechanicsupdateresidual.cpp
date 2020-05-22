/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/solid/solidmechanics.hpp>

//#include <feel/feelmodels/modelvf/solidmecstvenantkirchhoff.hpp>
#include <feel/feelmodels/modelvf/solidmecfirstpiolakirchhoff.hpp>
#include <feel/feelmodels/modelvf/solidmecincompressibility.hpp>
#include <feel/feelmodels/modelvf/solidmecgeomapeulerian.hpp>

namespace Feel
{
namespace FeelModels
{


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& X = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _buildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();


    std::string sc=(_buildCstPart)?" (cst part)":" (non cst part)";
    this->log("SolidMechanics","updateResidual", "start"+sc );
    this->timerTool("Solve").start();

    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;


    double timeSteppingScaling = 1.;
    bool timeSteppingEvaluateResidualWithoutTimeDerivative = false;
    if ( !this->isStationary() )
    {
        timeSteppingEvaluateResidualWithoutTimeDerivative = data.hasInfo( prefixvm(this->prefix(),"time-stepping.evaluate-residual-without-time-derivative") );
        if ( M_timeStepping == "Theta" )
        {
            if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                timeSteppingScaling = 1. - M_timeStepThetaValue;
            else
                timeSteppingScaling = M_timeStepThetaValue;
        }
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //--------------------------------------------------------------------------------------------------//

    mesh_ptrtype mesh = M_XhDisplacement->mesh();

    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormDisplacement = form1( _test=M_XhDisplacement, _vector=R,_rowstart=rowStartInVector );

    //--------------------------------------------------------------------------------------------------//
    auto const u = M_XhDisplacement->element(X, rowStartInVector);
    auto const& v = this->fieldDisplacement();

    auto se = this->symbolsExpr();
    //--------------------------------------------------------------------------------------------------//

    size_type blockIndexPressure = invalid_v<size_type>;
    std::decay_t<decltype(M_XhPressure->elementPtr(*X, blockIndexPressure))> p;
    if ( this->hasDisplacementPressureFormulation() )
    {
        // define pressure field
        blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
        p = M_XhPressure->elementPtr(*X, rowStartInVector+blockIndexPressure);
    }

    //#if defined(SOLIDMECHANICS_VOLUME_FORCE)
    //auto f = SOLIDMECHANICS_VOLUME_FORCE(this->shared_from_this());
    //#endif
    //auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());
#if 0
    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    //std::cout << "alpha_f " << alpha_f << "alpha_m " <<alpha_m <<"\n";
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);
#endif

    //--------------------------------------------------------------------------------------------------//

    // auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    // auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    // auto const& rho = this->mechanicalProperties()->fieldRho();
    //Identity Matrix
    auto Id = eye<nDim,nDim>();
    // auto Fv = Id + gradv(u);
    // auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    // auto Sv = idv(coeffLame1)*trace(Ev)*Id + 2*idv(coeffLame2)*Ev;
    //case elastic
    auto Ev_elastic = sym(gradv(u));

    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = this->materialsProperties()->materialProperties( matName );
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );

            //--------------------------------------------------------------------------------------------------//
            // stress tensor terms
            //--------------------------------------------------------------------------------------------------//
            this->timerTool("Solve").start();

            if ( physicSolidData->equation() == "Hyper-Elasticity" )
            {
                if (!BuildCstPart)
                {
                    auto FSv = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,*physicSolidData,matProperties,se);
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   //_expr= trace(val(Fv*Sv)*trans(grad(v))),
                                   //_expr=trace(Feel::FeelModels::stressStVenantKirchhoff(u,coeffLame1,coeffLame2)*trans(grad(v))),
                                   //_expr=Feel::FeelModels::stressStVenantKirchhoffResidual(u,coeffLame1,coeffLame2),// le dernier 
                                   _expr= timeSteppingScaling*inner(FSv,grad(v)),
                                   _geomap=this->geomap() );
                }
            }
            else if ( physicSolidData->equation() == "Elasticity" )
            {
                if (!BuildCstPart && !UseJacobianLinearTerms)
                {
                    auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                    auto lameSecondExpr = expr( matProperties.property( "Lame-second-parameter" ).exprScalar(), se );
                    auto Sv_elastic = lameFirstExpr*trace(Ev_elastic)*Id + 2*lameSecondExpr*Ev_elastic;
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner( Sv_elastic,grad(v) ),
                                   _geomap=this->geomap() );
                }
            }
            double timeElapsedBis = this->timerTool("Solve").stop();
            this->log("SolidMechanics","updateResidual",
                      "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

            //--------------------------------------------------------------------------------------------------//
            // constraint on compressibility
            //--------------------------------------------------------------------------------------------------//
            if ( physicSolidData->useDisplacementPressureFormulation() && !BuildCstPart )
            {
                if ( physicSolidData->equation() == "Hyper-Elasticity" )
                {
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   //_expr= trace(val(-idv(p)*trans(InvFv))*trans(grad(v))),
                                   //_expr=inner( -idv(p)*Feel::FeelModels::solidMecIncompressibilityPressure(u,p,*this->mechanicalProperties()),grad(v) ),
                                   _expr=inner( /*-idv(p)**/Feel::FeelModels::solidMecPressureFormulationMultiplier(u,*p,*physicSolidData),grad(v) ),
                                   _geomap=this->geomap() );
                }
                else if ( physicSolidData->equation() == "Elasticity" )
                {
                    linearFormDisplacement +=
                        integrate( _range= range,
                                   _expr= inner( -idv(p)*Id ,grad(v)),
                                   _geomap=this->geomap() );
                }

                auto linearFormPressure = form1( _test=M_XhPressure, _vector=R,_rowstart=rowStartInVector+blockIndexPressure );
                linearFormPressure +=
                    integrate( _range=range,
                               //_expr=val(Fav22+Fav11+Fav11*Fav22-Fav21*Fav12)*id(q),
                               //_expr= detFvM1*id(q),
                               _expr= Feel::FeelModels::solidMecPressureFormulationConstraint(u,*physicSolidData)*id(p),
                               _geomap=this->geomap() );


                if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
                {
                    auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                    linearFormPressure +=
                        integrate(_range=range,
                                  _expr= -(cst(1.)/lameFirstExpr)*idv(p)*id(p),
                                  _geomap=this->geomap() );
                }
                else
                {
                    auto K = expr( matProperties.property( "bulk-modulus" ).exprScalar(), se );
                    linearFormPressure +=
                        integrate( _range=range,
                                   _expr= -(cst(1.)/K)*idv(p)*id(p),
                                   _geomap=this->geomap() );
                }

            }

            //--------------------------------------------------------------------------------------------------//
            // discretisation acceleration term
            //--------------------------------------------------------------------------------------------------//
            if (!this->isStationary())
            {
                auto const& densityProp = this->materialsProperties()->density( matName );
                auto densityExpr = expr( densityProp.expr(), se );

                if ( M_timeStepping == "Newmark" )
                {
                    if (!BuildCstPart && !UseJacobianLinearTerms)
                    {
                        if ( !this->useMassMatrixLumped() )
                        {
                            linearFormDisplacement +=
                                integrate( _range=range,
                                           _expr= M_timeStepNewmark->polySecondDerivCoefficient()*densityExpr*inner(idv(u),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            if ( this->massMatrixLumped()->size1() == R->size() )
                            {
                                auto myvec = this->backend()->newVector(M_XhDisplacement);
                                *myvec = u;
                                myvec->scale(M_timeStepNewmark->polySecondDerivCoefficient());
                                R->close();
                                R->addVector( myvec, this->massMatrixLumped() );
                            }
                            else
                            {
                                CHECK( false ) << "TODO";
                            }
                        }
                    }
                    if (BuildCstPart)
                    {
                        auto polySecondDerivDisp = M_timeStepNewmark->polySecondDeriv();
                        if ( !this->useMassMatrixLumped() )
                        {
                            linearFormDisplacement +=
                                integrate( _range=range,
                                           _expr= -densityExpr*inner(idv(polySecondDerivDisp),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            if ( this->massMatrixLumped()->size1() == R->size() )
                            {
                                auto myvec = this->backend()->newVector(M_XhDisplacement);
                                *myvec = polySecondDerivDisp;
                                myvec->scale(-1.0);
                                R->close();
                                R->addVector( myvec, this->massMatrixLumped() );
                            }
                            else
                            {
                                R->close();
                                auto uAddResidual = M_XhDisplacement->element( R, rowStartInVector );
                                auto uDiagMassMatrixLumped = M_XhDisplacement->element( M_vecDiagMassMatrixLumped );
                                uAddResidual.add(-1.0, element_product( uDiagMassMatrixLumped, polySecondDerivDisp ) );
                            }
                        }
                    }
                } // Newmark
                else
                {
                    CHECK( this->hasStartSubBlockSpaceIndex( "velocity" ) ) << "no SubBlockSpaceIndex velocity";
                    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");
                    //std::cout << "RESIDUAL bdf\n";
                    if ( timeSteppingEvaluateResidualWithoutTimeDerivative )
                    {
                        if ( !BuildCstPart )
                        {
                            auto const curVel = M_XhDisplacement->element(X, rowStartInVector+startBlockIndexVelocity);
                            form1( _test=M_XhDisplacement, _vector=R,
                                   _rowstart=rowStartInVector+startBlockIndexVelocity ) +=
                                integrate( _range=range,
                                           _expr= -timeSteppingScaling*densityExpr*inner(idv(curVel),id(v)),
                                           _geomap=this->geomap() );
                        }
                    }
                    else if (!BuildCstPart && !UseJacobianLinearTerms)
                    {
                        auto const curVel = M_XhDisplacement->element(X, rowStartInVector+startBlockIndexVelocity);
                        if ( !this->useMassMatrixLumped() )
                        {
                            linearFormDisplacement +=
                                integrate( _range=range,
                                           _expr= M_timeStepBdfVelocity->polyDerivCoefficient(0)*densityExpr*inner(idv(curVel),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            CHECK( false ) << "TODO";
                        }

                        form1( _test=M_XhDisplacement, _vector=R,
                               _rowstart=rowStartInVector+startBlockIndexVelocity ) +=
                            integrate( _range=range,
                                       _expr= densityExpr*inner(M_timeStepBdfDisplacement->polyDerivCoefficient(0)*idv(u)-timeSteppingScaling*idv(curVel),id(v)),
                                       _geomap=this->geomap() );
                    }

                    if ( BuildCstPart && !timeSteppingEvaluateResidualWithoutTimeDerivative )
                    {
                        auto rhsTimeStepVelocity = M_timeStepBdfVelocity->polyDeriv();
                        if ( !this->useMassMatrixLumped() )
                        {
                            linearFormDisplacement +=
                                integrate( _range=range,
                                           _expr= -densityExpr*inner(idv(rhsTimeStepVelocity),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            R->close();
                            auto uAddRhs = M_XhDisplacement->element( R, rowStartInVector );
                            auto uDiagMassMatrixLumped = M_XhDisplacement->element( M_vecDiagMassMatrixLumped );
                            uAddRhs.add( -1., element_product( uDiagMassMatrixLumped, rhsTimeStepVelocity ) );
                        }
                        auto rhsTimeStepDisplacement = M_timeStepBdfDisplacement->polyDeriv();
                        form1( _test=M_XhDisplacement, _vector=R,
                               _rowstart=rowStartInVector+startBlockIndexVelocity ) +=
                            integrate( _range=range,
                                       _expr= -densityExpr*inner(idv(rhsTimeStepDisplacement),id(v)),
                                       _geomap=this->geomap() );
                    }
                } // BDF
            }

        } // matName
    } // physics
#if 0
    //--------------------------------------------------------------------------------------------------//
    // incompressibility terms
    if ( this->useDisplacementPressureFormulation() && !BuildCstPart)
    {
        // define pressure field
        size_type blockIndexPressure = rowStartInVector+this->startSubBlockSpaceIndex("pressure");
        auto const p = M_XhPressure->element(X, blockIndexPressure);
        // assemble
        this->updateResidualIncompressibilityTerms(u,p,R);
    }
#endif
    //--------------------------------------------------------------------------------------------------//
#if 0
    // viscoelastic terms
    // VERY OLD! : must be fix and updated
    this->updateResidualViscoElasticityTerms(u,R);
#endif
    //--------------------------------------------------------------------------------------------------//
    // source term
    if (BuildCstPart)
    {
        this->updateSourceTermResidual( R, timeSteppingScaling );
    }
    //--------------------------------------------------------------------------------------------------//
    // neumann bc
    if (BuildCstPart)
    {
        this->updateBCNeumannResidual( R, timeSteppingScaling );
    }
    // follower pressure bc
    if (!BuildCstPart)
    {
        this->updateBCFollowerPressureResidual( u,R, timeSteppingScaling );
    }
    // robin bc
    if ( !BuildCstPart )// && !UseJacobianLinearTerms)
    {
        this->updateBCRobinResidual( u, R, timeSteppingScaling );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateResidual",
              "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
#if 0
SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualIncompressibilityTerms( element_displacement_external_storage_type const& u,
                                                                              element_pressure_external_storage_type const& p, vector_ptrtype& R) const
{
    using namespace Feel::vf;

    boost::mpi::timer thetimer;
    this->log("SolidMechanics","updateResidualIncompressibilityTerms", "start" );

    auto mesh = M_XhDisplacement->mesh();
    auto const& v = this->fieldDisplacement();
    auto const& q = this->fieldPressure();
    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();

    /*double alpha_f=M_genAlpha_alpha_f;
     double alpha_m=M_genAlpha_alpha_m;
     double gamma=0.5+alpha_m-alpha_f;
     double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);*/

    size_type rowStartInVector = this->rowStartInVector();
    size_type blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
    auto linearFormDisplacement = form1( _test=M_XhDisplacement, _vector=R,
                                         _rowstart=rowStartInVector );
    auto linearFormPressure = form1( _test=M_XhPressure, _vector=R,
                                     _rowstart=rowStartInVector+blockIndexPressure );
#if 0
    if ( M_modelName == "Hyper-Elasticity" )
    {
        linearFormDisplacement +=
            integrate( _range=M_rangeMeshElements,
                       //_expr= trace(val(-idv(p)*trans(InvFv))*trans(grad(v))),
                       //_expr=inner( -idv(p)*Feel::FeelModels::solidMecIncompressibilityPressure(u,p,*this->mechanicalProperties()),grad(v) ),
                       _expr=inner( /*-idv(p)**/Feel::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()),grad(v) ),
                       _geomap=this->geomap() );
    }
    else if ( M_modelName == "Elasticity-Large-Deformation" || M_modelName == "Elasticity" )
    {
        //Identity Matrix
        auto Id = eye<nDim,nDim>();
        linearFormDisplacement +=
            integrate( _range= M_rangeMeshElements,
                       _expr= inner( -idv(p)*Id ,grad(v)),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    linearFormPressure +=
        integrate( _range=M_rangeMeshElements,
                   //_expr=val(Fav22+Fav11+Fav11*Fav22-Fav21*Fav12)*id(q),
                   //_expr= detFvM1*id(q),
                   _expr= Feel::FeelModels::solidMecPressureFormulationConstraint(u,*this->mechanicalProperties())*id(q),
                   _geomap=this->geomap() );


    if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
    {
        linearFormPressure +=
            integrate(_range=M_rangeMeshElements,
                      _expr= -(cst(1.)/idv(coeffLame1))*idv(p)*id(q),
                      _geomap=this->geomap() );
    }
    else
    {
        auto kappa = idv(this->mechanicalProperties()->fieldBulkModulus());
        linearFormPressure +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -(cst(1.)/kappa)*idv(p)*id(q),
                       _geomap=this->geomap() );
    }

#endif
    double timeElapsed=thetimer.elapsed();
    this->log("SolidMechanics","updateResidualIncompressibilityTerms",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );

}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualViscoElasticityTerms( element_displacement_external_storage_type const& u, vector_ptrtype& R) const
{
#if 0
    using namespace Feel::vf;

    if (this->verbose()) std::cout << "[SolidMechanics] : updateResidualViscoElasticityTerms start\n";

    auto mesh = M_XhDisplacement->mesh();

    auto u = U.element<0>();
    auto v = U.element<0>();

    auto Ev2 = 0.5*(gradv(u)+trans(gradv(u)) );// + 0.5*trans(gradv(u))*gradv(u);


    double gammav=0.01;

    auto Buzz1bis = M_bdf_displ_struct->polyDeriv();
    auto buzz1bis = Buzz1bis.element<0>();


    auto Evbis = 0.5*(gradv(buzz1bis)+trans(gradv(buzz1bis)) );

    form1( _test=M_XhDisplacement, _vector=R ) +=
        integrate( _range=M_rangeMeshElements,
                   _expr= gammav*M_bdf_displ_struct->polyDerivCoefficient(0)*trace( Ev2*trans(grad(v))),
                   _geomap=this->geomap() );

    form1( _test=M_XhDisplacement, _vector=R ) +=
        integrate( _range=M_rangeMeshElements,
                   _expr= -gammav*trace( Evbis*trans(grad(v))),
                   _geomap=this->geomap() );

    if (this->verbose()) std::cout << "[SolidMechanics] : updateResidualViscoElasticityTerms finish\n";
#endif
}
#endif
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( DataNewtonInitialGuess & data ) const
{
    if ( !this->hasDirichletBC() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "start",
                                               this->worldComm(),this->verboseAllProc());

    vector_ptrtype& U = data.initialGuess();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto u = Xh->element( U, this->rowStartInVector() );

    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=expression(d,this->symbolsExpr()) );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges),
                 _expr=expression(d,this->symbolsExpr()) );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(mesh,listMarkerPoints),
                 _expr=expression(d,this->symbolsExpr()) );
    }
    for ( auto const& bcDirComp : M_bcDirichletComponents )
    {
        ComponentType comp = bcDirComp.first;
        for( auto const& d : bcDirComp.second )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d), comp ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                u[comp].on(_range=markedfaces(mesh,listMarkerFaces ),
                           _expr=expression(d,this->symbolsExpr()) );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges),
                           _expr=expression(d,this->symbolsExpr()) );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(mesh,listMarkerPoints),
                           _expr=expression(d,this->symbolsExpr()) );
        }
    }

    // update info for synchronization
    this->updateDofEliminationIds( "displacement", data );

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual( vector_ptrtype& R, double timeSteppingScaling ) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                       _expr= -timeSteppingScaling*expression(d,this->symbolsExpr())*inner( N(),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                       _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                       _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr())*N(),id(v) ),
                       _geomap=this->geomap() );
}



SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinResidual(element_displacement_external_storage_type const& u, vector_ptrtype& R, double timeSteppingScaling ) const
{

    if ( this->M_bcRobin.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                       _expr= timeSteppingScaling*inner( expression1(d,this->symbolsExpr())(0,0)*idv(u) - expression2(d,this->symbolsExpr()) ,id(u) ),
                       _geomap=this->geomap() );

}




SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermResidual( vector_ptrtype& R, double timeSteppingScaling ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
    {
        auto rangeBodyForceUsed = markers(d).empty()? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
        myLinearForm +=
            integrate( _range=rangeBodyForceUsed,
                       _expr= -timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                       _geomap=this->geomap() );
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureResidual( element_displacement_external_storage_type const& u, vector_ptrtype& R, double timeSteppingScaling ) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -timeSteppingScaling*expression(d,this->symbolsExpr())*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d,this->symbolsExpr()),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d,this->symbolsExpr())*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidualDofElimination( DataUpdateResidual & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() )
        return;

    this->log("SolidMechanics","updateResidualDofElimination","start" );

    this->updateDofEliminationIds( "displacement", data );

    this->log("SolidMechanics","updateResidualDofElimination","finish" );
}


} // FeelModels
} // Feel



