/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobian( DataUpdateJacobian & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& X = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool BuildCstPart = data.buildCstPart();

    std::string sc=(BuildCstPart)?" (cst part)":" (non cst part)";
    this->log("SolidMechanics","updateJacobian", "start"+sc);
    this->timerTool("Solve").start();

    //--------------------------------------------------------------------------------------------------//

    mesh_ptrtype mesh = M_XhDisplacement->mesh();

    size_type rowStartInVector = this->rowStartInVector();
    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    auto bilinearForm_PatternDefault = form2( _test=M_XhDisplacement,_trial=M_XhDisplacement,_matrix=J,
                                              _pattern=size_type(Pattern::DEFAULT),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearForm_PatternCoupled = form2( _test=M_XhDisplacement,_trial=M_XhDisplacement,_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    auto const u = M_XhDisplacement->element(X, rowStartInVector);
    auto const& v = this->fieldDisplacement();

    //--------------------------------------------------------------------------------------------------//
#if 0
    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);
#endif

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }
    //--------------------------------------------------------------------------------------------------//

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    auto const& rho = this->mechanicalProperties()->fieldRho();
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    auto Fv = Id + gradv(u);
    auto dF = gradt(u);
    auto Ev = sym(gradv(u)) + 0.5*trans(gradv(u))*gradv(u);
    auto dE = sym(gradt(u)) + 0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
    auto Sv = idv(coeffLame1)*trace(Ev)*Id + 2*idv(coeffLame2)*Ev;
    auto dS = idv(coeffLame1)*trace(dE)*Id + 2*idv(coeffLame2)*dE;
    //case elastic
    auto dE_elastic = sym(gradt(u));
    auto dS_elastic = idv(coeffLame1)*trace(dE_elastic)*Id + 2*idv(coeffLame2)*dE_elastic;

    //--------------------------------------------------------------------------------------------------//
    // stress tensor terms
    //thetimerBis.restart();
    this->timerTool("Solve").start();

    if ( M_modelName == "Hyper-Elasticity" )
    {
        if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
        {
            if (!BuildCstPart)
            {
                auto const dFS_neohookean = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobian<3*(nOrderDisplacement-1)>(u,*this->mechanicalProperties());
                bilinearForm_PatternCoupled +=
                    integrate (_range=M_rangeMeshElements,
                               //_expr= trace( (dF*val(Sv) + val(Fv)*dS)*trans(grad(v)) ),
                               //_expr=Feel::vf::FSI::stressStVenantKirchhoffJacobian(u,coeffLame1,coeffLame2), //le dernier
                               _expr= timeSteppingScaling*inner( dFS_neohookean, grad(v) ),
                               _geomap=this->geomap() );
            }
        }
        else if (this->mechanicalProperties()->materialLaw() == "NeoHookean")
        {
            if ( !BuildCstPart )
            {
                auto const dFS_neohookean = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobian<2*nOrderDisplacement>(u,*this->mechanicalProperties());
                bilinearForm_PatternCoupled +=
                    integrate (_range=M_rangeMeshElements,
                               //_expr= trace(idv(coeffLame2)*dF*trans(grad(v)) ),
                               _expr= timeSteppingScaling*inner( dFS_neohookean, grad(v) ),
                               _quad=_Q<2*nOrderDisplacement+1>(),
                               _geomap=this->geomap() );
            }

        }
#if 0
        else if (this->mechanicalProperties()->materialLaw() == "MooneyRivlin")
        {
            if (!BuildCstPart)
            {
                bilinearForm_PatternCoupled +=
                    integrate (_range=M_rangeMeshElements,
                               _expr= trace( (dF*val(Sv_mooneyrivlin) + val(Fv)*dS_mooneyrivlin)*trans(grad(v)) ),
                               _geomap=this->geomap() );
            }
        }
#endif
    }
    else if ( M_modelName == "Elasticity-Large-Deformation" )
    {
        if (!BuildCstPart)
            bilinearForm_PatternCoupled +=
                integrate (_range=M_rangeMeshElements,
                           _expr= timeSteppingScaling*inner( dS, grad(v) ),
                           _geomap=this->geomap() );
    }
    else if ( M_modelName == "Elasticity" )
    {
        if (BuildCstPart)
            bilinearForm_PatternCoupled +=
                integrate (_range=M_rangeMeshElements,
                           _expr= timeSteppingScaling*inner( dS_elastic, grad(v) ),
                           _geomap=this->geomap() );
    }

    double timeElapsedBis = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateJacobian",
              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Newmark" )
        {
            if ( BuildCstPart )
            {
                if ( !this->useMassMatrixLumped() )
                {
                    bilinearForm_PatternDefault +=
                        integrate( _range=M_rangeMeshElements,
                                   _expr= M_timeStepNewmark->polyDerivCoefficient()*idv(rho)*inner( idt(u),id(v) ),
                                   _geomap=this->geomap() );
                }
                else
                {
                    J->close();
                    double thecoeff = M_timeStepNewmark->polyDerivCoefficient();
                    if ( this->massMatrixLumped()->size1() == J->size1() )
                        J->addMatrix( thecoeff, this->massMatrixLumped(), Feel::SUBSET_NONZERO_PATTERN );
                    else
                    {
                        auto vecAddDiagJ = this->backend()->newVector( J->mapRowPtr() );
                        auto uAddDiagJ = M_XhDisplacement->element( vecAddDiagJ, rowStartInVector );
                        uAddDiagJ = *M_vecDiagMassMatrixLumped;
                        uAddDiagJ.scale( thecoeff );
                        J->addDiagonal( vecAddDiagJ );
                    }
                }
            }
        } // Newmark
        else // bdf
        {
            if ( BuildCstPart )
            {
                CHECK( this->hasStartSubBlockSpaceIndex( "velocity" ) ) << "no SubBlockSpaceIndex velocity";
                size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");

                if ( !this->useMassMatrixLumped() )
                {
                    form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
                           _rowstart=rowStartInMatrix,
                           _colstart=colStartInMatrix+startBlockIndexVelocity ) +=
                        integrate( _range=M_rangeMeshElements,
                                   _expr= M_timeStepBdfVelocity->polyDerivCoefficient(0)*idv(rho)*inner(idt(u),id(v)),
                                   _geomap=this->geomap() );
                }
                else
                {
                    double thecoeff = M_timeStepBdfVelocity->polyDerivCoefficient(0);
                    for ( size_type i=0;i<M_XhDisplacement->nLocalDofWithoutGhost();++i)
                        J->add( J->mapRowPtr()->dofIdToContainerId(rowStartInMatrix)[i],
                                J->mapColPtr()->dofIdToContainerId(rowStartInMatrix+startBlockIndexVelocity)[i],
                                thecoeff*M_vecDiagMassMatrixLumped->operator()(i) );

                }
                form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
                       _rowstart=rowStartInMatrix+startBlockIndexVelocity,
                       _colstart=colStartInMatrix ) +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= M_timeStepBdfDisplacement->polyDerivCoefficient(0)*idv(rho)*inner(idt(u),id(v)),
                               _geomap=this->geomap() );
                form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
                       _rowstart=rowStartInMatrix+startBlockIndexVelocity,
                       _colstart=colStartInMatrix+startBlockIndexVelocity ) +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= -timeSteppingScaling*idv(rho)*inner(idt(u),id(v)),
                               _geomap=this->geomap() );
            }
        } // BDF
    }
    //--------------------------------------------------------------------------------------------------//
    // peusdo transient continuation
    if ( !BuildCstPart && data.hasInfo( "use-pseudo-transient-continuation" ) )
    {
        double pseudoTimeStepDelta = data.doubleInfo("pseudo-transient-continuation.delta");
        this->log("SolidMechanics","updateJacobian",(boost::format("pseudo-transient-continuation : delta=%1% s") %pseudoTimeStepDelta).str() );
        bilinearForm_PatternDefault +=
            integrate(_range=M_rangeMeshElements,
                      _expr=(1./pseudoTimeStepDelta)*inner(idt(u),id(u)),
                      _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    // incompressibility terms
    if (M_useDisplacementPressureFormulation && !BuildCstPart)
    {
        // define pressure field
        size_type blockIndexPressure = rowStartInVector+this->startSubBlockSpaceIndex("pressure");
        auto const p = M_XhPressure->element(X, blockIndexPressure);
        // assemble
        this->updateJacobianIncompressibilityTerms(u,p,J);
    }
    //--------------------------------------------------------------------------------------------------//
#if 0
    // viscoelastic terms
    // VERY OLD! : must be fix and updated
    this->updateJacobianViscoElasticityTerms(u,J);
#endif
    //--------------------------------------------------------------------------------------------------//
    // follower pressure bc
    if ( !BuildCstPart )
    {
        this->updateBCFollowerPressureJacobian( u, J, timeSteppingScaling );
    }
    //--------------------------------------------------------------------------------------------------//
    // robin bc
    if ( !BuildCstPart )
    {
        this->updateBCRobinJacobian( J, timeSteppingScaling );
    }
    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateJacobian","finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianIncompressibilityTerms( element_displacement_external_storage_type const& u,
                                                                              element_pressure_external_storage_type const& p,
                                                                              sparse_matrix_ptrtype& J) const
{
    using namespace Feel::vf;

    boost::mpi::timer thetimer;
    this->log("SolidMechanics","updateJacobianIncompressibilityTerms", "start " );

    auto mesh = M_XhDisplacement->mesh();
    auto const& v = this->fieldDisplacement();
    auto const& q = this->fieldPressure();

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    size_type startBlockIndexPressure = this->startSubBlockSpaceIndex("pressure");

    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();


    if ( M_modelName == "Hyper-Elasticity")
    {
        auto pFmtNLa = Feel::FeelModels::solidMecPressureFormulationMultiplierJacobianTrialPressure(u,p,*this->mechanicalProperties());
        form2( _test=M_XhDisplacement, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate ( _range=M_rangeMeshElements,
                        _expr= inner(pFmtNLa,grad(v) ),
                        _geomap=this->geomap() );

        // -dF*idv(p) or -d(F*C^{-1})*idv(p)
        auto pFmtNLb = /*-idv(p)**/Feel::FeelModels::solidMecPressureFormulationMultiplierJacobianTrialDisp(u,p,*this->mechanicalProperties());
        form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix ) +=
            integrate ( _range=M_rangeMeshElements,
                        _expr= inner(pFmtNLb,grad(v) ),
                        _geomap=this->geomap() );
    }
    else if ( M_modelName == "Elasticity-Large-Deformation" || M_modelName == "Elasticity")
    {
        form2( _test=M_XhDisplacement, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate ( _range=M_rangeMeshElements,
                        _expr= -trace(alpha_f*idt(p)*Id*trans(grad(v))),
                        _geomap=this->geomap() );
    }


    //--------------------------------------------------------------------------------------------//

    auto detJm1 = Feel::FeelModels::solidMecPressureFormulationConstraintJacobian(u,/*p,*/  *this->mechanicalProperties());

    form2( _test=M_XhPressure, _trial=M_XhDisplacement, _matrix=J,
           _rowstart=rowStartInMatrix+startBlockIndexPressure,
           _colstart=colStartInMatrix ) +=
        integrate( _range=M_rangeMeshElements,
                   _expr=detJm1*id(q),
                   _geomap=this->geomap() );


    if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
    {
        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix+startBlockIndexPressure,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate(_range=M_rangeMeshElements,
                      _expr= -(cst(1.)/idv(coeffLame1))*idt(p)*id(q),
                      _geomap=this->geomap() );
    }
    else
    {
        auto kappa = idv(this->mechanicalProperties()->fieldBulkModulus());
        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix+startBlockIndexPressure,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate( _range=M_rangeMeshElements,
                       _expr= -(cst(1.)/kappa)*idt(p)*id(q),
                       _geomap=this->geomap() );
    }


    //--------------------------------------------------------------------------------------------//

    double timeElapsed=thetimer.elapsed();
    this->log("SolidMechanics","updateJacobianIncompressibilityTerms",
              "finish in "+(boost::format("%1% s") % timeElapsed).str() );
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianViscoElasticityTerms( element_displacement_external_storage_type const& u, sparse_matrix_ptrtype& J) const
{
#if 0
    using namespace Feel::vf;

    if (this->verbose()) std::cout << "[SolidMechanics] : updateJacobianViscoElasticityTerms start\n";

    auto mesh = M_XhDisplacement->mesh();

    auto v=u;

    auto Et2 = 0.5*(gradt(u)+trans(gradt(u)) );// + 0.5*trans(gradv(u))*gradv(u);

    double gammav=0.01;

    form2( M_XhDisplacement, M_XhDisplacement, J)  +=
        integrate (_range=M_rangeMeshElements,
                   _expr=gammav*M_bdf_displ_struct->polyDerivCoefficient(0)*trace( Et2*trans(grad(v)) ),
                   _geomap=this->geomap() );


    if (this->verbose()) std::cout << "[SolidMechanics] : updateJacobianViscoElasticityTerms finish\n";

#endif
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateJacobianDofElimination( DataUpdateJacobian & data ) const
{
    if ( !this->hasMarkerDirichletBCelimination() ) return;

    this->log("SolidMechanics","updateJacobianDofElimination","start" );

    this->updateDofEliminationIds( "displacement", data );

    this->log("SolidMechanics","updateJacobianDofElimination","finish" );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinJacobian( sparse_matrix_ptrtype& J, double timeSteppingScaling ) const
{
    if ( this->M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                       _expr= timeSteppingScaling*expression1(d,this->symbolsExpr())(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureJacobian( element_displacement_external_storage_type const& u, sparse_matrix_ptrtype& J, double timeSteppingScaling ) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                       _expr= -timeSteppingScaling*expression(d,this->symbolsExpr())*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                       _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d,this->symbolsExpr()),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                       _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d,this->symbolsExpr())*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}

} // FeelModels

} // Feel



