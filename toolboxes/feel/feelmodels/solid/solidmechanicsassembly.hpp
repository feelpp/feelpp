/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */
#ifndef FEELPP_TOOLBOXES_SOLIDMECHANICS_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_SOLIDMECHANICS_ASSEMBLY_HPP 1

//#include <feel/feelmodels/modelcore/diffsymbolicexpr.hpp>

#include <feel/feelmodels/modelvf/solidmecincompressibility.hpp>
#include <feel/feelmodels/modelvf/solidmecgeomapeulerian.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, typename BasisDisplacementType >
template <typename ModelContextType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
    if ( !this->hasDirichletBC() ) return;

    this->log("SolidMechanics","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto u = this->functionSpaceDisplacement()->element( U, this->rowStartInVector() );
    auto const& se = mctx.symbolsExpr();

    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=expression(d,se) );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges),
                 _expr=expression(d,se) );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(mesh,listMarkerPoints),
                 _expr=expression(d,se) );
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
                           _expr=expression(d,se) );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges),
                           _expr=expression(d,se) );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(mesh,listMarkerPoints),
                           _expr=expression(d,se) );
        }
    }

    // update info for synchronization
    this->updateDofEliminationIds( "displacement", data );

    this->log("SolidMechanics","updateNewtonInitialGuess","finish" );
}

template< typename ConvexType, typename BasisDisplacementType >
template <typename ModelContextType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::updateJacobian( DataUpdateJacobian & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& X = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst part)":" (non cst part)";
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

    auto const& u = mctx.field( FieldTag::displacement(this), "displacement" );
    auto const& v = this->fieldDisplacement();
    auto const& se = mctx.symbolsExpr();

    size_type blockIndexPressure = invalid_v<size_type>;
    std::decay_t<decltype(M_XhPressure->elementPtr(*X, blockIndexPressure))> p;
    if ( this->hasDisplacementPressureFormulation() )
    {
        // define pressure field
        blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
        p = M_XhPressure->elementPtr(*X, rowStartInVector+blockIndexPressure);
    }

    //--------------------------------------------------------------------------------------------------//

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }
    //--------------------------------------------------------------------------------------------------//

    //Identity Matrix
    // auto const Id = eye<nDim,nDim>();

    //--------------------------------------------------------------------------------------------------//
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = this->materialsProperties()->materialProperties( matName );
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );

            //--------------------------------------------------------------------------------------------------//
            // first Piola-Kirchhoff tensor
            //--------------------------------------------------------------------------------------------------//
            if ( physicSolidData->equation() == "Hyper-Elasticity" || physicSolidData->equation() == "Elasticity" )
            {
                bool buildFPK = buildNonCstPart;
                if ( physicSolidData->equation() == "Elasticity" &&
                     matProperties.property( "Lame-first-parameter" ).isEvaluable() && matProperties.property( "Lame-second-parameter" ).isEvaluable() )
                    buildFPK = buildCstPart;

                if ( buildFPK )
                {
                    // stress tensor terms
                    this->timerTool("Solve").start();

                    auto dFPK_disp = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobianTrialDisplacement(u,p,*physicSolidData,matProperties,se);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner( dFPK_disp, grad(v) ),
                                   _geomap=this->geomap() );
                    if ( physicSolidData->useDisplacementPressureFormulation() )
                    {
                        auto dFPK_pressure = solidMecFirstPiolaKirchhoffTensorJacobianTrialPressure(u,p,*physicSolidData,matProperties,se);
                        form2( _test=M_XhDisplacement, _trial=M_XhPressure, _matrix=J,
                               _rowstart=rowStartInMatrix,
                               _colstart=colStartInMatrix+blockIndexPressure ) +=
                            integrate ( _range=range,
                                        _expr= inner(dFPK_pressure,grad(v) ),
                                        _geomap=this->geomap() );
                    }

                    double timeElapsedBis = this->timerTool("Solve").stop();
                    this->log("SolidMechanics","updateJacobian",
                              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );
                }
            }

            //--------------------------------------------------------------------------------------------------//
            // constraint on compressibility
            //--------------------------------------------------------------------------------------------------//
            if ( physicSolidData->useDisplacementPressureFormulation() )
            {
                if ( physicSolidData->equation() == "Elasticity" )
                {
                    if ( buildCstPart )
                    {
                        form2( _test=M_XhPressure, _trial=M_XhDisplacement, _matrix=J,
                               _rowstart=rowStartInMatrix+blockIndexPressure,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=range,
                                       _expr= id(p)*divt(u),
                                       _geomap=this->geomap() );
                    }
                }
                else
                {
                    if ( !buildCstPart )
                    {
                        auto detJm1 = Feel::FeelModels::solidMecPressureFormulationConstraintJacobian(u,/*p,*/  *physicSolidData);
                        form2( _test=M_XhPressure, _trial=M_XhDisplacement, _matrix=J,
                               _rowstart=rowStartInMatrix+blockIndexPressure,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=range,
                                       _expr=detJm1*id(p),
                                       _geomap=this->geomap() );
                    }
                }

                if ( !buildCstPart )
                {
                    if ( ( physicSolidData->equation() == "Hyper-Elasticity" && physicSolidData->materialModel() == "StVenantKirchhoff" ) ||  physicSolidData->equation() == "Elasticity" )
                    {
                        auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
                               _rowstart=rowStartInMatrix+blockIndexPressure,
                               _colstart=colStartInMatrix+blockIndexPressure ) +=
                            integrate(_range=range,
                                      _expr= -(cst(1.)/lameFirstExpr)*idt(p)*id(p),
                                      _geomap=this->geomap() );
                    }
                    else if ( physicSolidData->equation() == "Hyper-Elasticity" )
                    {
                        auto K = expr( matProperties.property( "bulk-modulus" ).exprScalar(), se );
                        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
                               _rowstart=rowStartInMatrix+blockIndexPressure,
                               _colstart=colStartInMatrix+blockIndexPressure ) +=
                            integrate( _range=range,
                                       _expr= -(cst(1.)/K)*idt(p)*id(p),
                                       _geomap=this->geomap() );
                    }
                }
            }
            //--------------------------------------------------------------------------------------------------//
            // discretisation acceleration term
            //--------------------------------------------------------------------------------------------------//
            if ( !this->isStationary() )
            {
                auto const& densityProp = this->materialsProperties()->density( matName );
                auto densityExpr = expr( densityProp.expr(), se );

                if ( M_timeStepping == "Newmark" )
                {
                    if ( buildCstPart )
                    {
                        if ( !this->useMassMatrixLumped() )
                        {
                            bilinearForm_PatternDefault +=
                                integrate( _range=range,
                                           _expr= M_timeStepNewmark->polyDerivCoefficient()*densityExpr*inner( idt(u),id(v) ),
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
                    if ( buildCstPart )
                    {
                        CHECK( this->hasStartSubBlockSpaceIndex( "velocity" ) ) << "no SubBlockSpaceIndex velocity";
                        size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");

                        if ( !this->useMassMatrixLumped() )
                        {
                            form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
                                   _rowstart=rowStartInMatrix,
                                   _colstart=colStartInMatrix+startBlockIndexVelocity ) +=
                                integrate( _range=range,
                                           _expr= M_timeStepBdfVelocity->polyDerivCoefficient(0)*densityExpr*inner(idt(u),id(v)),
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
                            integrate( _range=range,
                                       _expr= M_timeStepBdfDisplacement->polyDerivCoefficient(0)*densityExpr*inner(idt(u),id(v)),
                                       _geomap=this->geomap() );
                        form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
                               _rowstart=rowStartInMatrix+startBlockIndexVelocity,
                               _colstart=colStartInMatrix+startBlockIndexVelocity ) +=
                            integrate( _range=range,
                                       _expr= -timeSteppingScaling*densityExpr*inner(idt(u),id(v)),
                                       _geomap=this->geomap() );
                    }
                } // BDF
            } // if ( !this->isStationary() )

        } //matName
    } // physics

    //--------------------------------------------------------------------------------------------------//
    // peusdo transient continuation
    if ( !buildCstPart && data.hasInfo( "use-pseudo-transient-continuation" ) )
    {
        double pseudoTimeStepDelta = data.doubleInfo("pseudo-transient-continuation.delta");
        this->log("SolidMechanics","updateJacobian",(boost::format("pseudo-transient-continuation : delta=%1% s") %pseudoTimeStepDelta).str() );
        bilinearForm_PatternDefault +=
            integrate(_range=M_rangeMeshElements,
                      _expr=(1./pseudoTimeStepDelta)*inner(idt(u),id(u)),
                      _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//
#if 0
    // viscoelastic terms
    // VERY OLD! : must be fix and updated
    this->updateJacobianViscoElasticityTerms(u,J);
#endif
    //--------------------------------------------------------------------------------------------------//
    // follower pressure bc
    if ( !buildCstPart )
    {
        for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                           _expr= -timeSteppingScaling*expression(d,se)*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                           _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d,se),id(u) ),
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(this->mesh(),markers(d)) ,
                           _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d,se)*N(),id(u) ),
                           _geomap=this->geomap() );
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // robin bc
    if ( !buildCstPart )
    {
        // Warning : take only first component of expression1
        for( auto const& d : this->M_bcRobin )
            bilinearForm_PatternDefault +=
                integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                           _expr= timeSteppingScaling*expression1(d,se)(0,0)*inner( idt(u) ,id(u) ),
                           _geomap=this->geomap() );
    }
    //--------------------------------------------------------------------------------------------------//

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateJacobian","finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}

template< typename ConvexType, typename BasisDisplacementType >
template <typename ModelContextType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::updateResidual( DataUpdateResidual & data, ModelContextType const& mctx ) const
{
    const vector_ptrtype& X = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst part)":" (non cst part)";
    this->log("SolidMechanics","updateResidual", "start"+sc );
    this->timerTool("Solve").start();

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
    auto const& u = mctx.field( FieldTag::displacement(this), "displacement" );
    auto const& v = this->fieldDisplacement();
    auto const& se = mctx.symbolsExpr();

    //--------------------------------------------------------------------------------------------------//

    size_type blockIndexPressure = invalid_v<size_type>;
    std::decay_t<decltype(M_XhPressure->elementPtr(*X, blockIndexPressure))> p;
    if ( this->hasDisplacementPressureFormulation() )
    {
        // define pressure field
        blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
        p = M_XhPressure->elementPtr(*X, rowStartInVector+blockIndexPressure);
    }

    //Identity Matrix
    // auto Id = eye<nDim,nDim>();

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
            // first Piola-Kirchhoff tensor
            //--------------------------------------------------------------------------------------------------//
            if ( physicSolidData->equation() == "Hyper-Elasticity" || physicSolidData->equation() == "Elasticity" )
            {
                bool buildFPK = buildNonCstPart;
                if ( physicSolidData->equation() == "Elasticity" &&
                     matProperties.property( "Lame-first-parameter" ).isEvaluable() && matProperties.property( "Lame-second-parameter" ).isEvaluable() ) // TODO : to improve
                    buildFPK = buildNonCstPart && !UseJacobianLinearTerms;

                if ( buildFPK )
                {
                    this->timerTool("Solve").start();

                    auto evalFPK = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,p,*physicSolidData,matProperties,se,timeSteppingScaling);
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   _expr= inner(evalFPK,grad(v)),
                                   _geomap=this->geomap() );

                    double timeElapsedBis = this->timerTool("Solve").stop();
                    this->log("SolidMechanics","updateResidual",
                              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );
                }
            }

            //--------------------------------------------------------------------------------------------------//
            // constraint on compressibility
            //--------------------------------------------------------------------------------------------------//
            if ( physicSolidData->useDisplacementPressureFormulation() && buildNonCstPart )
            {
                auto linearFormPressure = form1( _test=M_XhPressure, _vector=R,_rowstart=rowStartInVector+blockIndexPressure );
                if ( physicSolidData->equation() == "Elasticity" )
                {
                    if ( buildNonCstPart && !UseJacobianLinearTerms )
                    {
                        linearFormPressure +=
                            integrate( _range=range,
                                       _expr= id(p)*divv(u),
                                       _geomap=this->geomap() );
                    }
                }
                else
                {
                    linearFormPressure +=
                        integrate( _range=range,
                                   _expr= Feel::FeelModels::solidMecPressureFormulationConstraint(u,*physicSolidData)*id(p),
                                   _geomap=this->geomap() );
                }

                if ( ( physicSolidData->equation() == "Hyper-Elasticity" && physicSolidData->materialModel() == "StVenantKirchhoff" ) ||  physicSolidData->equation() == "Elasticity" )
                {
                    auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                    linearFormPressure +=
                        integrate(_range=range,
                                  _expr= -(cst(1.)/lameFirstExpr)*idv(p)*id(p),
                                  _geomap=this->geomap() );
                }
                else if ( physicSolidData->equation() == "Hyper-Elasticity" )
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
                    if (buildNonCstPart && !UseJacobianLinearTerms)
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
                                *myvec = unwrap_ptr(u);
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
                    if (buildCstPart)
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
                        if ( buildNonCstPart )
                        {
                            auto const curVel = M_XhDisplacement->element(X, rowStartInVector+startBlockIndexVelocity);
                            form1( _test=M_XhDisplacement, _vector=R,
                                   _rowstart=rowStartInVector+startBlockIndexVelocity ) +=
                                integrate( _range=range,
                                           _expr= -timeSteppingScaling*densityExpr*inner(idv(curVel),id(v)),
                                           _geomap=this->geomap() );
                        }
                    }
                    else if (buildNonCstPart && !UseJacobianLinearTerms)
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

                    if ( buildCstPart && !timeSteppingEvaluateResidualWithoutTimeDerivative )
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
    //--------------------------------------------------------------------------------------------------//
#if 0
    // viscoelastic terms
    // VERY OLD! : must be fix and updated
    this->updateResidualViscoElasticityTerms(u,R);
#endif
    //--------------------------------------------------------------------------------------------------//
    // source term
    if (buildCstPart)
    {
        for( auto const& d : this->M_volumicForcesProperties )
        {
            auto rangeBodyForceUsed = markers(d).empty()? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
            linearFormDisplacement +=
                integrate( _range=rangeBodyForceUsed,
                           _expr= -timeSteppingScaling*inner( expression(d,se),id(v) ),
                           _geomap=this->geomap() );
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // neumann bc
    if (buildCstPart)
    {
        for( auto const& d : this->M_bcNeumannScalar )
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                           _expr= -timeSteppingScaling*expression(d,se)*inner( N(),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannVectorial )
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                           _expr= -timeSteppingScaling*inner( expression(d,se),id(v) ),
                           _geomap=this->geomap() );
        for( auto const& d : this->M_bcNeumannTensor2 )
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                           _expr= -timeSteppingScaling*inner( expression(d,se)*N(),id(v) ),
                           _geomap=this->geomap() );
    }
    // follower pressure bc
    if (buildNonCstPart)
    {
        for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
        {
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),markers(d)),
                           _expr= -timeSteppingScaling*expression(d,se)*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
        {
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),markers(d)),
                           _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d,se),id(u) ),
                           _geomap=this->geomap() );
        }
        for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
        {
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),markers(d)),
                           _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d,se)*N(),id(u) ),
                           _geomap=this->geomap() );
        }
    }
    // robin bc
    if ( buildNonCstPart )// && !UseJacobianLinearTerms)
    {
        // Warning : take only first component of expression1
        for( auto const& d : this->M_bcRobin )
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                           _expr= timeSteppingScaling*inner( expression1(d,se)(0,0)*idv(u) - expression2(d,se) ,id(u) ),
                           _geomap=this->geomap() );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateResidual",
              "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}


} // namespace FeelModels
} // namespace Feel

#endif
