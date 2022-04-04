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
SolidMechanics<ConvexType,BasisDisplacementType>::updateLinearPDE( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    if ( this->hasSolidEquation1dReduced() )
    {
        M_solid1dReduced->updateLinearPDE( data, mctx );
        return;
    }

    //const vector_ptrtype& X = data.
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _buildCstPart = data.buildCstPart();

    this->log("SolidMechanics","updateLinearPDE", "start" );
    this->timerTool("Solve").start();

    double timeSteppingScaling = 1.;
    if ( !this->isStationary() )
    {
        if ( M_timeStepping == "Theta" )
            timeSteppingScaling = M_timeStepThetaValue;
        data.addDoubleInfo( prefixvm(this->prefix(),"time-stepping.scaling"), timeSteppingScaling );
    }

    //---------------------------------------------------------------------------------------//

    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;
    bool BuildNonCstPart_TransientForm2Term = BuildNonCstPart;
    bool BuildNonCstPart_TransientForm1Term = BuildNonCstPart;
    bool BuildNonCstPart_SourceTerm = BuildNonCstPart;
    bool BuildNonCstPart_BoundaryNeumannTerm = BuildNonCstPart;
    if (this->useFSISemiImplicitScheme())
    {
        BuildNonCstPart_TransientForm2Term = BuildCstPart;
        BuildNonCstPart_TransientForm1Term=BuildCstPart;
        BuildNonCstPart_SourceTerm=BuildCstPart;
        BuildNonCstPart_BoundaryNeumannTerm=BuildCstPart;
    }
    if (M_timeStepNewmark->strategy()==TS_STRATEGY_DT_CONSTANT)
    {
        BuildNonCstPart_TransientForm2Term = BuildCstPart;
    }

    //---------------------------------------------------------------------------------------//

    auto mesh = M_XhDisplacement->mesh();
    auto Xh = M_XhDisplacement;

    size_type rowStartInMatrix = this->rowStartInMatrix();
    size_type colStartInMatrix = this->colStartInMatrix();
    size_type rowStartInVector = this->rowStartInVector();
#if 0
    auto u = Xh->element("u");//u = *X;
    auto v = Xh->element("v");
    for ( size_type k=0;k<M_Xh->nLocalDofWithGhost();++k )
        u(k) = X->operator()(rowStartInVector+k);
#else
    auto const& u = this->fieldDisplacement();
    auto const& v = this->fieldDisplacement();
#endif

    auto const& se = mctx.symbolsExpr();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A,
                                 _rowstart=this->rowStartInMatrix(),
                                 _colstart=this->colStartInMatrix() );
    auto linearFormDisp = form1( _test=Xh, _vector=F, _rowstart=rowStartInVector );

    //---------------------------------------------------------------------------------------//
    // Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // deformation tensor
    auto epst = sym(gradt(u));//0.5*(gradt(u)+trans(gradt(u)));
    auto eps = sym(grad(u));
    //---------------------------------------------------------------------------------------//
    // stress tensor
#if 0
    //#if (SOLIDMECHANICS_DIM==2) // cas plan
    double lll = 2*idv(coeffLame1)*idv(coeffLame2)/(idv(coeffLame1)+2*idv(coeffLame2));
    auto sigmat = lll*trace(epst)*Id + 2*idv(coeffLame2)*epst;
    auto sigmaold = lll*trace(epsold)*Id + 2*idv(coeffLame2)*epsold;
    //#endif
#endif
    //#if (SOLIDMECHANICS_DIM==3)  // cas 3d
    // auto sigmat = idv(coeffLame1)*trace(epst)*Id + 2*idv(coeffLame2)*epst;
    //auto sigmaold = idv(coeffLame1)*trace(epsold)*Id + 2*idv(coeffLame2)*epsold;
    //#endif
    //---------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------//

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = this->materialsProperties()->materialProperties( matName );
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );

            // internal force term
            if (BuildCstPart)
            {
                auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                auto lameSecondExpr = expr( matProperties.property( "Lame-second-parameter" ).exprScalar(), se );
                auto sigmat = lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst;

                if ( !physicSolidData->useDisplacementPressureFormulation() )
                {
                    bilinearFormDD +=
                        integrate (_range=range,
                                   _expr= timeSteppingScaling*inner(sigmat,grad(v)), // trace( sigmat*trans(grad(v))),
                                   _geomap=this->geomap() );
                }
                else
                {
                    bilinearFormDD +=
                        integrate (_range=range,
                                   //_expr= 2*lameSecondExpr*inner(epst,grad(v)),
                                   _expr= inner(2.0*lameSecondExpr*epst - (2.0/3.0)*lameSecondExpr*trace(epst)*Id ,grad(v)),
                                   _geomap=this->geomap() );

                    auto p = M_XhPressure->element();//*M_fieldPressure;
                    size_type startBlockIndexPressure = this->startSubBlockSpaceIndex("pressure");
                    form2( _test=Xh, _trial=M_XhPressure, _matrix=A,
                           _rowstart=rowStartInMatrix,
                           _colstart=colStartInMatrix+startBlockIndexPressure ) +=
                        integrate (_range=range,
                                   _expr= idt(p)*div(v),
                                   _geomap=this->geomap() );
                    form2( _test=M_XhPressure, _trial=Xh, _matrix=A,
                           _rowstart=rowStartInMatrix+startBlockIndexPressure,
                           _colstart=colStartInMatrix ) +=
                        integrate(_range=range,
                                  _expr= id(p)*divt(u),
                                  _geomap=this->geomap() );

                    auto K = expr( matProperties.property( "bulk-modulus" ).exprScalar(), se );
                    form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=A,
                           _rowstart=rowStartInMatrix+startBlockIndexPressure,
                           _colstart=colStartInMatrix+startBlockIndexPressure ) +=
                        integrate(_range=range,
                                  //_expr= -(cst(1.)/lameFirstExpr)*idt(p)*id(p),
                                  _expr= -(cst(1.)/K)*idt(p)*id(p),
                                  _geomap=this->geomap() );

                }
            }

            //--------------------------------------------------------------------------------------------------//
            // body forces
            //--------------------------------------------------------------------------------------------------//
            for ( auto const& bodyForce : physicSolidData->bodyForces() )
            {
                auto theExpr = bodyForce.expr( se );
                bool buildBodyForceTerm = BuildNonCstPart_SourceTerm;//theExpr.expression().isConstant()? buildCstPart : buildNonCstPart;
                if ( buildBodyForceTerm )
                {
                    linearFormDisp +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner( theExpr,id(v) ),
                                   _geomap=this->geomap() );
                }
            }

            //---------------------------------------------------------------------------------------//
            // discretisation acceleration term
            if ( !this->isStationary() )
            {
                auto const& densityProp = this->materialsProperties()->density( matName );
                auto densityExpr = expr( densityProp.expr(), se );

                if ( M_timeStepping == "Newmark" )
                {
                    if ( BuildNonCstPart_TransientForm2Term )
                    {
                        if ( !this->useMassMatrixLumped() )
                        {
                            form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                                integrate( _range=range,
                                           _expr= this->timeStepNewmark()->polySecondDerivCoefficient()*densityExpr*inner(idt(u),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            A->close();
                            double thecoeff = this->timeStepNewmark()->polyDerivCoefficient();
                            if ( this->massMatrixLumped()->size1() == A->size1() )
                                A->addMatrix( thecoeff, this->massMatrixLumped(), Feel::SUBSET_NONZERO_PATTERN );
                            else
                            {
                                auto vecAddDiagA = this->backend()->newVector( A->mapRowPtr() );
                                auto uAddDiagA = M_XhDisplacement->element( vecAddDiagA, rowStartInVector );
                                uAddDiagA = *M_vecDiagMassMatrixLumped;
                                uAddDiagA.scale( thecoeff );
                                A->addDiagonal( vecAddDiagA );
                            }
                        }
                    }
                    if ( BuildNonCstPart_TransientForm1Term )
                    {
                        auto polySecondDerivDisp = this->timeStepNewmark()->polySecondDeriv();
                        if ( !this->useMassMatrixLumped() )
                        {
                            form1( _test=Xh, _vector=F ) +=
                                integrate( _range=range,
                                           _expr= densityExpr*inner(idv(polySecondDerivDisp),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            if ( this->massMatrixLumped()->size1() == F->size() )
                            {
                                auto myvec = this->backend()->newVector(M_XhDisplacement);
                                *myvec = polySecondDerivDisp;
                                F->close();
                                F->addVector( myvec, this->massMatrixLumped() );
                            }
                            else
                            {
                                F->close();
                                auto uAddRhs = M_XhDisplacement->element( F, rowStartInVector );
                                auto uDiagMassMatrixLumped = M_XhDisplacement->element( M_vecDiagMassMatrixLumped );
                                uAddRhs.add( 1., element_product( uDiagMassMatrixLumped, polySecondDerivDisp ) );
                            }
                        }
                    }
                }
                else // if BDF
                {
                    CHECK( this->hasStartSubBlockSpaceIndex( "velocity" ) ) << "no SubBlockSpaceIndex velocity";
                    size_type startBlockIndexVelocity = this->startSubBlockSpaceIndex("velocity");

                    if ( BuildNonCstPart_TransientForm2Term )
                    {
                        if ( !this->useMassMatrixLumped() )
                        {
                            form2( _test=Xh, _trial=Xh, _matrix=A,
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
                                A->add( A->mapRowPtr()->dofIdToContainerId(rowStartInMatrix)[i],
                                        A->mapColPtr()->dofIdToContainerId(rowStartInMatrix+startBlockIndexVelocity)[i],
                                        thecoeff*M_vecDiagMassMatrixLumped->operator()(i) );
                        }
                        form2( _test=Xh, _trial=Xh, _matrix=A,
                               _rowstart=rowStartInMatrix+startBlockIndexVelocity,
                               _colstart=colStartInMatrix ) +=
                            integrate( _range=range,
                                       _expr= M_timeStepBdfDisplacement->polyDerivCoefficient(0)*densityExpr*inner(idt(u),id(v)),
                                       _geomap=this->geomap() );
                        form2( _test=Xh, _trial=Xh, _matrix=A,
                               _rowstart=rowStartInMatrix+startBlockIndexVelocity,
                               _colstart=colStartInMatrix+startBlockIndexVelocity ) +=
                            integrate( _range=range,
                                       _expr= -timeSteppingScaling*densityExpr*inner(idt(u),id(v)),
                                       _geomap=this->geomap() );
                    }
                    if ( BuildNonCstPart_TransientForm1Term )
                    {
                        auto rhsTimeStepVelocity = M_timeStepBdfVelocity->polyDeriv();
                        if ( !this->useMassMatrixLumped() )
                        {
                            form1( _test=Xh, _vector=F,
                                   _rowstart=rowStartInVector ) +=
                                integrate( _range=range,
                                           _expr= densityExpr*inner(idv(rhsTimeStepVelocity),id(v)),
                                           _geomap=this->geomap() );
                        }
                        else
                        {
                            F->close();
                            auto uAddRhs = M_XhDisplacement->element( F, rowStartInVector );
                            auto uDiagMassMatrixLumped = M_XhDisplacement->element( M_vecDiagMassMatrixLumped );
                            uAddRhs.add( 1., element_product( uDiagMassMatrixLumped, rhsTimeStepVelocity ) );
                        }
                        auto rhsTimeStepDisplacement = M_timeStepBdfDisplacement->polyDeriv();
                        form1( _test=Xh, _vector=F,
                               _rowstart=rowStartInVector+startBlockIndexVelocity ) +=
                            integrate( _range=range,
                                       _expr= densityExpr*inner(idv(rhsTimeStepDisplacement),id(v)),
                                       _geomap=this->geomap() );
                    }
                }
            }

        } // matName
    } // physics

    //---------------------------------------------------------------------------------------//

    for ( auto const& [bcName,bcData] : M_boundaryConditions->normalStress() )
    {
        if ( bcData->isFrameLagrangian() )
        {
            if ( bcData->isScalarExpr() )
            {
                if (BuildNonCstPart_BoundaryNeumannTerm)
                {
                    auto normalStessExpr = bcData->exprScalar( se );
                    linearFormDisp +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= timeSteppingScaling*normalStessExpr*inner( N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isVectorialExpr() )
            {
                if (BuildNonCstPart_BoundaryNeumannTerm)
                {
                    auto normalStessExpr = bcData->exprVectorial( se );
                    linearFormDisp +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= timeSteppingScaling*inner( normalStessExpr,id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isMatrixExpr() )
            {
                if (BuildNonCstPart_BoundaryNeumannTerm)
                {
                    auto normalStessExpr = bcData->exprMatrix( se );
                    linearFormDisp +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= timeSteppingScaling*inner( normalStessExpr*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
        }
        else
        {
            // TODO : eulerian frame
        }
    }

    for ( auto const& [bcName,bcData] : M_boundaryConditions->robin() )
    {
        if ( BuildNonCstPart )
        {
            auto exprRobin1 = bcData->expr1( se );
            bilinearFormDD +=
                integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                           _expr= timeSteppingScaling*exprRobin1*inner( idt(u) ,id(u) ),
                           _geomap=this->geomap() );
            auto exprRobin2 = bcData->expr2( se );
            linearFormDisp +=
                integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                           _expr= timeSteppingScaling*inner( exprRobin2 , id(u) ),
                           _geomap=this->geomap() );
        }
    }

    //---------------------------------------------------------------------------------------//
    double timeElapsed = this->timerTool("Solve").stop("createMesh");
    this->log("SolidMechanics","updateLinearPDE",
              (boost::format("finish in %1% s") % timeElapsed).str() );
}

template< typename ConvexType, typename BasisDisplacementType >
template <typename ModelContextType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::updateLinearPDEDofElimination( DataUpdateLinear & data, ModelContextType const& mctx ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    auto const& se = mctx.symbolsExpr();

    if ( this->hasSolidEquationStandard() )
    {
        if ( !M_boundaryConditions->hasTypeDofElimination() )
            return;

        auto Xh = this->functionSpace();
        auto mesh = Xh->mesh();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();

        M_boundaryConditions->applyDofEliminationLinear( bilinearForm, F, mesh, u, se );
    }

    if ( this->hasSolidEquation1dReduced() )
        M_solid1dReduced->updateLinearPDEDofElimination( data, mctx );
}


template< typename ConvexType, typename BasisDisplacementType >
template <typename ModelContextType>
void
SolidMechanics<ConvexType,BasisDisplacementType>::updateNewtonInitialGuess( DataNewtonInitialGuess & data, ModelContextType const& mctx ) const
{
    if ( !M_boundaryConditions->hasTypeDofElimination() )
        return;

    this->log("SolidMechanics","updateNewtonInitialGuess","start" );

    vector_ptrtype& U = data.initialGuess();
    auto Xh = this->functionSpace();
    auto mesh = this->mesh();
    auto u = this->functionSpaceDisplacement()->element( U, this->rowStartInVector() );
    auto const& se = mctx.symbolsExpr();

    M_boundaryConditions->applyNewtonInitialGuess( mesh, u, se );

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
    std::decay_t<decltype( mctx.field( FieldTag::pressure(this), "pressure" ) )> p;
    if ( this->hasDisplacementPressureFormulation() )
    {
        blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
        p = mctx.field( FieldTag::pressure(this), "pressure" );
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
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );

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

                    auto dFPK_disp = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobianTrialDisplacement(u,p,physicSolidData,matProperties,se);
                    bilinearForm_PatternCoupled +=
                        integrate( _range=range,
                                   _expr= timeSteppingScaling*inner( dFPK_disp, grad(v) ),
                                   _geomap=this->geomap() );
                    if ( physicSolidData->useDisplacementPressureFormulation() )
                    {
                        auto dFPK_pressure = solidMecFirstPiolaKirchhoffTensorJacobianTrialPressure(u,p,physicSolidData,matProperties,se);
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
                    if ( physicSolidData->equation() == "Hyper-Elasticity" && physicSolidData->materialModel() == "StVenantKirchhoff" )
                    {
                        auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
                               _rowstart=rowStartInMatrix+blockIndexPressure,
                               _colstart=colStartInMatrix+blockIndexPressure ) +=
                            integrate(_range=range,
                                      _expr= -(cst(1.)/lameFirstExpr)*idt(p)*id(p),
                                      _geomap=this->geomap() );
                    }
                    else if ( physicSolidData->equation() == "Hyper-Elasticity" || physicSolidData->equation() == "Elasticity" )
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

    for ( auto const& [bcName,bcData] : M_boundaryConditions->normalStress() )
    {
        if ( bcData->isFrameEulerian() )
        {
            if ( bcData->isScalarExpr() )
            {
                if ( buildNonCstPart )
                {
                    auto normalStessExpr = bcData->exprScalar( se );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*normalStessExpr*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isVectorialExpr() )
            {
                if ( buildNonCstPart )
                {
                    auto normalStessExpr = bcData->exprVectorial( se );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*normalStessExpr,id(u) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isMatrixExpr() )
            {
                if ( buildNonCstPart )
                {
                    auto normalStessExpr = bcData->exprMatrix( se );
                    bilinearForm_PatternCoupled +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*inner(Feel::FeelModels::solidMecGeomapEulerianJacobian(u)*normalStessExpr*N(),id(u) ),
                                   _geomap=this->geomap() );
                }
            }
        }
    }

    for ( auto const& [bcName,bcData] : M_boundaryConditions->robin() )
    {
        if ( buildNonCstPart )
        {
            auto exprRobin1 = bcData->expr1( se );
            bilinearForm_PatternDefault +=
                integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                           _expr= timeSteppingScaling*exprRobin1*inner( idt(u) ,id(u) ),
                           _geomap=this->geomap() );
        }
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
    std::decay_t<decltype( mctx.field( FieldTag::pressure(this), "pressure" ) )> p;
    if ( this->hasDisplacementPressureFormulation() )
    {
        // define pressure field
        blockIndexPressure = this->startSubBlockSpaceIndex("pressure");
        p = mctx.field( FieldTag::pressure(this), "pressure" );
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
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );

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

                    auto evalFPK = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor(u,p,physicSolidData,matProperties,se,timeSteppingScaling);
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

                if ( physicSolidData->equation() == "Hyper-Elasticity" && physicSolidData->materialModel() == "StVenantKirchhoff" )
                {
                    auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                    linearFormPressure +=
                        integrate(_range=range,
                                  _expr= -(cst(1.)/lameFirstExpr)*idv(p)*id(p),
                                  _geomap=this->geomap() );
                }
                else if ( physicSolidData->equation() == "Hyper-Elasticity" || physicSolidData->equation() == "Elasticity" )
                {
                    auto K = expr( matProperties.property( "bulk-modulus" ).exprScalar(), se );
                    linearFormPressure +=
                        integrate( _range=range,
                                   _expr= -(cst(1.)/K)*idv(p)*id(p),
                                   _geomap=this->geomap() );
                }

            }

            //--------------------------------------------------------------------------------------------------//
            // body forces
            //--------------------------------------------------------------------------------------------------//
            for ( auto const& bodyForce : physicSolidData->bodyForces() )
            {
                auto theExpr = bodyForce.expr( se );
                bool buildBodyForceTerm = buildCstPart;//theExpr.expression().isConstant()? buildCstPart : buildNonCstPart;
                if ( buildBodyForceTerm )
                {
                    linearFormDisplacement +=
                        integrate( _range=range,
                                   _expr= -timeSteppingScaling*inner( theExpr,id(v) ),
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

    for ( auto const& [bcName,bcData] : M_boundaryConditions->normalStress() )
    {
        if ( bcData->isFrameLagrangian() )
        {
            if ( bcData->isScalarExpr() )
            {
                if ( buildCstPart )
                {
                    auto normalStessExpr = bcData->exprScalar( se );
                    linearFormDisplacement +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*normalStessExpr*inner( N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isVectorialExpr() )
            {
                if ( buildCstPart )
                {
                    auto normalStessExpr = bcData->exprVectorial( se );
                    linearFormDisplacement +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*inner( normalStessExpr,id(v) ),
                                   _geomap=this->geomap() );
                }
            }
            else if ( bcData->isMatrixExpr() )
            {
                if ( buildCstPart )
                {
                    auto normalStessExpr = bcData->exprMatrix( se );
                    linearFormDisplacement +=
                        integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                                   _expr= -timeSteppingScaling*inner( normalStessExpr*N(),id(v) ),
                                   _geomap=this->geomap() );
                }
            }
        }
        else if ( bcData->isFrameEulerian() && buildNonCstPart )
        {
            if ( bcData->isScalarExpr() )
            {
                auto normalStessExpr = bcData->exprScalar( se );
                linearFormDisplacement +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*normalStessExpr*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                               _geomap=this->geomap() );
            }
            else if ( bcData->isVectorialExpr() )
            {
                auto normalStessExpr = bcData->exprVectorial( se );
                linearFormDisplacement +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*normalStessExpr,id(u) ),
                               _geomap=this->geomap() );
            }
            else if ( bcData->isMatrixExpr() )
            {
                auto normalStessExpr = bcData->exprMatrix( se );
                linearFormDisplacement +=
                    integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                               _expr= -timeSteppingScaling*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*normalStessExpr*N(),id(u) ),
                               _geomap=this->geomap() );
            }
        }
    }

    for ( auto const& [bcName,bcData] : M_boundaryConditions->robin() )
    {
        if ( buildNonCstPart )
        {
            auto exprRobin1 = bcData->expr1( se );
            auto exprRobin2 = bcData->expr2( se );
            linearFormDisplacement +=
                integrate( _range=markedfaces(this->mesh(),bcData->markers()),
                           _expr= timeSteppingScaling*inner( exprRobin1*idv(u) - exprRobin2 ,id(u) ),
                           _geomap=this->geomap() );
        }
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateResidual",
              "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}


} // namespace FeelModels
} // namespace Feel

#endif
