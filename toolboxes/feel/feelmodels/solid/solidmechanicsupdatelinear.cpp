/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */


#include <feel/feelmodels/solid/solidmechanics.hpp>

//#include <feel/feeldiscr/pchv.hpp>

namespace Feel
{
namespace FeelModels
{

//-------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDE( DataUpdateLinear & data ) const
{
    if ( M_modelName == "Generalised-String" )
    {
        this->updateLinearGeneralizedString( data );
        return;
    }

    //const vector_ptrtype& X = data.
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _buildCstPart = data.buildCstPart();

    this->log("SolidMechanics","updateLinearElasticityGeneralisedAlpha", "start" );
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

    auto se = this->symbolsExpr();
    //auto buzz1 = M_timeStepNewmark->previousUnknown();
    //---------------------------------------------------------------------------------------//
    // auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    // auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    // auto const& rho = this->mechanicalProperties()->fieldRho();
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
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( matName );

            // internal force term
            if (BuildCstPart)
            {
                auto lameFirstExpr = expr( matProperties.property( "Lame-first-parameter" ).exprScalar(), se );
                auto lameSecondExpr = expr( matProperties.property( "Lame-second-parameter" ).exprScalar(), se );
                auto sigmat = lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst;

                if ( !physicSolidData->useDisplacementPressureFormulation() )
                {
                    form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                        integrate (_range=range,
                                   _expr= timeSteppingScaling*inner(sigmat,grad(v)), // trace( sigmat*trans(grad(v))),
                                   _geomap=this->geomap() );
                }
                else
                {
                    form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                        integrate (_range=range,
                                   _expr= 2*lameSecondExpr*inner(epst,grad(v)),
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
                    form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=A,
                           _rowstart=rowStartInMatrix+startBlockIndexPressure,
                           _colstart=colStartInMatrix+startBlockIndexPressure ) +=
                        integrate(_range=range,
                                  _expr= -(cst(1.)/lameFirstExpr)*idt(p)*id(p),
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
    // source term
    if ( BuildNonCstPart_SourceTerm )
    {
        this->updateSourceTermLinearPDE( F, timeSteppingScaling );
    }
    // neumann boundary condition
    if (BuildNonCstPart_BoundaryNeumannTerm)
    {
        this->updateBCNeumannLinearPDE( F, timeSteppingScaling );
    }
    // robin bc
    if ( !BuildCstPart )
    {
        this->updateBCRobinLinearPDE( A, F, timeSteppingScaling );
    }

    //---------------------------------------------------------------------------------------//
    double timeElapsed = this->timerTool("Solve").stop("createMesh");
    this->log("SolidMechanics","updateLinearElasticityGeneralisedAlpha",
              (boost::format("finish in %1% s") % timeElapsed).str() );

} // updateLinearElasticityGeneralisedAlpha

//---------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearPDEDofElimination( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
 
    if ( this->isStandardModel() )
    {
        if ( !this->hasMarkerDirichletBCelimination() ) return;

        auto Xh = this->functionSpace();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();
        for( auto const& d : this->M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",name(d) ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u,_rhs=F,_expr=expression(d,this->symbolsExpr()),
                        _prefix=this->prefix() );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u,_rhs=F,_expr=expression(d,this->symbolsExpr()),
                        _prefix=this->prefix() );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u,_rhs=F,_expr=expression(d,this->symbolsExpr()),
                        _prefix=this->prefix() );
        }
        for ( auto const& bcDirComp : this->M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",name(d),comp ) );
                auto const& listMarkerFaces = std::get<0>( ret );
                auto const& listMarkerEdges = std::get<1>( ret );
                auto const& listMarkerPoints = std::get<2>( ret );
                if ( !listMarkerFaces.empty() )
                    bilinearForm +=
                        on( _range=markedfaces(this->mesh(), listMarkerFaces),
                            _element=u[comp],_rhs=F,_expr=expression(d,this->symbolsExpr()),
                            _prefix=this->prefix() );
                if ( !listMarkerEdges.empty() )
                    bilinearForm +=
                        on( _range=markededges(this->mesh(), listMarkerEdges),
                            _element=u[comp],_rhs=F,_expr=expression(d,this->symbolsExpr()),
                            _prefix=this->prefix() );
                if ( !listMarkerPoints.empty() )
                    bilinearForm +=
                        on( _range=markedpoints(this->mesh(), listMarkerPoints),
                            _element=u[comp],_rhs=F,_expr=expression(d,this->symbolsExpr()),
                            _prefix=this->prefix() );
            }
        }
    }
    else if ( this->is1dReducedModel() )
    {
        if ( this->M_bcDirichlet.empty() ) return;

        auto Xh = this->functionSpace1dReduced();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacementScal1dReduced();
        //WARNING : fixed at zero
        for( auto const& d : this->M_bcDirichlet )
            bilinearForm +=
                on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",name(d) ) ),
                    _element=u, _rhs=F, _expr=cst(0.),
                    _prefix=this->prefix() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F, double timeSteppingScaling ) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                       _expr= timeSteppingScaling*expression(d,this->symbolsExpr())*inner( N(),id(v) ),
                        _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                       _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                       _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr())*N(),id(v) ),
                       _geomap=this->geomap() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F, double timeSteppingScaling ) const
{
    if ( this->M_bcRobin.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );
    auto myLinearForm = form1( _test=Xh, _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& u = this->fieldDisplacement();

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                       _expr= timeSteppingScaling*expression1(d,this->symbolsExpr())(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                       _expr= timeSteppingScaling*inner( expression2(d,this->symbolsExpr()) , id(u) ),
                       _geomap=this->geomap() );

    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F, double timeSteppingScaling ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
    {
        auto rangeBodyForceUsed = markers(d).empty()? M_rangeMeshElements : markedelements(this->mesh(),markers(d));
        myLinearForm +=
            integrate( _range=rangeBodyForceUsed,
                       _expr= timeSteppingScaling*inner( expression(d,this->symbolsExpr()),id(v) ),
                       _geomap=this->geomap() );
    }
}


} // FeelModels

} // Feel


