/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/


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
    if ( M_pdeType=="Elasticity" )
    {
        this->updateLinearElasticityGeneralisedAlpha( data );
    }
    else if ( M_pdeType=="Generalised-String" )
    {
        this->updateLinearGeneralisedStringGeneralisedAlpha( data );
    }
}

//-------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateLinearElasticityGeneralisedAlpha( DataUpdateLinear & data ) const
{
    using namespace Feel::vf;

    //const vector_ptrtype& X = data.
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool _buildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();


    this->log("SolidMechanics","updateLinearElasticityGeneralisedAlpha", "start" );
    this->timerTool("Solve").start();
    //boost::timer thetimer;

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

    //auto buzz1 = M_timeStepNewmark->previousUnknown();
    //---------------------------------------------------------------------------------------//
    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    auto const& rho = this->mechanicalProperties()->fieldRho();
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
    auto sigmat = idv(coeffLame1)*trace(epst)*Id + 2*idv(coeffLame2)*epst;
    //auto sigmaold = idv(coeffLame1)*trace(epsold)*Id + 2*idv(coeffLame2)*epsold;
    //#endif
    //---------------------------------------------------------------------------------------//

    double rho_s=0.8;
    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);

    //---------------------------------------------------------------------------------------//
    // internal force term
    if (BuildCstPart)
    {
        if ( !this->useDisplacementPressureFormulation() )
        {
            form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                integrate (_range=elements(mesh),
                           _expr= alpha_f*trace( sigmat*trans(grad(v))),
                           _geomap=this->geomap() );
        }
        else
        {
            form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                integrate (_range=elements(mesh),
                           _expr= 2*idv(coeffLame2)*inner(epst,grad(v)),
                           _geomap=this->geomap() );

            auto p = M_XhPressure->element();//*M_fieldPressure;
            size_type startBlockIndexPressure = this->startBlockIndexFieldsInMatrix().find("pressure")->second;
            form2( _test=Xh, _trial=M_XhPressure, _matrix=A,
                   _rowstart=rowStartInMatrix,
                   _colstart=colStartInMatrix+startBlockIndexPressure ) +=
                integrate (_range=elements(mesh),
                           _expr= idt(p)*div(v),
                           _geomap=this->geomap() );
            form2( _test=M_XhPressure, _trial=Xh, _matrix=A,
                   _rowstart=rowStartInMatrix+startBlockIndexPressure,
                   _colstart=colStartInMatrix ) +=
                integrate(_range=elements(mesh),
                          _expr= id(p)*divt(u),
                          _geomap=this->geomap() );
            form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=A,
                   _rowstart=rowStartInMatrix+startBlockIndexPressure,
                   _colstart=colStartInMatrix+startBlockIndexPressure ) +=
                integrate(_range=elements(mesh),
                          _expr= -(cst(1.)/idv(coeffLame1))*idt(p)*id(p),
                          _geomap=this->geomap() );
        }
    }
    //---------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if (!this->isStationary() && BuildNonCstPart_TransientForm2Term)
    {
        form2( _test=Xh, _trial=Xh, _matrix=A )  +=
            integrate( _range=elements(mesh),
                       _expr= this->timeStepNewmark()->polySecondDerivCoefficient()*idv(rho)*inner(idt(u),id(v)),
                       _geomap=this->geomap() );
    }
    //---------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if (!this->isStationary() && BuildNonCstPart_TransientForm1Term)
    {
        form1( _test=Xh, _vector=F ) +=
            integrate( _range=elements(mesh),
                       _expr= idv(rho)*inner(idv(this->timeStepNewmark()->polyDeriv()),id(v)),
                       _geomap=this->geomap() );
    }
    //---------------------------------------------------------------------------------------//
    // source term
    if ( BuildNonCstPart_SourceTerm )
    {
        this->updateSourceTermLinearPDE( F );
    }
    //---------------------------------------------------------------------------------------//
    // neumann boundary condition
    if (BuildNonCstPart_BoundaryNeumannTerm)
    {
        this->updateBCNeumannLinearPDE( F );
    }
    //---------------------------------------------------------------------------------------//
#if 0
    if (this->markerNameFSI().size()>0)
    {
        // neumann boundary condition with normal stress (fsi boundary condition)
        if (BuildNonCstPart)
        {
            form1( _test=Xh, _vector=F) +=
                integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                           _expr= -alpha_f*trans(idv(*M_fieldNormalStressFromFluid))*id(v),
                           _geomap=this->geomap() );
        }

        if ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
             this->couplingFSIcondition() == "nitsche" )
        {

            double gammaRobinFSI = this->gammaNitschFSI();
            double muFluid = this->muFluidFSI();
#if 0
            // integrate on ref with variables change
            auto Fa = eye<nDim,nDim>() + gradv(this->timeStepNewmark()->previousUnknown());
            auto J = det(Fa);
            auto B = inv(Fa);
            auto variablechange = J*norm2( B*N() );
            if (BuildNonCstPart )
            {
                form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= variablechange*gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(v))/hFace(),
                               _geomap=this->geomap() );
                auto robinFSIRhs = idv(this->timeStepNewmark()->polyFirstDeriv() ) + idv(this->velocityInterfaceFromFluid());
                form1( _test=Xh, _vector=F) +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= variablechange*gammaRobinFSI*muFluid*inner( robinFSIRhs, id(v))/hFace(),
                               _geomap=this->geomap() );

            }
#else
            if (BuildNonCstPart )//BuildCstPart)
            {
                MeshMover<mesh_type> mymesh_mover;
                mesh_ptrtype mymesh = this->mesh();
                mymesh_mover.apply( mymesh, this->timeStepNewmark()->previousUnknown() );

                form2( _test=Xh, _trial=Xh, _matrix=A )  +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(v))/hFace(),
                               _geomap=this->geomap() );

                auto robinFSIRhs = idv(this->timeStepNewmark()->polyFirstDeriv() ) + idv(this->fieldVelocityInterfaceFromFluid());
                form1( _test=Xh, _vector=F) +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= gammaRobinFSI*muFluid*inner( robinFSIRhs, id(v))/hFace(),
                               _geomap=this->geomap() );

                auto dispInv = this->fieldDisplacement().functionSpace()->element(-idv(this->timeStepNewmark()->previousUnknown()));
                mymesh_mover.apply( mymesh, dispInv );
            }
#endif

        }
    }
#endif
    //---------------------------------------------------------------------------------------//

    // robin condition (used in fsi blood flow as external tissue)
    if ( !this->markerRobinBC().empty() && !BuildCstPart )
    {
        this->updateBCRobinLinearPDE( A, F );
    }

    //---------------------------------------------------------------------------------------//

    if ( this->hasMarkerDirichletBCelimination() && BuildNonCstPart && _doBCStrongDirichlet)
    {
        this->updateBCDirichletStrongLinearPDE( A,F );
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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongLinearPDE(sparse_matrix_ptrtype& A, vector_ptrtype& F) const
{
    if ( this->isStandardModel() )
    {
        if ( !this->hasDirichletBC() ) return;

        auto Xh = this->functionSpace();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _rowstart=this->rowStartInMatrix(),
                                   _colstart=this->colStartInMatrix() );
        auto const& u = this->fieldDisplacement();
        for( auto const& d : this->M_bcDirichlet )
        {
            auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) );
            auto const& listMarkerFaces = std::get<0>( ret );
            auto const& listMarkerEdges = std::get<1>( ret );
            auto const& listMarkerPoints = std::get<2>( ret );
            if ( !listMarkerFaces.empty() )
                bilinearForm +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
            if ( !listMarkerEdges.empty() )
                bilinearForm +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
            if ( !listMarkerPoints.empty() )
                bilinearForm +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u,_rhs=F,_expr=expression(d),
                        _prefix=this->prefix() );
        }
        for ( auto const& bcDirComp : this->M_bcDirichletComponents )
        {
            ComponentType comp = bcDirComp.first;
            for( auto const& d : bcDirComp.second )
            {
                auto ret = detail::distributeMarkerListOnSubEntity(this->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d),comp ) );
                auto const& listMarkerFaces = std::get<0>( ret );
                auto const& listMarkerEdges = std::get<1>( ret );
                auto const& listMarkerPoints = std::get<2>( ret );
                if ( !listMarkerFaces.empty() )
                    bilinearForm +=
                        on( _range=markedfaces(this->mesh(), listMarkerFaces),
                            _element=u[comp],_rhs=F,_expr=expression(d),
                            _prefix=this->prefix() );
                if ( !listMarkerEdges.empty() )
                    bilinearForm +=
                        on( _range=markededges(this->mesh(), listMarkerEdges),
                            _element=u[comp],_rhs=F,_expr=expression(d),
                            _prefix=this->prefix() );
                if ( !listMarkerPoints.empty() )
                    bilinearForm +=
                        on( _range=markedpoints(this->mesh(), listMarkerPoints),
                            _element=u[comp],_rhs=F,_expr=expression(d),
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
                on( _range=markedfaces(Xh->mesh(),this->markerDirichletBCByNameId( "elimination",marker(d) ) ),
                    _element=u, _rhs=F, _expr=cst(0.),
                    _prefix=this->prefix() );
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannLinearPDE( vector_ptrtype& F ) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,marker(d)) ),
                       _expr= expression(d)*inner( N(),id(v) ),
                        _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,marker(d)) ),
                       _expr= inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,marker(d)) ),
                       _expr= inner( expression(d)*N(),id(v) ),
                       _geomap=this->geomap() );
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinLinearPDE( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
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
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= expression1(d)(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= inner( expression2(d) , id(u) ),
                       _geomap=this->geomap() );

    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermLinearPDE( vector_ptrtype& F ) const
{
    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=F,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
    {
        if ( marker(d).empty() )
            myLinearForm +=
                integrate( _range=elements(this->mesh()),
                           _expr= inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
        else
            myLinearForm +=
                integrate( _range=markedelements(this->mesh(),marker(d)),
                           _expr= inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
    }
}


} // FeelModels

} // Feel


