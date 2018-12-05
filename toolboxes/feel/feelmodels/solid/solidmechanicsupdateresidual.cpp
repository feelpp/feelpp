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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateResidual( DataUpdateResidual & data ) const
{
    using namespace Feel::vf;

    const vector_ptrtype& X = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool _buildCstPart = data.buildCstPart();
    bool UseJacobianLinearTerms = data.useJacobianLinearTerms();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();


    std::string sc=(_buildCstPart)?" (cst part)":" (non cst part)";
    this->log("SolidMechanics","updateResidual", "start"+sc );
    this->timerTool("Solve").start();

    bool BuildNonCstPart = !_buildCstPart;
    bool BuildCstPart = _buildCstPart;
    bool BuildCstPart_BoundaryParoiMobile = BuildCstPart;
    if (this->useFSISemiImplicitScheme() )
    {
        BuildCstPart_BoundaryParoiMobile=BuildNonCstPart;
    }

    //--------------------------------------------------------------------------------------------------//

    mesh_ptrtype mesh = M_XhDisplacement->mesh();

    size_type rowStartInVector = this->rowStartInVector();
    auto linearFormDisplacement = form1( _test=M_XhDisplacement, _vector=R,_rowstart=rowStartInVector );

    //--------------------------------------------------------------------------------------------------//
    auto const u = M_XhDisplacement->element(X, rowStartInVector);
    auto const& v = this->fieldDisplacement();

    //--------------------------------------------------------------------------------------------------//

    //#if defined(SOLIDMECHANICS_VOLUME_FORCE)
    //auto f = SOLIDMECHANICS_VOLUME_FORCE(this->shared_from_this());
    //#endif
    //auto bcDef = SOLIDMECHANICS_BC(this->shared_from_this());

    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    //std::cout << "alpha_f " << alpha_f << "alpha_m " <<alpha_m <<"\n";
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);


    //--------------------------------------------------------------------------------------------------//

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    auto const& rho = this->mechanicalProperties()->fieldRho();
    //Identity Matrix
    auto Id = eye<nDim,nDim>();
    auto Fv = Id + alpha_f*gradv(u);
    auto Ev = alpha_f*sym(gradv(u)) + alpha_f*alpha_f*0.5*trans(gradv(u))*gradv(u);
    auto Sv = idv(coeffLame1)*trace(Ev)*Id + 2*idv(coeffLame2)*Ev;
    //case elastic
    auto Ev_elastic = alpha_f*sym(gradv(u));
    auto Sv_elastic = idv(coeffLame1)*trace(Ev_elastic)*Id + 2*idv(coeffLame2)*Ev_elastic;

    //--------------------------------------------------------------------------------------------------//
    //--------------------------------------------------------------------------------------------------//
    // stress tensor terms
    this->timerTool("Solve").start();

    if (M_pdeType=="Hyper-Elasticity")
    {
        if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
        {
            if (!BuildCstPart)
            {
                auto const FSv_neohookean = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor<3*(nOrderDisplacement-1)>(u,*this->mechanicalProperties());
                linearFormDisplacement +=
                    integrate( _range=M_rangeMeshElements,
                               //_expr= trace(val(Fv*Sv)*trans(grad(v))),
                               //_expr=trace(Feel::FeelModels::stressStVenantKirchhoff(u,coeffLame1,coeffLame2)*trans(grad(v))),
                               //_expr=Feel::FeelModels::stressStVenantKirchhoffResidual(u,coeffLame1,coeffLame2),// le dernier 
                               _expr= inner(FSv_neohookean,grad(v)),
                               _geomap=this->geomap() );
            }
        }
        else if (this->mechanicalProperties()->materialLaw() == "NeoHookean")
        {
            if (!BuildCstPart)
            {
                auto const FSv_neohookean = Feel::FeelModels::solidMecFirstPiolaKirchhoffTensor<2*nOrderDisplacement>(u,*this->mechanicalProperties());
                linearFormDisplacement +=
                    integrate( _range=M_rangeMeshElements,
                               //_expr= trace(val(FSv_neohookean)*trans(grad(v))),
                               _expr= inner(FSv_neohookean,grad(v)),
                               _quad=_Q<2*nOrderDisplacement+1>(),
                               _geomap=this->geomap() );
            }
        }
    } // if (M_pdeType=="Hyper-Elasticity")
    else if (M_pdeType=="Elasticity-Large-Deformation")
    {
        if (!BuildCstPart)
            linearFormDisplacement +=
                integrate( _range=M_rangeMeshElements,
                           _expr= trace(Sv*trans(grad(v))),
                           _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity")
    {
        if (!BuildCstPart && !UseJacobianLinearTerms)
            linearFormDisplacement +=
                integrate( _range=M_rangeMeshElements,
                           _expr= trace( Sv_elastic*trans(grad(v))),
                           _geomap=this->geomap() );
    }
    double timeElapsedBis = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateResidual",
              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

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
        this->updateSourceTermResidual( R );
    }
    //--------------------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if (!this->isStationary())
    {
        if (!BuildCstPart && !UseJacobianLinearTerms)
        {
            if ( !this->useMassMatrixLumped() )
            {
                linearFormDisplacement +=
                    integrate( _range=M_rangeMeshElements,
                               _expr= M_timeStepNewmark->polySecondDerivCoefficient()*idv(rho)*inner(idv(u),id(v)),
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
                    integrate( _range=M_rangeMeshElements,
                               _expr= -idv(rho)*inner(idv(polySecondDerivDisp),id(v)),
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
    }
    //--------------------------------------------------------------------------------------------------//
    // neumann bc
    if (BuildCstPart)
    {
        this->updateBCNeumannResidual( R );
    }
    // follower pressure bc
    if (!BuildCstPart)
    {
        this->updateBCFollowerPressureResidual( u,R );
    }

    //--------------------------------------------------------------------------------------------------//
#if 0
    // fsi bc
    if (this->markerNameFSI().size()>0)
    {
        // neumann boundary condition with normal stress (fsi boundary condition)
        if (BuildCstPart_BoundaryParoiMobile)
        {
            linearFormDisplacement +=
                integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                           _expr= alpha_f*trans(idv(*M_fieldNormalStressFromFluid))*id(v),
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
            auto Ja = det(Fa);
            auto Ba = inv(Fa);
            auto variablechange = Ja*norm2( Ba*N() );
            if (!BuildCstPart)
            {
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= variablechange*gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idv(u),id(v))/hFace(),
                               _geomap=this->geomap() );
            }
            if (!BuildCstPart )
            {
                auto robinFSIRhs = idv(this->timeStepNewmark()->polyFirstDeriv() ) + idv(this->velocityInterfaceFromFluid());
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= -variablechange*gammaRobinFSI*muFluid*inner( robinFSIRhs,id(v) )/hFace(),
                               _geomap=this->geomap() );
            }

#else
            MeshMover<mesh_type> mymesh_mover;
            mesh_ptrtype mymesh = this->mesh();
            mymesh_mover.apply( mymesh, this->timeStepNewmark()->previousUnknown() );

            if (!BuildCstPart)
            {
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idv(u),id(v))/hFace(),
                               _geomap=this->geomap() );
            }
            if (!BuildCstPart )
            {
                auto robinFSIRhs = idv(this->timeStepNewmark()->polyFirstDeriv() ) + idv(this->fieldVelocityInterfaceFromFluid());
                linearFormDisplacement +=
                    integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                               _expr= -gammaRobinFSI*muFluid*inner( robinFSIRhs,id(v) )/hFace(),
                               _geomap=this->geomap() );
            }

            auto dispInv = this->fieldDisplacement().functionSpace()->element(-idv(this->timeStepNewmark()->previousUnknown()));
            mymesh_mover.apply( mymesh, dispInv );
#endif
        } // robin-robin fsi

    }
#endif

    //--------------------------------------------------------------------------------------------------//
    // robin boundary condition (used in wavePressure3d as external tissue for arterial wall)
    //if ( this->markerRobinBC().size() > 0 && !BuildCstPart && !UseJacobianLinearTerms)
    if ( this->markerRobinBC().size() > 0 && !BuildCstPart )// && !UseJacobianLinearTerms)
    {
        this->updateBCRobinResidual( u, R );
    }
    //--------------------------------------------------------------------------------------------------//
    // dirichlet condition by elimination
    if (this->hasMarkerDirichletBCelimination() && !BuildCstPart)
    {
        //this->updateBCDirichletStrongResidual( R );
        R->close();

        auto resFeViewDisp = M_XhDisplacement->element(R,rowStartInVector);
        auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("displacement");
        auto const& dofsWithValueImposedDisp = ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )? itFindDofsWithValueImposed->second : std::set<size_type>();
        for ( size_type thedof : dofsWithValueImposedDisp )
            resFeViewDisp.set( thedof,0. );
        sync( resFeViewDisp, "=", dofsWithValueImposedDisp );
    }

    double timeElapsed = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateResidual",
              "finish"+sc+" in "+(boost::format("%1% s") % timeElapsed).str() );
}

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//

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

    if (M_pdeType=="Hyper-Elasticity")
    {
        linearFormDisplacement +=
            integrate( _range=M_rangeMeshElements,
                       //_expr= trace(val(-idv(p)*trans(InvFv))*trans(grad(v))),
                       //_expr=inner( -idv(p)*Feel::FeelModels::solidMecIncompressibilityPressure(u,p,*this->mechanicalProperties()),grad(v) ),
                       _expr=inner( /*-idv(p)**/Feel::FeelModels::solidMecPressureFormulationMultiplier(u,p,*this->mechanicalProperties()),grad(v) ),
                       _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity-Large-Deformation" || M_pdeType=="Elasticity")
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

//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateNewtonInitialGuess( vector_ptrtype& U ) const
{
    if ( !this->hasDirichletBC() ) return;

    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "start",
                                               this->worldComm(),this->verboseAllProc());

    auto Xh = this->functionSpace();
    auto mesh = this->mesh();

    if ( !Xh->worldsComm()[0]->isActive()) // only on Displacement Proc
        return;

    auto u = Xh->element( U, this->rowStartInVector() );

    for( auto const& d : M_bcDirichlet )
    {
        auto ret = detail::distributeMarkerListOnSubEntity(mesh,this->markerDirichletBCByNameId( "elimination",name(d) ) );
        auto const& listMarkerFaces = std::get<0>( ret );
        auto const& listMarkerEdges = std::get<1>( ret );
        auto const& listMarkerPoints = std::get<2>( ret );
        if ( !listMarkerFaces.empty() )
            u.on(_range=markedfaces(mesh,listMarkerFaces ),
                 _expr=expression(d) );
        if ( !listMarkerEdges.empty() )
            u.on(_range=markededges(mesh,listMarkerEdges),
                 _expr=expression(d) );
        if ( !listMarkerPoints.empty() )
            u.on(_range=markedpoints(mesh,listMarkerPoints),
                 _expr=expression(d) );
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
                           _expr=expression(d) );
            if ( !listMarkerEdges.empty() )
                u[comp].on(_range=markededges(mesh,listMarkerEdges),
                           _expr=expression(d) );
            if ( !listMarkerPoints.empty() )
                u[comp].on(_range=markedpoints(mesh,listMarkerPoints),
                           _expr=expression(d) );
        }
    }

    // synchronize velocity dof on interprocess
    auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("displacement");
    if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
        sync( u, "=", itFindDofsWithValueImposed->second );


    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".SolidMechanics","updateNewtonInitialGuess", "finish",
                                               this->worldComm(),this->verboseAllProc());
}

//--------------------------------------------------------------------------------------------------//

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCNeumannResidual(vector_ptrtype& R) const
{
    if ( this->M_bcNeumannScalar.empty() && this->M_bcNeumannVectorial.empty() && this->M_bcNeumannTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_bcNeumannScalar )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::SCALAR,name(d)) ),
                       _expr= -expression(d)*inner( N(),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannVectorial )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::VECTORIAL,name(d)) ),
                       _expr= -inner( expression(d),id(v) ),
                       _geomap=this->geomap() );
    for( auto const& d : this->M_bcNeumannTensor2 )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),this->markerNeumannBC(NeumannBCShape::TENSOR2,name(d)) ),
                       _expr= -inner( expression(d)*N(),id(v) ),
                       _geomap=this->geomap() );
}



SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinResidual(element_displacement_external_storage_type const& u, vector_ptrtype& R) const
{

    if ( this->M_bcRobin.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );

    // Warning : take only first component of expression1
    for( auto const& d : this->M_bcRobin )
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)/*this->markerRobinBC()*/),
                       _expr= inner( expression1(d)(0,0)*idv(u) - expression2(d) ,id(u) ),
                       _geomap=this->geomap() );

}




SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateSourceTermResidual( vector_ptrtype& R ) const
{

    if ( this->M_volumicForcesProperties.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    auto const& v = this->fieldDisplacement();

    for( auto const& d : this->M_volumicForcesProperties )
    {
        if ( markers(d).empty() )
            myLinearForm +=
                integrate( _range=M_rangeMeshElements,
                           _expr= -inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
        else
            myLinearForm +=
                integrate( _range=markedelements(this->mesh(),markers(d)),
                           _expr= -inner( expression(d),id(v) ),
                           _geomap=this->geomap() );
    }
}


SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureResidual( element_displacement_external_storage_type const& u, vector_ptrtype& R ) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto myLinearForm = form1( _test=this->functionSpaceDisplacement(), _vector=R,
                               _rowstart=this->rowStartInVector() );
    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -expression(d)*inner( Feel::FeelModels::solidMecGeomapEulerian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        myLinearForm +=
            integrate( _range=markedfaces(this->mesh(),markers(d)),
                       _expr= -inner( Feel::FeelModels::solidMecGeomapEulerian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}


} // FeelModels
} // Feel



