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
    vector_ptrtype& RBis = data.vectorUsedInStrongDirichlet();
    bool BuildCstPart = data.buildCstPart();
    bool _doBCStrongDirichlet = data.doBCStrongDirichlet();

    std::string sc=(BuildCstPart)?" (cst part)":" (non cst part)";
    this->log("SolidMechanics","updateJacobian", "start"+sc);
    this->timerTool("Solve").start();

    //boost::mpi::timer thetimer,thetimerBis;

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

    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);

    //--------------------------------------------------------------------------------------------------//

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    auto const& rho = this->mechanicalProperties()->fieldRho();
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    auto Fv = Id + alpha_f*gradv(u);
    auto dF = alpha_f*gradt(u);
    auto Ev = alpha_f*sym(gradv(u)) + alpha_f*alpha_f*0.5*trans(gradv(u))*gradv(u);
    auto dE = alpha_f*sym(gradt(u)) + alpha_f*alpha_f*0.5*(trans(gradv(u))*gradt(u) + trans(gradt(u))*gradv(u));
    auto Sv = idv(coeffLame1)*trace(Ev)*Id + 2*idv(coeffLame2)*Ev;
    auto dS = idv(coeffLame1)*trace(dE)*Id + 2*idv(coeffLame2)*dE;
    //case elastic
    auto dE_elastic = alpha_f*sym(gradt(u));
    auto dS_elastic = idv(coeffLame1)*trace(dE_elastic)*Id + 2*idv(coeffLame2)*dE_elastic;

    //--------------------------------------------------------------------------------------------------//
    // stress tensor terms
    //thetimerBis.restart();
    this->timerTool("Solve").start();

    if (M_pdeType=="Hyper-Elasticity")
    {
        if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
        {
            if (!BuildCstPart)
            {
                auto const dFS_neohookean = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobian<3*(nOrderDisplacement-1)>(u,*this->mechanicalProperties());
                bilinearForm_PatternCoupled +=
                    integrate (_range=elements(mesh),
                               //_expr= trace( (dF*val(Sv) + val(Fv)*dS)*trans(grad(v)) ),
                               //_expr=Feel::vf::FSI::stressStVenantKirchhoffJacobian(u,coeffLame1,coeffLame2), //le dernier
                               _expr= inner( dFS_neohookean, grad(v) ),
                               _geomap=this->geomap() );
            }
        }
        else if (this->mechanicalProperties()->materialLaw() == "NeoHookean")
        {
            if ( !BuildCstPart )
            {
                auto const dFS_neohookean = Feel::vf::FeelModels::solidMecFirstPiolaKirchhoffTensorJacobian<2*nOrderDisplacement>(u,*this->mechanicalProperties());
                bilinearForm_PatternCoupled +=
                    integrate (_range=elements(mesh),
                               //_expr= trace(idv(coeffLame2)*dF*trans(grad(v)) ),
                               _expr= inner( dFS_neohookean, grad(v) ),
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
                    integrate (_range=elements(mesh),
                               _expr= trace( (dF*val(Sv_mooneyrivlin) + val(Fv)*dS_mooneyrivlin)*trans(grad(v)) ),
                               _geomap=this->geomap() );
            }
        }
#endif
    }
    else if (M_pdeType=="Elasticity-Large-Deformation")
    {
        if (!BuildCstPart)
            bilinearForm_PatternCoupled +=
                integrate (_range=elements(mesh),
                           _expr= trace( dS*trans(grad(v)) ),
                           _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity")
    {
        if (BuildCstPart)
            bilinearForm_PatternCoupled +=
                integrate (_range=elements(mesh),
                           _expr= trace( dS_elastic*trans(grad(v)) ),
                           _geomap=this->geomap() );
    }

    double timeElapsedBis = this->timerTool("Solve").stop();
    this->log("SolidMechanics","updateJacobian",
              "build stresstensor term in "+(boost::format("%1% s") % timeElapsedBis ).str() );

    //--------------------------------------------------------------------------------------------------//
    // discretisation acceleration term
    if (BuildCstPart && !this->isStationary())
    {
        if ( !this->useMassMatrixLumped() )
        {
            bilinearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= M_timeStepNewmark->polyDerivCoefficient()*idv(rho)*inner( idt(u),id(v) ),
                           _geomap=this->geomap() );
        }
        else
        {
            J->close();
            double thecoeff = M_timeStepNewmark->polyDerivCoefficient()*this->mechanicalProperties()->cstRho();
            J->addMatrix( thecoeff, this->massMatrixLumped() );
        }
    }
    //--------------------------------------------------------------------------------------------------//
    // incompressibility terms
    if (M_useDisplacementPressureFormulation && !BuildCstPart)
    {
        // define pressure field
        size_type blockIndexPressure = rowStartInVector+this->startBlockIndexFieldsInMatrix().find("pressure")->second;
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
        this->updateBCFollowerPressureJacobian(u,J);
    }
    //--------------------------------------------------------------------------------------------------//
    // robin boundary condition (used in wavePressure3d as external tissue for arterial wall)
    if ( this->markerRobinBC().size() > 0 && !BuildCstPart )
    {
        this->updateBCRobinJacobian( J );
    }
    //--------------------------------------------------------------------------------------------------//
#if 0
    // fsi coupling using a robin boundary condition
    if (this->markerNameFSI().size()>0 && ( this->couplingFSIcondition() == "robin-robin" || this->couplingFSIcondition() == "robin-robin-genuine" ||
                                            this->couplingFSIcondition() == "nitsche" ) )
    {
        double gammaRobinFSI = this->gammaNitschFSI();
        double muFluid = this->muFluidFSI();
        if ( !BuildCstPart)
        {
#if 0
            // integrate on ref with variables change
            auto Fa = eye<nDim,nDim>() + gradv(this->timeStepNewmark()->previousUnknown());
            auto Ja = det(Fa);
            auto Ba = inv(Fa);
            auto variablechange = Ja*norm2( Ba*N() );
            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                           _expr= variablechange*gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(v))/hFace(),
                           _geomap=this->geomap() );

#else
            MeshMover<mesh_type> mymesh_mover;
            mesh_ptrtype mymesh = this->mesh();
            mymesh_mover.apply( mymesh, this->timeStepNewmark()->previousUnknown() );

            bilinearForm_PatternCoupled +=
                integrate( _range=markedfaces(mesh,this->markerNameFSI()),
                           _expr= gammaRobinFSI*muFluid*this->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(v))/hFace(),
                           _geomap=this->geomap() );

            auto dispInv = this->fieldDisplacement().functionSpace()->element(-idv(this->timeStepNewmark()->previousUnknown()));
            mymesh_mover.apply( mymesh, dispInv );
#endif
        }
    }
#endif
    //--------------------------------------------------------------------------------------------------//
    // strong Dirichlet bc
    if ( this->hasMarkerDirichletBCelimination() && !BuildCstPart && _doBCStrongDirichlet)
    {
        this->updateBCDirichletStrongJacobian( J, RBis );
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
    size_type startBlockIndexPressure = this->startBlockIndexFieldsInMatrix().find("pressure")->second;

    double alpha_f=M_genAlpha_alpha_f;
    double alpha_m=M_genAlpha_alpha_m;
    double gamma=0.5+alpha_m-alpha_f;
    double beta=0.25*(1+alpha_m-alpha_f)*(1+alpha_m-alpha_f);

    auto const& coeffLame1 = this->mechanicalProperties()->fieldCoeffLame1();
    auto const& coeffLame2 = this->mechanicalProperties()->fieldCoeffLame2();
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();


    if (M_pdeType=="Hyper-Elasticity")
    {
        auto pFmtNLa = Feel::vf::FeelModels::solidMecPressureFormulationMultiplierJacobianTrialPressure(u,p,*this->mechanicalProperties());
        form2( _test=M_XhDisplacement, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate ( _range=elements(mesh),
                        _expr= inner(pFmtNLa,grad(v) ),
                        _geomap=this->geomap() );

        // -dF*idv(p) or -d(F*C^{-1})*idv(p)
        auto pFmtNLb = /*-idv(p)**/Feel::vf::FeelModels::solidMecPressureFormulationMultiplierJacobianTrialDisp(u,p,*this->mechanicalProperties());
        form2( _test=M_XhDisplacement, _trial=M_XhDisplacement, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix ) +=
            integrate ( _range=elements(mesh),
                        _expr= inner(pFmtNLb,grad(v) ),
                        _geomap=this->geomap() );
    }
    else if (M_pdeType=="Elasticity-Large-Deformation" || M_pdeType=="Elasticity")
    {
        form2( _test=M_XhDisplacement, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate ( _range=elements(mesh),
                        _expr= -trace(alpha_f*idt(p)*Id*trans(grad(v))),
                        _geomap=this->geomap() );
    }


    //--------------------------------------------------------------------------------------------//

    auto detJm1 = Feel::vf::FeelModels::solidMecPressureFormulationConstraintJacobian(u,/*p,*/  *this->mechanicalProperties());

    form2( _test=M_XhPressure, _trial=M_XhDisplacement, _matrix=J,
           _rowstart=rowStartInMatrix+startBlockIndexPressure,
           _colstart=colStartInMatrix ) +=
        integrate( _range=elements(mesh),
                   _expr=detJm1*id(q),
                   _geomap=this->geomap() );


    if (this->mechanicalProperties()->materialLaw() == "StVenantKirchhoff")
    {
        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix+startBlockIndexPressure,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate(_range=elements(mesh),
                      _expr= -(cst(1.)/idv(coeffLame1))*idt(p)*id(q),
                      _geomap=this->geomap() );
    }
    else
    {
        auto kappa = idv(this->mechanicalProperties()->fieldBulkModulus());
        form2( _test=M_XhPressure, _trial=M_XhPressure, _matrix=J,
               _rowstart=rowStartInMatrix+startBlockIndexPressure,
               _colstart=colStartInMatrix+startBlockIndexPressure ) +=
            integrate( _range=elements(mesh),
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
        integrate (_range=elements(mesh),
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
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCDirichletStrongJacobian( sparse_matrix_ptrtype& J, vector_ptrtype& RBis ) const
{
    if ( !this->hasDirichletBC() ) return;

    //auto RBis = this->backend()->newVector( J->mapRowPtr() );
    auto bilinearForm_PatternCoupled = form2( _test=this->functionSpace(),_trial=this->functionSpace(),_matrix=J,
                                              _pattern=size_type(Pattern::COUPLED),
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
            bilinearForm_PatternCoupled +=
                on( _range=markedfaces(this->mesh(), listMarkerFaces),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
                    _prefix=this->prefix() );
        if ( !listMarkerEdges.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markededges(this->mesh(), listMarkerEdges),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
                    _prefix=this->prefix() );
        if ( !listMarkerPoints.empty() )
            bilinearForm_PatternCoupled +=
                on( _range=markedpoints(this->mesh(), listMarkerPoints),
                    _element=u,_rhs=RBis,_expr=0*one()/*Expression-idv(u)*/,
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
                bilinearForm_PatternCoupled +=
                    on( _range=markedfaces(this->mesh(), listMarkerFaces),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );
            if ( !listMarkerEdges.empty() )
                bilinearForm_PatternCoupled +=
                    on( _range=markededges(this->mesh(), listMarkerEdges),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );
            if ( !listMarkerPoints.empty() )
                bilinearForm_PatternCoupled +=
                    on( _range=markedpoints(this->mesh(), listMarkerPoints),
                        _element=u[comp],_rhs=RBis,_expr=cst(0.)/*Expression-idv(u)*/,
                        _prefix=this->prefix() );

        }
    }
}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCRobinJacobian( sparse_matrix_ptrtype& J) const
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
            integrate( _range=markedfaces(this->mesh(),marker(d)/*this->markerRobinBC()*/),
                       _expr= expression1(d)(0,0)*inner( idt(u) ,id(u) ),
                       _geomap=this->geomap() );

}

SOLIDMECHANICS_CLASS_TEMPLATE_DECLARATIONS
void
SOLIDMECHANICS_CLASS_TEMPLATE_TYPE::updateBCFollowerPressureJacobian( element_displacement_external_storage_type const& u, sparse_matrix_ptrtype& J) const
{
    if ( this->M_bcNeumannEulerianFrameScalar.empty() && this->M_bcNeumannEulerianFrameVectorial.empty() && this->M_bcNeumannEulerianFrameTensor2.empty() ) return;

    auto Xh = this->functionSpaceDisplacement();
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _rowstart=this->rowStartInMatrix(),
                               _colstart=this->colStartInMatrix() );

    for( auto const& d : this->M_bcNeumannEulerianFrameScalar )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -expression(d)*inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameVectorial )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d),id(u) ),
                       _geomap=this->geomap() );
    }
    for( auto const& d : this->M_bcNeumannEulerianFrameTensor2 )
    {
        bilinearForm +=
            integrate( _range=markedfaces(this->mesh(),marker(d)) ,
                       _expr= -inner(Feel::vf::FeelModels::solidMecGeomapEulerianJacobian(u)*expression(d)*N(),id(u) ),
                       _geomap=this->geomap() );
    }
}

} // FeelModels

} // Feel



