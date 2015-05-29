/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels2/fluid/fluidmecbase.hpp>

#include <feel/feelvf/vf.hpp>


namespace Feel
{
namespace FeelModels
{

FLUIDMECHANICSBASE_CLASS_TEMPLATE_DECLARATIONS
void
FLUIDMECHANICSBASE_CLASS_TEMPLATE_TYPE::updateOseen( sparse_matrix_ptrtype& A , vector_ptrtype& F, bool _BuildCstPart,
                                                     sparse_matrix_ptrtype& A_extended, bool _BuildExtendedPart,
                                                     bool _doClose, bool _doBCStrongDirichlet ) const
{
#if defined(FEELMODELS_FLUID_BUILD_LINEAR_CODE)
    using namespace Feel::vf;

    std::string sc=(_BuildCstPart)?" (build cst part)":" (build non cst part)";
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateOseen", "start"+sc,
                                               this->worldComm(),this->verboseAllProc());
    boost::mpi::timer thetimer;

    bool BuildNonCstPart = !_BuildCstPart;
    bool BuildCstPart = _BuildCstPart;
    bool BuildNonCstPart_ConvectiveTerm = BuildNonCstPart;
    bool BuildNonCstPart_Form2TransientTerm = BuildNonCstPart;
    bool BuildNonCstPart_Form1TransientTerm = BuildNonCstPart;
    //bool BuildNonCstPart_SourceTerm = BuildNonCstPart;
    bool BuildNonCstPart_BoundaryNeumannTerm = BuildNonCstPart;
    if ( this->timeStepBase()->strategy()==TS_STRATEGY_DT_CONSTANT )
    {
        BuildNonCstPart_Form2TransientTerm=BuildCstPart;
    }
    if (this->useFSISemiImplicitScheme())
    {
        BuildNonCstPart_ConvectiveTerm=BuildCstPart;
        BuildNonCstPart_Form2TransientTerm=BuildCstPart;
        BuildNonCstPart_Form1TransientTerm=BuildCstPart;
        //BuildNonCstPart_SourceTerm=BuildCstPart;
        BuildNonCstPart_BoundaryNeumannTerm=BuildCstPart;
    }

    auto mesh = this->mesh();
    auto Xh = this->functionSpace();

    auto const& U = this->fieldVelocityPressure();
    auto u = U.template element<0>();
    auto v = U.template element<0>();
    auto p = U.template element<1>();
    auto q = U.template element<1>();

    //Deformations tensor (trial)
    auto deft = sym(gradt(u));
    //Deformations tensor (test)
    //auto def = sym(gradv(u));
    //Identity Matrix
    auto const Id = eye<nDim,nDim>();
    // Strain tensor (trial)
    auto Sigmat = -idt(p)*Id + 2*idv(*M_P0Mu)*deft;

    //--------------------------------------------------------------------------------------------------//

    auto rowStartInMatrix = this->rowStartInMatrix();
    auto colStartInMatrix = this->colStartInMatrix();
    auto rowStartInVector = this->rowStartInVector();
    auto bilinearForm_PatternDefault = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::DEFAULT),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );
    auto bilinearForm_PatternCoupled = form2( _test=Xh,_trial=Xh,_matrix=A,
                                              _pattern=size_type(Pattern::COUPLED),
                                              _rowstart=rowStartInMatrix,
                                              _colstart=colStartInMatrix );

    auto myLinearForm =form1( _test=Xh, _vector=F,
                              _rowstart=rowStartInVector );

    //--------------------------------------------------------------------------------------------------//
    // sigma : grad(v) sur Omega
#if 1
    if (BuildCstPart)
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= trace(Sigmat*trans(grad(v))),
                       _geomap=this->geomap() );
#else
    form2( Xh, Xh, A ) +=
        integrate( _range=elements(mesh),
                   _expr= 2*idv(*M_P0Mu)*trace(trans(deft)*grad(v)),
                   _geomap=this->geomap() );
    form2( Xh, Xh, A ) +=
        integrate( _range=elements(mesh),
                   _expr= -div(v)*idt(p),
                   _geomap=this->geomap() );
#endif

    //--------------------------------------------------------------------------------------------------//
    // incompressibility term
    if (BuildCstPart)
        bilinearForm_PatternCoupled +=
            integrate( _range=elements(mesh),
                       _expr= divt(u)*id(q),
                       _geomap=this->geomap() );

    //--------------------------------------------------------------------------------------------------//
    // volume force
    //if (BuildNonCstPart_SourceTerm)
    this->updateSourceTermLinearPDE(F,BuildCstPart);

    if (M_haveSourceAdded && BuildNonCstPart)
    {
        myLinearForm +=
            integrate( _range=elements(mesh),
                       _expr= trans(idv(*M_SourceAdded))*id(v),
                       _geomap=this->geomap() );
    }

    //--------------------------------------------------------------------------------------------------//
    // define pressure cst
    if ( this->definePressureCst() )
    {
        if ( this->definePressureCstMethod() == "penalisation" && BuildCstPart  )
        {
            double beta = this->definePressureCstPenalisationBeta();
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(Xh->mesh()),
                           _expr=beta*idt(p)*id(q),
                           _geomap=this->geomap() );
        }
        if ( this->definePressureCstMethod() == "lagrange-multiplier" )
        {
            CHECK( this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm") != this->startDofIndexFieldsInMatrix().end() )
                << " start dof index for define-pressure-cst-lm is not present\n";
            size_type startDofIndexDefinePressureCstLM = this->startDofIndexFieldsInMatrix().find("define-pressure-cst-lm")->second;

#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_ONLY_ON_BOUNDARY)
            auto therange=boundaryfaces(mesh);
#else
            auto therange=elements(mesh);
#endif
            auto lambda = M_XhMeanPressureLM->element();
            if (BuildCstPart)
            {
                form2( _test=Xh, _trial=M_XhMeanPressureLM, _matrix=A,
                       _rowstart=this->rowStartInMatrix(),
                       _colstart=this->colStartInMatrix()+startDofIndexDefinePressureCstLM ) +=
                    integrate( _range=therange,//elements(mesh),
                               _expr= id(p)*idt(lambda) /*+ idt(p)*id(lambda)*/,
                               _geomap=this->geomap() );

                form2( _test=M_XhMeanPressureLM, _trial=Xh, _matrix=A,
                       _rowstart=this->rowStartInMatrix()+startDofIndexDefinePressureCstLM,
                       _colstart=this->colStartInMatrix() ) +=
                    integrate( _range=therange,//elements(mesh),
                               _expr= + idt(p)*id(lambda),
                               _geomap=this->geomap() );
            }

#if defined(FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE)
            if (BuildNonCstPart)
            {
                if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateOseen", "also add nonzero MEANPRESSURE",
                                                           this->worldComm(),this->verboseAllProc());
                form1( _test=M_XhMeanPressureLM, _vector=F,
                       _rowstart=this->rowStartInMatrix()+startDofIndexDefinePressureCstLM ) +=
                    integrate( _range=therange,//elements(mesh),
                               _expr= FLUIDMECHANICS_USE_LAGRANGEMULTIPLIER_MEANPRESSURE(this->shared_from_this())*id(lambda),
                               _geomap=this->geomap() );
            }
#endif
        } // if ( this->definePressureCstMethod() == "lagrange-multiplier" )
    } // if ( this->definePressureCst() )


    //--------------------------------------------------------------------------------------------------//
    // neumann condition
    if (BuildNonCstPart_BoundaryNeumannTerm)
    {
        this->updateBCNeumannLinearPDE( F );
    }

    //--------------------------------------------------------------------------------------------------//
    //pressure fix condition
    if (BuildCstPart && !this->markerPressureBC().empty() )
    {
        bilinearForm_PatternCoupled +=
            integrate( _range=markedfaces(mesh,this->markerPressureBC() ),
                       _expr= -trans(2*idv(*M_P0Mu)*deft*N())*id(v),
                       _geomap=this->geomap() );
    }
    if (BuildNonCstPart)
    {
        this->updateBCPressureLinearPDE( F );
    }

    //--------------------------------------------------------------------------------------------------//

#if 0
    // robin fsi
    if ( this->couplingFSIcondition()=="robin" && BuildNonCstPart )
    {
        //double alphaRobinFSI = 1.0;
        double alphaRobinFSI = 0.1*1.2/this->timeStepBDF()->timeStep() + (0.1*3e6/((1-0.3*0.3)*0.5*0.5) )*this->timeStepBDF()->timeStep() ;
        ForEachBC( bcDef,cl::paroi_mobile,
                   bilinearForm_PatternCoupled +=
                   /**/ integrate( _range=markedfaces(mesh,PhysicalName),
                                   _expr= alphaRobinFSI*trans(idt(u))*id(v),
                                   _geomap=this->geomap() );
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(mesh,PhysicalName),
                                   _expr= alphaRobinFSI*trans(idv(this->meshVelocity2()))*id(v),
                                   _geomap=this->geomap() );
                   myLinearForm +=
                   /**/ integrate( _range=markedfaces(mesh,PhysicalName),
                                   _expr= -trans(idv(this->normalStressFromStruct()))*id(v),
                                   _geomap=this->geomap() ); );
    }

#endif



    //--------------------------------------------------------------------------------------------------//
    // Todo : BuildNonCstPart ?
    if (M_pdeType == "Oseen" && BuildNonCstPart_ConvectiveTerm)
    {
        auto BetaU = M_bdf_fluid->poly();
        auto betaU = BetaU.template element<0>();
        if (M_isMoveDomain)
        {
#if defined( FEELPP_MODELS_HAS_MESHALE )
            bilinearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= idv(*M_P0Rho)*trans( gradt(u)*( idv(betaU) -idv( this->meshVelocity() )))*id(v),
                           _geomap=this->geomap() );
#endif
        }
        else
        {
            bilinearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= idv(*M_P0Rho)*trans( gradt(u)*idv(betaU) )*id(v),
                           _geomap=this->geomap() );
        }

        if ( this->doStabConvectionEnergy() )
        {
            bilinearForm_PatternCoupled +=
                integrate( _range=elements(mesh),
                           _expr= 0.5*idv(*M_P0Rho)*divt(u)*trans(idv(betaU))*id(v),
                           _geomap=this->geomap() );
        }
    }
    else if (M_pdeType == "Stokes" && BuildNonCstPart_ConvectiveTerm && M_isMoveDomain)
    {
#if defined( FEELPP_MODELS_HAS_MESHALE )
        bilinearForm_PatternDefault +=
            integrate( _range=elements(mesh),
                       _expr= -idv(*M_P0Rho)*trans( gradt(u)*(idv( this->meshVelocity() )))*id(v),
                       _geomap=this->geomap() );
#endif
    }

    //--------------------------------------------------------------------------------------------------//

    this->updateOseenStabilisation(A,F,_BuildCstPart,A_extended,_BuildExtendedPart);

    //--------------------------------------------------------------------------------------------------//

    this->updateOseenWeakBC(A,F,_BuildCstPart);

    //--------------------------------------------------------------------------------------------------//
    //transients terms
    if (!this->isStationary())
    {
        if (BuildNonCstPart_Form2TransientTerm)
        {
            bilinearForm_PatternDefault +=
                integrate( _range=elements(mesh),
                           _expr= idv(*M_P0Rho)*trans(idt(u))*id(v)*M_bdf_fluid->polyDerivCoefficient(0),
                           _geomap=this->geomap() );
        }

        if (BuildNonCstPart_Form1TransientTerm)
        {
            auto Buzz = M_bdf_fluid->polyDeriv();
            auto buzz = Buzz.template element<0>();
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= idv(*M_P0Rho)*trans(idv(buzz))*id(v),
                           _geomap=this->geomap() );
        }
    }


    //--------------------------------------------------------------------------------------------------//

    if (!this->velocityDivIsEqualToZero() && BuildNonCstPart)
    {
        myLinearForm +=
            integrate( _range=elements(mesh),
                       _expr= idv(this->velocityDiv())*id(q),
                       _geomap=this->geomap() );

        auto coeffDiv = (2./3.)*idv(*M_P0Mu); //(eps-2mu/3)
        myLinearForm +=
            integrate( _range=elements(mesh),
                       _expr= val(coeffDiv*gradv(this->velocityDiv()))*id(v),
                       _geomap=this->geomap() );

        if ( this->doStabConvectionEnergy() )
        {
            auto BetaU = M_bdf_fluid->poly();
            auto betaU = BetaU.template element<0>();
            myLinearForm +=
                integrate( _range=elements(mesh),
                           _expr= 0.5*idv(M_P0Rho)*idv(this->velocityDiv())*trans(idv(betaU))*id(v),
                           _geomap=this->geomap() );
        }
    }

    //--------------------------------------------------------------------------------------------------//

    /*if (_doClose)
     {
     A->close();
     F->close();
     }*/

    // strong formulation of the boundaries conditions
    if (this->hasMarkerDirichletBCelimination() && BuildNonCstPart && _doBCStrongDirichlet)
    {
        this->updateBCStrongDirichletLinearPDE(A,F);
    }

    double timeElapsed = thetimer.elapsed();
    if (this->verbose()) Feel::FeelModels::Log(this->prefix()+".FluidMechanics","updateOseen",
                                               "finish in "+(boost::format("%1% s") % timeElapsed).str(),
                                               this->worldComm(),this->verboseAllProc());

#endif // defined(FEELMODELS_FLUID_BUILD_LINEAR_CODE)
} // updateOseen

} // end namespace FeelModels
} // end namespace Feel


