
#include <feel/feelmodels/fsi/fsi.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{


template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDEDofElimination_Fluid( DataUpdateLinear & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        this->log("FSI","updateLinearPDEDofElimination_Fluid", "start" );

        sparse_matrix_ptrtype& A = data.matrix();
        vector_ptrtype& F = data.rhs();

        auto mesh = M_fluidModel->mesh();
        auto XhV = M_fluidModel->functionSpaceVelocity();
        auto bilinearForm = form2( _test=XhV,_trial=XhV,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=M_fluidModel->rowStartInMatrix(),
                                   _colstart=M_fluidModel->colStartInMatrix() );
        auto const& u = M_fluidModel->fieldVelocity();
        bilinearForm +=
            on( _range=M_rangeFSI_fluid,
                _element=u, _rhs=F,
                _expr=idv(this/*M_fluidModel*/->meshVelocity2()) );

        this->log("FSI","updateLinearPDEDofElimination_Fluid", "finish" );
    }
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateNewtonInitialGuess_Fluid( DataNewtonInitialGuess & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        this->log("FSI","updateNewtonInitialGuess_Fluid", "start" );

        vector_ptrtype& U = data.initialGuess();
        auto mesh = M_fluidModel->mesh();
        auto XhV = M_fluidModel->functionSpaceVelocity();
        auto u = XhV->element( U, M_fluidModel->rowStartInVector() );
        u.on(_range=M_rangeFSI_fluid,
             _expr=idv( this/*M_fluidModel*/->meshVelocity2() ) );
        // update info for synchronization
        M_fluidModel->updateDofEliminationIds( "velocity", this->dofEliminationIds( "fluid.velocity" ), data );

        this->log("FSI","updateNewtonInitialGuess_Fluid", "finish" );
    }
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateJacobianDofElimination_Fluid( DataUpdateJacobian & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "dirichlet-neumann" )
        return;

    this->log("FSI","updateJacobianDofElimination_Fluid", "start" );

    M_fluidModel->updateDofEliminationIds( "velocity", this->dofEliminationIds( "fluid.velocity" ), data );

    this->log("FSI","updateJacobianDofElimination_Fluid", "finish" );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateResidualDofElimination_Fluid( DataUpdateResidual & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "dirichlet-neumann" )
        return;

    this->log("FSI","updateResidualDofElimination_Fluid", "start" );

    M_fluidModel->updateDofEliminationIds( "velocity", this->dofEliminationIds( "fluid.velocity" ), data );

    this->log("FSI","updateResidualDofElimination_Fluid", "finish" );
}


template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDE_Fluid( DataUpdateLinear & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "robin-neumann" && this->fsiCouplingBoundaryCondition() != "robin-neumann-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-robin" && this->fsiCouplingBoundaryCondition() != "robin-robin-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-neumann-generalized" && this->fsiCouplingBoundaryCondition() != "nitsche" )
        return;

    const vector_ptrtype& vecCurrentPicardSolution = data.currentSolution();
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool buildNonCstPart_robinFSI = buildNonCstPart;
    if ( this->useFSISemiImplicitScheme() )
        buildNonCstPart_robinFSI = buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateLinearPDE_Fluid", "start"+sc );
    M_fluidModel->timerTool("Solve").start();

    auto mesh = M_fluidModel->mesh();
    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto XhP = M_fluidModel->functionSpacePressure();

    auto const& u = M_fluidModel->fieldVelocity();
    auto const& p = M_fluidModel->fieldPressure();
    auto const& uEval = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0);
    auto const& pEval = M_fluidModel->fieldPressure();

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmav = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(uEval),idv(pEval),*M_fluidModel->materialProperties(),matName,true);
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(uEval),*M_fluidModel->materialProperties(),matName);
    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto linearForm = form1( _test=XhV, _vector=F,
                             _rowstart=M_fluidModel->rowStartInVector() );
    auto bilinearForm = form2( _test=XhV,_trial=XhV,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() );

    auto rangeFSI = M_rangeFSI_fluid;

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
        //---------------------------------------------------------------------------//
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;
        if ( buildCstPart )
        {
            bilinearForm +=
                integrate( _range=rangeFSI,
                           _expr= ( gammaRobinFSI*muExpr/hFace() )*inner(idt(u),id(u)),
                           _geomap=this->geomap() );
        }
        if ( buildNonCstPart )
        {
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= inner( sigmav*N(),id(u)),
                           _geomap=this->geomap() );

            //M_fluidModel->meshALE()->revertReferenceMesh();
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= ( gammaRobinFSI*muExpr/hFace() )*inner(idv(this/*M_fluidModel*/->meshVelocity2()),id(u)),
                           _geomap=this->geomap() );
            //M_fluidModel->meshALE()->revertMovingMesh();
        }
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
             this->fsiCouplingBoundaryCondition() == "nitsche" )
        {
            double alpha = M_couplingNitscheFamily_alpha;
            double gamma0RobinFSI = M_couplingNitscheFamily_gamma0;
            auto mysigma_p = id(p)*Id;
            auto mysigma_u = 2*alpha*muExpr*sym(grad(u));
            if ( buildCstPart )
            {
                auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=A,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=M_fluidModel->rowStartInMatrix()+1,
                                             _colstart=M_fluidModel->colStartInMatrix()+1 );
                auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=A,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=M_fluidModel->rowStartInMatrix()+1,
                                             _colstart=M_fluidModel->colStartInMatrix() );
                bilinearFormPP +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idt(p)*id(p),
                               _geomap=this->geomap() );
                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    bilinearFormPV +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), vf::N() )*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma_u*vf::N() ),
                                   _geomap=this->geomap() );
                    bilinearFormPV +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma_p*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
            if ( buildNonCstPart )
            {
                auto linearFormP = form1( _test=XhP, _vector=F,
                                          _rowstart=M_fluidModel->rowStartInVector()+1 );
                linearFormP +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(pEval)*id(p),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    //M_fluidModel->meshALE()->revertReferenceMesh();
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner(idv(this/*M_fluidModel*/->meshVelocity2()),vf::N())*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    //M_fluidModel->meshALE()->revertReferenceMesh();
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner(idv(this/*M_fluidModel*/->meshVelocity2()),mysigma_u*vf::N()),
                                   _geomap=this->geomap() );
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner(idv(this/*M_fluidModel*/->meshVelocity2()),mysigma_p*vf::N()),
                                   _geomap=this->geomap() );
                }
                //M_fluidModel->meshALE()->revertMovingMesh();
            }
        }
    }
    //---------------------------------------------------------------------------//
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
        //---------------------------------------------------------------------------//
    {
        if ( !M_coulingRNG_usePrecomputeBC )
        {
            auto myB = this->couplingRNG_operatorExpr( mpl::int_<fluid_type::nDim>() );
            if ( buildCstPart )
            {
                M_fluidModel->meshALE()->revertReferenceMesh();
                bilinearForm +=
                    integrate( _range=rangeFSI,
                               //_expr=this->couplingRNG_coeffForm2()*inner(idt(u),id(u)),
                               _expr=this->couplingRNG_coeffForm2()*inner(myB*idt(u),id(u)),
                               _geomap=this->geomap() );
            }
            if ( buildNonCstPart )
            {
                M_fluidModel->meshALE()->revertReferenceMesh();
                linearForm +=
                    integrate( _range=rangeFSI,
                               //_expr= -inner(idv(this->couplingRNG_evalForm1()),id(u)),
                               _expr= -inner(myB*idv(this->couplingRNG_evalForm1()),id(u)),
                               _geomap=this->geomap() );
                auto sigmaSolidN = idv( M_fieldNormalStressRefMesh_fluid );
                //auto sigmaSolidN = -idv( this->fluidModel()->normalStressFromStruct() );
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= inner( sigmaSolidN,id(u)),
                               _geomap=this->geomap() );
            }
            M_fluidModel->meshALE()->revertMovingMesh();
        }
        else
        {
            if ( buildCstPart )
            {
                A->close();
                A->addMatrix( this->couplingRNG_coeffForm2(), M_coulingRNG_matrixTimeDerivative, Feel::SUBSET_NONZERO_PATTERN );
            }
            if ( buildNonCstPart )
            {
                auto uWrap = XhV->element( M_coulingRNG_vectorTimeDerivative, 0 );
                uWrap.zero();
                uWrap.add( -1., *this->couplingRNG_evalForm1() );
                F->close();
                F->addVector( M_coulingRNG_vectorTimeDerivative, M_coulingRNG_matrixTimeDerivative );

                *M_coulingRNG_vectorStress = *M_fieldNormalStressRefMesh_fluid;
                F->addVector( M_coulingRNG_vectorStress, M_coulingRNG_matrixStress );
            }
        }
    }

    double timeElapsed = M_fluidModel->timerTool("Solve").stop();
    this->log("FSI","updateLinearPDE_Fluid","finish in "+(boost::format("%1% s") %timeElapsed ).str() );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateJacobian_Fluid( DataUpdateJacobian & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "robin-neumann" && this->fsiCouplingBoundaryCondition() != "robin-neumann-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-robin" && this->fsiCouplingBoundaryCondition() != "robin-robin-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-neumann-generalized" && this->fsiCouplingBoundaryCondition() != "nitsche" )
        return;

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateJacobian_Fluid", "start"+sc );

    auto mesh = M_fluidModel->mesh();
    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto XhP = M_fluidModel->functionSpacePressure();

    auto u = XhV->element(XVec, M_fluidModel->rowStartInVector());
    auto const& q =  M_fluidModel->fieldPressure();
    auto const& uPrevious = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0);
    auto const& pPrevious = M_fluidModel->fieldPressure();

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto bilinearForm = form2( _test=XhV,_trial=XhV,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() );

    auto rangeFSI = M_rangeFSI_fluid;

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(uPrevious),*M_fluidModel->materialProperties(),matName);


    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;
        if ( buildCstPart )
        {
            bilinearForm +=
                integrate( _range=rangeFSI,
                           _expr= ( gammaRobinFSI*muExpr/hFace() )*inner(idt(u),id(u)),
                           _geomap=this->geomap() );

            if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
                 this->fsiCouplingBoundaryCondition() == "nitsche" )
            {
                double alpha = M_couplingNitscheFamily_alpha;
                double gamma0RobinFSI = M_couplingNitscheFamily_gamma0;
                auto mysigma_u = 2*alpha*muExpr*sym(grad(u));
                auto mysigma_p = id(q)*Id;
                auto bilinearFormPV = form2( _test=XhP,_trial=XhV,_matrix=J,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=M_fluidModel->rowStartInMatrix()+1,
                                             _colstart=M_fluidModel->colStartInMatrix() );
                auto bilinearFormPP = form2( _test=XhP,_trial=XhP,_matrix=J,
                                             _pattern=size_type(Pattern::COUPLED),
                                             _rowstart=M_fluidModel->rowStartInMatrix()+1,
                                             _colstart=M_fluidModel->colStartInMatrix()+1 );

                bilinearFormPP +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idt(q)*id(q),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    bilinearFormPV +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), vf::N() )*id(q),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma_u*vf::N() ),
                                   _geomap=this->geomap() );
                    bilinearFormPV +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma_p*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
        }
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        if ( !M_coulingRNG_usePrecomputeBC )
        {
            if ( buildCstPart )
            {
                auto myB = this->couplingRNG_operatorExpr( mpl::int_<fluid_type::nDim>() );
                M_fluidModel->meshALE()->revertReferenceMesh();
                bilinearForm +=
                    integrate( _range=rangeFSI,
                               _expr=this->couplingRNG_coeffForm2()*inner(myB*idt(u),id(u)),
                               _geomap=this->geomap() );
                M_fluidModel->meshALE()->revertMovingMesh();
            }
        }
        else
        {
            if ( buildCstPart )
            {
                J->close();
                J->addMatrix( this->couplingRNG_coeffForm2(), M_coulingRNG_matrixTimeDerivative, Feel::SUBSET_NONZERO_PATTERN );
            }
        }
    }

    this->log("FSI","updateJacobian_Fluid", "finish"+sc );

}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateResidual_Fluid( DataUpdateResidual & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "robin-neumann" && this->fsiCouplingBoundaryCondition() != "robin-neumann-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-robin" && this->fsiCouplingBoundaryCondition() != "robin-robin-genuine" &&
         this->fsiCouplingBoundaryCondition() != "robin-neumann-generalized" && this->fsiCouplingBoundaryCondition() != "nitsche" )
        return;

    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateResidual_Fluid", "start"+sc );

    auto mesh = M_fluidModel->mesh();
    auto XhV = M_fluidModel->functionSpaceVelocity();
    auto XhP = M_fluidModel->functionSpacePressure();
    auto u = XhV->element(XVec, M_fluidModel->rowStartInVector());
    auto p = XhP->element(XVec, M_fluidModel->rowStartInVector() +1);
    auto const& q = M_fluidModel->fieldPressure();

    auto linearForm = form1( _test=XhV, _vector=R,
                             _rowstart=M_fluidModel->rowStartInVector() );
    auto linearFormP = form1( _test=XhP, _vector=R,
                             _rowstart=M_fluidModel->rowStartInVector()+1 );

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto rangeFSI = M_rangeFSI_fluid;

    auto const& uPrevious = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0);
    auto const& pPrevious = M_fluidModel->fieldPressure();

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmavPrevious = Feel::FeelModels::fluidMecNewtonianStressTensor(gradv(uPrevious),idv(pPrevious),*M_fluidModel->materialProperties(),matName,true);
    auto muExpr = Feel::FeelModels::fluidMecViscosity(gradv(uPrevious),*M_fluidModel->materialProperties(),matName);

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "robin-neumann" || this->fsiCouplingBoundaryCondition() == "robin-neumann-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;
        if ( buildNonCstPart && !useJacobianLinearTerms )
        {
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= (gammaRobinFSI*muExpr/hFace())*inner(idv(u),id(u)),
                           _geomap=this->geomap() );
        }

        if ( buildCstPart )
        {
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= -inner( sigmavPrevious*N(),id(u)),
                           _geomap=this->geomap() );
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= -(gammaRobinFSI*muExpr/hFace())*inner(idv(this/*M_fluidModel*/->meshVelocity2()),id(u)),
                           _geomap=this->geomap() );
        }
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
             this->fsiCouplingBoundaryCondition() == "nitsche" )
        {
            double alpha = M_couplingNitscheFamily_alpha;
            double gamma0RobinFSI = M_couplingNitscheFamily_gamma0;
            auto mysigma_u = 2*alpha*muExpr*sym(grad(u));
            auto mysigma_p = id(q)*Id;
            if ( buildNonCstPart && !useJacobianLinearTerms )
            {
                auto p = XhP->element(XVec, M_fluidModel->rowStartInVector()+1);
                linearFormP +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(p)*id(p),
                               _geomap=this->geomap() );
                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idv(u), vf::N() )*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idv(u), mysigma_u*vf::N() ),
                                   _geomap=this->geomap() );
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idv(u), mysigma_p*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
            if ( buildCstPart )
            {
                linearFormP +=
                    integrate( _range=rangeFSI,
                               _expr= -( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(pPrevious)*id(p),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= inner(idv(this/*M_fluidModel*/->meshVelocity2()),vf::N())*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= inner(idv(this/*M_fluidModel*/->meshVelocity2()),mysigma_u*vf::N()),
                                   _geomap=this->geomap() );
                    linearFormP +=
                        integrate( _range=rangeFSI,
                                   _expr= inner(idv(this/*M_fluidModel*/->meshVelocity2()),mysigma_p*vf::N()),
                                   _geomap=this->geomap() );
                }
            }
        }
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        if ( !M_coulingRNG_usePrecomputeBC )
        {
            auto myB = this->couplingRNG_operatorExpr( mpl::int_<fluid_type::nDim>() );
            if ( buildNonCstPart && !useJacobianLinearTerms )
            {
                M_fluidModel->meshALE()->revertReferenceMesh();
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr=this->couplingRNG_coeffForm2()*inner(myB*idv(u),id(u)),
                               _geomap=this->geomap() );
            }
            if ( buildCstPart )
            {
                M_fluidModel->meshALE()->revertReferenceMesh();
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= inner(myB*idv(this->couplingRNG_evalForm1()),id(u)),
                               _geomap=this->geomap() );
                auto sigmaSolidN = idv( M_fieldNormalStressRefMesh_fluid );
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= -inner( sigmaSolidN,id(u)),
                               _geomap=this->geomap() );
            }
            M_fluidModel->meshALE()->revertMovingMesh();
        }
        else
        {
            if ( buildNonCstPart && !useJacobianLinearTerms )
            {
                CHECK(false) << "Not implemented";
            }
            if ( buildCstPart )
            {
                auto uWrap = XhV->element( M_coulingRNG_vectorTimeDerivative, 0 );
                uWrap.zero();
                uWrap.add( 1., *this->couplingRNG_evalForm1() );
                R->close();
                R->addVector( M_coulingRNG_vectorTimeDerivative, M_coulingRNG_matrixTimeDerivative );

                *M_coulingRNG_vectorStress = *M_fieldNormalStressRefMesh_fluid;
                M_coulingRNG_vectorStress->scale(-1.);
                R->addVector( M_coulingRNG_vectorStress, M_coulingRNG_matrixStress );
            }
        }
    }

    this->log("FSI","updateResidual_Fluid", "finish"+sc );
}


template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDE_Solid( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

     std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateLinearPDE_Solid", "start"+sc );

    auto mesh = M_solidModel->mesh();
    auto Xh = M_solidModel->functionSpaceDisplacement();
    auto const& u = M_solidModel->fieldDisplacement();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_solidModel->rowStartInMatrix(),
                               _colstart=M_solidModel->colStartInMatrix() );
    auto linearForm = form1( _test=Xh, _vector=F,
                             _rowstart=M_solidModel->rowStartInVector() );

    auto rangeFSI = M_rangeFSI_solid;

    double timeSteppingScaling = 1.;
    if ( !this->solidModel()->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->solidModel()->prefix(),"time-stepping.scaling") );

    // neumann boundary condition with normal stress (fsi boundary condition)
    if ( buildNonCstPart)
    {
        linearForm +=
            integrate( _range=rangeFSI,
                       _expr= -timeSteppingScaling*inner( idv(this->fieldNormalStressFromFluidPtr_solid()),id(u) ),
                       _geomap=this->geomap() );
    }

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;

        CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
        std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
#if 0
        auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
        CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
        double muFluid = dynamicViscosity.newtonian().value();
#endif

        auto gradVelocityExpr = gradVelocityExpr_fluid2solid( hana::int_<fluid_type::nDim>() );
        auto muFluid = Feel::FeelModels::fluidMecViscosity( gradVelocityExpr, *M_fluidModel->materialProperties(), matName, invalid_uint16_type_value, true );

#if 0
        MeshMover<mesh_type> mymesh_mover;
        mesh_ptrtype mymesh = this->mesh();
        mymesh_mover.apply( mymesh, this->timeStepNewmark()->previousUnknown() );
#endif

        if ( this->solidModel()->timeStepping() == "Newmark" )
        {
            if ( buildCstPart )
            {
                bilinearForm +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(u))/hFace(),
                               _geomap=this->geomap() );
            }
            if ( buildNonCstPart )
            {
                auto robinFSIRhs = idv(M_solidModel->timeStepNewmark()->polyFirstDeriv() ) + idv(this->fieldVelocityInterfaceFromFluid_solid());
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*inner( robinFSIRhs, id(u))/hFace(),
                               _geomap=this->geomap() );
            }
        }
        else
        {
            size_type startBlockIndexVelocity = this->solidModel()->startSubBlockSpaceIndex("velocity");
            if ( buildCstPart )
            {
                form2( _test=Xh,_trial=Xh,_matrix=A,
                       _pattern=size_type(Pattern::COUPLED),
                       _rowstart=this->solidModel()->rowStartInMatrix(),
                       _colstart=this->solidModel()->colStartInMatrix() + startBlockIndexVelocity ) +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*inner(idt(u),id(u))/hFace(),
                               _geomap=this->geomap() );
            }
            if ( buildNonCstPart )
            {
                auto robinFSIRhs = idv(this->fieldVelocityInterfaceFromFluid_solid());
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*inner( robinFSIRhs, id(u))/hFace(),
                               _geomap=this->geomap() );
            }
        }
#if 0
        auto dispInv = this->fieldDisplacement().functionSpace()->element(-idv(this->timeStepNewmark()->previousUnknown()));
        mymesh_mover.apply( mymesh, dispInv );
#endif
    }

     this->log("FSI","updateLinearPDE_Solid", "finish"+sc );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateJacobian_Solid( DataUpdateJacobian & data ) const
{
    if ( this->fsiCouplingBoundaryCondition() != "robin-robin" && this->fsiCouplingBoundaryCondition() != "robin-robin-genuine" &&
         this->fsiCouplingBoundaryCondition() != "nitsche" )
        return;

    const vector_ptrtype& XVec = data.currentSolution();
    sparse_matrix_ptrtype& J = data.jacobian();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateJacobian_Solid", "start"+sc );

    auto mesh = M_solidModel->mesh();
    auto Xh = M_solidModel->functionSpaceDisplacement();

    auto u = Xh->element(XVec, M_solidModel->rowStartInVector());

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_solidModel->rowStartInMatrix(),
                               _colstart=M_solidModel->colStartInMatrix() );


    double timeSteppingScaling = 1.;
    if ( !this->solidModel()->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->solidModel()->prefix(),"time-stepping.scaling") );


    double gammaRobinFSI = M_couplingNitscheFamily_gamma;

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
#if 0
    auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
    CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
    double muFluid = dynamicViscosity.newtonian().value();
#endif

    if ( buildCstPart )
    {
        auto gradVelocityExpr = gradVelocityExpr_fluid2solid( hana::int_<fluid_type::nDim>() );
        auto muFluid = Feel::FeelModels::fluidMecViscosity( gradVelocityExpr, *M_fluidModel->materialProperties(),matName, invalid_uint16_type_value, true );
        auto rangeFSI = M_rangeFSI_solid;

        if ( this->solidModel()->timeStepping() == "Newmark" )
        {
            bilinearForm +=
                integrate( _range=rangeFSI,
                           _expr= timeSteppingScaling*gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(u))/hFace(),
                           _geomap=this->geomap() );
        }
        else
        {
            size_type startBlockIndexVelocity = this->solidModel()->startSubBlockSpaceIndex("velocity");
            form2( _test=Xh,_trial=Xh,_matrix=J,
                   _pattern=size_type(Pattern::DEFAULT),
                   _rowstart=this->solidModel()->rowStartInMatrix(),
                   _colstart=this->solidModel()->colStartInMatrix() + startBlockIndexVelocity ) +=
                integrate( _range=rangeFSI,
                           _expr= timeSteppingScaling*gammaRobinFSI*muFluid*inner(idt(u),id(u))/hFace(),
                           _geomap=this->geomap() );
        }
    }

    this->log("FSI","updateJacobian_Solid", "finish"+sc );
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateResidual_Solid( DataUpdateResidual & data ) const
{
    const vector_ptrtype& XVec = data.currentSolution();
    vector_ptrtype& R = data.residual();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    bool useJacobianLinearTerms = data.useJacobianLinearTerms();

    bool buildCstPart_BoundaryParoiMobile = buildCstPart;
    if ( this->useFSISemiImplicitScheme() )
        buildCstPart_BoundaryParoiMobile = buildNonCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateResidual_Solid", "start"+sc );

    double timeSteppingScaling = 1.;
    if ( !this->solidModel()->isStationary() )
        timeSteppingScaling = data.doubleInfo( prefixvm(this->solidModel()->prefix(),"time-stepping.scaling") );

    auto mesh = M_solidModel->mesh();
    auto Xh = M_solidModel->functionSpaceDisplacement();
    auto linearForm = form1( _test=Xh, _vector=R,
                             _rowstart=M_solidModel->rowStartInVector() );

    auto u = Xh->element(XVec, M_solidModel->rowStartInVector());

    auto rangeFSI = M_rangeFSI_solid;

    // neumann boundary condition with normal stress (fsi boundary condition)
    if ( buildCstPart )
    {
        linearForm +=
            integrate( _range=rangeFSI,
                       _expr= timeSteppingScaling*inner(idv(this->fieldNormalStressFromFluidPtr_solid()),id(u)),
                       _geomap=this->geomap() );
    }

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;

        CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
        std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
#if 0
        auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
        CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
        double muFluid = dynamicViscosity.newtonian().value();
#endif
        auto gradVelocityExpr = gradVelocityExpr_fluid2solid( hana::int_<fluid_type::nDim>() );
        auto muFluid = Feel::FeelModels::fluidMecViscosity( gradVelocityExpr, *M_fluidModel->materialProperties(), matName, invalid_uint16_type_value, true );

        if ( buildNonCstPart && !useJacobianLinearTerms )
        {
            if ( this->solidModel()->timeStepping() == "Newmark" )
            {
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idv(u),id(u))/hFace(),
                               _geomap=this->geomap() );
            }
            else
            {
                size_type startBlockIndexVelocity = this->solidModel()->startSubBlockSpaceIndex("velocity");
                auto const curVel = Xh->element(XVec, this->solidModel()->rowStartInVector()+startBlockIndexVelocity);
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= timeSteppingScaling*gammaRobinFSI*muFluid*inner(idv(curVel),id(u))/hFace(),
                               _geomap=this->geomap() );
            }
        }
        if ( buildCstPart )
        {
            if ( this->solidModel()->timeStepping() == "Newmark" )
            {
                auto robinFSIRhs = idv(M_solidModel->timeStepNewmark()->polyFirstDeriv() ) + idv(this->fieldVelocityInterfaceFromFluid_solid());
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= -timeSteppingScaling*gammaRobinFSI*muFluid*inner( robinFSIRhs,id(u) )/hFace(),
                               _geomap=this->geomap() );
            }
            else
            {
                auto robinFSIRhs = idv(this->fieldVelocityInterfaceFromFluid_solid());
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= -timeSteppingScaling*gammaRobinFSI*muFluid*inner( robinFSIRhs,id(u) )/hFace(),
                               _geomap=this->geomap() );

            }
        }

    } // robin-robin fsi
    this->log("FSI","updateResidual_Solid", "finish"+sc );

}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDE_Solid1dReduced( DataUpdateLinear & data ) const
{
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    std::string sc=(buildCstPart)?" (cst)":" (non cst)";
    this->log("FSI","updateLinearPDE_Solid1dReduced", "start"+sc );

    if ( buildNonCstPart )
    {
        auto mesh = M_solidModel->mesh();
        auto Xh = M_solidModel->solid1dReduced()->functionSpace1dReduced();
        auto const& v = M_solidModel->solid1dReduced()->fieldDisplacementScal1dReduced();
        auto linearForm = form1( _test=Xh, _vector=F,
                                 _rowstart=M_solidModel->rowStartInVector() );
        auto rangeMeshElements1dReduced = elements(M_solidModel->solid1dReduced()->mesh());

        linearForm +=
            integrate( _range=rangeMeshElements1dReduced,
                       _expr=idv(*M_fieldNormalStressFromFluidScalar_solid1dReduced)*id(v) );
    }

    this->log("FSI","updateLinearPDE_Solid1dReduced", "finish"+sc );
}


} // namespace FeelModels
} // namespace Feel
