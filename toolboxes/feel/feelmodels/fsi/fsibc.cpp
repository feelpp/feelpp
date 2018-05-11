
#include <feel/feelmodels/fsi/fsi.hpp>
#include <feel/feelmodels/modelvf/fluidmecstresstensor.hpp>

namespace Feel
{
namespace FeelModels
{

#if 0
template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateLinearPDEStrongDirichletBC_Fluid( sparse_matrix_ptrtype& A, vector_ptrtype& F ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        this->log("FSI","updateLinearPDEStrongDirichletBC_Fluid", "start" );

        auto mesh = M_fluidModel->mesh();
        auto Xh = M_fluidModel->spaceVelocityPressure();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=M_fluidModel->rowStartInMatrix(),
                                   _colstart=M_fluidModel->colStartInMatrix() );
        auto const& u = M_fluidModel->fieldVelocity();
        bilinearForm +=
            on( _range=markedfaces(mesh, M_fluidModel->markersNameMovingBoundary()),
                _element=u, _rhs=F,
                _expr=idv(M_fluidModel->meshVelocity2()) );

        this->log("FSI","updateLinearPDEStrongDirichletBC_Fluid", "finish" );
    }
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateNewtonInitialGuess_Fluid( vector_ptrtype& U ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        auto mesh = M_fluidModel->mesh();
        auto Xh = M_fluidModel->spaceVelocityPressure();
        auto up = Xh->element( U, M_fluidModel->rowStartInVector() );
        auto u = up.template element<0>();
        u.on(_range=markedfaces(mesh, M_fluidModel->markersNameMovingBoundary()),
             _expr=idv( M_fluidModel->meshVelocity2() ) );
#if 0
        // synchronize velocity dof on interprocess
        auto itFindDofsWithValueImposed = M_dofsWithValueImposed.find("velocity");
        if ( itFindDofsWithValueImposed != M_dofsWithValueImposed.end() )
            sync( u, "=", itFindDofsWithValueImposed->second );
#endif
    }
}

template< class FluidType, class SolidType >
void
FSI<FluidType,SolidType>::updateJacobianStrongDirichletBC_Fluid( sparse_matrix_ptrtype& J,vector_ptrtype& RBis ) const
{
    if ( this->fsiCouplingBoundaryCondition() == "dirichlet-neumann" )
    {
        this->log("FSI","updateJacobianStrongDirichletBC_Fluid", "start" );

        auto mesh = M_fluidModel->mesh();
        auto Xh = M_fluidModel->spaceVelocityPressure();
        auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                                   _pattern=size_type(Pattern::COUPLED),
                                   _rowstart=M_fluidModel->rowStartInMatrix(),
                                   _colstart=M_fluidModel->colStartInMatrix() );
        auto const& u = M_fluidModel->fieldVelocity();
        bilinearForm +=
            on( _range=markedfaces(mesh, M_fluidModel->markersNameMovingBoundary()),
                _element=u, _rhs=RBis,
                _expr= vf::zero<nDim,1>() );

        this->log("FSI","updateJacobianStrongDirichletBC_Fluid", "finish" );
    }
}
#endif

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
    auto Xh = M_fluidModel->spaceVelocityPressure();

    auto const& u = M_fluidModel->fieldVelocity();
    auto const& p = M_fluidModel->fieldPressure();
    auto const& uEval = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0).template element<0>();
    auto const& pEval = (true)? M_fluidModel->fieldPressure() : M_fluidModel->timeStepBDF()->unknown(0).template element<1>();

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmav = Feel::vf::FeelModels::fluidMecNewtonianStressTensor<2*fluid_type::nOrderVelocity>(uEval,pEval,*M_fluidModel->materialProperties(),matName,true);
    auto muExpr = Feel::vf::FeelModels::fluidMecViscosity<2*fluid_type::nOrderVelocity>(uEval,pEval,*M_fluidModel->materialProperties(),matName);
    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto linearForm = form1( _test=Xh, _vector=F,
                             _rowstart=M_fluidModel->rowStartInVector() );
    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() );

    auto rangeFSI = markedfaces(mesh,M_fluidModel->markersNameMovingBoundary());

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
                           _expr= ( gammaRobinFSI*muExpr/hFace() )*inner(idv(M_fluidModel->meshVelocity2()),id(u)),
                           _geomap=this->geomap() );
            //M_fluidModel->meshALE()->revertMovingMesh();
        }
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
             this->fsiCouplingBoundaryCondition() == "nitsche" )
        {
            double alpha = M_couplingNitscheFamily_alpha;
            double gamma0RobinFSI = M_couplingNitscheFamily_gamma0;
            auto mysigma = id(p)*Id+2*alpha*muExpr*sym(grad(u));
            if ( buildCstPart )
            {
                bilinearForm +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idt(p)*id(p),
                               _geomap=this->geomap() );
                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), vf::N() )*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
            if ( buildNonCstPart )
            {
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(pEval)*id(p),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    //M_fluidModel->meshALE()->revertReferenceMesh();
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner(idv(M_fluidModel->meshVelocity2()),vf::N())*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    //M_fluidModel->meshALE()->revertReferenceMesh();
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner(idv(M_fluidModel->meshVelocity2()),mysigma*vf::N()),
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
        if ( true )
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
                auto sigmaSolidN = idv( this->fluidModel()->fieldNormalStressRefMeshPtr() );
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
            //CHECK( false ) << "Not implemented";
            if ( buildCstPart )
            {
                A->close();
                A->addMatrix( this->couplingRNG_coeffForm2(), M_coulingRNG_matrixTimeDerivative );
            }
            if ( buildNonCstPart )
            {
                auto myvec = this->fluidModel()->backend()->newVector(this->fluidModel()->spaceVelocityPressure() );
                *myvec = *this->couplingRNG_evalForm1();
                myvec->scale( -1. );
                F->close();
                F->addVector( myvec, M_coulingRNG_matrixTimeDerivative );

                auto myvec2 = this->fluidModel()->backend()->newVector( this->fluidModel()->fieldNormalStressRefMeshPtr()->functionSpace() );
                *myvec2 = *this->fluidModel()->fieldNormalStressRefMeshPtr();
                F->addVector( myvec2, M_coulingRNG_matrixStress );
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
    auto Xh = M_fluidModel->spaceVelocityPressure();

    auto U = Xh->element(XVec, M_fluidModel->rowStartInVector());
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto const& uPrevious = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0).template element<0>();
    auto const& pPrevious = (true)? M_fluidModel->fieldPressure() : M_fluidModel->timeStepBDF()->unknown(0).template element<1>();

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=J,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_fluidModel->rowStartInMatrix(),
                               _colstart=M_fluidModel->colStartInMatrix() );

    auto rangeFSI = markedfaces(mesh,M_fluidModel->markersNameMovingBoundary());

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto muExpr = Feel::vf::FeelModels::fluidMecViscosity<2*fluid_type::nOrderVelocity>(uPrevious,pPrevious,*M_fluidModel->materialProperties(),matName);


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
                auto mysigma = id(p)*Id+2*alpha*muExpr*sym(grad(u));

                bilinearForm +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idt(p)*id(p),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), vf::N() )*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    bilinearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idt(u), mysigma*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
        }
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        if ( true )
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
            CHECK( false ) << "Not implemented";
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
    auto Xh = M_fluidModel->spaceVelocityPressure();

    auto U = Xh->element(XVec, M_fluidModel->rowStartInVector());
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto linearForm = form1( _test=Xh, _vector=R,
                             _rowstart=M_fluidModel->rowStartInVector() );

    auto const Id = eye<fluid_type::nDim,fluid_type::nDim>();

    auto rangeFSI = markedfaces(mesh,M_fluidModel->markersNameMovingBoundary());

    auto const& uPrevious = (true)? M_fluidModel->fieldVelocity() : M_fluidModel->timeStepBDF()->unknown(0).template element<0>();
    auto const& pPrevious = (true)? M_fluidModel->fieldPressure() : M_fluidModel->timeStepBDF()->unknown(0).template element<1>();

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto sigmavPrevious = Feel::vf::FeelModels::fluidMecNewtonianStressTensor<2*fluid_type::nOrderVelocity>(uPrevious,pPrevious,*M_fluidModel->materialProperties(),matName,true);
    auto muExpr = Feel::vf::FeelModels::fluidMecViscosity<2*fluid_type::nOrderVelocity>(uPrevious,pPrevious,*M_fluidModel->materialProperties(),matName);

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
                           _expr= -(gammaRobinFSI*muExpr/hFace())*inner(idv(M_fluidModel->meshVelocity2()),id(u)),
                           _geomap=this->geomap() );
        }
        if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" ||
             this->fsiCouplingBoundaryCondition() == "nitsche" )
        {
            double alpha = M_couplingNitscheFamily_alpha;
            double gamma0RobinFSI = M_couplingNitscheFamily_gamma0;
            auto mysigma = id(p)*Id+2*alpha*muExpr*sym(grad(u));
            if ( buildNonCstPart && !useJacobianLinearTerms )
            {
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= ( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(p)*id(p),
                               _geomap=this->geomap() );
                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idv(u), vf::N() )*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= -inner( idv(u), mysigma*vf::N() ),
                                   _geomap=this->geomap() );
                }
            }
            if ( buildCstPart )
            {
                linearForm +=
                    integrate( _range=rangeFSI,
                               _expr= -( gamma0RobinFSI*hFace()/(gammaRobinFSI*muExpr) )*idv(pPrevious)*id(p),
                               _geomap=this->geomap() );

                if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-neumann" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= inner(idv(M_fluidModel->meshVelocity2()),vf::N())*id(p),
                                   _geomap=this->geomap() );
                }
                else if ( this->fsiCouplingBoundaryCondition() == "nitsche" )
                {
                    linearForm +=
                        integrate( _range=rangeFSI,
                                   _expr= inner(idv(M_fluidModel->meshVelocity2()),mysigma*vf::N()),
                                   _geomap=this->geomap() );
                }
            }
        }
    }
    else if ( this->fsiCouplingBoundaryCondition() == "robin-neumann-generalized" )
    {
        if ( buildCstPart )
        {
            M_fluidModel->meshALE()->revertReferenceMesh();
            auto sigmaSolidN = idv( this->fluidModel()->fieldNormalStressRefMeshPtr() );
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= -inner( sigmaSolidN,id(u)),
                           _geomap=this->geomap() );
        }

        if ( true )
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
            }
        }
        else
        {
            CHECK(false) << "Not implemented";
        }
        M_fluidModel->meshALE()->revertMovingMesh();
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

    auto mesh = M_solidModel->mesh();
    auto Xh = M_solidModel->functionSpaceDisplacement();
    auto const& u = M_solidModel->fieldDisplacement();

    auto bilinearForm = form2( _test=Xh,_trial=Xh,_matrix=A,
                               _pattern=size_type(Pattern::COUPLED),
                               _rowstart=M_solidModel->rowStartInMatrix(),
                               _colstart=M_solidModel->colStartInMatrix() );
    auto linearForm = form1( _test=Xh, _vector=F,
                             _rowstart=M_solidModel->rowStartInVector() );

    auto rangeFSI = markedfaces(mesh,M_solidModel->markerNameFSI());

    // neumann boundary condition with normal stress (fsi boundary condition)
    if ( buildNonCstPart)
    {
        form1( _test=Xh, _vector=F) +=
            integrate( _range=rangeFSI,
                       _expr= -inner( idv(M_solidModel->fieldNormalStressFromFluidPtr()),id(u) ),
                       _geomap=this->geomap() );
    }

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;

        CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
        std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
        auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
        CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
        double muFluid = dynamicViscosity.newtonian().value();

#if 0
            MeshMover<mesh_type> mymesh_mover;
            mesh_ptrtype mymesh = this->mesh();
            mymesh_mover.apply( mymesh, this->timeStepNewmark()->previousUnknown() );
#endif

        if ( buildCstPart )
        {
            bilinearForm +=
                integrate( _range=rangeFSI,
                           _expr= gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(u))/hFace(),
                           _geomap=this->geomap() );
        }
        if ( buildNonCstPart )
        {
            auto robinFSIRhs = idv(M_solidModel->timeStepNewmark()->polyFirstDeriv() ) + idv(M_solidModel->fieldVelocityInterfaceFromFluid());
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= gammaRobinFSI*muFluid*inner( robinFSIRhs, id(u))/hFace(),
                           _geomap=this->geomap() );
        }
#if 0
            auto dispInv = this->fieldDisplacement().functionSpace()->element(-idv(this->timeStepNewmark()->previousUnknown()));
            mymesh_mover.apply( mymesh, dispInv );
#endif
    }
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



    double gammaRobinFSI = M_couplingNitscheFamily_gamma;

    CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
    std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
    auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
    CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
    double muFluid = dynamicViscosity.newtonian().value();

    if ( buildCstPart )
    {
#if 0
        MeshMover<typename solid_type::mesh_type> mymesh_mover;
        mymesh_mover.apply( mesh, M_solidModel->timeStepNewmark()->previousUnknown() );
#endif
        bilinearForm +=
            integrate( _range=markedfaces(mesh,M_solidModel->markerNameFSI()),
                       _expr= gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idt(u),id(u))/hFace(),
                       _geomap=this->geomap() );

#if 0
        auto dispInv = M_solidModel->fieldDisplacement().functionSpace()->element(-idv(M_solidModel->timeStepNewmark()->previousUnknown()));
        mymesh_mover.apply( mesh, dispInv );
#endif
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


    auto mesh = M_solidModel->mesh();
    auto Xh = M_solidModel->functionSpaceDisplacement();
    auto linearForm = form1( _test=Xh, _vector=R,
                             _rowstart=M_solidModel->rowStartInVector() );

    auto u = Xh->element(XVec, M_solidModel->rowStartInVector());

    auto rangeFSI = markedfaces(mesh,M_solidModel->markerNameFSI());

    // neumann boundary condition with normal stress (fsi boundary condition)
    if ( buildCstPart )
    {
        linearForm +=
            integrate( _range=rangeFSI,
                       _expr= trans(idv(M_solidModel->fieldNormalStressFromFluidPtr()))*id(u),
                       _geomap=this->geomap() );
    }

    if ( this->fsiCouplingBoundaryCondition() == "robin-robin" || this->fsiCouplingBoundaryCondition() == "robin-robin-genuine" ||
         this->fsiCouplingBoundaryCondition() == "nitsche" )
    {
        double gammaRobinFSI = M_couplingNitscheFamily_gamma;

        CHECK( M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().size() == 1 ) << "support only one";
        std::string matName = M_fluidModel->materialProperties()->rangeMeshElementsByMaterial().begin()->first;
        auto const& dynamicViscosity = M_fluidModel->materialProperties()->dynamicViscosity(matName);
        CHECK( dynamicViscosity.isNewtonianLaw() && dynamicViscosity.newtonian().isConstant() ) << "TODO";
        double muFluid = dynamicViscosity.newtonian().value();

#if 0
        MeshMover<typename solid_type::mesh_type> mymesh_mover;
        mymesh_mover.apply( mesh, M_solidModel->timeStepNewmark()->previousUnknown() );
#endif
        if ( buildNonCstPart && !useJacobianLinearTerms )
        {
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= gammaRobinFSI*muFluid*M_solidModel->timeStepNewmark()->polyFirstDerivCoefficient()*inner(idv(u),id(u))/hFace(),
                           _geomap=this->geomap() );
        }
        if ( buildCstPart )
        {
            auto robinFSIRhs = idv(M_solidModel->timeStepNewmark()->polyFirstDeriv() ) + idv(M_solidModel->fieldVelocityInterfaceFromFluid());
            linearForm +=
                integrate( _range=rangeFSI,
                           _expr= -gammaRobinFSI*muFluid*inner( robinFSIRhs,id(u) )/hFace(),
                           _geomap=this->geomap() );
        }
#if 0
        auto dispInv = M_solidModel->fieldDisplacement().functionSpace()->element(-idv(M_solidModel->timeStepNewmark()->previousUnknown()));
        mymesh_mover.apply( mesh, dispInv );
#endif
    } // robin-robin fsi
    this->log("FSI","updateResidual_Solid", "finish"+sc );

}

} // namespace FeelModels
} // namespace Feel
