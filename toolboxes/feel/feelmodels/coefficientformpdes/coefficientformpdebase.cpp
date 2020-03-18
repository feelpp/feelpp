/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdebase.hpp>

#include <feel/feelmodels/modelmesh/createmesh.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType>
CoefficientFormPDEBase<ConvexType>::CoefficientFormPDEBase( std::string const& prefix,
                                                            std::string const& keyword,
                                                            worldcomm_ptr_t const& worldComm,
                                                            std::string const& subPrefix,
                                                            ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep )
{
    this->log("CoefficientFormPDE","constructor", "start" );

    std::string nameFileConstructor = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".CoefficientFormPDEConstructor.data";
    std::string nameFileSolve = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".CoefficientFormPDESolve.data";
    std::string nameFilePostProcessing = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".CoefficientFormPDEPostProcessing.data";
    std::string nameFileTimeStepping = this->scalabilityPath() + "/" + this->scalabilityFilename() + ".CoefficientFormPDETimeStepping.data";
    this->addTimerTool("Constructor",nameFileConstructor);
    this->addTimerTool("Solve",nameFileSolve);
    this->addTimerTool("PostProcessing",nameFilePostProcessing);
    this->addTimerTool("TimeStepping",nameFileTimeStepping);

    //-----------------------------------------------------------------------------//
    // option in cfg files
    this->loadParameterFromOptionsVm();
    //-----------------------------------------------------------------------------//
    this->log("CoefficientFormPDE","constructor", "finish");
}


template< typename ConvexType>
void
CoefficientFormPDEBase<ConvexType>::loadParameterFromOptionsVm()
{
#if 0
    // M_fieldVelocityConvectionIsUsed = boption(_name="use_velocity-convection",_prefix=this->prefix()) ||
    //     Environment::vm().count(prefixvm(this->prefix(),"velocity-convection").c_str());
    // M_fieldVelocityConvectionIsIncompressible = boption(_name="velocity-convection_is_incompressible",_prefix=this->prefix());

    M_stabilizationGLS = boption(_name="stabilization-gls",_prefix=this->prefix());
    M_stabilizationGLSType = soption(_name="stabilization-gls.type",_prefix=this->prefix());

    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix());
#endif
}

template< typename ConvexType>
void
CoefficientFormPDEBase<ConvexType>::initMesh()
{
    this->log("CoefficientFormPDE","initMesh", "start");
    this->timerTool("Constructor").start();

    createMeshModel<mesh_type>(*this,M_mesh,this->fileNameMeshPath());
    CHECK( M_mesh ) << "mesh generation fail";

    double tElpased = this->timerTool("Constructor").stop("initMesh");
    this->log("CoefficientFormPDE","initMesh",(boost::format("finish in %1% s")%tElpased).str() );

} // createMesh()

template< typename ConvexType>
void
CoefficientFormPDEBase<ConvexType>::initMaterialProperties()
{
    this->log("CoefficientFormPDE","initMaterialProperties", "start" );
    this->timerTool("Constructor").start();

    if ( !M_materialsProperties )
    {
        auto paramValues = this->modelProperties().parameters().toParameterValues();
        this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->prefix(), this->repository().expr() ) );
        M_materialsProperties->updateForUse( M_mesh, this->modelProperties().materials(), *this );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("CoefficientFormPDE","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

template class CoefficientFormPDEBase< Simplex<2,1> >;
template class CoefficientFormPDEBase< Simplex<3,1> >;

} // namespace Feel
} // namespace FeelModels
