/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#include <feel/feelmodels/coefficientformpdes/coefficientformpdebase.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType>
CoefficientFormPDEBase<ConvexType>::CoefficientFormPDEBase( typename super2_type::infos_type const& infosPDE,
                                                            std::string const& prefix,
                                                            std::string const& keyword,
                                                            worldcomm_ptr_t const& worldComm,
                                                            std::string const& subPrefix,
                                                            ModelBaseRepository const& modelRep )
    :
    super_type( prefix, keyword, worldComm, subPrefix, modelRep, ModelBaseCommandLineOptions( coefficientformpde_options( prefix ) ) ),
    super2_type( infosPDE )
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
    M_applyStabilization = boption(_name="stabilization",_prefix=this->prefix(),_vm=this->clovm());
    M_stabilizationType = soption(_name="stabilization.type",_prefix=this->prefix(),_vm=this->clovm());
    M_stabilizationGLS_applyShockCapturing = boption(_name="stabilization.gls.shock-capturing",_prefix=this->prefix(),_vm=this->clovm());

    M_stabilizationDoAssemblyWithGradDiffusionCoeff = boption(_name="stabilization.do-assembly-with-grad-diffusion-coeff",_prefix=this->prefix(),_vm=this->clovm());

    // time stepping
    M_timeStepping = soption(_name="time-stepping",_prefix=this->prefix(),_vm=this->clovm());
    M_timeStepThetaValue = doption(_name="time-stepping.theta.value",_prefix=this->prefix(),_vm=this->clovm());
}

template< typename ConvexType>
void
CoefficientFormPDEBase<ConvexType>::initMesh()
{
    this->log("CoefficientFormPDE","initMesh", "start");
    this->timerTool("Constructor").start();

    // std::string fileNameMeshPath = prefixvm(this->prefix(),"CoefficientFormPDEMesh.path");
    // createMeshModel<mesh_type>(*this,M_mesh,fileNameMeshPath);
    // CHECK( M_mesh ) << "mesh generation fail";

    if ( this->doRestart() )
        super_type::super_model_meshes_type::setupRestart( this->keyword() );

    //super_type::super_model_meshes_type::setMesh( this->keyword(), M_mesh );
    super_type::super_model_meshes_type::updateForUse<mesh_type>( this->keyword() );

    CHECK( this->mesh() ) << "mesh generation fail";

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
        // auto paramValues = this->modelProperties().parameters().toParameterValues();
        // this->modelProperties().materials().setParameterValues( paramValues );
        M_materialsProperties.reset( new materialsproperties_type( this->shared_from_this_cfpdebase() ) );
        M_materialsProperties->updateForUse( this->modelProperties().materials() );
    }

    double tElpased = this->timerTool("Constructor").stop("initMaterialProperties");
    this->log("CoefficientFormPDE","initMaterialProperties",(boost::format("finish in %1% s")%tElpased).str() );
}

template< typename ConvexType>
void
CoefficientFormPDEBase<ConvexType>::initBasePostProcess()
{
    this->log("Heat","initPostProcess", "start");
    this->timerTool("Constructor").start();

    this->setPostProcessExportsAllFieldsAvailable( { this->unknownName() } );
    this->addPostProcessExportsAllFieldsAvailable( this->materialsProperties()->postProcessExportsAllFieldsAvailable( this->mesh(),this->physicsAvailable() ) );
    this->setPostProcessExportsPidName( "pid" );
    this->setPostProcessSaveAllFieldsAvailable( { this->unknownName() } );
    super_type::initPostProcess();

    if ( !this->postProcessExportsFields().empty() || this->hasPostProcessExportsExpr() )
    {
        std::string geoExportType="static";//change_coords_only, change, static
        M_exporter = exporter( _mesh=this->mesh(),
                               _name="Export",
                               _geo=geoExportType,
                               _path=this->exporterPath() );

        // restart exporter
        if ( M_exporter->doExport() && this->doRestart() && this->restartPath().empty() )
            M_exporter->restart(this->timeInitial());
    }

    // start or restart the export of measures
    if ( !this->isStationary() )
    {
        if ( this->doRestart() )
            this->postProcessMeasuresIO().restart( "time", this->timeInitial() );
        else
            this->postProcessMeasuresIO().setMeasure( "time", this->timeInitial() ); //just for have time in first column
    }

}

template class CoefficientFormPDEBase< Simplex<2,1> >;
template class CoefficientFormPDEBase< Simplex<3,1> >;
template class CoefficientFormPDEBase< Simplex<2,2> >;
template class CoefficientFormPDEBase< Simplex<3,2> >;

} // namespace Feel
} // namespace FeelModels
