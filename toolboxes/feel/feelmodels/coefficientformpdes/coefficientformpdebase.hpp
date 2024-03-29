/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP 1

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/tsbase.hpp>

#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
//#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/modelcore/modelgenericpde.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{
template< typename ConvexType>
class CoefficientFormPDEBase : public ModelNumerical,
                               public ModelGenericPDE<ConvexType::nDim>
    //public ModelPhysics<ConvexType::nRealDim>
{
protected :
    using super_type = ModelNumerical;
    //using super_physics_type = ModelPhysics<ConvexType::nRealDim>;
    using super_physics_type = ModelGenericPDE<ConvexType::nDim>;
    using super2_type = super_physics_type;
    using Coefficient = typename super_physics_type::Coefficient;
public :
    using self_type = CoefficientFormPDEBase<ConvexType>;
    using size_type = typename super_type::size_type;

    // mesh
    typedef ConvexType convex_type;
    static inline const uint16_type nDim = convex_type::nDim;
    static inline const uint16_type nOrderGeo = convex_type::nOrder;
    static inline const uint16_type nRealDim = convex_type::nDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // stabilization
    typedef StabilizationGLSParameterBase<mesh_type> stab_gls_parameter_type;
    typedef std::shared_ptr<stab_gls_parameter_type> stab_gls_parameter_ptrtype;

    // exporter
    typedef Exporter<mesh_type,nOrderGeo> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

public :
    CoefficientFormPDEBase( typename super_physics_type::infos_ptrtype const& infosPDE,
                            std::string const& prefix,
                            std::string const& keyword,
                            worldcomm_ptr_t const& worldComm,
                            std::string const& subPrefix,
                            ModelBaseRepository const& modelRep );

    //! return current shared_ptr of type CoefficientFormPDEBase
    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }

    //! return true is the unknown is scalar
    virtual bool unknownIsScalar() const = 0;

    //! return tabulate informations
    virtual tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const = 0;

    //___________________________________________________________________________________//
    // mesh
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    Range<mesh_type,MESH_ELEMENTS> const& rangeMeshElements() const { return M_rangeMeshElements; }
    void setMesh( mesh_ptrtype const& mesh ) { super_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    //___________________________________________________________________________________//
    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    //___________________________________________________________________________________//
    // stabilization
    bool applyStabilization() const { return M_applyStabilization; }
    std::string const& stabilizationType() const { return M_stabilizationType; }
    stab_gls_parameter_ptrtype const& stabilizationGLSParameter() const { return M_stabilizationGLSParameter; }
    bool stabilizationGLS_applyShockCapturing() const { return M_stabilizationGLS_applyShockCapturing; }
    bool stabilizationDoAssemblyWithGradDiffusionCoeff() const { return M_stabilizationDoAssemblyWithGradDiffusionCoeff; }

    //___________________________________________________________________________________//
    // time discretisation
    virtual std::shared_ptr<TSBase> timeStepBase() const = 0;
    virtual void startTimeStep() = 0;
    virtual void updateTimeStep() = 0;
    std::string const& timeStepping() const { return M_timeStepping; }

    //___________________________________________________________________________________//

    virtual void setParameterValues( std::map<std::string,double> const& paramValues ) = 0;

    template <typename SymbolsExprType>
    bool hasSymbolDependencyInCoefficients( std::set<std::string> const& symbs, SymbolsExprType const& se ) const;

protected :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initMaterialProperties();
    void initBasePostProcess();

protected :

    Range<mesh_type,MESH_ELEMENTS> M_rangeMeshElements;

    // physical parameters
    materialsproperties_ptrtype M_materialsProperties;

    // time discretisation
    std::string M_timeStepping;
    double M_timeStepThetaValue;
    vector_ptrtype M_timeStepThetaSchemePreviousContrib;

    // stabilization
    bool M_applyStabilization;
    std::string M_stabilizationType;
    stab_gls_parameter_ptrtype M_stabilizationGLSParameter;
    bool M_stabilizationGLS_applyShockCapturing;
    bool M_stabilizationDoAssemblyWithGradDiffusionCoeff;

    // post-process
    export_ptrtype M_exporter;
};


template< typename ConvexType>
template <typename SymbolsExprType>
bool
CoefficientFormPDEBase<ConvexType>::hasSymbolDependencyInCoefficients( std::set<std::string> const& symbs, SymbolsExprType const& se ) const
{
    for ( auto const& [physicId,physicData] : this->physicsFromCurrentType() )
    {
        auto physicCFPDEData = std::static_pointer_cast<ModelPhysicCoefficientFormPDE<nDim>>(physicData);
        //for ( Coefficient c : { convection, diffusion, reaction, firstTimeDerivative, secondTimeDerivative, source, conservativeFluxConvection, conservativeFluxSource, curlCurl } )
        for ( auto const& [c,props] : physicCFPDEData->infos()->coefficientProperties() )
            if ( physicCFPDEData->hasCoefficient( c ) && physicCFPDEData->coefficient( c ).hasSymbolDependency( symbs, se ) )
                return true;
    }
    return false;
}



} // namespace Feel
} // namespace FeelModels

#endif
