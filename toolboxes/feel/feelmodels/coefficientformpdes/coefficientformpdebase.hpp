/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP 1

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelts/tsbase.hpp>

#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelgenericpde.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>
#include <feel/feelmodels/modelcore/stabilizationglsparameterbase.hpp>

namespace Feel
{
namespace FeelModels
{
template< typename ConvexType>
class CoefficientFormPDEBase : public ModelNumerical,
                               public ModelGenericPDE<ConvexType::nDim>
{
public :
    using super_type = ModelNumerical;
    using super2_type = ModelGenericPDE<ConvexType::nDim>;
    using size_type = typename super_type::size_type;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    static const uint16_type nRealDim = convex_type::nDim;
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
    CoefficientFormPDEBase( typename super2_type::infos_type const& infosPDE,
                            std::string const& prefix,
                            std::string const& keyword,
                            worldcomm_ptr_t const& worldComm,
                            std::string const& subPrefix,
                            ModelBaseRepository const& modelRep );
#if 0
    CoefficientFormPDEBase( std::string const& prefix,
                            std::string const& keyword,
                            worldcomm_ptr_t const& worldComm,
                            std::string const& subPrefix,
                            ModelBaseRepository const& modelRep )
        :
        CoefficientFormPDEBase( super2_type(), prefix, keyword, worldComm, subPrefix, modelRep )
        {}
#endif

    //! return current shared_ptr of type CoefficientFormPDEBase
    virtual std::shared_ptr<CoefficientFormPDEBase<ConvexType>> shared_from_this_cfpdebase() = 0;

    //! return current shared_ptr of type CoefficientFormPDEBase
    virtual std::shared_ptr<const CoefficientFormPDEBase<ConvexType>> shared_from_this_cfpdebase() const = 0;

    //! return true is the unknown is scalar
    virtual bool unknownIsScalar() const = 0;

    //___________________________________________________________________________________//
    // mesh
    mesh_ptrtype mesh() const { return super_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }
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
    // algebraic data and solver
    backend_ptrtype const& backend() const { return M_backend; }
    BlocksBaseVector<double> const& blockVectorSolution() const { return M_blockVectorSolution; }
    BlocksBaseVector<double> & blockVectorSolution() { return M_blockVectorSolution; }
    size_type nLocalDof() const;
    model_algebraic_factory_ptrtype const& algebraicFactory() const { return M_algebraicFactory; }
    model_algebraic_factory_ptrtype & algebraicFactory() { return M_algebraicFactory; }

    //int nBlockMatrixGraph() const { return 1; }

    virtual void setParameterValues( std::map<std::string,double> const& paramValues ) = 0;

    template <typename SymbolsExprType>
    bool hasSymbolDependencyInCoefficients( std::set<std::string> const& symbs, SymbolsExprType const& se ) const;

protected :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initMaterialProperties();
    void initBasePostProcess();

protected :

    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

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

    // algebraic data/tools
    backend_ptrtype M_backend;
    model_algebraic_factory_ptrtype M_algebraicFactory;
    BlocksBaseVector<double> M_blockVectorSolution;
};


template< typename ConvexType>
template <typename SymbolsExprType>
bool
CoefficientFormPDEBase<ConvexType>::hasSymbolDependencyInCoefficients( std::set<std::string> const& symbs, SymbolsExprType const& se ) const
{
    bool hasDependency = false;
    for ( std::string const& matName : this->materialsProperties()->physicToMaterials( this->physicDefault() ) )
    {
        if ( ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->reactionCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->reactionCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->secondTimeDerivativeCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->secondTimeDerivativeCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxConvectionCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->conservativeFluxConvectionCoefficientName() ).hasSymbolDependency( symbs, se ) ) ||
             ( this->materialsProperties()->hasProperty( matName, this->conservativeFluxSourceCoefficientName() ) &&
               this->materialsProperties()->materialProperty( matName, this->conservativeFluxSourceCoefficientName() ).hasSymbolDependency( symbs, se ) )

             )
        {
            hasDependency = true;
            break;
        }
    }
    return hasDependency;
}



} // namespace Feel
} // namespace FeelModels

#endif
