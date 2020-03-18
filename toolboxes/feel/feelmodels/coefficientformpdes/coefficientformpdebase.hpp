/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 */

#ifndef FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP
#define FEELPP_TOOLBOXES_COEFFICIENTFORMPDEBASE_HPP 1

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelmodels/modelcore/options.hpp>
#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelgenericpde.hpp>
#include <feel/feelmodels/modelcore/markermanagement.hpp>
#include <feel/feelmodels/modelalg/modelalgebraicfactory.hpp>
#include <feel/feelmodels/modelmaterials/materialsproperties.hpp>

namespace Feel
{
namespace FeelModels
{
template< typename ConvexType>
class CoefficientFormPDEBase : public ModelNumerical,
                               public ModelGenericPDE<ConvexType::nDim>
{
public :
    typedef ModelNumerical super_type;
    using size_type = typename super_type::size_type;

    // mesh
    typedef ConvexType convex_type;
    static const uint16_type nDim = convex_type::nDim;
    static const uint16_type nOrderGeo = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // materials properties
    typedef MaterialsProperties<mesh_type> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // algebraic solver
    typedef ModelAlgebraicFactory model_algebraic_factory_type;
    typedef std::shared_ptr< model_algebraic_factory_type > model_algebraic_factory_ptrtype;

public :
    CoefficientFormPDEBase( std::string const& prefix,
                            std::string const& keyword,
                            worldcomm_ptr_t const& worldComm,
                            std::string const& subPrefix,
                            ModelBaseRepository const& modelRep );

    std::string fileNameMeshPath() const { return prefixvm(this->prefix(),"CoefficientFormPDEMesh.path"); }

    //___________________________________________________________________________________//
    // mesh
    mesh_ptrtype const& mesh() const { return M_mesh; }
    elements_reference_wrapper_t<mesh_type> const& rangeMeshElements() const { return M_rangeMeshElements; }
    void setMesh( mesh_ptrtype const& mesh ) { M_mesh = mesh; }

    //___________________________________________________________________________________//
    // physical parameters
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }


protected :
    void loadParameterFromOptionsVm();
    void initMesh();
    void initMaterialProperties();

protected :

    mesh_ptrtype M_mesh;
    elements_reference_wrapper_t<mesh_type> M_rangeMeshElements;

    // physical parameters
    materialsproperties_ptrtype M_materialsProperties;

};


} // namespace Feel
} // namespace FeelModels

#endif
