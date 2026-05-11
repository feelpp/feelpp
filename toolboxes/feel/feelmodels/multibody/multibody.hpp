/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_BODYMOTION_HPP
#define FEELPP_TOOLBOXES_BODYMOTION_HPP 1

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>
#include <feel/feelmodels/multibody/body.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType>
class Multibody : public ModelNumerical,
                  public ModelPhysics<ConvexType::nDim>
{
    using super_type = ModelNumerical;
    using super_numerical_type = super_type;
    using super_physics_type = ModelPhysics<ConvexType::nDim>;
public:

    using size_type = typename super_type::size_type;
    typedef Multibody<ConvexType> self_type;
    //typedef std::shared_ptr<self_type> self_ptrtype;
    //___________________________________________________________________________________//
    // mesh
    typedef ConvexType convex_type;
    static inline const uint16_type nDim = convex_type::nDim;
    static inline const uint16_type nOrderGeo = convex_type::nOrder;
    static inline const uint16_type nRealDim = convex_type::nRealDim;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    // materials properties
    typedef MaterialsProperties<nRealDim> materialsproperties_type;
    typedef std::shared_ptr<materialsproperties_type> materialsproperties_ptrtype;

    // body
    using body_type = Body<convex_type>;

    Multibody( std::string const& prefix,
               std::string const& keyword = "multibody",
               worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
               ModelBaseRepository const& modelRep = ModelBaseRepository() );

    std::shared_ptr<self_type> shared_from_this() { return std::dynamic_pointer_cast<self_type>( super_type::shared_from_this() ); }


    void init();

    void applyRemesh( mesh_ptrtype oldMesh, mesh_ptrtype newMesh, std::shared_ptr<RemeshInterpolation> remeshInterp = std::make_shared<RemeshInterpolation>() );

    // information
    void updateInformationObject( nl::json & p ) const override;
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override;

    // materials properties
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    //materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

    // mesh
    mesh_ptrtype mesh() const { return super_numerical_type::super_model_meshes_type::mesh<mesh_type>( this->keyword() ); }
    void setMesh( mesh_ptrtype const& mesh ) { super_numerical_type::super_model_meshes_type::setMesh( this->keyword(), mesh ); }

    // bodies
    std::map<std::string, std::unique_ptr<body_type>> const& bodies() const noexcept { return M_bodies; }
    bool hasBody( std::string const& name ) const { return M_bodies.find( name ) != M_bodies.end(); }
    body_type const& body( std::string const& name ) const { return *M_bodies.at( name ); }
    body_type & body( std::string const& name ) { return *M_bodies.at( name ); }


private:
    materialsproperties_ptrtype M_materialsProperties;

    std::map<std::string, std::unique_ptr<body_type>> M_bodies;
};



} // namespace FeelModels
} // namespace Feel

#endif
