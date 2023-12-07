/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_BODYMOTION_HPP
#define FEELPP_TOOLBOXES_BODYMOTION_HPP 1

#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType>
class BodyMotion : public ModelNumerical,
                   public ModelPhysics<ConvexType::nDim>,
                   public std::enable_shared_from_this< BodyMotion<ConvexType> >
{
    using super_numerical_type = ModelNumerical;
    using super_physics_type = ModelPhysics<ConvexType::nDim>;
public:

    using size_type = typename super_type::size_type;
    typedef BodyMotion<ConvexType> self_type;
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

    BodyMotion( std::string const& prefix,
                std::string const& keyword = "body",
                worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                ModelBaseRepository const& modelRep = ModelBaseRepository() );

    void init();

    void updateInformationObject( nl::json & p ) const override {}
    tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp ) const override
        {
            auto tabInfo = TabulateInformationsSections::New();
            return tabInfo;
        }


    // materials properties
    materialsproperties_ptrtype const& materialsProperties() const { return M_materialsProperties; }
    materialsproperties_ptrtype & materialsProperties() { return M_materialsProperties; }
    void setMaterialsProperties( materialsproperties_ptrtype mp ) { M_materialsProperties = mp; }

private:
    materialsproperties_ptrtype M_materialsProperties;
};


template< typename ConvexType>
BodyMotion<ConvexType>::BodyMotion( std::string const& prefix,
                                    std::string const& keyword,
                                    worldcomm_ptr_t const& worldComm,
                                    ModelBaseRepository const& modelRep )
    :
    super_numerical_type( prefix, keyword, worldComm, "", modelRep ),
    super_physics_type( "body" ),
    ModelBase( prefix, keyword, worldComm, "", modelRep )
{}


template< typename ConvexType>
void
BodyMotion<ConvexType>::init()
{

}


} // namespace FeelModels
} // namespace Feel

#endif
