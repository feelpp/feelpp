/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_MODELGENERICPDE_H
#define FEELPP_TOOLBOXES_MODELCORE_MODELGENERICPDE_H 1

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class ModelGenericPDE : public ModelPhysics<Dim>
{
    using super_type = ModelPhysics<Dim>;
    using material_property_shape_dim_type = typename super_type::material_property_shape_dim_type;
    inline static const uint16_type nDim = Dim;
public :
    using self_type = ModelGenericPDE<Dim>;
    struct Infos
    {
        Infos() = default;
        Infos( std::string const& name, pt::ptree const& p );
        Infos( std::string const& name, std::string const& unknownName, std::string const& unknownSymbol, std::string const& unknownBasis )
            :
            M_equationName( name ),
            M_unknownName( unknownName ),
            M_unknownSymbol( unknownSymbol ),
            M_unknownBasis( unknownBasis )
            {}
        Infos( Infos const& ) = default;
        Infos( Infos && ) = default;
        std::string const& equationName() const { return M_equationName; }
        std::string const& unknownName() const { return M_unknownName; }
        std::string const& unknownSymbol() const { return M_unknownSymbol; }
        std::string const& unknownBasis() const { return M_unknownBasis; }
    private :
        std::string M_equationName, M_unknownName, M_unknownSymbol, M_unknownBasis;
    };
    using infos_type = self_type::Infos;

    //ModelGenericPDE();
    ModelGenericPDE( infos_type const& infos );
    ModelGenericPDE( ModelGenericPDE const& ) = default;
    ModelGenericPDE( ModelGenericPDE && ) = default;

    std::string const& physic() const { return this->physicDefault(); }
    std::string const& equationName() const { return M_infos.equationName(); }
    std::string const& unknownName() const { return M_infos.unknownName(); }
    std::string const& unknownSymbol() const { return M_infos.unknownSymbol(); }
    std::string const& unknownBasis() const { return M_infos.unknownBasis(); }

    std::string convectionCoefficientName() const { return prefixvm( this->physic(), "beta", "_" ); }
    std::string diffusionCoefficientName() const { return prefixvm( this->physic(), "c", "_" ); }
    std::string reactionCoefficientName() const { return prefixvm( this->physic(), "a", "_" ); }
    std::string firstTimeDerivativeCoefficientName() const { return prefixvm( this->physic(), "d", "_" ); }
    std::string secondTimeDerivativeCoefficientName() const { return prefixvm( this->physic(), "m", "_" ); }
    std::string sourceCoefficientName() const { return prefixvm( this->physic(), "f", "_" ); }
    std::string conservativeFluxConvectionCoefficientName() const { return prefixvm( this->physic(), "alpha", "_" ); }
    std::string conservativeFluxSourceCoefficientName() const { return prefixvm( this->physic(), "gamma", "_" ); }
protected :
    void setupGenericPDE();

private :
    infos_type M_infos;
};

template <uint16_type Dim>
class ModelGenericPDEs : public ModelPhysics<Dim>
{
    using super_type = ModelPhysics<Dim>;
    using material_property_shape_dim_type = typename super_type::material_property_shape_dim_type;
    inline static const uint16_type nDim = Dim;
public :
    ModelGenericPDEs();
    ModelGenericPDEs( ModelGenericPDEs const& ) = default;
    ModelGenericPDEs( ModelGenericPDEs && ) = default;

    auto const& pdes() const { return M_pdes; }
    auto & pdes() { return M_pdes; }

    void addGenericPDE( typename ModelGenericPDE<nDim>::infos_type const& infos );
protected :
    void initGenericPDEs( std::string const& name );
    void setupGenericPDEs( pt::ptree const& p );
    void updateForUseGenericPDEs();
private :
    std::vector<std::tuple<typename ModelGenericPDE<nDim>::infos_type,std::shared_ptr<ModelGenericPDE<nDim>>>> M_pdes;
};

} // namespace FeelModels
} // namespace Feel

#endif
