/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4*/

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
    ModelGenericPDE();
    ModelGenericPDE( std::string const& name, pt::ptree const& p );
    ModelGenericPDE( ModelGenericPDE const& ) = default;
    ModelGenericPDE( ModelGenericPDE && ) = default;

    std::string const& physic() const { return this->physicDefault(); }

    std::string const& unknownName() const { return M_unknownName; }
    std::string const& unknownSymbol() const { return M_unknownSymbol; }
    std::string const& unknownBasis() const { return M_unknownBasis; }

    std::string convectionCoefficientName() const { return prefixvm( this->physic(), "beta", "_" ); }
    std::string diffusionCoefficientName() const { return prefixvm( this->physic(), "c", "_" ); }
    std::string reactionCoefficientName() const { return prefixvm( this->physic(), "a", "_" ); }
    std::string firstTimeDerivativeCoefficientName() const { return prefixvm( this->physic(), "d", "_" ); }
    std::string secondTimeDerivativeCoefficientName() const { return prefixvm( this->physic(), "m", "_" ); }
    std::string sourceCoefficientName() const { return prefixvm( this->physic(), "f", "_" ); }

protected :
    void setupGenericPDE( std::string const& name, pt::ptree const& p );

private :
    std::string M_unknownName, M_unknownSymbol, M_unknownBasis;

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

    std::vector<ModelGenericPDE<nDim>> const& pdes() const { return M_pdes; }
protected :
    void setupGenericPDEs( std::string const& name, pt::ptree const& p );
private :
    std::vector<ModelGenericPDE<nDim>> M_pdes;
};

} // namespace FeelModels
} // namespace Feel

#endif
