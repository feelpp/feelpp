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
    using infos_type = typename ModelPhysicCoefficientFormPDE<Dim>::infos_type;
    using Coefficient = typename ModelPhysicCoefficientFormPDE<Dim>::Coefficient;

    using infos_ptrtype = std::shared_ptr<infos_type>;

    ModelGenericPDE( infos_ptrtype const& infos );
    ModelGenericPDE( ModelGenericPDE const& ) = default;
    ModelGenericPDE( ModelGenericPDE && ) = default;

    std::string const& physic() const { return this->equationName(); }
    std::string const& equationName() const { return M_infos->equationName(); }
    std::string const& unknownName() const { return M_infos->unknownName(); }
    std::string const& unknownSymbol() const { return M_infos->unknownSymbol(); }
    std::string const& unknownBasis() const { return M_infos->unknownBasis(); }

private :
    friend class ModelPhysicCoefficientFormPDE<Dim>;

    infos_ptrtype M_infos;
};

template <uint16_type Dim>
class ModelGenericPDEs : public ModelPhysics<Dim>
{
    using super_type = ModelPhysics<Dim>;
    inline static const uint16_type nDim = Dim;
public :
    ModelGenericPDEs();
    ModelGenericPDEs( ModelGenericPDEs const& ) = default;
    ModelGenericPDEs( ModelGenericPDEs && ) = default;

    //auto const& pdes() const { return M_pdes; }
    //auto & pdes() { return M_pdes; }

private :
    //std::vector<std::tuple<typename ModelGenericPDE<nDim>::infos_type,std::shared_ptr<ModelGenericPDE<nDim>>>> M_pdes;
};

} // namespace FeelModels
} // namespace Feel

#endif
