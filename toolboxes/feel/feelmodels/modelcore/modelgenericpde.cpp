/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelgenericpde.hpp>


namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelGenericPDE<Dim>::ModelGenericPDE( infos_ptrtype const& infos )
    :
    super_type( "GenericPDE" ),
    ModelBase(""),
    M_infos( infos )
{}

template class ModelGenericPDE<2>;
template class ModelGenericPDE<3>;

} // namespace FeelModels
} // namespace Feel

