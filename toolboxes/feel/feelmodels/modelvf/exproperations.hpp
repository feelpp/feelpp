/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_VF_EXPROPERATIONS_HPP
#define FEELPP_TOOLBOXES_VF_EXPROPERATIONS_HPP 1

namespace Feel
{
namespace FeelModels
{
namespace vfdetail
{

template<typename... Dummy>
auto
addExpr( hana::tuple<> const& /**/ )
{
    CHECK( false ) << "not allow";
    return 0*one();
}
template<typename T1>
auto
addExpr( hana::tuple<T1> const& t )
{
    return hana::at_c<0>( t );
}
template<typename T1,typename T2,typename... TOther>
auto
addExpr( hana::tuple<T1,T2,TOther...> const& t )
{
    return hana::at_c<0>( t ) + addExpr( hana::remove_at_c<0>( t ) );
}

} // namespace vfdetail
} // namespace FeelModels
} // namespace Feel

#endif
