/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim1,uint16_type Dim2>
void
GenericDirichletBoundaryCondition<Dim1,Dim2>::setup( nl::json const& jarg, ModelIndexes const& indexes )
{
    auto tbParent = M_toolboxParent.lock();
    if ( jarg.contains( "expr" ) )
        M_mexpr.setExpr( jarg.at( "expr" ), tbParent->worldComm(), tbParent->repository().expr(), indexes );
     if ( jarg.contains( "markers" ) )
     {
         ModelMarkers markers;
         markers.setup( jarg.at("markers"), indexes );
         M_markers = markers;
     }
     else
         M_markers = { M_name };

     if ( Dim1 > 1 )
     {
         if ( jarg.contains( "component" ) )
         {
             auto const& j_comp = jarg.at( "component" );
             if ( j_comp.is_string() )
             {
                 std::string compStr = j_comp.template get<std::string>();
                 if ( compStr == "x" )
                     M_comp = ComponentType::X;
                 else if ( compStr == "y" )
                     M_comp = ComponentType::Y;
                 else if ( compStr == "z" )
                     M_comp = ComponentType::Z;
             }
         }
     }

     if ( jarg.contains( "method" ) )
     {
         std::string const& methodStr = jarg.at( "method" ).template get<std::string>();
         if ( methodStr == "elimination" )
             M_method = Method::elimination;
         else if ( methodStr == "nitsche" )
             M_method = Method::nitsche;
         else if ( methodStr == "lagrange_multiplier" ||  methodStr == "lm" )
             M_method = Method::lagrange_multiplier;
     }
}

template <uint16_type Dim1,uint16_type Dim2>
void
GenericDirichletBoundaryCondition<Dim1,Dim2>::updateInformationObject( nl::json & p ) const
{
    auto [exprStr,compInfo] = M_mexpr.exprInformations();
    p["expr"] = exprStr;
    p["markers"] = M_markers;
    if ( M_comp != ComponentType::NO_COMPONENT )
        p["components"] = M_comp;

    std::map<Method,std::string> methodToStr = { { Method::elimination, "elimination" }, { Method::nitsche, "Nitsche" }, { Method::lagrange_multiplier, "lagrange_multiplier" } };
    p["method"] = methodToStr[M_method];
}

template <uint16_type Dim1,uint16_type Dim2>
tabulate_informations_ptr_t
GenericDirichletBoundaryCondition<Dim1,Dim2>::tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp )
{
    Feel::Table tabInfo;
    TabulateInformationTools::FromJSON::addKeyToValues( tabInfo, jsonInfo, tabInfoProp, { /*"name","type",*/"method","expr" } );
    tabInfo.format()
        .setShowAllBorders( false )
        .setColumnSeparator(":")
        .setHasRowSeparator( false );
    if ( jsonInfo.contains("components") )
    {
        auto const& j_comps = jsonInfo.at("components");
        ComponentType comp = j_comps.template get<ComponentType>();
        Feel::Table tabInfoComp;
        tabInfoComp.format()
            .setShowAllBorders( false )
            .setColumnSeparator(":")
            .setHasRowSeparator( false );
        switch ( comp )
        {
        case ComponentType::X: tabInfoComp.add_row( {"x"} ); break;
        case ComponentType::Y: tabInfoComp.add_row( {"y"} ); break;
        case ComponentType::Z: tabInfoComp.add_row( {"z"} ); break;
        default: break;
        }
        CHECK( tabInfoComp.nRow() != 0 ) << "something wrong in json info";
        tabInfo.add_row( { "components", tabInfoComp } );
    }
    if ( jsonInfo.contains("markers") )
    {
        Feel::Table tabInfoMarkers = TabulateInformationTools::FromJSON::createTableFromArray( jsonInfo.at("markers") , true );
        if ( tabInfoMarkers.nRow() > 0 )
            tabInfo.add_row( { "markers", tabInfoMarkers } );
    }
    return TabulateInformations::New( tabInfo, tabInfoProp );
}




template class GenericDirichletBoundaryCondition<1,1>;
template class GenericDirichletBoundaryCondition<2,1>;
template class GenericDirichletBoundaryCondition<3,1>;


} // namespace FeelModels
} // namespace Feel
