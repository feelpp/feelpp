/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP

#include <feel/feelcore/json.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>

namespace Feel
{
namespace FeelModels
{

class HeatBoundaryConditions
{
public:
    class TemperatureImposed
    {
    public:
        TemperatureImposed( std::string const& name ) : M_name( name ) {}
        TemperatureImposed( TemperatureImposed const& ) = default;
        TemperatureImposed( TemperatureImposed && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<1,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

        //! update informations
        void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    private:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
    };

    class HeatFlux
    {
    public:
        HeatFlux( std::string const& name ) : M_name( name ) {}
        HeatFlux( HeatFlux const& ) = default;
        HeatFlux( HeatFlux && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<1,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

        //! update informations
        void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    private:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
    };

    class ConvectiveHeatFlux
    {
    public:
        ConvectiveHeatFlux( std::string const& name ) : M_name( name ) {}
        ConvectiveHeatFlux( ConvectiveHeatFlux const& ) = default;
        ConvectiveHeatFlux( ConvectiveHeatFlux && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return h expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr_h( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr_h.template expr<1,1>(), se );
        }
        //! return T exterior expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr_Text( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr_Text.template expr<1,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues )
            {
                M_mexpr_h.setParameterValues( paramValues );
                M_mexpr_Text.setParameterValues( paramValues );
            }

        //! update informations
        void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    private:
        std::string M_name;
        ModelExpression M_mexpr_h, M_mexpr_Text;
        std::set<std::string> M_markers;
    };

    HeatBoundaryConditions() = default;
    HeatBoundaryConditions( HeatBoundaryConditions const& ) = default;
    HeatBoundaryConditions( HeatBoundaryConditions && ) = default;

    //! return temperature imposed
    std::map<std::string,TemperatureImposed> const& temperatureImposed() const { return M_temperatureImposed; }
    //! return heat flux
    std::map<std::string,HeatFlux> const& heatFlux() const { return M_heatFlux; }
    //! return convective heat flux
    std::map<std::string,ConvectiveHeatFlux> const& convectiveHeatFlux() const { return M_convectiveHeatFlux; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_temperatureImposed.empty(); }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( ModelBase const& mparent, nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::string,TemperatureImposed> M_temperatureImposed;
    std::map<std::string,HeatFlux> M_heatFlux;
    std::map<std::string,ConvectiveHeatFlux> M_convectiveHeatFlux;

}; // HeatBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
