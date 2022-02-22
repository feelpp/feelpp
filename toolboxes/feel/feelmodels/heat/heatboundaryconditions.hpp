/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

class HeatBoundaryConditions
{
public:
    class TemperatureImposed : public GenericDirichletBoundaryCondition<1,1>
    {
        using super_type = GenericDirichletBoundaryCondition<1,1>;
    public:
        TemperatureImposed( std::string const& name ) : super_type( name ) {}
        TemperatureImposed( TemperatureImposed const& ) = default;
        TemperatureImposed( TemperatureImposed && ) = default;
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
    std::map<std::string,std::shared_ptr<TemperatureImposed>> const& temperatureImposed() const { return M_temperatureImposed; }
    //! return heat flux
    std::map<std::string,std::shared_ptr<HeatFlux>> const& heatFlux() const { return M_heatFlux; }
    //! return convective heat flux
    std::map<std::string,std::shared_ptr<ConvectiveHeatFlux>> const& convectiveHeatFlux() const { return M_convectiveHeatFlux; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_temperatureImposed.empty(); }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_temperatureImposed, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_temperatureImposed, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( ModelBase const& mparent, nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::string,std::shared_ptr<TemperatureImposed>> M_temperatureImposed;
    std::map<std::string,std::shared_ptr<HeatFlux>> M_heatFlux;
    std::map<std::string,std::shared_ptr<ConvectiveHeatFlux>> M_convectiveHeatFlux;

}; // HeatBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
