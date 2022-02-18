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

class ElectricBoundaryConditions
{
public:
    using self_type = ElectricBoundaryConditions;

    enum class Type { ElectricPotentialImposed=0, Ground, SurfaceChargeDensity };

    class ElectricPotentialImposed
    {
    public:
        ElectricPotentialImposed( std::string const& name ) : M_name( name ) {}
        ElectricPotentialImposed( ElectricPotentialImposed const& ) = default;
        ElectricPotentialImposed( ElectricPotentialImposed && ) = default;

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

        virtual self_type::Type type() const { return self_type::Type::ElectricPotentialImposed; }

        //! update informations
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
    };

    class Ground : public ElectricPotentialImposed
    {
        using super_type = ElectricPotentialImposed;
    public:

        Ground( std::string const& name ) : super_type( name ) {}
        Ground( Ground const& ) = default;
        Ground( Ground && ) = default;
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        self_type::Type type() const override { return self_type::Type::ElectricPotentialImposed; }

        //! update informations
        void updateInformationObject( nl::json & p ) const override;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
    };

    class SurfaceChargeDensity
    {
    public:
        SurfaceChargeDensity( std::string const& name ) : M_name( name ) {}
        SurfaceChargeDensity( SurfaceChargeDensity const& ) = default;
        SurfaceChargeDensity( SurfaceChargeDensity && ) = default;

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


    ElectricBoundaryConditions() = default;
    ElectricBoundaryConditions( ElectricBoundaryConditions const& ) = default;
    ElectricBoundaryConditions( ElectricBoundaryConditions && ) = default;

    //! return electric potential imposed
    std::map<std::pair<Type,std::string>,std::shared_ptr<ElectricPotentialImposed>> const& electricPotentialImposed() const { return M_electricPotentialImposed; }
    //! return charge density
    std::map<std::string,std::shared_ptr<SurfaceChargeDensity>> const& surfaceChargeDensity() const { return M_surfaceChargeDensity; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_electricPotentialImposed.empty(); }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( ModelBase const& mparent, nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::pair<Type,std::string>,std::shared_ptr<ElectricPotentialImposed>> M_electricPotentialImposed;
    std::map<std::string,std::shared_ptr<SurfaceChargeDensity>> M_surfaceChargeDensity;

}; // ElectricBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
