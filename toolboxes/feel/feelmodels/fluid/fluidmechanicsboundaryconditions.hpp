/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_FLUID_FLUIDMECHANICSBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_FLUID_FLUIDMECHANICSBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class FluidMechanicsBoundaryConditions
{
public:
    using self_type = FluidMechanicsBoundaryConditions<Dim>;
    enum class Type { VelocityImposed=0, MeshVelocityImposed };

    class VelocityImposed : public GenericDirichletBoundaryCondition<Dim,1>
    {
        using super_type = GenericDirichletBoundaryCondition<Dim,1>;
    public:
        VelocityImposed( std::string const& name ) : super_type( name ) {}
        VelocityImposed( VelocityImposed const& ) = default;
        VelocityImposed( VelocityImposed && ) = default;

        virtual self_type::Type type() const { return self_type::Type::VelocityImposed; }
    };

    class MeshVelocityImposed : public VelocityImposed
    {
        using super_type = VelocityImposed;
    public:
        MeshVelocityImposed( std::string const& name ) : super_type( name ) {}
        MeshVelocityImposed( MeshVelocityImposed const& ) = default;
        MeshVelocityImposed( MeshVelocityImposed && ) = default;
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes ) override;

        self_type::Type type() const override { return self_type::Type::MeshVelocityImposed; }

        //! update informations
        void updateInformationObject( nl::json & p ) const override;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );
    };


#if 0
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
#endif
    FluidMechanicsBoundaryConditions() = default;
    FluidMechanicsBoundaryConditions( FluidMechanicsBoundaryConditions const& ) = default;
    FluidMechanicsBoundaryConditions( FluidMechanicsBoundaryConditions && ) = default;

    //! return temperature imposed
    std::map<std::pair<Type,std::string>,std::shared_ptr<VelocityImposed>> const& velocityImposed() const { return M_velocityImposed; }
    //! return heat flux
    //std::map<std::string,std::shared_ptr<HeatFlux>> const& heatFlux() const { return M_heatFlux; }
    //! return convective heat flux
    //std::map<std::string,std::shared_ptr<ConvectiveHeatFlux>> const& convectiveHeatFlux() const { return M_convectiveHeatFlux; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_velocityImposed.empty(); }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_velocityImposed, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_velocityImposed, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( ModelBase const& mparent, nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::pair<Type,std::string>,std::shared_ptr<VelocityImposed>> M_velocityImposed;
    //std::map<std::string,std::shared_ptr<HeatFlux>> M_heatFlux;
    //std::map<std::string,std::shared_ptr<ConvectiveHeatFlux>> M_convectiveHeatFlux;

}; // HeatBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
