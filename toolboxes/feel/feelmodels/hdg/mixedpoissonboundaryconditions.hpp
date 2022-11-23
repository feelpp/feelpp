/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_HDG_MIXEDPOISSONBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_HDG_MIXEDPOISSONBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim,uint8_type EquationRank>
class HDGMixedPoissonBoundaryConditions : public BoundaryConditionsBase
{
    using super_type = BoundaryConditionsBase;
public:
    using self_type = HDGMixedPoissonBoundaryConditions<Dim,EquationRank>;

    class Dirichlet : public GenericDirichletBoundaryCondition<(EquationRank>=1)?Dim:1,(EquationRank>=2)?Dim:1>
    {
        using super_type = GenericDirichletBoundaryCondition<(EquationRank>=1)?Dim:1,(EquationRank>=2)?Dim:1>;
    public:
        Dirichlet( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name,tbParent ) {}
        Dirichlet( Dirichlet const& ) = default;
        Dirichlet( Dirichlet && ) = default;
    };

    class Neumann
    {
    public:
        Neumann( std::string const& name ) : M_name( name ) {}
        Neumann( Neumann const& ) = default;
        Neumann( Neumann && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            static constexpr uint16_type dim1 = (EquationRank>=1)?Dim:1;
            return Feel::vf::expr( M_mexpr.template expr<dim1,1>(), se );
        }

        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

        //! update information
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
    };

    class Robin
    {
    public:
        Robin( std::string const& name ) : M_name( name ) {}
        Robin( Robin const& ) = default;
        Robin( Robin && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr1( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr1.template expr<1,1>(), se );
        }

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr2( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            static constexpr uint16_type dim1 = (EquationRank>=1)?Dim:1;
            return Feel::vf::expr( M_mexpr2.template expr<dim1,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr1.setParameterValues( paramValues ); M_mexpr2.setParameterValues( paramValues ); }

        //! update information
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr1, M_mexpr2;
        std::set<std::string> M_markers;
    };

    class Integral
    {
    public:
        Integral( std::string const& name ) : M_name( name ) {}
        Integral( Integral const& ) = default;
        Integral( Integral && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            static constexpr uint16_type dim1 = (EquationRank>=1)?Dim:1;
            return Feel::vf::expr( M_mexpr.template expr<dim1,1>(), se );
        }

        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

        //! update information
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
    };

    class CouplingODEs
    {
    public:
        CouplingODEs( std::string const& name ) : M_name( name ) {}
        CouplingODEs( CouplingODEs const& ) = default;
        CouplingODEs( CouplingODEs && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        std::string const& circuit() const { return M_circuit; }
        std::string const& capacitor() const { return M_capacitor; }
        std::string const& resistor() const { return M_resistor; }
        std::string const& buffer() const { return M_buffer; }

        //! update information
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        std::set<std::string> M_markers;

        std::string M_circuit;
        std::string M_capacitor;
        std::string M_resistor;
        std::string M_buffer;
    };


    HDGMixedPoissonBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    HDGMixedPoissonBoundaryConditions( HDGMixedPoissonBoundaryConditions const& ) = default;
    HDGMixedPoissonBoundaryConditions( HDGMixedPoissonBoundaryConditions && ) = default;

    //! return Dirichlet boundary conditions
    std::map<std::string,std::shared_ptr<Dirichlet>> const& dirichlet() const { return M_dirichlet; }

    //! return Neumann boundary conditions
    std::map<std::string,std::shared_ptr<Neumann>> const& neumann() const { return M_neumann; }

    //! return Neumann boundary conditions
    std::map<std::string,std::shared_ptr<Robin>> const& robin() const { return M_robin; }

    //! return Integral boundary conditions
    std::map<std::string,std::shared_ptr<Integral>> const& integral() const { return M_integral; }

    std::map<std::string,std::shared_ptr<CouplingODEs>> const& couplingODEs() const { return M_couplingODEs; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_dirichlet.empty(); }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_dirichlet, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_dirichlet, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( nl::json const& jarg );

    //! update information
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::string,std::shared_ptr<Dirichlet>> M_dirichlet;
    std::map<std::string,std::shared_ptr<Neumann>> M_neumann;
    std::map<std::string,std::shared_ptr<Robin>> M_robin;
    std::map<std::string,std::shared_ptr<Integral>> M_integral;
    std::map<std::string,std::shared_ptr<CouplingODEs>> M_couplingODEs;
}; // HDGMixedPoissonBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
