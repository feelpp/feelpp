/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_SOLID_SOLIDMECHANICSBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_SOLID_SOLIDMECHANICSBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class SolidMechanicsBoundaryConditions : public BoundaryConditionsBase
{
    using super_type = BoundaryConditionsBase;
public:
    using self_type = SolidMechanicsBoundaryConditions<Dim>;
    enum class Type { DisplacementImposed=0 };

    class DisplacementImposed : public GenericDirichletBoundaryCondition<Dim,1>
    {
        using super_type = GenericDirichletBoundaryCondition<Dim,1>;
    public:
        DisplacementImposed( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name, tbParent ) {}
        DisplacementImposed( DisplacementImposed const& ) = default;
        DisplacementImposed( DisplacementImposed && ) = default;

        virtual self_type::Type type() const { return self_type::Type::DisplacementImposed; }
    };

    class NormalStress
    {
    public:
        enum class Frame { Lagrangian=0, Eulerian };

        NormalStress( std::string const& name ) : M_name( name ), M_frame( Frame::Lagrangian ) {}
        NormalStress( NormalStress const& ) = default;
        NormalStress( NormalStress && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return frame
        Frame frame() const { return M_frame; }

        bool isFrameLagrangian() const { return M_frame == Frame::Lagrangian; }
        bool isFrameEulerian() const { return M_frame == Frame::Eulerian; }

        //! return true if expr is scalar
        bool isScalarExpr() const { return M_mexpr.template hasExpr<1,1>(); }

        //! return true if expr is vectorial
        bool isVectorialExpr() const { return M_mexpr.template hasExpr<Dim,1>(); }

        //! return true if expr is matrix
        bool isMatrixExpr() const { return M_mexpr.template hasExpr<Dim,Dim>(); }

        //! return scalar expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto exprScalar( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<1,1>(), se );
        }

        //! return vectorial expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto exprVectorial( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<Dim,1>(), se );
        }

        //! return matrix expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto exprMatrix( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<Dim,Dim>(), se );
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
        Frame M_frame;
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
            return Feel::vf::expr( M_mexpr2.template expr<Dim,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr1.setParameterValues( paramValues ); M_mexpr2.setParameterValues( paramValues ); }

        //! update informations
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr1, M_mexpr2;
        std::set<std::string> M_markers;
    };


    SolidMechanicsBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    SolidMechanicsBoundaryConditions( SolidMechanicsBoundaryConditions const& ) = default;
    SolidMechanicsBoundaryConditions( SolidMechanicsBoundaryConditions && ) = default;

    //! return velocity imposed
    std::map<std::pair<Type,std::string>,std::shared_ptr<DisplacementImposed>> const& displacementImposed() const { return M_displacementImposed; }
    //! return normal stress
    std::map<std::string,std::shared_ptr<NormalStress>> const& normalStress() const { return M_normalStress; }
    //! return robin
    std::map<std::string,std::shared_ptr<Robin>> const& robin() const { return M_robin; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const
        {
            for ( auto const& [bcId,bcData] : M_displacementImposed )
                if ( bcData->isMethodElimination() )
                    return true;
            return false;
        }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_displacementImposed, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_displacementImposed, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::pair<Type,std::string>,std::shared_ptr<DisplacementImposed>> M_displacementImposed;
    std::map<std::string,std::shared_ptr<NormalStress>> M_normalStress;
    std::map<std::string,std::shared_ptr<Robin>> M_robin;

}; // SolidMechanicsBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
