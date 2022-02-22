/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_HEAT_HEATBOUNDARYCONDITIONS_HPP

#include <feel/feelcore/json.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feeldiscr/enums.hpp>

#include <feel/feelmodels/modelcore/markermanagement.hpp>

namespace Feel
{
namespace FeelModels
{

namespace detail
{

template <ElementsType ET>
constexpr int indexInDistributeMarker()
{
    if constexpr ( ET == MESH_FACES )
        return 0;
    else if constexpr ( ET == MESH_EDGES )
        return 1;
    else if constexpr ( ET == MESH_POINTS )
        return 2;
    else //if constexpr ( ET == MESH_ELEMENTS )
        return 3;
}

template <ElementsType ET, typename MeshType, typename MarkerType>
auto rangeOfMarkedEntity( MeshType const& mesh, MarkerType const& markers )
{
    if constexpr ( ET == MESH_FACES )
        return markedfaces(mesh,markers );
    else if constexpr ( ET == MESH_EDGES )
        return markededges(mesh,markers );
    else if constexpr ( ET == MESH_POINTS )
        return markedpoints(mesh,markers );
    else //if constexpr ( ET == MESH_ELEMENTS )
        return markedelements(mesh,markers );
}

}

template <uint16_type Dim,uint8_type EquationRank>
class CoefficientFormPDEBoundaryConditions
{
public:
    using self_type = CoefficientFormPDEBoundaryConditions<Dim,EquationRank>;

    class Dirichlet
    {
    public:
        Dirichlet( std::string const& name ) : M_name( name ), M_comp( ComponentType::NO_COMPONENT ) {}
        Dirichlet( Dirichlet const& ) = default;
        Dirichlet( Dirichlet && ) = default;

        //! setup bc from json
        void setup( ModelBase const& mparent, nl::json const& jarg, ModelIndexes const& indexes );

        //! return expression
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            static constexpr uint16_type dim1 = (EquationRank>=1)?Dim:1;
            static constexpr uint16_type dim2 = (EquationRank>=2)?Dim:1;
            return Feel::vf::expr( M_mexpr.template expr<dim1,dim2>(), se );
        }
        template <typename SymbolsExprType = symbols_expression_empty_t>
        auto exprComponent( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<1,1>(), se );
        }
        //! return markers
        std::set<std::string> const& markers() const { return M_markers; }

        void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

        //! update informations
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );


        template <typename ToolboxType,typename SpaceType>
        void updateDofEliminationIds( ToolboxType & tb, std::string const& fieldName, std::shared_ptr<SpaceType> const& space ) const
            {
                auto mesh = space->mesh();
                using mesh_type = typename SpaceType::mesh_type;

                auto meshMarkersByEntities = Feel::FeelModels::detail::distributeMarkerListOnSubEntity( mesh, M_markers );
                ComponentType comp = M_comp;

                // on element
                auto const& listMarkedElementUnknown = std::get<3>( meshMarkersByEntities );
                if ( !listMarkedElementUnknown.empty() )
                    tb.updateDofEliminationIds( fieldName, space, markedelements( mesh,listMarkedElementUnknown ), comp );
                // on topological faces
                auto const& listMarkedFacesUnknown = std::get<0>( meshMarkersByEntities );
                if ( !listMarkedFacesUnknown.empty() )
                    tb.updateDofEliminationIds( fieldName, space, markedfaces( mesh,listMarkedFacesUnknown ), comp );
                // on marked edges (only 3d)
                if constexpr ( mesh_type::nDim == 3)
                {
                    auto const& listMarkedEdgesUnknown = std::get<1>( meshMarkersByEntities );
                    if ( !listMarkedEdgesUnknown.empty() )
                        tb.updateDofEliminationIds( fieldName, space, markededges( mesh,listMarkedEdgesUnknown ), comp );
                }
                // on marked points
                auto const& listMarkedPointsUnknown = std::get<2>( meshMarkersByEntities );
                if ( !listMarkedPointsUnknown.empty() )
                    tb.updateDofEliminationIds( fieldName, space, markedpoints( mesh,listMarkedPointsUnknown ), comp );
            }

        template <ElementsType ET, typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
        void
        applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
            {
                auto meshMarkersByEntities = Feel::FeelModels::detail::distributeMarkerListOnSubEntity( mesh, M_markers );
                ComponentType comp = M_comp;
                static const int indexDistrib = Feel::FeelModels::detail::indexInDistributeMarker<ET>();
                auto const& listMarkedEntities = std::get</*0*/indexDistrib >( meshMarkersByEntities );
                if ( listMarkedEntities.empty() )
                    return;
                if ( comp == ComponentType::NO_COMPONENT )
                {
                    auto theExpr = this->expr(se);
                    bilinearForm +=
                        on( _range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh,listMarkedEntities),
                            _element=u,_rhs=F,_expr=theExpr );
                }
                else if constexpr ( EquationRank == 1 )
                {
                    auto theExpr = this->exprComponent(se);
                    bilinearForm +=
                        on( _range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                            _element=u.comp( comp ),_rhs=F,_expr=theExpr );
                }
            }

        template <ElementsType ET, typename MeshType, typename EltType, typename SymbolsExprType>
        void
        applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
            {
                auto meshMarkersByEntities = Feel::FeelModels::detail::distributeMarkerListOnSubEntity( mesh, M_markers );
                ComponentType comp = ComponentType::NO_COMPONENT;
                static const int indexDistrib = Feel::FeelModels::detail::indexInDistributeMarker<ET>();
                auto const& listMarkedEntities = std::get</*0*/indexDistrib >( meshMarkersByEntities );
                if ( listMarkedEntities.empty() )
                    return;
                if ( comp == ComponentType::NO_COMPONENT )
                {
                    auto theExpr = this->expr(se);
                    u.on(_range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                         _expr=theExpr );
                }
                else if constexpr ( EquationRank == 1 )
                {
                    auto theExpr = this->exprComponent(se);
                    u.comp( comp ).on(_range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                                      _expr=theExpr );
                }
            }

    protected:
        std::string M_name;
        ModelExpression M_mexpr;
        std::set<std::string> M_markers;
        ComponentType M_comp;
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

        //! update informations
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

        //! update informations
        virtual void updateInformationObject( nl::json & p ) const;
        //! return tabulate information from json info
        static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    protected:
        std::string M_name;
        ModelExpression M_mexpr1, M_mexpr2;
        std::set<std::string> M_markers;
    };


    CoefficientFormPDEBoundaryConditions() = default;
    CoefficientFormPDEBoundaryConditions( CoefficientFormPDEBoundaryConditions const& ) = default;
    CoefficientFormPDEBoundaryConditions( CoefficientFormPDEBoundaryConditions && ) = default;

    //! return Dirichlet boundary conditions
    std::map<std::string,std::shared_ptr<Dirichlet>> const& dirichlet() const { return M_dirichlet; }

    //! return Neumann boundary conditions
    std::map<std::string,std::shared_ptr<Neumann>> const& neumann() const { return M_neumann; }

    //! return Neumann boundary conditions
    std::map<std::string,std::shared_ptr<Robin>> const& robin() const { return M_robin; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_dirichlet.empty(); }

    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            using space_type = typename unwrap_ptr_t<EltType>::functionspace_type;
            //using trial_space_type = typename BfType::trial_space_type;
            using mesh_type = unwrap_ptr_t<MeshType>;
            for ( auto const& [bcName,bcData] : M_dirichlet )
                bcData->template applyDofEliminationLinear<MESH_ELEMENTS>( bilinearForm,F,mesh,u,se );
            for ( auto const& [bcName,bcData] : M_dirichlet )
                bcData->template applyDofEliminationLinear<MESH_FACES>( bilinearForm,F,mesh,u,se );
            if constexpr ( !is_hcurl_conforming_v<typename space_type::fe_type> )
            {
                if constexpr ( mesh_type::nDim == 3 )
                    for ( auto const& [bcName,bcData] : M_dirichlet )
                        bcData->template applyDofEliminationLinear<MESH_EDGES>( bilinearForm, F, mesh, u, se );
                for ( auto const& [bcName,bcData] : M_dirichlet )
                    bcData->template applyDofEliminationLinear<MESH_POINTS>( bilinearForm, F, mesh, u, se );
            }
        }
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            using space_type = typename unwrap_ptr_t<EltType>::functionspace_type;
            using mesh_type = unwrap_ptr_t<MeshType>;
            for ( auto const& [bcName,bcData] : M_dirichlet )
                bcData->template applyNewtonInitialGuess<MESH_ELEMENTS>( mesh,u,se );
            for ( auto const& [bcName,bcData] : M_dirichlet )
                bcData->template applyNewtonInitialGuess<MESH_FACES>( mesh,u,se );
            if constexpr ( !is_hcurl_conforming_v<typename space_type::fe_type> )
            {
                if constexpr ( mesh_type::nDim == 3 )
                    for ( auto const& [bcName,bcData] : M_dirichlet )
                        bcData->template applyNewtonInitialGuess<MESH_EDGES>( mesh, u, se );
                for ( auto const& [bcName,bcData] : M_dirichlet )
                    bcData->template applyNewtonInitialGuess<MESH_POINTS>( mesh, u, se );
            }
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( ModelBase const& mparent, nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::string,std::shared_ptr<Dirichlet>> M_dirichlet;
    std::map<std::string,std::shared_ptr<Neumann>> M_neumann;
    std::map<std::string,std::shared_ptr<Robin>> M_robin;

}; // CoefficientFormPDEBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
