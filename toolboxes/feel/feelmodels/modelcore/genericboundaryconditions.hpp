/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_CORE_GENERICBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_CORE_GENERICBOUNDARYCONDITIONS_HPP

#include <feel/feelcore/json.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelcore/modelbase.hpp>
#include <feel/feeldiscr/enums.hpp>

#include <feel/feelmodels/modelcore/markermanagement.hpp>


namespace Feel
{
namespace FeelModels
{

class BoundaryConditionsBase
{
public :
    BoundaryConditionsBase( std::shared_ptr<ModelBase> const& tbParent ) : M_toolboxParent( tbParent ) {}
    BoundaryConditionsBase( BoundaryConditionsBase const& ) = default;
    BoundaryConditionsBase( BoundaryConditionsBase && ) = default;

    std::shared_ptr<ModelBase> toolboxParent() const { return M_toolboxParent.lock(); }

protected:
    std::weak_ptr<ModelBase> M_toolboxParent;
};

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

//! apply dof elimination in linear context
template <typename BoundaryConditionsType, typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
void
applyDofEliminationLinearOnBoundaryConditions( BoundaryConditionsType const& bcs, BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se )
{
    using space_type = typename unwrap_ptr_t<EltType>::functionspace_type;
    //using trial_space_type = typename BfType::trial_space_type;
    using mesh_type = unwrap_ptr_t<MeshType>;
    for ( auto const& [bcId,bcData] : bcs )
        bcData->template applyDofEliminationLinear<MESH_ELEMENTS>( bilinearForm,F,mesh,u,se );
    for ( auto const& [bcId,bcData] : bcs )
        bcData->template applyDofEliminationLinear<MESH_FACES>( bilinearForm,F,mesh,u,se );
    if constexpr ( !is_hcurl_conforming_v<typename space_type::fe_type> )
    {
        if constexpr ( mesh_type::nDim == 3 )
            for ( auto const& [bcId,bcData] : bcs )
                bcData->template applyDofEliminationLinear<MESH_EDGES>( bilinearForm, F, mesh, u, se );
        for ( auto const& [bcId,bcData] : bcs )
            bcData->template applyDofEliminationLinear<MESH_POINTS>( bilinearForm, F, mesh, u, se );
    }
}

//! apply Newton initial guess (on dof elimination context)
template <typename BoundaryConditionsType, typename MeshType, typename EltType, typename SymbolsExprType>
void
applyNewtonInitialGuessOnBoundaryConditions( BoundaryConditionsType const& bcs, MeshType const& mesh, EltType & u, SymbolsExprType const& se )
{
    using space_type = typename unwrap_ptr_t<EltType>::functionspace_type;
    using mesh_type = unwrap_ptr_t<MeshType>;
    for ( auto const& [bcId,bcData] : bcs )
        bcData->template applyNewtonInitialGuess<MESH_ELEMENTS>( mesh,u,se );
    for ( auto const& [bcId,bcData] : bcs )
        bcData->template applyNewtonInitialGuess<MESH_FACES>( mesh,u,se );
    if constexpr ( !is_hcurl_conforming_v<typename space_type::fe_type> )
    {
        if constexpr ( mesh_type::nDim == 3 )
            for ( auto const& [bcId,bcData] : bcs )
                bcData->template applyNewtonInitialGuess<MESH_EDGES>( mesh, u, se );
        for ( auto const& [bcId,bcData] : bcs )
            bcData->template applyNewtonInitialGuess<MESH_POINTS>( mesh, u, se );
    }
}

} // detail

template <uint16_type Dim1,uint16_type Dim2>
class GenericDirichletBoundaryCondition
{
public:
    enum class Method { elimination=0, nitsche, lagrange_multiplier };

    GenericDirichletBoundaryCondition( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : M_name( name ), M_toolboxParent( tbParent ), M_comp( ComponentType::NO_COMPONENT ), M_method( Method::elimination ) {}
    GenericDirichletBoundaryCondition( GenericDirichletBoundaryCondition const& ) = default;
    GenericDirichletBoundaryCondition( GenericDirichletBoundaryCondition && ) = default;

    //! setup bc from json
    virtual void setup( nl::json const& jarg, ModelIndexes const& indexes );

    //! return expression
    template <typename SymbolsExprType = symbols_expression_empty_t>
    auto expr( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<Dim1,Dim2>(), se );
        }
    template <typename SymbolsExprType = symbols_expression_empty_t>
    auto exprComponent( SymbolsExprType const& se = symbols_expression_empty_t{} ) const
        {
            return Feel::vf::expr( M_mexpr.template expr<1,1>(), se );
        }

    //! return markers
    std::set<std::string> const& markers() const { return M_markers; }

    void setParameterValues( std::map<std::string,double> const& paramValues ) { M_mexpr.setParameterValues( paramValues ); }

    //! return method used for apply Dirichlet treatment
    Method method() const { return M_method; }
    //! return true if this method is used
    bool isMethod( Method method ) const { return M_method == method; }
    //! return true if method is elimination
    bool isMethodElimination() const { return this->isMethod( Method::elimination ); }
    //! return true if method is nitsche
    bool isMethodNitsche() const { return this->isMethod( Method::nitsche ); }
    //! return true if method is lagrange_multiplier
    bool isMethodLagrangeMultiplier() const { return this->isMethod( Method::lagrange_multiplier ); }

    //! update information
    virtual void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );


    template <typename ToolboxType,typename SpaceType>
    void updateDofEliminationIds( ToolboxType & tb, std::string const& fieldName, std::shared_ptr<SpaceType> const& space ) const
        {
            if ( !this->isMethodElimination() )
                return;
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
            if ( !this->isMethodElimination() )
                return;
            auto tbParent = M_toolboxParent.lock();
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
                        _element=u,_rhs=F,_expr=theExpr,
                        _vm=tbParent->clovm(),_prefix=tbParent->prefix() );
            }
            else if constexpr ( Dim1 > 1 && Dim2 == 1 )
            {
                auto theExpr = this->exprComponent(se);
                bilinearForm +=
                    on( _range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                        _element=u.comp( comp ),_rhs=F,_expr=theExpr,
                        _vm=tbParent->clovm(),_prefix=tbParent->prefix() );
            }
        }

    template <ElementsType ET, typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            if ( !this->isMethodElimination() )
                return;
            auto meshMarkersByEntities = Feel::FeelModels::detail::distributeMarkerListOnSubEntity( mesh, M_markers );
            ComponentType comp = M_comp;
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
            else if constexpr ( Dim1 > 1 && Dim2 == 1 )
            {
                auto theExpr = this->exprComponent(se);
                u.comp( comp ).on(_range=Feel::FeelModels::detail::rangeOfMarkedEntity<ET>(mesh, listMarkedEntities),
                                  _expr=theExpr );
            }
        }

protected:
    std::string M_name;
    std::weak_ptr<ModelBase> M_toolboxParent;
    ModelExpression M_mexpr;
    std::set<std::string> M_markers;
    ComponentType M_comp;
    Method M_method;
};

} // namespace FeelModels
} // namespace Feel


#endif
