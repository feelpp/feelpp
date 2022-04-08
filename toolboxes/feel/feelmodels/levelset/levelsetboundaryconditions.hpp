/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_LEVELSET_LEVELSETBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_LEVELSET_LEVELSETBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

class LevelSetBoundaryConditions : public BoundaryConditionsBase
{
    using super_type = BoundaryConditionsBase;
public:

    LevelSetBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    LevelSetBoundaryConditions( LevelSetBoundaryConditions const& ) = default;
    LevelSetBoundaryConditions( LevelSetBoundaryConditions && ) = default;

    //! setup bc from json
    void setup( nl::json const& jarg );

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! return true if a bc is type of dof eliminitation
    //bool hasTypeDofElimination() const { //TODO }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
    {
        //TODO
    }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
    {
        //TODO
    }

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:

}; // LevelSetBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
