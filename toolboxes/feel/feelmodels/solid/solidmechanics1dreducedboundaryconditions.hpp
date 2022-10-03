/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_SOLID_SOLIDMECHANICS1DREDUCEDBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_SOLID_SOLIDMECHANICS1DREDUCEDBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class SolidMechanics1dReducedBoundaryConditions : public BoundaryConditionsBase
{
    using super_type = BoundaryConditionsBase;
public:
    using self_type = SolidMechanics1dReducedBoundaryConditions<Dim>;
    enum class Type { DisplacementImposed=0 };

    class DisplacementImposed : public GenericDirichletBoundaryCondition<1,1>
    {
        using super_type = GenericDirichletBoundaryCondition<1,1>;
    public:
        DisplacementImposed( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name, tbParent ) {}
        DisplacementImposed( DisplacementImposed const& ) = default;
        DisplacementImposed( DisplacementImposed && ) = default;

        virtual self_type::Type type() const { return self_type::Type::DisplacementImposed; }
    };


    SolidMechanics1dReducedBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    SolidMechanics1dReducedBoundaryConditions( SolidMechanics1dReducedBoundaryConditions const& ) = default;
    SolidMechanics1dReducedBoundaryConditions( SolidMechanics1dReducedBoundaryConditions && ) = default;

    //! return velocity imposed
    std::map<std::pair<Type,std::string>,std::shared_ptr<DisplacementImposed>> const& displacementImposed() const { return M_displacementImposed; }

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

}; // SolidMechanics1dReducedBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
