/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#ifndef FEELPP_TOOLBOXES_MAGNETIC_MAGNETICBOUNDARYCONDITIONS_HPP
#define FEELPP_TOOLBOXES_MAGNETIC_MAGNETICBOUNDARYCONDITIONS_HPP

#include <feel/feelmodels/modelcore/genericboundaryconditions.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
class MagneticBoundaryConditions : public BoundaryConditionsBase
{
    using super_type = BoundaryConditionsBase;
    using self_type = MagneticBoundaryConditions<Dim>;
    static constexpr uint16_type nRealDim = Dim;
public:
    enum class Type { MagneticPotentialImposed=0 };

    class MagneticPotentialImposed : public GenericDirichletBoundaryCondition<nRealDim,1>
    {
        using super_type = GenericDirichletBoundaryCondition<nRealDim,1>;
    public:
        MagneticPotentialImposed( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name, tbParent ) {}
        MagneticPotentialImposed( MagneticPotentialImposed const& ) = default;
        MagneticPotentialImposed( MagneticPotentialImposed && ) = default;
    };

    MagneticBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    MagneticBoundaryConditions( MagneticBoundaryConditions const& ) = default;
    MagneticBoundaryConditions( MagneticBoundaryConditions && ) = default;

    //! return magnetic potential imposed
    std::map<std::string,std::shared_ptr<MagneticPotentialImposed>> const& magneticPotentialImposed() const { return M_magneticPotentialImposed; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_magneticPotentialImposed.empty(); }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_magneticPotentialImposed, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_magneticPotentialImposed, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

private:
    std::map<std::string,std::shared_ptr<MagneticPotentialImposed>> M_magneticPotentialImposed;

}; // MagneticBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
