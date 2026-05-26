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

    using vector_potential_imposed_base_type = GenericDirichletBoundaryCondition<nRealDim,1>;
    //! n x A = g
    class MagneticPotentialImposed : public GenericDirichletBoundaryCondition<nRealDim,1>
    {
        using super_type = GenericDirichletBoundaryCondition<nRealDim,1>;
    public:
        MagneticPotentialImposed( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name, tbParent ) {}
        MagneticPotentialImposed( MagneticPotentialImposed const& ) = default;
        MagneticPotentialImposed( MagneticPotentialImposed && ) = default;
    };

    //! n x A = 0
    class MagneticInsulation : public GenericDirichletBoundaryCondition<nRealDim,1>
    {
        using super_type = GenericDirichletBoundaryCondition<nRealDim,1>;
    public:
        MagneticInsulation( std::string const& name, std::shared_ptr<ModelBase> const& tbParent ) : super_type( name, tbParent ) {}
        MagneticInsulation( MagneticInsulation const& ) = default;
        MagneticInsulation( MagneticInsulation && ) = default;
    };

    MagneticBoundaryConditions( std::shared_ptr<ModelBase> const& tbParent ) : super_type( tbParent ) {}
    MagneticBoundaryConditions( MagneticBoundaryConditions const& ) = default;
    MagneticBoundaryConditions( MagneticBoundaryConditions && ) = default;

    //! return magnetic potential imposed
    std::map<std::string,std::shared_ptr<MagneticPotentialImposed>> const& magneticPotentialImposed() const { return M_magneticPotentialImposed; }
    //! return magnetic insulation
    std::map<std::string,std::shared_ptr<MagneticInsulation>> const& magneticInsulation() const { return M_magneticInsulation; }

    //! return true if a bc is type of dof eliminitation
    bool hasTypeDofElimination() const { return !M_magneticPotentialImposed.empty() || !M_magneticInsulation.empty(); }

    //! apply dof elimination in linear context
    template <typename BfType, typename RhsType,typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyDofEliminationLinear( BfType& bilinearForm, RhsType& F, MeshType const& mesh, EltType const& u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_magneticPotentialImposed, bilinearForm, F, mesh, u, se );
            Feel::FeelModels::detail::applyDofEliminationLinearOnBoundaryConditions( M_magneticInsulation, bilinearForm, F, mesh, u, se );
        }

    //! apply Newton initial guess (on dof elimination context)
    template <typename MeshType, typename EltType, typename SymbolsExprType>
    void
    applyNewtonInitialGuess( MeshType const& mesh, EltType & u, SymbolsExprType const& se ) const
        {
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_magneticPotentialImposed, mesh, u, se );
            Feel::FeelModels::detail::applyNewtonInitialGuessOnBoundaryConditions( M_magneticInsulation, mesh, u, se );
        }

    void setParameterValues( std::map<std::string,double> const& paramValues );

    //! setup bc from json
    void setup( nl::json const& jarg );

    //! update informations
    void updateInformationObject( nl::json & p ) const;
    //! return tabulate information from json info
    static tabulate_informations_ptr_t tabulateInformations( nl::json const& jsonInfo, TabulateInformationProperties const& tabInfoProp );

    // helper functions to get vector potential bc by method
    std::vector<std::tuple<std::string,std::shared_ptr<vector_potential_imposed_base_type>>> anyBcWithVectorPotentialImposed( typename vector_potential_imposed_base_type::Method method ) const
        {
            std::vector<std::tuple<std::string,std::shared_ptr<vector_potential_imposed_base_type>>> ret;
            for ( auto const& [bcId,bcData] : M_magneticPotentialImposed )
                if ( bcData->isMethod( method ) )
                    ret.push_back( std::make_tuple( bcId, bcData ) );
            for ( auto const& [bcId,bcData] : M_magneticInsulation )
                if ( bcData->isMethod( method ) )
                    ret.push_back( std::make_tuple( bcId, bcData ) );
            return ret;
        }


private:
    std::map<std::string,std::shared_ptr<MagneticPotentialImposed>> M_magneticPotentialImposed;
    std::map<std::string,std::shared_ptr<MagneticInsulation>> M_magneticInsulation;

}; // MagneticBoundaryConditions

} // namespace FeelModels
} // namespace Feel

#endif
