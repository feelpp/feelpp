#ifndef FEELPP_TOOLBOXES_MIXEDPOISSON_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_MIXEDPOISSON_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template< typename ConvexType, int Dim, int E_Order>
template< typename ModelConvexType>
void
MixedPoisson<ConvexType, Dim, E_Order>::updateLinearPDE( DataUpdateHDG & data, ModelConvexType const& mctx ) const
{
    condensed_matrix_ptr_t<value_type>& A = data.matrix();
    condensed_vector_ptr_t<value_type>& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    auto mesh = this->mesh();
    auto ps = this->spaceProduct();
    auto bbf = blockform2( ps, A );
    auto blf = blockform1( ps, F );
    auto u = this->fieldFlux();
    auto p = this->fieldPotential();
    auto phat = this->fieldTrace();
    auto tau_constant = cst(M_tauCst);

    auto const& symbolsExpr = mctx.symbolsExpr();

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& condExpr = this->materialsProperties()->materialProperty( matName, this->M_physicMap.at("condK"));
            auto cond = expr( condExpr.expr(), symbolsExpr);
            bbf(0_c, 0_c) += integrate(_range=range, _expr=inner(idt(u),id(u))/cond);
        }
    }

    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "VolumicForces") )
    {
        auto range = bc.emptyMarkers() ? elements(support(M_Wh)) : markedelements(support(M_Wh), bc.markers());
        auto g = expr(bc.expr(), symbolsExpr);
        blf(1_c) += integrate(_range=range, _expr=inner(id(p),g));
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Dirichlet") )
    {
        bbf(2_c, 2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                   _expr=inner(idt(phat),id(phat)) );
        auto g = expr(bc.expr(), symbolsExpr);
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g));
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Neumann") )
    {
        // <j.n,mu>_Gamma_N
        bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=id(phat)*normalt(u) );
        // <tau p, mu>_Gamma_N
        bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=tau_constant*inner(id(phat), idt(p)) );
        // <-tau phat, mu>_Gamma_N
        bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=-tau_constant*inner(idt(phat), id(phat)) );
        auto g = expr(bc.expr(), symbolsExpr);
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g));
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Robin") )
    {
        // <j.n,mu>_Gamma_R
        bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=id(phat)*normalt(u) );
        // <tau p, mu>_Gamma_R
        bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=tau_constant*inner(id(phat), idt(p)) );
        // <-tau phat, mu>_Gamma_R
        bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=-tau_constant*inner(idt(phat), id(phat)) );

        auto g1 = expr(bc.expr1(), symbolsExpr);
        bbf(2_c, 2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                   _expr=g1*inner(idt(phat),id(phat)) );
        auto g2 = expr(bc.expr2(), symbolsExpr);
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g2));
    }

}

} // namespace FeelModels

} // namespace Feel

#endif
