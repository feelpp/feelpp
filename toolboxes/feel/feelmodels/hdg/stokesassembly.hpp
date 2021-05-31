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
    auto l = this->M_Ch->element();
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

    if( !this->isStationary() )
    {
        auto polyDeriv = this->timeStepBdfPotential()->polyDeriv();
        blf(1_c) += integrate(_range=M_rangeMeshElements, _expr=inner(id(p), idv(polyDeriv)) );
    }

    // needed for static condensation
    blf(1_c) += integrate(_range=M_rangeMeshElements, _expr=cst(0.));
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
    int i = 0;
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_fluxKey, "Integral") )
    {
        // <lambda, v.n>_Gamma_I
        bbf( 0_c, 3_c, 0, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr= idt(l) * normal(u) );

        // -<lambda, tau w>_Gamma_I
        bbf( 1_c, 3_c, 0, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=-tau_constant*inner(idt(l), id(p)) );

        // <j.n, m>_Gamma_I
        bbf( 3_c, 0_c, i, 0 ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=normalt(u) * id(l) );

        // <tau p, m>_Gamma_I
        bbf( 3_c, 1_c, i, 0 ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=tau_constant*inner(idt(p), id(l)) );

        // -<lambda2, m>_Gamma_I
        bbf( 3_c, 3_c, i, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=-tau_constant*inner(id(l), idt(l)) );

        auto g = expr(bc.expr(), symbolsExpr);
        double meas = integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                 _expr=cst(1.)).evaluate()(0,0);
        blf(3_c, i) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                  _expr=inner(g,id(l))/meas);
        i++;
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_fluxKey, "Interface") )
    {
        auto g = expr(bc.expr(), symbolsExpr);
        blf(2_c) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                               _expr=inner(id(phat),g) );
    }
}

template< typename ConvexType, int Dim, int E_Order>
template< typename ModelConvexType>
void
MixedPoisson<ConvexType, Dim, E_Order>::updatePostPDE( DataUpdateHDG & data, ModelConvexType const& mctx ) const
{
    condensed_matrix_ptr_t<value_type>& A = data.matrix();
    condensed_vector_ptr_t<value_type>& F = data.rhs();
    bool buildCstPart = data.buildCstPart();
    auto mesh = this->mesh();
    auto ps = product(M_Whp);
    auto bbf = blockform2( ps, A );
    auto blf = blockform1( ps, F );
    auto pp = this->fieldPostPotential();
    auto u = this->fieldFlux();

    auto const& symbolsExpr = mctx.symbolsExpr();

    if( buildCstPart )
    {
        bbf(0_c, 0_c) = integrate( _range=elements(support(M_Wh)),
                                   _expr=inner(gradt(pp), grad(pp)) );
        M_postMatrixInit = true;
    }
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            auto const& condExpr = this->materialsProperties()->materialProperty( matName, this->M_physicMap.at("condK"));
            auto cond = expr( condExpr.expr(), symbolsExpr);
            blf(0_c) = integrate(_range=range, _expr=-grad(pp)*idv(u)/cond );
        }
    }

}

} // namespace FeelModels

} // namespace Feel

#endif
