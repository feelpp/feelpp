#ifndef FEELPP_TOOLBOXES_MIXEDPOISSON_ASSEMBLY_HPP
#define FEELPP_TOOLBOXES_MIXEDPOISSON_ASSEMBLY_HPP 1

namespace Feel
{
namespace FeelModels
{

template<typename ConvexType, int Order, template<uint16_type> class PolySetType, int E_Order>
template< typename ModelConvexType>
void
MixedPoisson<ConvexType, Order, PolySetType, E_Order>::updateLinearPDE( DataUpdateLinear & data, ModelConvexType const& mctx ) const
{
    auto A = std::dynamic_pointer_cast<condensed_matrix_t<value_type>>(data.matrix());
    auto F = std::dynamic_pointer_cast<condensed_vector_t<value_type>>(data.rhs());
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;
    if( buildNonCstPart )
        return;

    auto mesh = this->mesh();
    auto ps = this->spaceProduct();
    auto bbf = blockform2( ps, A );
    auto blf = blockform1( ps, F );
    auto u = this->fieldFlux();
    auto p = this->fieldPotential();
    auto phat = this->fieldTrace();
    auto l = this->M_Ch->element();
    auto tau_constant = cst(M_tauCst);
    auto sc_param = M_useSC ? 0.5 : 1.0;

    auto const& symbolsExpr = mctx.symbolsExpr();

    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            if( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
            {
                auto coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
                auto coeff_c_expr = expr( coeff_c.expr(), symbolsExpr );
                bbf(0_c, 0_c) += integrate(_range=range, _expr=inner(idt(u),id(u))/coeff_c_expr);
            }
            if( this->materialsProperties()->hasProperty( matName, this->lameLambdaCoefficientName() ) && this->materialsProperties()->hasProperty( matName, this->lameMuCoefficientName() ) )
            {
                auto coeff_lambda = this->materialsProperties()->materialProperty( matName, this->lameLambdaCoefficientName() );
                auto coeff_lambda_expr = expr( coeff_lambda.expr(), symbolsExpr );
                auto coeff_mu = this->materialsProperties()->materialProperty( matName, this->lameMuCoefficientName() );
                auto coeff_mu_expr = expr( coeff_mu.expr(), symbolsExpr );
                auto c1 = cst(0.5)/coeff_mu_expr;
                auto c2 = -coeff_lambda_expr/(cst(2.) * coeff_mu_expr * (nDim*coeff_lambda_expr + cst(2.)*coeff_mu_expr));
                bbf( 0_c, 0_c ) += integrate(_range=range, _expr=-c1*inner(idt(u),id(u)) );
                if constexpr( is_tensor2symm ) {
                    bbf( 0_c, 0_c ) += integrate(_range=range, _expr=-c2*trace(idt(u))*trace(id(u)) );
                }
            }
        }
    }

    // -(p,div(v))_Omega
    bbf( 0_c, 1_c ) += integrate(_range=elements(support(M_Wh)), _expr=-inner(idt(p),div(u)) );

    // <phat,v.n>_Gamma\Gamma_I
    bbf( 0_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=inner( idt(phat), leftface(normal(u))+rightface(normal(u)) ) );
    bbf( 0_c, 2_c ) += integrate(_range=M_gammaMinusIntegral,
                                 _expr=inner(idt(phat), normal(u)) );

    // (div(j),q)_Omega
    bbf( 1_c, 0_c ) += integrate(_range=elements(support(M_Wh)), _expr=inner(id(p), divt(u)) );


    // <tau p, w>_Gamma
    bbf( 1_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=tau_constant *
                                 ( inner(leftfacet( idt(p)), leftface(id(p))) +
                                   inner(rightfacet( idt(p)), rightface(id(p))) ) );
    bbf( 1_c, 1_c ) += integrate(_range=boundaryfaces(support(M_Wh)),
                                 _expr=tau_constant * inner(id(p), idt(p)) );
    if( !this->isStationary() )
    {
        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                if( this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
                {
                    auto coeff = this->timeStepBdfPotential()->polyDerivCoefficient(0);
                    auto coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
                    auto coeff_d_expr = expr(coeff_d.expr(), symbolsExpr);
                    // (1/delta_t p, w)_Omega  [only if it is not stationary]
                    bbf( 1_c, 1_c) += integrate(_range=range,
                                                _expr=coeff_d_expr*coeff*inner(idt(p), id(p)) );
                }
                if( this->materialsProperties()->hasProperty( matName, this->secondTimeDerivativeCoefficientName() ) )
                {
                    auto dt = this->timeStep();
                    auto coeff_d2 = this->materialsProperties()->materialProperty( matName, this->secondTimeDerivativeCoefficientName() );
                    auto coeff_d2_expr = expr(coeff_d2.expr(), symbolsExpr);
                    bbf( 1_c, 1_c) += integrate(_range=range,
                                                _expr=coeff_d2_expr*inner(idt(p), id(p))/(dt*dt) );
                }
            }
        }
    }

    // <-tau phat, w>_Gamma\Gamma_I
    bbf( 1_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=-tau_constant * inner(idt(phat),
                                                             leftface( id(p) )+
                                                             rightface( id(p) ) ) );
    bbf( 1_c, 2_c ) += integrate(_range=M_gammaMinusIntegral,
                                 _expr=-tau_constant * inner(idt(phat), id(p)) );


    // <j.n,mu>_Omega/Gamma
    bbf( 2_c, 0_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=inner( id(phat), leftfacet(normalt(u))+rightfacet(normalt(u)) ) );

    // <tau p, mu>_Omega/Gamma
    bbf( 2_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=tau_constant * inner(id(phat),
                                                            leftfacet( idt(p) )+
                                                            rightfacet( idt(p) )) );

    // <-tau phat, mu>_Omega/Gamma
    bbf( 2_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=-sc_param*tau_constant * inner(idt(phat), id(phat) ) );

    if( !this->isStationary() )
    {
        for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                if( this->materialsProperties()->hasProperty( matName, this->firstTimeDerivativeCoefficientName() ) )
                {
                    auto polyDeriv = this->timeStepBdfPotential()->polyDeriv();
                    auto coeff_d = this->materialsProperties()->materialProperty( matName, this->firstTimeDerivativeCoefficientName() );
                    auto coeff_d_expr = expr(coeff_d.expr(), symbolsExpr);
                    blf(1_c) += integrate(_range=range,
                                          _expr=coeff_d_expr*inner(id(p), idv(polyDeriv)) );
                }
                if( this->materialsProperties()->hasProperty( matName, this->secondTimeDerivativeCoefficientName() ) )
                {
                    auto dt = this->timeStep();
                    auto coeff_d2 = this->materialsProperties()->materialProperty( matName, this->secondTimeDerivativeCoefficientName() );
                    auto coeff_d2_expr = expr(coeff_d2.expr(), symbolsExpr);
                    // blf(1_c) += integrate(_range=range,
                    //                       _expr=coeff_d2_expr*inner(idt(p), id(p))/(dt*dt) );
                }
            }
        }
    }

    // needed for static condensation
    blf(1_c) += integrate(_range=M_rangeMeshElements, _expr=cst(0.));
    for ( auto const& [physicName,physicData] : this->physicsFromCurrentType() )
    {
        for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
            if( this->materialsProperties()->hasProperty( matName, this->sourceCoefficientName() ) )
            {
                auto coeff_f = this->materialsProperties()->materialProperty( matName, this->sourceCoefficientName() );
                auto coeff_f_expr = [&coeff_f,&symbolsExpr]() -> decltype(auto) {
                                        if constexpr( is_scalar ) {
                                            return expr(coeff_f.expr(), symbolsExpr);
                                        } else {
                                            return expr(coeff_f.template expr<nDim>(), symbolsExpr);
                                        }
                                    }();
                blf(1_c) += integrate(_range=range, _expr=inner(id(p), coeff_f_expr));
            }
        }
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Dirichlet") )
    {
        bbf(2_c, 2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                   _expr=inner(idt(phat),id(phat)) );
        auto g = [&bc = bc,&symbolsExpr]() -> decltype(auto) {
                     if constexpr( is_scalar ) {
                         return expr(bc.expr(), symbolsExpr);
                     } else {
                         return expr(bc.template expr<nDim>(), symbolsExpr);
                     }
                 }();
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g));
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Neumann") )
    {
        // <j.n,mu>_Gamma_N
        bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=inner(id(phat), normalt(u)) );
        // <tau p, mu>_Gamma_N
        bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=tau_constant*inner(id(phat), idt(p)) );
        // <-tau phat, mu>_Gamma_N
        bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=-tau_constant*inner(idt(phat), id(phat)) );
        auto g = [&bc = bc,&symbolsExpr]() -> decltype(auto) {
                     if constexpr( is_scalar ) {
                         return expr(bc.expr(), symbolsExpr);
                     } else {
                         return expr(bc.template expr<nDim>(), symbolsExpr);
                     }
                 }();
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g));
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_potentialKey, "Robin") )
    {
        // <j.n,mu>_Gamma_R
        bbf( 2_c, 0_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=inner(id(phat), normalt(u)) );
        // <tau p, mu>_Gamma_R
        bbf( 2_c, 1_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=tau_constant*inner(id(phat), idt(p)) );
        // <-tau phat, mu>_Gamma_R
        bbf( 2_c, 2_c ) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                     _expr=-tau_constant*inner(idt(phat), id(phat)) );

        auto g1 = expr(bc.expr1(), symbolsExpr);
        bbf(2_c, 2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                                   _expr=g1*inner(idt(phat),id(phat)) );
        auto g2 = [&bc = bc,&symbolsExpr]() -> decltype(auto) {
                      if constexpr( is_scalar ) {
                          return expr(bc.expr2(), symbolsExpr);
                      } else {
                          return expr(bc.template expr2<nDim>(), symbolsExpr);
                      }
                 }();
        blf(2_c) += integrate(_range=markedfaces(support(M_Wh), bc.markers()),
                              _expr=inner(id(phat),g2));
    }
    int i = 0;
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_fluxKey, "Integral") )
    {
        // <lambda, v.n>_Gamma_I
        bbf( 0_c, 3_c, 0, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr= inner(idt(l), normal(u)) );

        // -<lambda, tau w>_Gamma_I
        bbf( 1_c, 3_c, 0, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=-tau_constant*inner(idt(l), id(p)) );

        // <j.n, m>_Gamma_I
        bbf( 3_c, 0_c, i, 0 ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=inner(id(l), normalt(u)) );

        // <tau p, m>_Gamma_I
        bbf( 3_c, 1_c, i, 0 ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=tau_constant*inner(idt(p), id(l)) );

        // -<lambda2, m>_Gamma_I
        bbf( 3_c, 3_c, i, i ) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                            _expr=-tau_constant*inner(id(l), idt(l)) );

        double meas = integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                 _expr=cst(1.)).evaluate()(0,0);
        auto g = [&bc = bc,&symbolsExpr]() -> decltype(auto) {
                     if constexpr( is_scalar ) {
                         return expr(bc.expr(), symbolsExpr);
                     } else {
                         return expr(bc.template expr<nDim>(), symbolsExpr);
                     }
                 }();
        blf(3_c, i) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                                  _expr=inner(g,id(l))/meas);
        i++;
    }
    for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_fluxKey, "Interface") )
    {
        auto g = [&bc = bc,&symbolsExpr]() -> decltype(auto) {
                     if constexpr( is_scalar ) {
                         return expr(bc.expr(), symbolsExpr);
                     } else {
                         return expr(bc.template expr<nDim>(), symbolsExpr);
                     }
                 }();
        blf(2_c) += integrate( _range=markedfaces(support(M_Wh), bc.markers()),
                               _expr=inner(id(phat),g) );
    }
}

template<typename ConvexType, int Order, template<uint16_type> class PolySetType, int E_Order>
template< typename ModelConvexType>
void
MixedPoisson<ConvexType, Order, PolySetType, E_Order>::updatePostPDE( DataUpdateLinear & data, ModelConvexType const& mctx ) const
{
    auto A = std::dynamic_pointer_cast<condensed_matrix_t<value_type>>(data.matrix());
    auto F = std::dynamic_pointer_cast<condensed_vector_t<value_type>>(data.rhs());
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
            if( this->materialsProperties()->hasProperty( matName, this->diffusionCoefficientName() ) )
            {
                auto coeff_c = this->materialsProperties()->materialProperty( matName, this->diffusionCoefficientName() );
                auto coeff_c_expr = expr( coeff_c.expr(), symbolsExpr );
                if constexpr( is_scalar ) {
                    blf(0_c) = integrate(_range=range, _expr=-grad(pp)*idv(u)/coeff_c_expr );
                } else {
                    blf(0_c) = integrate(_range=range, _expr=-inner(grad(pp), idv(u))/coeff_c_expr );
                }
            }
        }
    }

}

} // namespace FeelModels

} // namespace Feel

#endif
