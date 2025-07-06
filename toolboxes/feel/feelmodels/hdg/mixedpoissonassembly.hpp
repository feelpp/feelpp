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
    auto el_param = is_tensor2symm ? -1.0 : 1.0;

    auto const& symbolsExpr = mctx.symbolsExpr();

    forEachMaterialWithCoefficientExpr<1,1>(
        *this, this->diffusionCoefficientName(), symbolsExpr,
        [&](std::string const& matName, auto const& range, auto const& coeff_c_expr) {
            LOG(INFO) << "diffusion term with coefficient " << " in material " << matName;
            //auto coeff_c_expr = expr<1,1>( coeff_c.expr(), symbolsExpr );

            // 1/c*(j,v)_Omega
            bbf(0_c, 0_c) += integrate(_range=range, _expr=inner(idt(u), id(u))/coeff_c_expr );

            // -(p,div(v))_Omega
            bbf(0_c, 1_c) += integrate(_range=range, _expr=inner(-el_param*idt(p), div(u)) );

            if ( !hasMaterialWithProperty( *this, matName, this->conservativeFluxConvectionCoefficientName() ) )
            {
                // integration by parts
                // <phat,v.n>_Gamma\Gamma_I
                bbf( 0_c, 2_c ) += integrate(_range=faces(support(M_Wh),range, isNotIbcFace ),
                                             _expr=el_param*inner( idt(phat), leftface(normal(u))+rightface(normal(u)) ) );
            } // else see conservative flux convection term below
        });
#if 0  
    forEachMaterialWithCoefficient(
        this->materialsProperties(), symbolsExpr, {this->lameLambdaCoefficientName(), this->lameMuCoefficientName()},
        [&](std::string const& matName, auto const& range, auto const& coeff_lambda_expr, auto const& coeff_mu_expr) {
            auto c1 = cst(0.5) / coeff_mu_expr;
            auto c2 = -coeff_lambda_expr / (cst(2.) * coeff_mu_expr * (nDim * coeff_lambda_expr + cst(2.) * coeff_mu_expr));
            bbf(0_c, 0_c) += integrate(_range=range, _expr=-el_param * c1 * inner(idt(u), id(u)));
            if constexpr (is_tensor2symm) {
                bbf(0_c, 0_c) += integrate(_range=range, _expr=-el_param * c2 * trace(idt(u)) * trace(id(u)));
            }
        });
#endif        
    // convection term
    forEachMaterialWithCoefficientExpr<nDim,1>(
        *this, this->conservativeFluxConvectionCoefficientName(), symbolsExpr,
        [&](std::string const& matName, auto const& range, auto const& coeff_alpha_expr) {
            if constexpr ( is_scalar )
            {
                LOG(INFO) << "convection term with scalar coefficient " << " in material " << matName;
                //auto coeff_alpha_expr = expr( coeff_alpha.template expr<nDim, 1>(), symbolsExpr );

                if ( hasMaterialWithProperty( *this, matName, this->diffusionCoefficientName() ) )
                {
                    auto const& prop_c = this->materialsProperties()->materialProperty(matName, this->diffusionCoefficientName());
                    auto coeff_c_expr = expr(prop_c.expr(), symbolsExpr);

                    bbf(0_c, 1_c) += integrate(_range=range, _expr=inner(coeff_alpha_expr * idt(p), id(u))/coeff_c_expr );

                    auto alpha_N  = inner( coeff_alpha_expr, N());
                    auto inflow = chi(alpha_N <= 0.);
                    auto outflow = chi(alpha_N > 0.);
                    auto inmaterial = chi( emarker() == this->mesh()->markerId( matName ) );
                    // integration by parts
                    // <phat,v.n>_Gamma\Gamma_I

                    bbf( 0_c, 1_c ) += integrate(_range=faces(support(M_Wh),range, isNotIbcFace ),
                                                _expr=el_param*inner( leftfacet(outflow*idt(p))+rightfacet(outflow*idt(p)), leftface(normal(u))+rightface(normal(u)) ) );

                    bbf( 0_c, 2_c ) += integrate(_range=faces(support(M_Wh),range, isNotIbcFace ),
                                                _expr=el_param*(leftface(normal(u))+rightface(normal(u))*(leftfacev(inflow)+rightfacev(inflow))*idt(phat) ) );

                }
            }
        });

#if 1
    

    // (div(j),q)_Omega
    bbf( 1_c, 0_c ) += integrate(_range=elements(support(M_Wh)), _expr=el_param*inner(id(p), divt(u)) );


    // <tau p, w>_Gamma
    bbf( 1_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=tau_constant *
                                 ( inner(leftfacet( idt(p)), leftface(id(p))) +
                                   inner(rightfacet( idt(p)), rightface(id(p))) ) );
    bbf( 1_c, 1_c ) += integrate(_range=boundaryfaces(support(M_Wh)),
                                 _expr=tau_constant * inner(id(p), idt(p)) );

    forEachMaterialWithCoefficientExpr<1,1>(
        *this, this->reactionCoefficientName(), symbolsExpr,
        [&](std::string const& matName, auto const& range, auto const& coeff_a_expr)
        {
            bbf(1_c, 1_c) += integrate(_range = range, _expr = inner(coeff_a_expr * idt(p), id(p)));
        }
    );

    if (!this->isStationary())
    {
        // First-order time derivative term
        forEachMaterialWithCoefficientExpr<1,1>(
            *this, this->firstTimeDerivativeCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_d_expr)
            {
                LOG(INFO) << "First-order time derivative term in material " << matName;
                auto coeff = this->timeStepBdfPotential()->polyDerivCoefficient(0);
                bbf(1_c, 1_c) += integrate(_range=range, _expr=coeff_d_expr * coeff * inner(idt(p), id(p)));
            });

        // Second-order time derivative term
        forEachMaterialWithCoefficientExpr<1,1>(
            *this, this->secondTimeDerivativeCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_d2_expr)
            {
                LOG(INFO) << "Second-order time derivative term in material " << matName;
                auto dt = this->timeStep();
                bbf(1_c, 1_c) += integrate(_range=range, _expr=coeff_d2_expr * inner(idt(p), id(p)) / (dt * dt));
            });
    }
#endif
    // <-tau phat, w>_Gamma\Gamma_I
    bbf( 1_c, 2_c ) += integrate(_range=faces(support(M_Wh), isNotIbcFace ),
                                 _expr=-tau_constant * inner(idt(phat),
                                                             leftface( id(p) )+
                                                             rightface( id(p) ) ) );
    // <j.n,mu>_Omega/Gamma
    bbf( 2_c, 0_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=inner( id(phat), leftfacet(normalt(u))+rightfacet(normalt(u)) ) );

    auto tau_D = tau_constant/h();
    // <tau p, mu>_Omega/Gamma
    bbf( 2_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=el_param*tau_D * inner(id(phat),
                                                            leftfacet( idt(p) )+
                                                            rightfacet( idt(p) )) );

    // <-tau phat, mu>_Omega/Gamma
    bbf( 2_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                 _expr=-sc_param*el_param*tau_D * inner(idt(phat), id(phat) ) );

    forEachMaterialWithCoefficientExpr<nDim,1>(
        *this, this->conservativeFluxConvectionCoefficientName(), symbolsExpr,
        [&](std::string const& matName, auto const& range, auto const& coeff_alpha_expr)
        {
            if constexpr (is_scalar)
            {
                LOG(INFO) << "convection term stabilization in material " << matName;

                auto alpha_N = inner(coeff_alpha_expr, N());
                auto tau_C   = max(alpha_N, cst(0.));
                auto inmaterial = chi( emarker() == this->mesh()->markerId( matName ) );
                // <tau * p, phat> on interfaces
                bbf(2_c, 1_c) += integrate(
                    _range = faces(support(M_Wh),range),
                    _expr  = el_param * inner(id(phat),
                            leftfacet((inmaterial *tau_C * idt(p)) + rightfacet( inmaterial*tau_C * idt(p))) ) );

                // <-tau * phat, mu> on interfaces (symmetric or skew-symmetric form)
                bbf(2_c, 2_c) += integrate(
                    _range = faces(support(M_Wh),range),
                    _expr  = -sc_param * el_param * inner(idt(phat), tau_C * id(phat)) );
            }
        }
    );
#if 0
    if constexpr ( is_scalar )
    {
        for ( auto const& [physicName, physicData] : this->physicsFromCurrentType() )
        {
            for ( std::string const& matName : this->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto const& range = this->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(), matName );
                if ( this->materialsProperties()->hasProperty( matName, this->convectionCoefficientName() ) )
                {
                    auto coeff_alpha = this->materialsProperties()->materialProperty( matName, this->convectionCoefficientName() );
                    auto coeff_alpha_expr = expr( coeff_alpha.template expr<nDim, 1>(), symbolsExpr );
                    auto tau_C = max(trans(coeff_alpha_expr)*N(), cst(0.));
                    // <tau p, mu>_Omega/Gamma
                    bbf( 2_c, 1_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                                 _expr=el_param* tau_C * inner(id(phat),
                                                                            leftfacet( idt(p) )+
                                                                            rightfacet( idt(p) )) );
                    // <-tau phat, mu>_Omega/Gamma
                    bbf( 2_c, 2_c ) += integrate(_range=internalfaces(support(M_Wh)),
                                                 _expr=-sc_param*el_param* tau_C * inner(idt(phat), id(phat) ) );
                }
            }
        }
    }
#endif    

    if (!this->isStationary())
    {
        // First-order time derivative RHS term
        forEachMaterialWithCoefficientExpr<1,1>(
            *this, this->firstTimeDerivativeCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_d_expr)
            {
                auto polyDeriv = this->timeStepBdfPotential()->polyDeriv();
                blf(1_c) += integrate(_range=range,
                                    _expr=coeff_d_expr * inner(id(p), idv(polyDeriv)));
            });

        // Second-order time derivative RHS term (if enabled later)
        forEachMaterialWithCoefficientExpr<1,1>(
            *this, this->secondTimeDerivativeCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_d2_expr)
            {
                auto dt = this->timeStep();
                // Uncomment when needed
                // blf(1_c) += integrate(_range=range,
                //                       _expr=coeff_d2_expr * inner(id(p), idv(...)) / (dt * dt));
            });
    }

    // needed for static condensation
    blf(1_c) += integrate(_range=M_rangeMeshElements, _expr=cst(0.));

    forEachMaterialWithCoefficientExpr<is_scalar ? 1 : nDim, 1>(
        *this, this->sourceCoefficientName(), symbolsExpr,
        [&](auto const& matName, auto const& range, auto const& coeff_f_expr)
        {
            LOG(INFO) << "source term in material " << matName;
            blf(1_c) += integrate(_range = range,
                                _expr  = el_param * inner(id(p), coeff_f_expr));
        }
    );


    for ( auto const& [bcName,bcData] : M_boundaryConditions->dirichlet() )
    {
        auto bcRangeFaces = markedfaces(support(M_Wh), bcData->markers());
        bbf(2_c, 2_c) += integrate(_range=bcRangeFaces,
                                   _expr=inner(idt(phat),id(phat)) );
        auto g = bcData->expr( symbolsExpr );
        blf(2_c) += integrate(_range=bcRangeFaces,
                              _expr=inner(id(phat),g));
    }
    for (auto const& [bcName, bcData] : M_boundaryConditions->neumann())
    {
        LOG(INFO) << "neumann boundary condition " << bcName << "\n";
        auto tau_D = tau_constant / h();
        auto bcRangeFaces = markedfaces(support(M_Wh), bcData->markers());

        // Classical Neumann flux terms
        bbf(2_c, 0_c) += integrate(_range = bcRangeFaces,
                                _expr  = inner(id(phat), normalt(u)));

        bbf(2_c, 1_c) += integrate(_range = bcRangeFaces,
                                _expr  = el_param * tau_D * inner(id(phat), idt(p)));

        bbf(2_c, 2_c) += integrate(_range = bcRangeFaces,
                                _expr  = -el_param * tau_D * inner(idt(phat), id(phat)));

        // Optional convection stabilization on Gamma_N
        forEachMaterialWithCoefficientExpr<nDim,1>(
            *this, this->conservativeFluxConvectionCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_alpha_expr)
            {
                if constexpr (is_scalar)
                {
                    if (!intersectionIsEmpty(range, bcRangeFaces))
                    {
                        LOG(INFO) << "[hdg] convection Neumann stabilization for material " << matName;

                        auto alpha_N = inner(coeff_alpha_expr, N());
                        auto tau_C = max(alpha_N, cst(0.));

                        bbf(2_c, 1_c) += integrate(_range = bcRangeFaces,
                                                _expr = el_param * tau_C * inner(id(phat), idt(p)));

                        bbf(2_c, 2_c) += integrate(_range = bcRangeFaces,
                                                _expr = -el_param * tau_C * inner(idt(phat), id(phat)));
                    }
                }
            });

        // Neumann source term
        auto g = bcData->expr(symbolsExpr);
        blf(2_c) += integrate(_range = bcRangeFaces,
                              _expr  = inner(id(phat), g));
    }
   for (auto const& [bcName, bcData] : M_boundaryConditions->robin())
    {
        auto tau_D = tau_constant / h();
        auto bcRangeFaces = markedfaces(support(M_Wh), bcData->markers());

        // Classical Robin flux terms
        bbf(2_c, 0_c) += integrate(_range = bcRangeFaces,
                                _expr  = inner(id(phat), normalt(u)));

        bbf(2_c, 1_c) += integrate(_range = bcRangeFaces,
                                _expr  = el_param * tau_D * inner(id(phat), idt(p)));

        bbf(2_c, 2_c) += integrate(_range = bcRangeFaces,
                                _expr  = -el_param * tau_D * inner(idt(phat), id(phat)));

        // Convection stabilization on Robin boundary (filtered by intersection)
        forEachMaterialWithCoefficientExpr<nDim,1>(
            *this, this->conservativeFluxConvectionCoefficientName(), symbolsExpr,
            [&](std::string const& matName, auto const& range, auto const& coeff_alpha_expr)
            {
                if constexpr (is_scalar)
                {
                    if (!intersectionIsEmpty(range, bcRangeFaces))
                    {
                        LOG(INFO) << "[hdg] convection Robin stabilization for material " << matName;

                        auto alpha_N = inner(coeff_alpha_expr, N());
                        auto tau_C   = max(alpha_N, cst(0.));

                        bbf(2_c, 1_c) += integrate(_range = bcRangeFaces,
                                                _expr  = el_param * tau_C * inner(id(phat), idt(p)));

                        bbf(2_c, 2_c) += integrate(_range = bcRangeFaces,
                                                _expr  = -el_param * tau_C * inner(idt(phat), id(phat)));
                    }
                }
            });

        // Robin coefficient terms
        auto g1 = bcData->expr1(symbolsExpr);
        bbf(2_c, 2_c) += integrate(_range = bcRangeFaces,
                                _expr  = g1 * inner(idt(phat), id(phat)));

        auto g2 = bcData->expr2(symbolsExpr);
        blf(2_c) += integrate(_range = bcRangeFaces,
                            _expr  = inner(id(phat), g2));
    }
    int i = 0;
    for ( auto const& [bcName,bcData] : M_boundaryConditions->integral() )
    {
        auto bcRangeFaces = markedfaces(support(M_Wh), bcData->markers());
        // <lambda, v.n>_Gamma_I
        bbf( 0_c, 3_c, 0, i ) += integrate( _range=bcRangeFaces,
                                            _expr= inner(idt(l), normal(u)) );

        // -<lambda, tau w>_Gamma_I
        bbf( 1_c, 3_c, 0, i ) += integrate( _range=bcRangeFaces,
                                            _expr=-tau_constant*inner(idt(l), id(p)) );

        // <j.n, m>_Gamma_I
        bbf( 3_c, 0_c, i, 0 ) += integrate( _range=bcRangeFaces,
                                            _expr=inner(id(l), normalt(u)) );

        // <tau p, m>_Gamma_I
        bbf( 3_c, 1_c, i, 0 ) += integrate( _range=bcRangeFaces,
                                            _expr=tau_constant*inner(idt(p), id(l)) );

        // -<lambda2, m>_Gamma_I
        bbf( 3_c, 3_c, i, i ) += integrate( _range=bcRangeFaces,
                                            _expr=-tau_constant*inner(id(l), idt(l)) );

        double meas = integrate( _range=bcRangeFaces,
                                 _expr=cst(1.)).evaluate()(0,0);
        auto g = bcData->expr( symbolsExpr );
        blf(3_c, i) += integrate( _range=bcRangeFaces,
                                  _expr=inner(g,id(l))/meas);
        i++;
    }
#if 0
    //for( auto& [name, bc] : this->modelProperties().boundaryConditions2().byFieldType( M_fluxKey, "Interface") )
    for ( auto const& [bcName,bcData] : M_boundaryConditions->interface() )
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
#endif
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
    forEachMaterialWithCoefficientExpr<1,1>(
        *this, this->diffusionCoefficientName(), symbolsExpr,
        [&](std::string const& matName, auto const& range, auto const& coeff_c_expr)
        {
            if constexpr (is_scalar)
            {
                blf(0_c) += integrate(
                    _range = range,
                    _expr  = -grad(pp) * idv(u) / coeff_c_expr );
            }
            else
            {
                blf(0_c) += integrate(
                    _range = range,
                    _expr  = -inner(grad(pp), idv(u)) / coeff_c_expr );
            }
        }
    );

}

} // namespace FeelModels

} // namespace Feel

#endif
