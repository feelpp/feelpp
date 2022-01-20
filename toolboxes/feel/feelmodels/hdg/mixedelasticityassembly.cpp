#include <feel/feelmodels/hdg/mixedelasticity.hpp>

namespace Feel
{
namespace FeelModels
{

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assembleMatrixIBC( int i , std::string markerOpt)
{


    auto bbf = blockform2( *M_ps, M_A_cst);

    auto v = M_Vh->element( "v" );
    auto u = M_Wh->element( "u" );
    auto w = M_Wh->element( "w" );
    auto nu = M_Ch->element( "nu" );
    auto uI = M_Ch->element( "uI" );


    auto H = M_M0h->element( "H" );

    if ( M_hFace == 0 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( M_hFace == 1 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( M_hFace == 2 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
        H.on( _range=elements(M_M0h->mesh()), _expr=h() );

    // stabilisation parameter
    auto tau_constant = cst(M_tauCst);

    std::string marker;
    if ( !markerOpt.empty())
    {
        marker = markerOpt;
    }
    else
    {
        auto exAtMarker = M_IBCList[i];
        marker = exAtMarker.marker();
        Feel::cout << "Integral on: " << marker << std::endl;
    }


    // <lambda, v.n>_Gamma_I
    bbf( 0_c, 3_c, 0, i) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh, marker), _expr=-trans(idt(uI))*(id(v)*N()) );

    // <lambda, tau w>_Gamma_I
    bbf( 1_c, 3_c, 1, i ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr= tau_constant * trans(idt(uI)) * pow(idv(H),M_tauOrder)*id(w) );

    // <sigma.n, m>_Gamma_I
    bbf( 3_c, 0_c, i, 0 ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr= inner(idt(v)*N(),id(nu)) );


    // <tau u, m>_Gamma_I
    bbf( 3_c, 1_c, i, 1 ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                       _expr= tau_constant * pow(idv(H),M_tauOrder)* inner(idt(u),id(nu)) ),

    // -<lambda2, m>_Gamma_I
    bbf( 3_c, 3_c, i, i ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                       _expr=-tau_constant * pow(idv(H),M_tauOrder) * inner(idt(uI),id(nu)) );



}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assembleRhsIBC( int i, std::string markerOpt, double intjn )
{
    auto blf = blockform1( *M_ps, M_F );
    auto nu = M_Ch->element( "nu" );

    auto exAtMarker = M_IBCList[i];
    auto marker = exAtMarker.marker();

    auto g = expr<Dim,1,expr_order>(exAtMarker.expression());
    if ( !this->isStationary() )
        g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );

    Feel::cout << "IBC condition: " << g << std::endl;

    double meas = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=cst(1.0)).evaluate()(0,0);
    Feel::cout << "Measure of the ibc: " << meas << std::endl;

    // <F_target,m>_Gamma_I
    blf(3_c,i) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=inner(g,id(nu))/meas);

}




MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::initExporter( mesh_ptrtype meshVisu )
{
    std::string geoExportType="static"; //change_coords_only, change, static
    M_exporter = exporter ( _mesh=meshVisu?meshVisu:this->mesh(),
                            _name="Export",
                            _geo=geoExportType,
                            _path=this->exporterPath() );
}


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assemble()
{

    tic();
    solve::strategy s = M_useSC ? solve::strategy::static_condensation : solve::strategy::monolithic;

    auto U = M_ps -> element();
    M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps), *M_backend ); //M_backend->newBlockMatrix(_block=csrGraphBlocks(ps));
    //M_A_cst = makeSharedMatrixCondensed<value_type>(s, csrGraphBlocks(*M_ps, (s==solve::strategy::static_condensation)?Pattern::COUPLED:pattern), *M_backend, (s==solve::strategy::static_condensation)?false:true);

    M_F = makeSharedVectorCondensed<value_type>(s, blockVector(*M_ps), *M_backend, false);//M_backend->newBlockVector(_block=blockVector(ps), _copy_values=false);
    //    M_A_cst = M_backend->newBlockMatrix(_block=csrGraphBlocks(*M_ps));
    //M_F = M_backend->newBlockVector(_block=blockVector(*M_ps), _copy_values=false);
    toc("creating matrices and vectors", this->verbose() || FLAGS_v > 0);


}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assembleCst()
{
    // Assembling standard matrix
    tic();
    M_A_cst->zero();
    this->assembleSTD();
    M_timers["asbStd"].push_back(toc("assembleStandardMatrix", this->verbose() || FLAGS_v > 0));

    // Assembling ibc part
    tic();
    for ( int i = 0; i < M_IBCList.size(); i++ )
        this->assembleMatrixIBC( i );
    M_timers["asbIbc"].push_back(toc("assembleIbcMatrix", this->verbose() || FLAGS_v > 0));

    M_A_cst->close();
}


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assembleNonCst()
{
    tic();
    M_F->zero();
    this->assembleF( );

    for ( int i = 0; i < M_IBCList.size(); i++ )
        this->assembleRhsIBC( i );
    M_F->close();
    M_timers["asbRHS"].push_back(toc("assembleRHS", this->verbose() || FLAGS_v > 0));

}


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::solve()
{

    auto U = M_ps -> element();
    auto bbf = blockform2(*M_ps, M_A_cst);

    auto blf = blockform1(*M_ps, M_F);


    std::shared_ptr<NullSpace<double> > myNullSpace( new NullSpace<double>(get_backend(),hdgNullSpace(M_Wh,mpl::int_<Dim>())) );
    get_backend()->attachNearNullSpace( myNullSpace );
    if ( M_nullspace )
        get_backend()->attachNearNullSpace( myNullSpace );


    std::string solver_string = "MixedElasticity : ";
    if( M_useSC )
        solver_string += "static condensation";
    else
        solver_string += "monolithic";

    tic();
    tic();
    bbf.solve(_solution=U, _rhs=blf, _condense=M_useSC, _name= this->prefix());
    M_timers["solver"].push_back(toc("solver", this->verbose() || FLAGS_v > 0));
    toc(solver_string, this->verbose() || FLAGS_v > 0);

    M_up = U(0_c);
    M_pp = U(1_c);

    if ( VLOG_IS_ON(2) )
    {
        Feel::cout << "u_hat=" << U(2_c) << std::endl;
        Feel::cout << "u=" << U(1_c) << std::endl;
        Feel::cout << "sigma=" << U(0_c) << std::endl;
    }

    for( int i = 0; i < M_integralCondition; i++ )
        M_mup.push_back(U(3_c,i));

}


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::assembleSTD()
{
    auto tau_constant = cst(M_tauCst);

    auto sigma = M_Vh->element( "sigma" );
    auto v     = M_Vh->element( "v" );
    auto u     = M_Wh->element( "u" );
    auto w     = M_Wh->element( "w" );
    auto uhat  = M_Mh->element( "uhat" );
    auto m     = M_Mh->element( "m" );
    auto H     = M_M0h->element( "H" );

    auto gammaMinusIntegral = complement(boundaryfaces(M_mesh),[this]( auto const& e ) {
                                                                   for( auto exAtMarker : this->M_IBCList)
                                                                   {
                                                                       if ( e.marker().value() == this->M_mesh->markerName( exAtMarker.marker() ) )
                                                                           return true;
                                                                   }
                                                                   return false; });



    if ( M_hFace == 0 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMax()) );
    else if ( M_hFace == 1 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hMin()) );
    else if ( M_hFace == 2 )
        H.on( _range=elements(M_M0h->mesh()), _expr=cst(M_Vh->mesh()->hAverage()) );
    else
        H.on( _range=elements(M_M0h->mesh()), _expr=h() );

    auto sc_param = M_useSC ? 0.5 : 1.0;

    auto bbf = blockform2 ( *M_ps, M_A_cst );


    for( auto const& pairMat : modelProperties().materials() )
    {
        auto material = pairMat.second;
        auto lambda = material.property("lambda").exprScalar();
        Feel::cout << "Lambda: " << lambda << std::endl;
        auto mu = material.property("mu").exprScalar();
        Feel::cout << "Mu: " << mu << std::endl;
        auto c1 = cst(0.5)/mu;
        auto c2 = -lambda/(cst(2.) * mu * (cst(Dim)*lambda + cst(2.)*mu));
        Feel::cout << "c1: " << mean(_range=elements(M_mesh),_expr=c1) << std::endl;
        Feel::cout << "c2: " << mean(_range=elements(M_mesh),_expr=c2) << std::endl;

        bbf( 0_c, 0_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(c1*inner(idt(sigma),id(v))) );
        bbf( 0_c, 0_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(c2*trace(idt(sigma))*trace(id(v))) );
    }


    bbf( 0_c, 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),_expr=(trans(idt(u))*div(v)));

    bbf( 0_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
                                _expr=-( trans(idt(uhat))*leftface(id(v)*N())+
                                         trans(idt(uhat))*rightface(id(v)*N())) );
    bbf( 0_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=gammaMinusIntegral,
                                _expr=-trans(idt(uhat))*(id(v)*N()) );

    bbf( 1_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
                                _expr=(trans(id(w))*divt(sigma)));

    // ( d^2u/dt^2, w)_Omega  [only if it is not stationary]
    if ( !this->isStationary() ) {
        auto dt = this->timeStep();
        for( auto const& pairMat : modelProperties().materials() )
        {
            auto material = pairMat.second;
            auto rho = material.property("rho").exprScalar();
            // bbf( 1_c, 1_c ) += integrate(_range=elements(M_mesh),
            //                              _expr = this->timeStepNM()->polySecondDerivCoefficient()*rho*inner(idt(u),id(w)) );
            bbf( 1_c, 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
                                         _expr = -rho*inner(idt(u),id(w))/(dt*dt) );
        }
    }

    // begin dp: here we need to put the projection of u on the faces
    bbf( 1_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),_expr=-tau_constant *
                                ( leftfacet( pow(idv(H),M_tauOrder)*trans(idt(u)))*leftface(id(w)) +
                                  rightfacet( pow(idv(H),M_tauOrder)*trans(idt(u)))*rightface(id(w) )));

    bbf( 1_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=boundaryfaces(M_mesh),
                                _expr=-(tau_constant * pow(idv(H),M_tauOrder)*trans(idt(u))*id(w)));

    bbf( 1_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
                                _expr=tau_constant *
                                ( leftfacet(trans(idt(uhat)))*leftface( pow(idv(H),M_tauOrder)*id(w))+
                                  rightfacet(trans(idt(uhat)))*rightface( pow(idv(H),M_tauOrder)*id(w) )));

    bbf( 1_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=gammaMinusIntegral,
                                _expr=tau_constant * trans(idt(uhat)) * pow(idv(H),M_tauOrder)*id(w) );


    bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
                                _expr=( trans(id(m))*(leftfacet(idt(sigma)*N())+
                                                      rightfacet(idt(sigma)*N())) ) );


    // BC
    bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
                                _expr=-tau_constant * trans(id(m)) * (leftfacet( pow(idv(H),M_tauOrder)*idt(u) )+
                                                                      rightfacet( pow(idv(H),M_tauOrder)*idt(u) )));

    bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=internalfaces(M_mesh),
                                _expr=sc_param*tau_constant * trans(idt(uhat)) * id(m) * ( leftface( pow(idv(H),M_tauOrder) )+
                                                                                           rightface( pow(idv(H),M_tauOrder) )));

    auto itField = modelProperties().boundaryConditions().find( "displacement");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Dirichlet" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                Feel::cout << "Dirichlet on " << marker << std::endl;
                bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=trans(idt(uhat)) * id(m) );
            }
        }

    }

    itField = modelProperties().boundaryConditions().find( "stress");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find( "Neumann" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                Feel::cout << "Neumann on " << marker << std::endl;
                bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
                                            _expr=( trans(id(m))*(idt(sigma)*N()) ));

                bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tauOrder)*idt(u) ) );

                bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tauOrder) ) );
            }
        }
        itType = mapField.find( "Neumann_scalar" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                Feel::cout << "Neumann on " << marker << std::endl;
                bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
                                            _expr=( trans(id(m))*(idt(sigma)*N()) ));

                bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tauOrder)*idt(u) ) );

                bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tauOrder) ) );
            }
        }
        itType = mapField.find( "Neumann_exact" );
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                std::string marker = exAtMarker.marker();
                Feel::cout << "Neumann on " << marker << std::endl;
                bbf( 2_c, 0_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker ),
                                            _expr=( trans(id(m))*(idt(sigma)*N()) ));

                bbf( 2_c, 1_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=-tau_constant * trans(id(m)) * ( pow(idv(H),M_tauOrder)*idt(u) ) );

                bbf( 2_c, 2_c) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                            _expr=tau_constant * trans(idt(uhat)) * id(m) * ( pow(idv(H),M_tauOrder) ) );
            }
        }
    }

} // end assemble STD


MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MixedElasticity<Dim, Order, G_Order,E_Order>::assembleF()
{

    auto blf = blockform1( *M_ps, M_F );

    auto w     = M_Wh->element( "w" );
    auto m     = M_Mh->element( "m" );

    // Building the RHS

    auto itField = modelProperties().boundaryConditions().find("stress");
    if (itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("SourceTerm");

        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto g = expr<Dim,1,expr_order> (exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
                blf( 1_c ) += integrate(_quad=_Q<expr_order>(),_range=elements(M_mesh),
                                        _expr=trans(g)*id(w));
            }
        }
        itType = mapField.find("Neumann");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto marker = exAtMarker.marker();
                auto g = expr<Dim,1,expr_order> (exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
                Feel::cout << "Neumann condition on " << marker << ": " << g << std::endl;
                blf( 2_c ) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                        _expr=trans(id(m))* g );
            }
        }

        itType = mapField.find("Neumann_scalar");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto marker = exAtMarker.marker();
                auto g = expr<expr_order> (exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });
                Feel::cout << "Neumann condition on " << marker << ": " << g << std::endl;
                blf( 2_c ) += integrate(_quad=_Q<expr_order>(),_range=markedfaces(M_mesh,marker),
                                        _expr= inner(g*N(), id(m)) );
            }
        }

        itType = mapField.find("Neumann_exact");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto marker = exAtMarker.marker();
                auto g = expr<Dim,1, expr_order>(exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues({ {"t", M_nm_mixedelasticity->time()} });

                for( auto const& pairMat : modelProperties().materials() )
                {
                    auto gradu_exact = grad<Dim>( g );
                    auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
                    auto material = pairMat.second;
                    auto lambda = material.property("lambda").exprScalar();
                    auto mu = material.property("mu").exprScalar();
                    auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;

                    Feel::cout << "Neumann condition computed from displacement on " << marker << std::endl;
                    blf( 2_c ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(id(m)) * sigma_exact *N() );
                }
            }
        }
    }

    itField = modelProperties().boundaryConditions().find("displacement");
    if (itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = (*itField).second;
        auto itType = mapField.find("Dirichlet");
        if ( itType != mapField.end() )
        {
            for ( auto const& exAtMarker : (*itType).second )
            {
                auto marker = exAtMarker.marker();
                auto g = expr<Dim,1, expr_order>(exAtMarker.expression());
                if ( !this->isStationary() )
                    g.setParameterValues( { {"t", M_nm_mixedelasticity->time()} } );
                Feel::cout << "Dirichlet condition on " << marker << ": " << g << std::endl;
                blf( 2_c ) += integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(id(m))*g);
            }
        }
    }


    // (u_old,w)_Omega
    if ( !this->isStationary() )
    {
        for( auto const& pairMat : modelProperties().materials() )
        {
            auto material = pairMat.second;
            auto rho = material.property("rho").exprScalar();
            auto u = this->timeStepNM()->previousUnknown(0);
            auto u1 = this->timeStepNM()->previousUnknown(1);
            auto dt = this-> timeStep();
            // blf(1_c) += integrate( _range=elements(M_mesh),
            //                         _expr= rho*inner(idv(this->timeStepNM()->polyDeriv()),id(w)) );
            blf(1_c) += integrate(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr= -rho*inner( 2*idv(u)-idv(u1) ,id(w))/(dt*dt) );
        }
    }

} // end assembleF



MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::exportResults( double time, mesh_ptrtype mesh , op_interp_ptrtype Idh , opv_interp_ptrtype Idhv )
{
    this->log("MixedElasticity","exportResults", "start");
    this->timerTool("PostProcessing").start();


    if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC && mesh  )
    {
        LOG(INFO) << "exporting on visualisation mesh at time " << time;
        M_exporter->step( time )->setMesh( mesh );
    }
    else if ( M_exporter->exporterGeometry() != EXPORTER_GEOMETRY_STATIC )
    {
        LOG(INFO) << "exporting on computational mesh at time " << time;
        M_exporter->step( time )->setMesh( M_mesh );
    }

    // Export computed solutions
    {
        for ( auto const& field : modelProperties().postProcess().exports().fields() )
        {
            if ( field == "stress" )
            {
                LOG(INFO) << "exporting stress at time " << time;

                // Exporting the stress component by component
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::X)) );
                M_exporter->step(time)->add(prefixvm(prefix(),"sigmaXX"), SXX );

                if (Dim > 1)
                {
                    auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Y, Component::Y)) );
                    auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::Y)) );
                    M_exporter->step(time)->add(prefixvm(prefix(),"sigmaYY"), SYY );
                    M_exporter->step(time)->add(prefixvm(prefix(),"sigmaXY"), SXY );
                }
                if (Dim > 2)
                {
                    auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Z, Component::Z)) );
                    auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::Y, Component::Z)) );
                    auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_up.comp(Component::X, Component::Z)) );
                    M_exporter->step(time)->add(prefixvm(prefix(),"sigmaZZ"), SZZ );
                    M_exporter->step(time)->add(prefixvm(prefix(),"sigmaYZ"), SYZ );
                    M_exporter->step(time)->add(prefixvm(prefix(),"sigmaXZ"), SXZ );
                }

                // M_exporter->step(time)->add(prefixvm(prefix(), "stress"), Idhv?(*Idhv)( M_up):M_up );

                if (M_integralCondition)
                {

                    for( auto exAtMarker : this->M_IBCList)
                    {
                        std::vector<double> force_integral(Dim);
                        auto marker = exAtMarker.marker();
                        LOG(INFO) << "exporting integral flux at time "
                                  << time << " on marker " << marker;
                        auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                                    _expr=trans(idv(M_up))*N());

                        Feel::cout << "Force computed: " << std::endl;
                        for( auto i=0;i < Dim;i++ )
                        {
                            std::string stringForce_help = (boost::format("integralForce_%1%")%i).str();
                            force_integral[i] = j_integral.evaluate()(i,0);
                            Feel::cout << force_integral[i] << std::endl;
                            M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                        }
                    }

                }

            }
            else if ( M_mesh->hasFaceMarker(field) )
            {
                auto marker = field;
                LOG(INFO) << "exporting computed force on " << marker << " at time " << time;
                std::vector<double> force_integral(Dim);

                auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker),
                                            _expr=trans(idv(M_up))*N());

                Feel::cout << "Force computed: " << std::endl;
                for( auto i=0;i < Dim;i++ )
                {
                    std::string stringForce_help = (boost::format("integralForce_%1%")%i).str();
                    force_integral[i] = j_integral.evaluate()(i,0);
                    Feel::cout << force_integral[i] << std::endl;
                    M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                }
            }
            else if ( field == "scaled_displacement" )
            {
                auto scaled_displ = M_Wh->element("scaled_displacement");
                for( auto const& pairMat : modelProperties().materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.property( "scale_displacement" ).exprScalar();

                    scaled_displ.on( _range=markedelements(M_mesh,marker) , _expr= kk*idv(M_pp));
                }

                M_exporter->step(time)->add(prefixvm(prefix(), "displacement"),Idh?(*Idh)( scaled_displ):scaled_displ ) ;



            }
            else if ( field == "scaled_stress" )
            {
                auto scaled_stress = M_Vh->element("scaled_stress");
                for( auto const& pairMat : modelProperties().materials() )
                {
                    auto marker = pairMat.first;
                    auto material = pairMat.second;
                    auto kk = material.property( "scale_stress" ).exprScalar();

                    scaled_stress.on( _range=markedelements(M_mesh,marker) , _expr= kk*idv(M_up));
                }

                // Exporting the scaled stress component by component
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::X)) );
                M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaXX"), SXX );

                if (Dim > 1)
                {
                    auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Y, Component::Y)) );
                    auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::Y)) );
                    M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaYY"), SYY );
                    M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaXY"), SXY );
                }
                if (Dim > 2)
                {
                    auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Z, Component::Z)) );
                    auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::Y, Component::Z)) );
                    auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (scaled_stress.comp(Component::X, Component::Z)) );
                    M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaZZ"), SZZ );
                    M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaYZ"), SYZ );
                    M_exporter->step(time)->add(prefixvm(prefix(),"scaled_sigmaXZ"), SXZ );
                }

                // Exporting scaled stress integral
                if (M_integralCondition)
                {
                    for( auto exAtMarker : this->M_IBCList)
                    {
                        std::vector<double> force_integral(Dim);
                        auto marker = exAtMarker.marker();
                        LOG(INFO) << "exporting scaled integral flux at time "
                                  << time << " on marker " << marker;
                        auto j_integral = integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,marker), _expr=trans(idv(scaled_stress))*N());

                        Feel::cout << "Force computed: " << std::endl;
                        for( auto i=0;i < Dim;i++ )
                        {
                            std::string stringForce_help = (boost::format("scaled_integralForce_%1%")%i).str();
                            force_integral[i] = j_integral.evaluate()(i,0);
                            Feel::cout << force_integral[i] << std::endl;
                            M_exporter->step( time )->add(prefixvm(prefix(), stringForce_help),force_integral[i]);
                        }
                    }

                }


            }
            else if ( field == "displacement" )
            {
                LOG(INFO) << "exporting displacement at time " << time;
                M_exporter->step(time)->add(prefixvm(prefix(), "displacement"),Idh?(*Idh)( M_pp):M_pp ) ;
                // Projecting on L2 space for continuity.
                auto Sh = Pch<Order> (M_mesh);
                auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                auto UX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::X]) );
                M_exporter->step(time)->add(prefixvm(prefix(), "UX"),UX ) ;

                if (Dim > 1)
                {
                    auto UY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::Y]) );
                    M_exporter->step(time)->add(prefixvm(prefix(), "UY"),UY ) ;
                }
                if (Dim > 2)
                {
                    auto UZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (M_pp[Component::Z]) );
                    M_exporter->step(time)->add(prefixvm(prefix(), "UZ"),UZ ) ;
                }

                auto itField = modelProperties().boundaryConditions().find("ExactSolution");
                if ( itField != modelProperties().boundaryConditions().end() )
                {
                    auto mapField = (*itField).second;
                    auto itType = mapField.find( "u_exact" );
                    if (itType != mapField.end() )
                    {
                        for (auto const& exAtMarker : (*itType).second )
                        {
                            if (exAtMarker.isExpression() )
                            {
                                auto u_exact = expr<Dim,1,expr_order> (exAtMarker.expression());
                                if ( !this->isStationary() )
                                    u_exact.setParameterValues( { {"t", time } } );

                                auto export_uEX = project(_quad=_Q<expr_order>(), _space=M_Wh, _range=elements( M_mesh ), _expr=u_exact);
                                M_exporter->step(time)->add(prefixvm(prefix(), "u_exact"), Idh?(*Idh)( export_uEX): export_uEX );

                                auto l2err_u = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=u_exact - idv(M_pp) );

                                auto l2norm_uex = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=u_exact );
                                if (l2norm_uex < 1)
                                    l2norm_uex = 1.0;

                                Feel::cout << "----- Computed Errors -----" << std::endl;
                                // Feel::cout << "||u-u_ex||_L2=\t" << l2err_u/l2norm_uex << std::endl;
                                Feel::cout << "||u-u_ex||_L2=\t" << l2err_u << std::endl;
                                // Export the errors
                                M_exporter -> step( time )->add(prefixvm(prefix(), "u_error_L2"), l2err_u/l2norm_uex );
                                //------ Sigma  ------//
                                auto gradu_exact = grad<Dim>(u_exact);
                                auto eps_exact   = cst(0.5) * ( gradu_exact + trans(gradu_exact) );
                                for( auto const& pairMat : modelProperties().materials() )
                                {
                                    auto material = pairMat.second;
                                    auto lambda = material.property("lambda").exprScalar();
                                    auto mu = material.property("mu").exprScalar();
                                    auto sigma_exact = lambda * trace(eps_exact) * eye<Dim>() + cst(2.) * mu * eps_exact;

                                    // EXPORT SIGMA EXACT
                                    auto export_sigmaEX = project(_quad=_Q<expr_order>(), _space=M_Vh, _range=elements(M_mesh), _expr=sigma_exact);

                                    auto Sh = Pch<Order> (M_mesh);
                                    auto l2p = opProjection(_domainSpace=Sh, _imageSpace=Sh);

                                    auto SXX = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::X)) );
                                    M_exporter->step(time)->add(prefixvm(prefix(),"s_exactXX"), SXX );

                                    if (Dim > 1)
                                    {
                                        auto SYY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Y, Component::Y)) );
                                        auto SXY = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::Y)) );
                                        M_exporter->step(time)->add(prefixvm(prefix(),"s_exactYY"), SYY );
                                        M_exporter->step(time)->add(prefixvm(prefix(),"s_exactXY"), SXY );
                                    }
                                    if (Dim > 2)
                                    {
                                        auto SYZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Y, Component::Z)) );
                                        auto SXZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::X, Component::Z)) );
                                        auto SZZ = l2p -> project (_quad=_Q<expr_order>(),_expr = idv (export_sigmaEX.comp(Component::Z, Component::Z)) );
                                        M_exporter->step(time)->add(prefixvm(prefix(),"s_exactZZ"), SZZ );
                                        M_exporter->step(time)->add(prefixvm(prefix(),"s_exactYZ"), SYZ );
                                        M_exporter->step(time)->add(prefixvm(prefix(),"s_exactXZ"), SXZ );
                                    }

                                    // M_exporter->add(prefixvm(prefix(), "sigma_exact"), Idhv?(*Idhv)( export_sigmaEX): export_sigmaEX );

                                    auto l2err_sigma = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=sigma_exact - idv(M_up) );
                                    auto l2norm_sigmaex = normL2(_quad=_Q<expr_order>(), _range=elements(M_mesh), _expr=sigma_exact );
                                    if (l2norm_sigmaex < 1)
                                        l2norm_sigmaex = 1.0;
                                    Feel::cout << "||sigma-sigma_ex||_L2=\t" << l2err_sigma/l2norm_sigmaex << std::endl;
                                    Feel::cout << "---------------------------" << std::endl;
                                    // Export the errors
                                    M_exporter -> step( time )->add(prefixvm(prefix(), "sigma_error_L2"), l2err_sigma/l2norm_sigmaex );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    this->timerTool("PostProcessing").stop("exportResults");
    if ( this->scalabilitySave() )
    {
        if ( !this->isStationary() )
            this->timerTool("PostProcessing").setAdditionalParameter("time",this->currentTime());
        this->timerTool("PostProcessing").save();
    }
    this->log("MixedElasticity","exportResults", "finish");
}

MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::geometricTest( ){

    auto itField = modelProperties().boundaryConditions().find("GeometricalTest");
    if ( itField != modelProperties().boundaryConditions().end() )
    {
        auto mapField = itField -> second;
        auto itType = mapField.find( "force_F" );
        if (itType != mapField.end() )
        {
            for (auto const& exAtMarker : itType->second )
            {
                auto curvedForce = (integrate(_quad=_Q<expr_order>(), _range=markedfaces(M_mesh,exAtMarker.marker()), _expr = idv(M_up)*N() )).evaluate();
                auto forceF = (expr<Dim,1,expr_order> (exAtMarker.expression() )).evaluate();
                /*
                 Feel::cout << "Force F uploaded:\t" << forceF << std::endl;
                 Feel::cout << "Force F computed from M_up:\t" << curvedForce << std::endl;
                 */
                auto curveError = (curvedForce - forceF).cwiseAbs();

                Feel::cout << "Error for geometrical order:\t" << curveError << std::endl;
            }
        }
        /*
         itType = mapField.find( "force_F_2" );
         if (itType != mapField.end() )
         {
         for (auto const& exAtMarker : itType->second )
         {
         auto forceF_2 = expr<Dim,1,expr_order> (exAtMarker.expression() );
         auto forceIntegral = (integrate( _range=markedfaces(M_mesh,exAtMarker.marker()), _expr = forceF_2 * N() )).evaluate();
         Feel::cout << "Force F computed from input:\t" << forceIntegral << std::endl;

         }
         }
         */
    }
}


// Time exporter
MIXEDELASTICITY_CLASS_TEMPLATE_DECLARATIONS
void
MIXEDELASTICITY_CLASS_TEMPLATE_TYPE::exportTimers()
{
    if( Environment::isMasterRank() )
    {
        std::ofstream timers( "timers.dat", std::ios::out | std::ios::trunc);
        std::string fmtS = "";
        for( int i = 1; i <= M_timers.size(); ++i )
            fmtS += "%" + std::to_string(i) + "% %|" + std::to_string(14*i) + "t|";
        boost::format fmt(fmtS);
        for( auto const& pair : M_timers )
            fmt % pair.first;
        timers << fmt << std::endl;
        if (isStationary())
        {
            for( auto const& pair : M_timers )
                fmt % pair.second;
            timers << fmt << std::endl;
        }

        /*
         //( for( int i = 0; this->timeStepBase()->isFinished() ; ++i )
         {
         for( auto const& pair : M_timers )
         fmt % pair.second[i];
         timers << fmt << std::endl;
         }*/

        timers.close();
    }
}

} // namespace Feel
} // namespace FeelModels
