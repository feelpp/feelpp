#include <feel/feelmodels/hdg/mixedpoisson.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

inline
po::options_description
makeThermoElectricHDGOptions( std::string prefix = "thermoelectric",
                              std::string prefixThermo = "thermo",
                              std::string prefixElectro = "electro" )
{
    po::options_description teOptions( "ThermoElectric HDG options");
    teOptions.add_options()
        ( prefixvm( prefix, "tolerance").c_str(), po::value<double>()->default_value( 1e-5 ),
          "picard tolerance" )
        ( prefixvm( prefix, "itmax").c_str(), po::value<int>()->default_value( 20 ),
          "picard max iteration")
        ( prefixvm( prefix, "input-marker").c_str(), po::value<std::string>()->default_value("V1"),
          "marker for input (where V=0)")
        ( prefixvm( prefix, "output-marker").c_str(), po::value<std::string>()->default_value("V0"),
          "marker for output (where V!=0)")
        ( prefixvm( prefix, "continuation").c_str(), po::value<bool>()->default_value(false),
          "use continuation for Picard" )
        ( prefixvm( prefix, "continuation.steps").c_str(), po::value<int>()->default_value(1),
          "number of steps to use for continuation for Picard" )
        ( prefixvm( prefix, "load-initial-guess").c_str(), po::value<bool>()->default_value(false), "load initial guess" )
        ( prefixvm( prefix, "load-initial-guess.path").c_str(), po::value<std::string>()->default_value(""), "path to load initial guess" )
        ( prefixvm( prefix, "use-current-joule").c_str(), po::value<int>()->default_value(0), "use idv(current) or gradv(potential)" )
        ( "quad", po::value<int>()->default_value(2), "quadrature order" )
        ( prefixvm( prefix, "current-expr").c_str(), po::value<std::string>()->default_value("{0,0,0}"), "expression for current density" )
        ;
    teOptions.add( makeMixedPoissonOptions( prefixThermo ) );
    teOptions.add( makeMixedPoissonOptions( prefixElectro ) );
    return teOptions;
}

template<int Dim, int OrderT, int OrderV, int OrderG>
class ThermoElectricHDG
{
public:
    using thermo_type = FeelModels::MixedPoisson<Dim, OrderT, OrderG>;
    using thermo_ptrtype = std::shared_ptr<thermo_type>;
    using temp_type = typename thermo_type::Wh_element_t;
    using temp_ptrtype = typename thermo_type::Wh_element_ptr_t;
    using tempflux_type = typename thermo_type::Vh_element_t;
    using tempflux_ptrtype = typename thermo_type::Vh_element_ptr_t;

    using electro_type = FeelModels::MixedPoisson<Dim, OrderV, OrderG>;
    using electro_ptrtype = std::shared_ptr<electro_type>;
    using potential_type = typename electro_type::Wh_element_t;
    using potential_ptrtype = typename electro_type::Wh_element_ptr_t;
    using current_type = typename electro_type::Vh_element_t;
    using current_ptrtype = typename electro_type::Vh_element_ptr_t;

    using mesh_type = typename thermo_type::mesh_type;
    using mesh_ptrtype = typename thermo_type::mesh_ptrtype;

    using init_guess_space_type = FunctionSpace<mesh_type, bases<Lagrange<OrderT, Scalar, Continuous, PointSetFekete> > >;
    using init_guess_space_ptrtype = std::shared_ptr<init_guess_space_type>;
    using init_guess_type = typename init_guess_space_type::element_type;

private:
    std::string M_prefix;
    mesh_ptrtype M_mesh;

    std::string M_prefixThermo;
    thermo_ptrtype M_thermo;
    temp_type M_temperature;
    tempflux_type M_tempflux;

    std::string M_prefixElectro;
    electro_ptrtype M_electro;
    potential_type M_potential;
    current_type M_current;

    init_guess_space_ptrtype M_initGuessSpace;
    init_guess_type M_initGuess;

    potential_type M_joule;
    potential_type M_sigma;
    potential_type M_k;


    int M_itMax;
    double M_tolerance;
    bool M_useContinuation;
    int M_continuationSteps;
    std::string M_inputMarker;
    std::string M_outputMarker;

public:
    ThermoElectricHDG( std::string prefix = "thermoelectric",
                       std::string prefixThermo = "thermo",
                       std::string prefixElectro = "electro" );
    void init( mesh_ptrtype const& mesh = nullptr );
    void solve();
    void exportResults();

    mesh_ptrtype mesh() const { return M_mesh; }
    thermo_ptrtype thermo() const { return M_thermo; }
    electro_ptrtype electro() const { return M_electro; }
    ModelProperties modelProperties() const { return M_electro->modelProperties(); }
};

template<int Dim, int OrderT, int OrderV, int OrderG>
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::ThermoElectricHDG( std::string prefix,
                                                           std::string prefixThermo,
                                                           std::string prefixElectro )
    :
    M_prefix(prefix),
    M_prefixThermo("hdg.poisson."+prefixThermo),
    M_prefixElectro("hdg.poisson."+prefixElectro),
    M_itMax(ioption(prefixvm( M_prefix, "itmax"))),
    M_tolerance(doption(prefixvm( M_prefix, "tolerance"))),
    M_useContinuation(boption(prefixvm( M_prefix, "continuation"))),
    M_continuationSteps(ioption(prefixvm( M_prefix, "continuation.steps"))),
    M_inputMarker(soption(prefixvm( M_prefix, "input-marker"))),
    M_outputMarker(soption(prefixvm( M_prefix, "output-marker")))
{
    tic();
    M_thermo = thermo_type::New(M_prefixThermo);
    M_electro = electro_type::New(M_prefixElectro);
    toc("construct");
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::init( mesh_ptrtype const& mesh )
{
    tic();
    if( mesh )
        M_mesh = mesh;
    else
        M_mesh = loadMesh( new mesh_type );
    M_thermo->init( M_mesh );
    M_electro->init( M_mesh );
    toc("init");
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::solve()
{
    auto electroProp = M_electro->modelProperties();
    auto electroMat = electroProp.materials();
    auto thermoProp = M_thermo->modelProperties();
    auto thermoMat = M_thermo->modelProperties().materials();
    double V = 0, Vinit = 0, I = 0, Iinit = 0;
    bool conditionIbc = false;
    if( M_useContinuation )
    {
        auto itField = electroProp.boundaryConditions().find( "potential");
        if ( itField != electroProp.boundaryConditions().end() )
        {
            auto mapField = (*itField).second;
            auto itType = mapField.find( "Dirichlet" );
            if ( itType != mapField.end() )
            {
                for ( auto const& exAtMarker : (*itType).second )
                {
                    std::string marker = exAtMarker.marker();
                    if(marker == M_outputMarker )
                        Vinit = std::stod(exAtMarker.expression());
                }
            }
        }
        if( Vinit == 0 )
        {
            itField = electroProp.boundaryConditions().find( "flux");
            if ( itField != electroProp.boundaryConditions().end() )
            {
                auto mapField = (*itField).second;
                auto itType = mapField.find( "Integral" );
                if ( itType != mapField.end() )
                {
                    for ( auto const& exAtMarker : (*itType).second )
                    {
                        std::string marker = exAtMarker.marker();
                        if(marker == M_outputMarker )
                            Iinit = std::stod(exAtMarker.expression());
                    }
                }
            }
            if( Iinit != 0 )
            {
                conditionIbc = true;
                I = Iinit/M_continuationSteps;
            }
            else
                Feel::cout << tc::red << "marker for continuation " << M_outputMarker
                           << " not found in Dirichlet or Ibc. Continuation not used !"
                           << tc::reset << std::endl;
        }
        else
        {
            V = Vinit/M_continuationSteps;
        }
    }

    M_joule = M_electro->potentialSpace()->element();
    M_k = M_electro->potentialSpace()->element();
    M_sigma = M_electro->potentialSpace()->element();

    M_electro->setCstMatrixToZero();
    M_electro->setVectorToZero();
    M_thermo->setCstMatrixToZero();
    M_thermo->setVectorToZero();

    tic();
    M_electro->assembleCstPart();
    toc("assembleCstElectro");

    // first iteration
#ifndef USE_SAME_MAT
    Feel::cout << tc::red << "WARNING copy cst part called !!!" << tc::reset << std::endl;
    M_electro->copyCstPart();
#endif
    if( boption("thermoelectric.load-initial-guess") )
    {
        M_initGuessSpace = std::make_shared<init_guess_space_type>(M_mesh);
        M_initGuess = M_initGuessSpace->element();
        M_initGuess.load(_path=soption("thermoelectric.load-initial-guess.path"));
        for( auto const& pairMat : electroMat )
        {
            std::string marker = pairMat.first;
            auto mat = pairMat.second;
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            auto sigma = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                                       {"T"}, {idv(M_initGuess)},
                                       {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            M_electro->updateConductivityTerm( sigma, marker);
            M_sigma += vf::project( _space=M_electro->potentialSpace(),
                                       _range=markedelements(M_mesh,marker),
                                       _expr=sigma );
        }
    }
    else
    {
        M_electro->updateConductivityTerm( false );
    }
    M_electro->assembleRhsBoundaryCond();
    if( M_useContinuation )
    {
        if( conditionIbc )
        {
            M_electro->assembleRhsIBC(0, M_outputMarker, I-Iinit);
        }
        else
        {
            auto g = expr(std::to_string(V-Vinit));
            M_electro->assembleRhsDirichlet(g,M_outputMarker);
        }
    }
    M_electro->assemblePotentialRHS( cst(0.), "");
    M_electro->solve();
    M_potential = M_electro->potentialField();
    M_current = M_electro->fluxField();


    tic();
    M_thermo->assembleCstPart();
    toc("assembleCstThermo");

#ifndef USE_SAME_MAT
    Feel::cout << tc::red << "WARNING copy cst part called !!!" << tc::reset << std::endl;
    M_thermo->copyCstPart();
#endif
    for( auto const& pairMat : thermoMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        double alpha = mat.getDouble("alpha");
        double T0 = mat.getDouble("T0");
        auto sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
        double k0 = 0, L = 0;
        if( mat.hasProperty("Lorentz") )
            L = mat.getDouble("Lorentz");
        if( mat.hasProperty(soption(prefixvm(M_prefixThermo,"conductivity_json"))) )
            k0 = mat.getDouble(soption(prefixvm(M_prefixThermo,"conductivity_json")));
        auto k = mat.getScalar(soption(prefixvm(M_prefixThermo,"conductivityNL_json")),
                               {{"T",T0},{"k0",k0},{"T0",T0},{"alpha",alpha},{"Lorentz",L},{"sigma0",sigma0}});
        M_thermo->updateConductivityTerm( k, marker);
        M_k += vf::project( _space=M_electro->potentialSpace(),
                               _range=markedelements(M_mesh,marker),
                               _expr=k );
    }
    M_thermo->assembleRhsBoundaryCond();
    for( auto const& pairMat : electroMat )
    {
        std::string marker = pairMat.first;
        auto mat = pairMat.second;
        if( boption("thermoelectric.load-initial-guess") )
        {
            double alpha = mat.getDouble("alpha");
            double T0 = mat.getDouble("T0");
            double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            auto cond = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                                      {"T"}, {idv(M_initGuess)},
                                      {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
            auto rhs = inner(idv(M_current),idv(M_current))/cond;
            M_joule += vf::project( _space=M_electro->potentialSpace(),
                                  _range=markedelements(M_mesh,marker),
                                  _expr=rhs );
            M_thermo->assemblePotentialRHS( rhs, marker);
        }
        else
        {
            auto cond = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            auto sigma = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
            if( ioption("thermoelectric.use-current-joule") == 1 )
            {
                auto rhs = inner(idv(M_current),idv(M_current))/sigma;
                M_joule += vf::project( _space=M_electro->potentialSpace(),
                                      _range=markedelements(M_mesh,marker),
                                      _expr=rhs );
                double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                Feel::cout << "norm J     = " << n << std::endl
                           << "integral J = " << nn << std::endl
                           << "min J      = " << mima.min() << std::endl
                           << "max J      = " << mima.max() << std::endl;
                M_thermo->assemblePotentialRHS( rhs, marker);
            }
            else if( ioption("thermoelectric.use-current-joule") == 0 )
            {
                auto rhs = sigma*inner(gradv(M_potential));
                M_joule += vf::project( _space=M_electro->potentialSpace(),
                                      _range=markedelements(M_mesh,marker),
                                      _expr=rhs );
                double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                Feel::cout << "norm J     = " << n << std::endl
                           << "integral J = " << nn << std::endl
                           << "min J      = " << mima.min() << std::endl
                           << "max J      = " << mima.max() << std::endl;
                M_thermo->assemblePotentialRHS( rhs, marker);
            }
            else if( ioption("thermoelectric.use-current-joule") == 2 )
            {
                auto j = expr<3,1>(soption("thermoelectric.current-expr"));
                auto rhs = sigma*inner(j);
                M_joule += vf::project( _space=M_electro->potentialSpace(),
                                      _range=markedelements(M_mesh,marker),
                                      _expr=rhs,
                                      _quad=ioption("quad") );
                double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                Feel::cout << "use current expr" << std::endl
                           << "\tnorm J     = " << n << std::endl
                           << "\tintegral J = " << nn << std::endl
                           << "\tmin J      = " << mima.min() << std::endl
                           << "\tmax J      = " << mima.max() << std::endl;
                M_thermo->assemblePotentialRHS( rhs, marker);
            }
        }
    }
    M_thermo->solve();
    M_temperature = M_thermo->potentialField();
    M_tempflux = M_thermo->fluxField();

    potential_type oldPotential = M_electro->potentialSpace()->element();
    current_type oldCurrent = M_electro->fluxSpace()->element();
    temp_type oldTemperature = M_thermo->potentialSpace()->element();
    tempflux_type oldTempflux = M_thermo->fluxSpace()->element();
    // auto rangeV = M_electro->potentialSpace()->dof()->meshSupport()->rangeElements();
    // auto rangeT = M_thermo->potentialSpace()->dof()->meshSupport()->rangeElements();
    // double normV = normL2( rangeV, idv(M_potential) );
    // double normT = normL2( rangeT, idv(M_temperature) );
    // double incrV = normL2( rangeV, idv(M_potential) - idv(oldPotential) );
    // double incrT = normL2( rangeT, idv(M_temperature) - idv(oldTemperature) );
    double normV = normL2( elements(M_mesh), idv(M_potential) );
    double normT = normL2( elements(M_mesh), idv(M_temperature) );
    double incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
    double incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );
    double relV = incrV/normV;
    double relT = incrT/normT;

    double tol = doption( prefixvm( M_prefix, "tolerance") );
    int itmax = ioption( prefixvm( M_prefix, "itmax") );


    int steps = M_useContinuation ? M_continuationSteps : 1;
    for( int c = 1; c <= steps; ++c)
    {
        if( M_useContinuation )
        {
            if( conditionIbc )
            {
                I = c*Iinit/M_continuationSteps;
                Feel::cout << "Start Picard with I =  " << I << " (" << itmax
                           << " iterations max and tolerance = " << tol << ")" << std::endl;
            }
            else
            {
                V = c*Vinit/M_continuationSteps;
                Feel::cout << "Start Picard with V =  " << V << " (" << itmax
                           << " iterations max and tolerance = " << tol << ")" << std::endl;
            }
        }
        else
            Feel::cout << "Start Picard  (" << itmax
                       << " iterations max and tolerance = " << tol << ")" << std::endl;

        int i = 0;
        relV = 1;
        relT = 1;

        // Picard loop
        while( i++ < itmax && ( relV > tol || relT > tol ) )
        {
            tic();

            auto meanT = mean(_range=elements(M_mesh), _expr=idv(M_temperature));
            auto mima = minmax(_range=elements(M_mesh), _pset=_Q<>(), _expr=idv(M_temperature));
            auto mi = mima.min();
            auto ma = mima.max();
            auto power = integrate(_range=elements(M_mesh), _expr=idv(M_joule)).evaluate()(0,0);
            auto int0 = integrate(_range=markedfaces(M_mesh,M_outputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);
            auto int1 = integrate(_range=markedfaces(M_mesh,M_inputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);

            Feel::cout << std::setw(20) << "tempMax";
            Feel::cout << std::setw(20) << "tempMoy";
            Feel::cout << std::setw(20) << "tempMin";
            Feel::cout << std::setw(20) << "intensite0";
            Feel::cout << std::setw(20) << "intensite1";
            Feel::cout << std::setw(20) << "power";
            Feel::cout << std::endl;
            Feel::cout << std::setw(20) << ma;
            Feel::cout << std::setw(20) << meanT;
            Feel::cout << std::setw(20) << mi;
            Feel::cout << std::setw(20) << int0;
            Feel::cout << std::setw(20) << int1;
            Feel::cout << std::setw(20) << power;
            Feel::cout << std::endl;

            oldPotential = M_potential;
            oldCurrent = M_current;
            oldTemperature = M_temperature;
            oldTempflux = M_tempflux;
            M_sigma.zero();
            M_k.zero();
            M_joule.zero();

            tic();
#ifdef USE_SAME_MAT
            M_electro->setCstMatrixToZero();
            M_electro->setVectorToZero();
            M_electro->assembleCstPart();
#else
            M_electro->copyCstPart();
#endif
            for( auto const& pairMat : electroMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
                auto sigma = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                                           {"T"}, {idv(M_temperature)},
                                           {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                M_electro->updateConductivityTerm( sigma, marker);
                M_sigma += vf::project( _space=M_electro->potentialSpace(),
                                           _range=markedelements(M_mesh,marker),
                                           _expr=sigma );
            }
            M_electro->assembleRhsBoundaryCond();
            if( M_useContinuation )
            {
                if( conditionIbc )
                {
                    M_electro->assembleRhsIBC(0, M_outputMarker, I-Iinit);
                }
                else
                {
                    auto g = expr(std::to_string(V-Vinit));
                    M_electro->assembleRhsDirichlet(g,M_outputMarker);
                }
            }
            M_electro->assemblePotentialRHS( cst(0.), "");
            toc("assembleElectro");
            tic();
            M_electro->solve();
            toc("solveElectro");
            M_potential = M_electro->potentialField();
            M_current = M_electro->fluxField();

            tic();
#ifdef USE_SAME_MAT
            M_thermo->setCstMatrixToZero();
            M_thermo->setVectorToZero();
            M_thermo->assembleCstPart();
#else
            M_thermo->copyCstPart();
#endif
            for( auto const& pairMat : thermoMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                auto sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
                double k0 = 0, L = 0;
                if( mat.hasProperty("Lorentz") )
                    L = mat.getDouble("Lorentz");
                if( mat.hasProperty(soption(prefixvm(M_prefixThermo,"conductivity_json"))) )
                    k0 = mat.getDouble(soption(prefixvm(M_prefixThermo,"conductivity_json")));
                auto k = mat.getScalar(soption(prefixvm(M_prefixThermo,"conductivityNL_json")),
                                       {"T"}, {idv(M_temperature)},
                                       {{"k0",k0},{"T0",T0},{"alpha",alpha},
                                        {"Lorentz",L},{"sigma0",sigma0}});
                M_thermo->updateConductivityTerm( k, marker);
                M_k += vf::project( _space=M_electro->potentialSpace(),
                                       _range=markedelements(M_mesh,marker),
                                       _expr=k );
            }
            M_thermo->assembleRhsBoundaryCond();
            for( auto const& pairMat : electroMat )
            {
                std::string marker = pairMat.first;
                auto mat = pairMat.second;
                double alpha = mat.getDouble("alpha");
                double T0 = mat.getDouble("T0");
                double sigma0 = mat.getDouble(soption(prefixvm(M_prefixElectro,"conductivity_json")));
                // auto sigma = mat.getScalar(soption(prefixvm(M_prefixElectro,"conductivityNL_json")),
                //                            {"T"}, {idv(M_temperature)},
                //                            {{"sigma0",sigma0},{"alpha",alpha},{"T0",T0}});
                auto sigma = cst(sigma0)/(1+cst(alpha)*(idv(M_temperature)-cst(T0)));
                if( ioption("thermoelectric.use-current-joule") == 1 )
                {
                    auto rhs = inner(idv(M_current),idv(M_current))/sigma;
                    M_joule += vf::project( _space=M_electro->potentialSpace(),
                                          _range=markedelements(M_mesh,marker),
                                          _expr=rhs,
                                          _quad=ioption("quad"));
                    double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                    double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                    auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                    Feel::cout << "use inner(current)" << std::endl
                               << "\tnorm J     = " << n << std::endl
                               << "\tintegral J = " << nn << std::endl
                               << "\tmin J      = " << mima.min() << std::endl
                               << "\tmax J      = " << mima.max() << std::endl;
                    M_thermo->assemblePotentialRHS( rhs, marker);
                }
                else if( ioption("thermoelectric.use-current-joule") == 0 )
                {
                    auto rhs = sigma*inner(gradv(M_potential));
                    M_joule += vf::project( _space=M_electro->potentialSpace(),
                                          _range=markedelements(M_mesh,marker),
                                          _expr=rhs,
                                          _quad=ioption("quad") );
                    double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                    double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                    auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                    Feel::cout << "use inner(grad(v))" << std::endl
                               << "\tnorm J     = " << n << std::endl
                               << "\tintegral J = " << nn << std::endl
                               << "\tmin J      = " << mima.min() << std::endl
                               << "\tmax J      = " << mima.max() << std::endl;
                    M_thermo->assemblePotentialRHS( rhs, marker);
                }
                else if( ioption("thermoelectric.use-current-joule") == 2 )
                {
                    auto j = expr<3,1>(soption("thermoelectric.current-expr"));
                    auto rhs = sigma*inner(j);
                    M_joule += vf::project( _space=M_electro->potentialSpace(),
                                          _range=markedelements(M_mesh,marker),
                                          _expr=rhs,
                                          _quad=ioption("quad") );
                    double n = normL2(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad"));
                    double nn = integrate(_range=markedelements(M_mesh,marker),_expr=rhs,_quad=ioption("quad")).evaluate()(0,0);
                    auto mima = minmax(_range=markedelements(M_mesh,marker), _pset=_Q<>(ioption("quad")), _expr=rhs);
                    Feel::cout << "use current expr" << std::endl
                               << "\tnorm J     = " << n << std::endl
                               << "\tintegral J = " << nn << std::endl
                               << "\tmin J      = " << mima.min() << std::endl
                               << "\tmax J      = " << mima.max() << std::endl;
                    M_thermo->assemblePotentialRHS( rhs, marker);
                }
            }
            toc("assembleThermo");
            tic();
            M_thermo->solve();
            toc("solveThermo");
            M_temperature = M_thermo->potentialField();
            M_tempflux = M_thermo->fluxField();

            incrV = normL2( elements(M_mesh), idv(M_potential) - idv(oldPotential) );
            incrT = normL2( elements(M_mesh), idv(M_temperature) - idv(oldTemperature) );
            normV = normL2( elements(M_mesh), idv(oldPotential) );
            normT = normL2( elements(M_mesh), idv(oldTemperature) );
            relV = incrV/normV;
            relT = incrT/normT;

            Feel::cout << "iteration " << i
                       << "\n\tincrement on V = " << relV << " (absolute = " << incrV << ")"
                       << "\n\tincrement on T = " << relT << " (absolute = " << incrT << ")" << std::endl;

            toc("loop");
        } // Picard loop
    }
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::exportResults()
{
    M_electro->exportResults();
    M_thermo->exportResults();
    auto e = exporter( M_mesh );
    e->add("potential", M_potential);
    e->add("current", M_current);
    e->add("temperature", M_temperature);
    e->add("tempflux", M_tempflux);
    e->add("joule", M_joule);
    e->add("sigma", M_sigma);
    e->add("k", M_k);
    e->save();

    auto meanT = mean(_range=elements(M_mesh), _expr=idv(M_temperature));
    auto mima = minmax(_range=elements(M_mesh), _pset=_Q<>(), _expr=idv(M_temperature));
    auto mi = mima.min();
    auto ma = mima.max();
    auto power = integrate(_range=elements(M_mesh), _expr=idv(M_joule)).evaluate()(0,0);
    auto int0 = integrate(_range=markedfaces(M_mesh,M_outputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);
    auto int1 = integrate(_range=markedfaces(M_mesh,M_inputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);

    Feel::cout << std::setw(20) << "tempMax";
    Feel::cout << std::setw(20) << "tempMoy";
    Feel::cout << std::setw(20) << "tempMin";
    Feel::cout << std::setw(20) << "intensite0";
    Feel::cout << std::setw(20) << "intensite1";
    Feel::cout << std::setw(20) << "power";
    Feel::cout << std::endl;
    Feel::cout << std::setw(20) << ma;
    Feel::cout << std::setw(20) << meanT;
    Feel::cout << std::setw(20) << mi;
    Feel::cout << std::setw(20) << int0;
    Feel::cout << std::setw(20) << int1;
    Feel::cout << std::setw(20) << power;
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        std::ofstream file ( "measures.csv" );
        if( file )
        {
            file << std::setw(20) << "tempMax";
            file << std::setw(20) << "tempMoy";
            file << std::setw(20) << "tempMin";
            file << std::setw(20) << "intensite0";
            file << std::setw(20) << "intensite1";
            file << std::setw(20) << "power";
            file << std::endl;
            file << std::setw(20) << ma;
            file << std::setw(20) << meanT;
            file << std::setw(20) << mi;
            file << std::setw(20) << int0;
            file << std::setw(20) << int1;
            file << std::setw(20) << power;
            file << std::endl;
            file.close();
        }
    }
}
