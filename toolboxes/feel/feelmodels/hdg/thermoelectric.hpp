#include <feel/feelmodels/hdg/mixedpoisson.hpp>
#include <feel/feelfilters/loadmesh.hpp>

using namespace Feel;
using namespace Feel::FeelModels;

inline
po::options_description
makeThermoElectricHDGOptions( std::string prefix = "thermoelectric",
                              std::string prefixHeat = "heat",
                              std::string prefixElectric = "electric" )
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
        ;
    teOptions.add( mixedpoisson_options( prefixHeat ) );
    teOptions.add( mixedpoisson_options( prefixElectric ) );
    return teOptions;
}

template<int Dim, int OrderT, int OrderV, int OrderG>
class ThermoElectricHDG
{
public:
    using value_type = double;
    using self_type = ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>;
    using heat_type = FeelModels::MixedPoisson<Simplex<Dim, OrderG>, OrderT>;
    using heat_ptrtype = std::shared_ptr<heat_type>;
    using temp_type = typename heat_type::element_potential_type;
    using temp_ptrtype = typename heat_type::element_potential_ptrtype;
    using tempflux_type = typename heat_type::element_flux_type;
    using tempflux_ptrtype = typename heat_type::element_flux_ptrtype;

    using electric_type = FeelModels::MixedPoisson<Simplex<Dim, OrderG>, OrderV>;
    using electric_ptrtype = std::shared_ptr<electric_type>;
    using potential_type = typename electric_type::element_potential_type;
    using potential_ptrtype = typename electric_type::element_potential_ptrtype;
    using current_type = typename electric_type::element_flux_type;
    using current_ptrtype = typename electric_type::element_flux_ptrtype;

    using property_type = typename electric_type::element_potential_type;

    using mesh_type = typename heat_type::mesh_type;
    using mesh_ptrtype = typename heat_type::mesh_ptrtype;

    using init_guess_space_type = FunctionSpace<mesh_type, bases<Lagrange<OrderT, Scalar, Continuous, PointSetFekete> > >;
    using init_guess_space_ptrtype = std::shared_ptr<init_guess_space_type>;
    using init_guess_type = typename init_guess_space_type::element_type;

private:
    std::string M_prefix;
    mesh_ptrtype M_mesh;

    std::string M_prefixHeat;
    heat_ptrtype M_heatModel;
    temp_type M_temperature;
    tempflux_type M_tempFlux;

    std::string M_prefixElectric;
    electric_ptrtype M_electricModel;
    potential_type M_potential;
    current_type M_current;

    int M_itMax;
    double M_tolerance;
    bool M_useContinuation;
    int M_continuationSteps;
    std::string M_inputMarker;
    std::string M_outputMarker;

public:
    ThermoElectricHDG( std::string prefix = "thermoelectric",
                       std::string prefixHeat = "heat",
                       std::string prefixElectric = "electric" );
    void init( mesh_ptrtype const& mesh = nullptr );
    void printAndSaveInfo();
    void solve();
    void exportResults();
    bool checkResults() const {
        bool ch = M_heatModel->checkResults();
        bool ce = M_electricModel->checkResults();
        return ch && ce;
    }

    mesh_ptrtype mesh() const { return M_mesh; }
    heat_ptrtype heatModel() const { return M_heatModel; }
    electric_ptrtype electricModel() const { return M_electricModel; }
    ModelProperties modelProperties() const { return M_electricModel->modelProperties(); }
    void updateLinear_Heat( Feel::FeelModels::ModelAlgebraic::DataUpdateLinear & data ) const;
    void updateLinear_Electric( Feel::FeelModels::ModelAlgebraic::DataUpdateLinear & data ) const;

    auto modelFields( std::string const& prefix = "" ) const
        {
            return Feel::FeelModels::modelFields( this->heatModel()->modelFields( prefixvm( prefix, this->heatModel()->keyword() ) ),
                                                  this->electricModel()->modelFields( prefixvm( prefix, this->electricModel()->keyword() ) ) );
        }
    auto symbolsExpr( std::string const& prefix = "") const
        {
            auto seHeat = this->heatModel()->symbolsExpr( prefixvm( prefix, this->heatModel()->keyword() ) );
            auto seElectric = this->electricModel()->symbolsExpr( prefixvm( prefix, this->electricModel()->keyword() ) );
            return Feel::vf::symbolsExpr( seHeat,seElectric );
        }
    auto modelContext( std::string const& prefix = "" ) const
        {
            auto mfields = this->modelFields(prefix);
            auto se = this->symbolsExpr(prefix).template createTensorContext<mesh_type>();
            return Feel::FeelModels::modelContext( std::move( mfields ), std::move( se ) );
        }

};

template<int Dim, int OrderT, int OrderV, int OrderG>
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::ThermoElectricHDG( std::string prefix,
                                                                   std::string prefixHeat,
                                                                   std::string prefixElectric )
    :
    M_prefix(prefix),
    M_prefixHeat(prefixHeat),
    M_prefixElectric(prefixElectric),
    M_itMax(ioption(prefixvm( M_prefix, "itmax"))),
    M_tolerance(doption(prefixvm( M_prefix, "tolerance"))),
    M_useContinuation(boption(prefixvm( M_prefix, "continuation"))),
    M_continuationSteps(ioption(prefixvm( M_prefix, "continuation.steps"))),
    M_inputMarker(soption(prefixvm( M_prefix, "input-marker"))),
    M_outputMarker(soption(prefixvm( M_prefix, "output-marker")))
{
    tic();
    M_heatModel = heat_type::New(M_prefixHeat, MixedPoissonPhysics::Heat);
    M_electricModel = electric_type::New(M_prefixElectric, MixedPoissonPhysics::Electric);
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

    M_heatModel->setMesh( M_mesh );
    M_heatModel->init();
    M_heatModel->algebraicFactory()->setFunctionLinearAssembly(std::bind( &self_type::updateLinear_Heat, std::ref(*this), std::placeholders::_1 ));

    M_electricModel->setMesh( M_mesh );
    M_electricModel->init();
    M_electricModel->algebraicFactory()->setFunctionLinearAssembly(std::bind( &self_type::updateLinear_Electric, std::ref(*this), std::placeholders::_1 ));
    toc("init");
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::printAndSaveInfo()
{
    M_electricModel->printAndSaveInfo();
    M_heatModel->printAndSaveInfo();
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::updateLinear_Heat( Feel::FeelModels::ModelAlgebraic::DataUpdateLinear & data ) const
{
    bool buildCstPart = data.buildCstPart();
    bool buildNonCstPart = !buildCstPart;

    auto mctx = this->modelContext();
    M_heatModel->updateLinearPDE( data, mctx );
    auto const& symbolsExpr = mctx.symbolsExpr();

    auto F = std::dynamic_pointer_cast<condensed_vector_t<value_type>>(data.rhs());
    auto mesh = M_heatModel->mesh();
    auto ps = M_heatModel->spaceProduct();
    auto blf = blockform1( ps, F );
    auto v = M_electricModel->fieldPotential();
    auto j = M_electricModel->fieldFlux();
    auto phi = M_heatModel->fieldPotential();
    if ( buildNonCstPart )
    {
        for ( auto const& [physicName,physicData] : M_electricModel->physicsFromCurrentType() )
        {
            for ( std::string const& matName : M_electricModel->materialsProperties()->physicToMaterials( physicName ) )
            {
                auto coeff_c = M_electricModel->materialsProperties()->materialProperty( matName, M_electricModel->diffusionCoefficientName() );
                auto coeff_c_expr = expr( coeff_c.expr(), symbolsExpr );
                auto const& range = M_electricModel->materialsProperties()->rangeMeshElementsByMaterial( this->mesh(),matName );
                blf(1_c) += integrate( _range=range,
                                       _expr=inner(idv(M_current),idv(M_current))/coeff_c_expr*id(phi) );
            }
        }
    }
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::updateLinear_Electric( Feel::FeelModels::ModelAlgebraic::DataUpdateLinear & data ) const
{
    auto mctx = this->modelContext();
    M_electricModel->updateLinearPDE( data, mctx );
}


template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::solve()
{
    potential_type oldPotential;
    current_type oldCurrent;
    temp_type oldTemperature;
    tempflux_type oldTempFlux;

    int i = 0;
    double relV = 1, relT = 1, incrV, incrT, normV, normT;

    M_potential = M_electricModel->fieldPotential();
    M_current = M_electricModel->fieldFlux();
    M_temperature = M_heatModel->fieldPotential();
    M_tempFlux = M_heatModel->fieldFlux();

    // Picard loop
    while( i++ < M_itMax && ( relV > M_tolerance || relT > M_tolerance ) )
    {
        tic();
        oldPotential = M_electricModel->fieldPotential();
        oldCurrent = M_electricModel->fieldFlux();
        oldTemperature = M_heatModel->fieldPotential();
        oldTempFlux = M_heatModel->fieldFlux();

        tic();
        M_electricModel->solve();
        toc("solve electric model");
        M_potential = M_electricModel->fieldPotential();
        M_current = M_electricModel->fieldFlux();

        tic();
        M_heatModel->solve();
        toc("solve heat model");
        M_temperature = M_heatModel->fieldPotential();
        M_tempFlux = M_heatModel->fieldFlux();

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
    }
}

template<int Dim, int OrderT, int OrderV, int OrderG>
void
ThermoElectricHDG<Dim, OrderT, OrderV, OrderG>::exportResults()
{
    M_electricModel->exportResults();
    M_heatModel->exportResults();

    double meanT = mean(_range=elements(M_mesh), _expr=idv(M_temperature))(0,0);
    auto mima = minmax(_range=elements(M_mesh), _pset=_Q<>(), _expr=idv(M_temperature));
    double miT = mima.min();
    double maT = mima.max();
    // auto power = integrate(_range=elements(M_mesh), _expr=idv(M_joule)).evaluate()(0,0);
    auto int0 = integrate(_range=markedfaces(M_mesh,M_outputMarker),_expr=inner(idv(M_current),N())).evaluate()(0,0);
    auto int1 = integrate(_range=markedfaces(M_mesh,M_inputMarker),_expr=-inner(idv(M_current),N())).evaluate()(0,0);
    double meanP = mean(_range=markedfaces(M_mesh,M_outputMarker), _expr=idv(M_potential))(0,0);
    auto mimaP = minmax(_range=markedfaces(M_mesh,M_outputMarker), _pset=_Q<>(), _expr=idv(M_potential));
    double miP = mimaP.min();
    double maP = mimaP.max();

    Feel::cout << std::setw(25) << "tempMax";
    Feel::cout << std::setw(25) << "tempMoy";
    Feel::cout << std::setw(25) << "tempMin";
    Feel::cout << std::setw(25) << "intensite0";
    Feel::cout << std::setw(25) << "intensite1";
    // Feel::cout << std::setw(25) << "power";
    Feel::cout << std::setw(25) << "potMax";
    Feel::cout << std::setw(25) << "potMoy";
    Feel::cout << std::setw(25) << "potMin";
    Feel::cout << std::endl;
    Feel::cout << std::setw(25) << maT;
    Feel::cout << std::setw(25) << meanT;
    Feel::cout << std::setw(25) << miT;
    Feel::cout << std::setw(25) << int0;
    Feel::cout << std::setw(25) << int1;
    // Feel::cout << std::setw(25) << power;
    Feel::cout << std::setw(25) << maP;
    Feel::cout << std::setw(25) << meanP;
    Feel::cout << std::setw(25) << miP;
    Feel::cout << std::endl;
    if( Environment::isMasterRank() )
    {
        std::ofstream file ( "measures.csv" );
        if( file )
        {
            file << std::setw(25) << "tempMax";
            file << std::setw(25) << "tempMoy";
            file << std::setw(25) << "tempMin";
            file << std::setw(25) << "intensite0";
            file << std::setw(25) << "intensite1";
            // file << std::setw(25) << "power";
            file << std::setw(25) << "potMax";
            file << std::setw(25) << "potMoy";
            file << std::setw(25) << "potMin";
            file << std::endl;
            file << std::setw(25) << maT;
            file << std::setw(25) << meanT;
            file << std::setw(25) << miT;
            file << std::setw(25) << int0;
            file << std::setw(25) << int1;
            // file << std::setw(25) << power;
            file << std::setw(25) << maP;
            file << std::setw(25) << meanP;
            file << std::setw(25) << miP;
            file << std::endl;
            file.close();
        }
    }
}
