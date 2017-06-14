#include <feel/feel.hpp>
#include <feel/feelopt/nlopt.hpp>

#include "biotsavart.hpp"
#include "thermoelectric.hpp"

int iter=0;

int main(int argc, char**argv )
{
    using namespace Feel;

    po::options_description nloptoptions( "NLOpt options" );
    nloptoptions.add_options()
        ( "Tcritic", po::value<double>()->default_value( 400 ), "critical value for the temperature")
        ( "nlopt.algo", po::value<std::string>()->default_value( "LN_NEWUOA" ), "NLOPT algorithm [LN_NEWUOA,LD_LBFGS]" )
        ;

    nloptoptions.add(biotsavartOptions());
    nloptoptions.add(crbOptions());
    nloptoptions.add(crbSEROptions());
    nloptoptions.add(eimOptions());
    nloptoptions.add(podOptions());
    nloptoptions.add(backend_options("backend-primal"));
    nloptoptions.add(backend_options("backend-dual"));
    nloptoptions.add(backend_options("backend-l2"));
    nloptoptions.add(bdf_options("ThermoElectricCRB"));
    nloptoptions.add(makeOptions());

    Environment env( _argc=argc, _argv=argv,
                     _desc=nloptoptions,
                     _about=about(_name="biotsavart_nlopt",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    // List of NLOPT algorithm
    const boost::unordered_map< const std::string, ::nlopt::algorithm >& authAlgo = boost::assign::map_list_of
        ("LN_NEWUOA", ::nlopt::LN_NEWUOA )
        ("LN_COBYLA", ::nlopt::LN_COBYLA )
        ("LN_BOBYQA", ::nlopt::LN_BOBYQA )
        ("LD_LBFGS", ::nlopt::LD_LBFGS )
        ("LD_MMA", ::nlopt::LD_MMA )
        ("LD_SLSQP", ::nlopt::LD_SLSQP );

    BiotSavartCRB<Thermoelectric> BS = BiotSavartCRB<Thermoelectric>();
    BS.offline();
    int N = BS.nbParameters();

    auto salgo = soption("nlopt.algo");
    Feel::cout << "NLOP algorithm: " << salgo << "\n";
    auto algo = authAlgo.at(salgo);
    ::nlopt::opt opt( algo, N );

    // lower and upper bounds
    auto muMin = BS.parameterSpace()->min();
    auto muMax = BS.parameterSpace()->max();
    std::vector<double> x(N);
    for( int i = 0; i < N; ++i ) x[i] = muMin(i);

    std::vector<double> lb(muMin.data(), muMin.data()+muMin.size());
    std::vector<double> ub(muMax.data(), muMax.data()+muMax.size());

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    // stopping criteria
    opt.set_maxeval( ioption("nlopt.maxeval") );
    opt.set_xtol_rel( doption("nlopt.xtol_rel") );
    opt.set_ftol_rel( doption("nlopt.ftol_rel") );
    opt.set_xtol_abs( doption("nlopt.xtol_abs") );
    opt.set_ftol_abs( doption("nlopt.ftol_abs") );

    // Objective function.
    // We could use the vector version here (see ::NLOPT::vfunc).
    auto myfunc = []( unsigned n, const double *x, double *grad, void *my_func_data )->double
    {
        // iter++;
        BiotSavartCRB<Thermoelectric>* BS = reinterpret_cast<BiotSavartCRB<Thermoelectric>*>(my_func_data);
        int N = BS->nbParameters();
        auto mu = BS->newParameter();
        for( int i = 0; i < N; ++i ) mu(i) = x[i];
        BS->online(mu);
        auto B = BS->magneticFlux();
        return B.max();
    };

    opt.set_max_objective( myfunc, &BS );

    // inequality constraints
    auto myconstraint = []( unsigned n, const double *x, double *grad, void *data)->double
    {
        BiotSavartCRB<Thermoelectric>* BS = reinterpret_cast<BiotSavartCRB<Thermoelectric>*>(data);
        auto VT = BS->potentialTemperature();
        auto T = VT.template element<1>();
        return T.max() - doption("Tcritic");
    };

    opt.add_inequality_constraint(myconstraint, &BS, 1e-8);

    // optimization
    double minf;

    tic();
    ::nlopt::result result = opt.optimize(x, minf);
    double optiTime = toc("optimization", false);
    Feel::cout << iter << " iterations in " << optiTime << " (" << optiTime/iter << "/iter)" << std::endl;

    // export
    auto mu = BS.newParameter();
    for( int i = 0; i < N; ++i ) mu(i) = x[i];

    BS.exportResults();

    switch( result )
    {
        case ::nlopt::FAILURE:
            Feel::cerr << "NLOPT Generic Failure!" << "\n";
            break;
        case ::nlopt::INVALID_ARGS:
            Feel::cerr << "NLOPT Invalid arguments!" << "\n";
            break;
        case ::nlopt::OUT_OF_MEMORY:
            Feel::cerr << "NLOPT Out of memory!" << "\n";
            break;
        case ::nlopt::ROUNDOFF_LIMITED:
            Feel::cerr << "NLOPT Roundoff limited!" << "\n";
            break;
        case ::nlopt::FORCED_STOP:
            Feel::cerr << "NLOPT Forced stop!" << "\n";
            break;
        case::nlopt::SUCCESS:
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
        case ::nlopt::STOPVAL_REACHED:
            Feel::cout << "NLOPT Stop value reached!" << "\n";
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
        case ::nlopt::FTOL_REACHED:
            Feel::cout << "NLOPT ftol reached!" << "\n";
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
        case ::nlopt::XTOL_REACHED:
            Feel::cout << "NLOPT xtol reached!" << "\n";
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
        case ::nlopt::MAXEVAL_REACHED:
            Feel::cout << "NLOPT Maximum number of evaluation reached!" << "\n";
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
        case ::nlopt::MAXTIME_REACHED:
            Feel::cout << "NLOPT Maximum time reached" << "\n";
            Feel::cout << "NLOPT coefficient found! (status " << result << ") \n"
                       << mu << "\n"
                       << "Evaluation number: " << iter << "\n";
            break;
    }

    return 0;
}
