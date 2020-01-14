#include "feel/feelmodels/hdg/mixedpoisson.hpp"
#include "feel/feelcore/checker.hpp"

namespace Feel {

// namespace FeelModels {

inline
po::options_description
makeConvOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "cvg.refine-nb", po::value<int>()->default_value(3), "number of refinement" )
        ( "cvg.refine-factor", po::value<double>()->default_value(2), "factor of refinement" )
        ( "cvg.p_exact", po::value<std::string>()->default_value( "0" ), "exact potential" )
        ( "cvg.cond", po::value<double>()->default_value( -1.), "conductivity" )
        ( "cvg.normal-compute", po::value<bool>()->default_value( true ), "compute the error on normal" )
        ( "cvg.normal-exact", po::value<double>()->default_value( 0 ), "exact value of int(u.n) on cvg.normal-marker" )
        ( "cvg.normal-marker", po::value<std::string>()->default_value( "ext"), "marker on which compute the error on the normal" )
        ( "cvg.use-dirichlet", po::value<bool>()->default_value( true ), "use dirichlet conditions" )
        ( "cvg.bc-dirichlet", po::value<std::vector<std::string> >(), "list of markers for Dirichlet condition ")
        ( "cvg.use-neumann", po::value<bool>()->default_value( false ), "use neumann conditions" )
        ( "cvg.bc-neumann", po::value<std::vector<std::string> >()->multitoken(), "list of markers for Neumann condition ")
        ( "cvg.use-ibc", po::value<bool>()->default_value( false ), "use ibc conditions" )
        ( "cvg.bc-ibc", po::value<std::vector<std::string> >()->multitoken(), "list of markers for ibc condition ")
        ;
    options.add( FeelModels::makeMixedPoissonOptions() );
    return options;
}

inline
po::options_description
makeConvLibOptions()
{
    po::options_description options ( "Electro-Thermal lib options");
    options.add( FeelModels::makeMixedPoissonLibOptions() );
    return options ;
}

template<int Dim, int Order,int G_Order = 1,int E_Order = 2>
class ConvergenceTest
{
public:
    using convex_type = Simplex<Dim,G_Order>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    static const uint16_type expr_order = Order+E_Order;
    using expr_scalar_type = scalar_field_expression<expr_order>;
    using expr_vectorial_type = vector_field_expression<Dim,G_Order,expr_order>;

    using mixed_poisson_type = FeelModels::MixedPoisson<Dim,Order,G_Order,E_Order>;
    using mixed_poisson_ptrtype = std::shared_ptr<mixed_poisson_type> ;

    using flux_type = typename mixed_poisson_type::Vh_element_t;
    using potential_type = typename mixed_poisson_type::Wh_element_t;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef std::shared_ptr<export_type> export_ptrtype;

    using markers_type = std::vector<std::string>;
private:
    mesh_ptrtype M_mesh;

    mixed_poisson_ptrtype M_model;

    flux_type M_u;
    potential_type M_p;
    expr_scalar_type M_p_exact;

    markers_type M_markersDirichlet;
    markers_type M_markersNeumann;
    markers_type M_markersIBC;

    std::map<std::string, std::vector<double> > M_timers;

public:
    ConvergenceTest();
    void assembleExact();
    int run();
    void exportTimers();
};

template<int Dim, int Order, int G_Order, int E_Order>
ConvergenceTest<Dim,Order,G_Order,E_Order>::ConvergenceTest()
{
    M_model = mixed_poisson_type::New();
    auto solution = expr<expr_order>( checker().solution(), "solution");
    M_p_exact = checker().check() ? solution : expr<expr_order>(soption("cvg.p_exact"));
    Feel::cout << "using conductivity: " << doption("cvg.cond") << std::endl;
    Feel::cout << "using p exact     : " << M_p_exact << std::endl;
    Feel::cout << "grad_p_exact      : " << grad<Dim,expr_order>(M_p_exact) << std::endl;
    Feel::cout << "laplacian(p_exact): " << laplacian(M_p_exact) << std::endl;
    if( boption("cvg.normal-compute") )
    {
        Feel::cout << "compute normale on " << soption("cvg.normal-marker") << std::endl;
        Feel::cout << "int(u.n)          : " << doption("cvg.normal-exact") << std::endl;
    }

    if ( boption("cvg.use-dirichlet") )
        M_markersDirichlet = vsoption("cvg.bc-dirichlet");
    if ( boption("cvg.use-neumann") )
        M_markersNeumann = vsoption("cvg.bc-neumann");
    if ( boption("cvg.use-ibc") )
        M_markersIBC = vsoption("cvg.bc-ibc");
    M_model->setIBCList(M_markersIBC);

    auto repo = boost::format( "conv_mixedpoisson/%1%D/P%2%/N%3%/E%4%/Dir%5%Neu%6%Ibc%7%/" )
        % Dim
        % Order
        % G_Order
        % E_Order
        % M_markersDirichlet.size()
        % M_markersNeumann.size()
        % M_markersIBC.size();
    Environment::changeRepository( repo );
    Feel::cout << "Results are now stored in " << tc::red << repo << tc::reset << std::endl;
}

// assemble the matrices with the exact solutions
// using options p_exact, cond and bc-<type>
template<int Dim, int Order, int G_Order, int E_Order>
void
ConvergenceTest<Dim,Order,G_Order,E_Order>::assembleExact()
{
    tic();
    tic();
    M_model->assembleCstPart();
    M_timers["asbCstPart"].push_back(toc("assembleCstPart"));

    auto cond = cst(doption("cvg.cond"));
    tic();
    M_model->updateConductivityTerm(cond, "");
    M_timers["asbCndPart"].push_back(toc("assembleConductivityPart"));

    auto grad_p_exact = grad<Dim,expr_order>(M_p_exact);
    auto u_exact = -cond*trans(grad_p_exact);
    auto f = -cond*laplacian(M_p_exact);
    tic();
    M_model->assemblePotentialRHS(f, "");
    M_timers["asbPotRHS"].push_back(toc("assemblePotentialRHS"));

    tic();
    for( auto const& marker : M_markersDirichlet )
    {
        if( M_mesh->hasFaceMarker(marker))
        {
            tic();
            M_model->assembleDirichlet(marker);
            M_model->assembleRhsDirichlet(M_p_exact, marker);
            std::string s = (boost::format("asbDir_%1%") % marker).str();
            M_timers[s].push_back(toc(s));
            Feel::cout << "add Dirichlet on " << marker << std::endl;
        }
        else
        {
            LOG(FATAL) << "marker for Dirichlet condition " << marker << "does not exist";
        }
    }
    for( auto const& marker : M_markersNeumann )
    {
        if( M_mesh->hasFaceMarker(marker))
        {
            tic();
            M_model->assembleNeumann( marker);
            M_model->assembleRhsNeumann( -cond*grad_p_exact*N(), marker);
            std::string s = (boost::format("asbNeu_%1%") % marker).str();
            M_timers[s].push_back(toc(s));
            Feel::cout << "add Neumann on " << marker << std::endl;
        }
        else
        {
            LOG(FATAL) << "marker for Neumann condition " << marker << "does not exist";
        }
    }
    int i = 0;
    for( auto const& marker : M_markersIBC )
    {
        if( M_mesh->hasFaceMarker(marker))
        {
            double intjn = integrate(_range=markedfaces(M_mesh,marker), _expr=inner(u_exact,N())).evaluate()(0,0);
            // double intjn = -std::log(2)/(2.*boost::math::constants::pi<double>());
            tic();
            M_model->assembleRhsIBC(i++, marker, intjn);
            std::string s = (boost::format("asbIBC_%1%") % marker).str();
            M_timers[s].push_back(toc(s));
            Feel::cout << "add ibc on " << marker << " with value of " << intjn << std::endl;
        }
        else
        {
            LOG(FATAL) << "marker for IBC condition " << marker << "does not exist";
        }
    }
    M_timers["asbBC"].push_back(toc("assembleBC"));
    M_timers["cvgAsb"].push_back(toc("convergence assembly"));
}

template<int Dim, int Order, int G_Order, int E_Order>
int
ConvergenceTest<Dim,Order,G_Order,E_Order>::run()
{
    double h = doption("gmsh.hsize");

    auto cond = cst(doption("cvg.cond"));
    auto grad_p_exact = grad<Dim,expr_order>(M_p_exact);
    auto u_exact = -cond*trans(grad_p_exact);

    std::ofstream cvg_u, cvg_p, cvg_n;
    if( Environment::isMasterRank() )
    {
        cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
        cvg_p.open( "convergence_p.dat", std::ios::out | std::ios::trunc);
        if( boption("cvg.normal-compute") )
            cvg_n.open( "../../../../../convergence_n.dat", std::ios::out | std::ios::app);
    }
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    boost::format fmterOut;
    if( boption("cvg.normal-compute") )
        fmterOut = boost::format("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5% %|70t|%6%\n");
    else
        fmterOut = boost::format("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    if( Environment::isMasterRank() )
    {
        cvg_u << fmter % "#h" % "nDof" % "l2err";
        cvg_p << fmter % "#h" % "nDof" % "l2err";
        if( boption("cvg.normal-compute") )
            cvg_n << fmter % "#P" % "G" % "err";
    }


    export_ptrtype e( export_type::New( "convergence") );

    for ( int i = 0; i < ioption("cvg.refine-nb"); i++)
    {
        M_mesh = loadMesh( _mesh=new mesh_type, _h=h);
        tic();
        M_model->init(M_mesh);
        auto nDofU = M_model->fluxSpace()->nDof();
        auto nDofP = M_model->potentialSpace()->nDof();
        M_timers["nDofU"].push_back(nDofU);
        M_timers["nDofP"].push_back(nDofP);
        M_timers["initModel"].push_back(toc("initModel"));

        auto u_ex = M_model->fluxSpace()->element();
        auto p_ex = M_model->potentialSpace()->element();
        auto p_ex_mean = M_model->potentialSpace()->element();
        auto p_mean = M_model->potentialSpace()->element();

        this->assembleExact();
        tic();
        M_model->solve();
        M_timers["solve"].push_back(toc("solve"));

        M_u = M_model->fluxField();
        M_p = M_model->potentialField();

        double errU = normL2(_range=elements(M_mesh), _expr=idv(M_u)-u_exact, _quad=_Q<expr_order>());
        double mean_p_exact = mean( elements(M_mesh), M_p_exact)(0,0);
        double mean_p = mean( elements(M_mesh), idv(M_p))(0,0);
        double errP = normL2(_range=elements(M_mesh), _expr=idv(M_p)-cst(mean_p)-M_p_exact+cst(mean_p_exact), _quad=_Q<expr_order>());

        if (boption("cvg.normal-compute") )
        {
            double un_exact = doption("cvg.normal-exact");
            double un_approx = integrate(_range=markedfaces(M_mesh, soption("cvg.normal-marker")),
                                         _expr=inner(idv(M_u),N()),
                                         _quad=_Q<expr_order>() ).evaluate()(0,0);
            double errN = math::abs(un_exact - un_approx);
            cout << std::setprecision(14) << "geo : un_exact=" << un_exact << " un_approx=" << un_approx << std::endl;
            // double length = integrate(_range=markedfaces(M_mesh, soption("cvg.normal-marker")),
            //                           _expr=expr(soption("functions.f")) ).evaluate()(0,0);
            // double errN = math::abs(length - doption("cvg.normal-exact") );
            if( Environment::isMasterRank() )
            {
                cout << fmterOut % "h" % "nDofU" % "errU" % "nDofP" % "errP" % "errN";
                cout << fmterOut % h % nDofU % errU % nDofP % errP % errN;
                cvg_n << fmter % Order % G_Order % errN;
            }
        }
        else
        {
            if( Environment::isMasterRank() )
            {
                cout << fmterOut % "h" % "nDofU" % "errU" % "nDofP" % "errP";
                cout << fmterOut % h % nDofU % errU % nDofP % errP;
            }
        }
        if( Environment::isMasterRank() )
        {
            cvg_u << fmter % h % nDofU % errU;
            cvg_p << fmter % h % nDofP % errP;
        }

        u_ex.on( _range=elements(M_mesh), _expr=u_exact);
        p_ex.on( _range=elements(M_mesh), _expr=M_p_exact);
        p_ex_mean.on( _range=elements(M_mesh), _expr=M_p_exact-cst(mean_p_exact));
        p_mean.on( _range=elements(M_mesh), _expr=idv(M_p)-cst(mean_p));

        e->step(i)->setMesh(M_mesh);
        e->step(i)->add("u", M_u);
        e->step(i)->add("p", M_p);
        e->step(i)->add("p_mean", p_mean);
        e->step(i)->add("p_ex_mean", p_ex_mean);
        e->step(i)->add("u_ex", u_ex);
        e->step(i)->add("p_ex", p_ex);
        e->save();

        h /= doption("cvg.refine-factor");
        toc("convergence loop");
    }

    this->exportTimers();

    if( Environment::isMasterRank() )
    {
        cvg_u.close();
        cvg_p.close();
        if( boption("cvg.normal-compute") )
            cvg_n.close();
    }

    auto norms = [=]( std::string const& solution ) ->std::map<std::string,double>
        {
            tic();
            double errU = normL2(_range=elements(M_mesh), _expr=idv(M_u)-u_exact, _quad=_Q<expr_order>());
            toc("u error norm");
            tic();
            double mean_p_exact = mean( elements(M_mesh), M_p_exact)(0,0);
            double mean_p = mean( elements(M_mesh), idv(M_p))(0,0);
            double errP = normL2(_range=elements(M_mesh), _expr=idv(M_p)-cst(mean_p)-M_p_exact+cst(mean_p_exact), _quad=_Q<expr_order>());
            toc("p error norm");
            return { { "u", errU }, {  "p", errP } };
        };
    int status = checker().runOnce( norms, rate::hp( M_mesh->hMax(), M_u.functionSpace()->fe()->order() ) );

    return !status;
}

template<int Dim, int Order, int G_Order, int E_Order>
void
ConvergenceTest<Dim,Order,G_Order,E_Order>::exportTimers()
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
        for( int i = 0; i < ioption("cvg.refine-nb"); ++i )
        {
            for( auto const& pair : M_timers )
                fmt % pair.second[i];
            timers << fmt << std::endl;
        }
        timers.close();
    }
}

// } // end namespace FeelModels
} // end namespace Feel
