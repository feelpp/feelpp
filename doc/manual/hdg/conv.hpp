#include "mixedpoisson2.hpp"

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
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    static const uint16_type expr_order = Order+E_Order;
    using expr_scalar_type = scalar_field_expression<expr_order>;
    using expr_vectorial_type = vector_field_expression<Dim,G_Order,expr_order>;

    using mixed_poisson_type = FeelModels::MixedPoisson<Dim,Order,G_Order,E_Order>;
    using mixed_poisson_ptrtype = boost::shared_ptr<mixed_poisson_type> ;

    using flux_type = typename mixed_poisson_type::Vh_element_t;
    using potential_type = typename mixed_poisson_type::Wh_element_t;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

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

public:
    ConvergenceTest();
    void assembleExact();
    void run();
};

template<int Dim, int Order, int G_Order, int E_Order>
ConvergenceTest<Dim,Order,G_Order,E_Order>::ConvergenceTest()
{
    M_model = mixed_poisson_type::New("mixedpoisson");
    M_p_exact = expr<expr_order>(soption("cvg.p_exact"));
    Feel::cout << "using conductivity: " << doption("cvg.cond") << std::endl;
    Feel::cout << "using p exact     : " << M_p_exact << std::endl;
    Feel::cout << "grad_p_exact      : " << grad<Dim,expr_order>(M_p_exact) << std::endl;
    Feel::cout << "laplacian(p_exact): " << laplacian(M_p_exact) << std::endl;

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
    M_model->assembleCstPart();

    auto cond = cst(doption("cvg.cond"));
    M_model->updateConductivityTerm(cond, "");

    auto grad_p_exact = grad<Dim,expr_order>(M_p_exact);
    auto u_exact = -cond*trans(grad_p_exact);
    auto f = -cond*laplacian(M_p_exact);
    M_model->assemblePotentialRHS(f, "");

    for( auto const& marker : M_markersDirichlet )
    {
        if( M_mesh->hasFaceMarker(marker))
        {
            M_model->assembleDirichlet(M_p_exact, marker);
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
            M_model->assembleNeumann( -cond*grad_p_exact*N(), marker);
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
            M_model->assembleIBC(i++, marker, intjn);
            Feel::cout << "add ibc on " << marker << " with value of " << intjn << std::endl;
        }
        else
        {
            LOG(FATAL) << "marker for IBC condition " << marker << "does not exist";
        }
    }
}

template<int Dim, int Order, int G_Order, int E_Order>
void
ConvergenceTest<Dim,Order,G_Order,E_Order>::run()
{
    double h = doption("gmsh.hsize");

    auto cond = cst(doption("cvg.cond"));
    auto grad_p_exact = grad<Dim,expr_order>(M_p_exact);
    auto u_exact = -cond*trans(grad_p_exact);

    std::ofstream cvg_u, cvg_p;
    cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
    cvg_p.open( "convergence_p.dat", std::ios::out | std::ios::trunc);
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    boost::format fmterOut("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    cvg_u << fmter % "#h" % "nDof" % "l2err";
    cvg_p << fmter % "#h" % "nDof" % "l2err";

    export_ptrtype e( export_type::New( "convergence") );

    for ( int i = 0; i < ioption("cvg.refine-nb"); i++)
    {
        M_mesh = loadMesh( _mesh=new mesh_type, _h=h);
        M_model -> init(M_mesh);

        auto nDofU = M_model->fluxSpace()->nDof();
        auto nDofP = M_model->potentialSpace()->nDof();
        auto u_ex = M_model->fluxSpace()->element();
        auto p_ex = M_model->potentialSpace()->element();
        auto p_ex_mean = M_model->potentialSpace()->element();
        auto p_mean = M_model->potentialSpace()->element();

        this->assembleExact();
        M_model->solve();

        M_u = M_model->fluxField();
        M_p = M_model->potentialField();

        double errU = normL2(_range=elements(M_mesh), _expr=idv(M_u)-u_exact, _quad=_Q<expr_order>());
        double mean_p_exact = mean( elements(M_mesh), M_p_exact)(0,0);
        double mean_p = mean( elements(M_mesh), idv(M_p))(0,0);
        double errP = normL2(_range=elements(M_mesh), _expr=idv(M_p)-cst(mean_p)-M_p_exact+cst(mean_p_exact), _quad=_Q<expr_order>());

        cout << fmterOut % "h" % "nDofU" % "errU" % "nDofP" % "errP";
        cout << fmterOut % h % nDofU % errU % nDofP % errP;
        cvg_u << fmter % h % nDofU % errU;
        cvg_p << fmter % h % nDofP % errP;

        u_ex.on( elements(M_mesh), u_exact);
        p_ex.on( elements(M_mesh), M_p_exact);
        p_ex_mean.on( elements(M_mesh), M_p_exact-cst(mean_p_exact));
        p_mean.on( elements(M_mesh), idv(M_p)-cst(mean_p));

        e->step(i)->setMesh(M_mesh);
        e->step(i)->add("u", M_u);
        e->step(i)->add("p", M_p);
        e->step(i)->add("p_mean", p_mean);
        e->step(i)->add("u_ex", u_ex);
        e->step(i)->add("p_ex", p_ex);
        e->save();

        h /= doption("cvg.refine-factor");
    }

    cvg_u.close();
    cvg_p.close();
}

// } // end namespace FeelModels
} // end namespace Feel
