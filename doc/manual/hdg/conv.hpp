#include "mixedpoisson2.hpp"

namespace Feel {

// namespace FeelModels {

inline
po::options_description
makeConvOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "refine-nb", po::value<int>()->default_value(3), "number of refinement" )
        ( "refine-factor", po::value<double>()->default_value(2), "factor of refinement" )
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

template<int Dim, int Order,int G_Order = 1>
class ConvergenceTest
{
public:
    using convex_type = Simplex<Dim,G_Order>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    static const uint16_type expr_order = Order+2;
    using expr_scalar_type = scalar_field_expression<expr_order>;
    using expr_vectorial_type = vector_field_expression<Dim,G_Order,expr_order>;

    using mixed_poisson_type = FeelModels::MixedPoisson<Dim,Order,G_Order>;
    using mixed_poisson_ptrtype = boost::shared_ptr<mixed_poisson_type> ;

    using flux_type = typename mixed_poisson_type::Vh_element_t;
    using potential_type = typename mixed_poisson_type::Wh_element_t;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

private:
    mesh_ptrtype M_mesh;

    mixed_poisson_ptrtype M_model;

    flux_type M_u;
    potential_type M_p;
    expr_vectorial_type M_u_exact;
    expr_scalar_type M_p_exact;

public:
    ConvergenceTest();
    void run();
};

template<int Dim, int Order, int G_Order>
ConvergenceTest<Dim,Order,G_Order>::ConvergenceTest()
{
    M_model = mixed_poisson_type::New("mixedpoisson");
    auto u_expr = M_model->modelProperties().functions()["u"].expressionString();
    auto p_expr = M_model->modelProperties().functions()["p"].expressionString();
    M_u_exact = expr<FEELPP_DIM,1,expr_order>(u_expr);
    M_p_exact = expr<expr_order>(p_expr);
// #if FEELPP_DIM == 2
//     M_u_exact = M_model -> modelProperties().functions()["u"].expressionVectorial2();
// #else
//     M_u_exact = M_model -> modelProperties().functions()["u"].expressionVectorial3();
// #endif
//     M_p_exact = M_model -> modelProperties().functions()["p"].expressionScalar();
}

template<int Dim, int Order, int G_Order>
void
ConvergenceTest<Dim,Order,G_Order>::run()
{
    double h = doption("gmsh.hsize");

    std::ofstream cvg_u, cvg_p;
    std::string fileU =  (boost::format( "convergence_u_%1%DP%2%N%3%.dat" ) % Dim % Order % G_Order ).str();
    std::string fileP =  (boost::format( "convergence_p_%1%DP%2%N%3%.dat" ) % Dim % Order % G_Order ).str();
    boost::format fmter("%1% %|14t|%2% %|28t|%3%\n");
    boost::format fmterOut("%1% %|14t|%2% %|28t|%3% %|42t|%4% %|56t|%5%\n");
    cvg_u.open( fileU, std::ios::out | std::ios::trunc);
    cvg_p.open( fileP, std::ios::out | std::ios::trunc);
    cvg_u << fmter % "#h" % "nDof" % "l2err";
    cvg_p << fmter % "#h" % "nDof" % "l2err";

    export_ptrtype e( export_type::New( "convergence") );

    for ( int i = 0; i < ioption("refine-nb"); i++)
    {
        M_mesh = loadMesh( _mesh=new mesh_type, _h=h);
        M_model -> init(M_mesh);

        auto nDofU = M_model->fluxSpace()->nDof();
        auto nDofP = M_model->potentialSpace()->nDof();
        auto u_ex = M_model->fluxSpace()->element();
        auto p_ex = M_model->potentialSpace()->element();
        auto p_mean = M_model->potentialSpace()->element();

        M_model->solve();

        M_u = M_model->fluxField();
        M_p = M_model->potentialField();

        double errU = normL2(_range=elements(M_mesh), _expr=idv(M_u)-M_u_exact, _quad=_Q<expr_order>());
        double mean_p_exact = mean( elements(M_mesh), M_p_exact)(0,0);
        double mean_p = mean( elements(M_mesh), idv(M_p))(0,0);
        double errP = normL2(_range=elements(M_mesh), _expr=idv(M_p)-cst(mean_p)-M_p_exact+cst(mean_p_exact), _quad=_Q<expr_order>());

        cout << fmterOut % "h" % "nDofU" % "errU" % "nDofP" % "errP";
        cout << fmterOut % h % nDofU % errU % nDofP % errP;
        cvg_u << fmter % h % nDofU % errU;
        cvg_p << fmter % h % nDofP % errP;

        u_ex.on( elements(M_mesh), M_u_exact);
        p_ex.on( elements(M_mesh), M_p_exact-cst(mean_p_exact));
        p_mean.on( elements(M_mesh), idv(M_p)-cst(mean_p));

        e->step(i)->setMesh(M_mesh);
        e->step(i)->add("u", M_u);
        e->step(i)->add("p", M_p);
        e->step(i)->add("p_mean", p_mean);
        e->step(i)->add("u_ex", u_ex);
        e->step(i)->add("p_ex", p_ex);
        e->save();

        h /= doption("refine-factor");
    }

    cvg_u.close();
    cvg_p.close();
}

// } // end namespace FeelModels
} // end namespace Feel
