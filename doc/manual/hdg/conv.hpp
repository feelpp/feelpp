#include "mixedpoisson.hpp"

namespace Feel {

inline
po::options_description
makeConvOptions()
{
    po::options_description options ( "Electro-Thermal options");
    options.add_options()
        ( "nb_refine", po::value<int>()->default_value(3), "number of refinement" )
        ;
    options.add( makeMixedPoissonOptions());
    return options;
}

inline
po::options_description
makeConvLibOptions()
{
    po::options_description options ( "Electro-Thermal lib options");
    options.add( makeMixedPoissonLibOptions());
    return options.add( feel_options());
}

template<int Dim, int Order>
class ConvergenceTest
{
public:
    using convex_type = Simplex<Dim,1>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

    using expr_scalar_type = scalar_field_expression<2>;
    using expr_vectorial_type = vector_field_expression<Dim,1,2>;

    using mixed_poisson_type = MixedPoisson<Dim,Order>;

    using flux_ptrtype = typename mixed_poisson_type::Vh_element_ptr_t;
    using potential_ptrtype = typename mixed_poisson_type::Wh_element_ptr_t;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

private:
    mesh_ptrtype M_mesh;

    mixed_poisson_type M_model;

    flux_ptrtype M_u;
    potential_ptrtype M_p;
    expr_vectorial_type M_u_exact;
    expr_scalar_type M_p_exact;

public:
    ConvergenceTest();
    void run();
};

template<int Dim, int Order>
ConvergenceTest<Dim,Order>::ConvergenceTest()
{
    M_model = mixed_poisson_type();
#if FEELPP_DIM == 2
    M_u_exact = M_model.modelProperties().functions()["u"].expressionVectorial2();
#else
    M_u_exact = M_model.modelProperties().functions()["u"].expressionVectorial3();
#endif
    M_p_exact = M_model.modelProperties().functions()["p"].expressionScalar();
}

template<int Dim, int Order>
void
ConvergenceTest<Dim,Order>::run()
{
    double h = doption("gmsh.hsize");

    std::ofstream cvg_u, cvg_p;
    cvg_u.open( "convergence_u.dat", std::ios::out | std::ios::trunc);
    cvg_p.open( "convergence_p.dat", std::ios::out | std::ios::trunc);
    cvg_u << "h" << "\t" << "nDof" << "\t" << "l2err" << std::endl;
    cvg_p << "h" << "\t" << "nDof" << "\t" << "l2err" << std::endl;

    export_ptrtype e( export_type::New( "convergence") );

    for ( int i = 0; i < ioption("nb_refine"); i++)
    {
        M_mesh = loadMesh( _mesh=new mesh_type, _h=h);
        M_model.init(M_mesh);

        auto nDofU = M_model.fluxSpace()->nDof();
        auto nDofP = M_model.potentialSpace()->nDof();
        auto u_ex = M_model.fluxSpace()->element();
        auto p_ex = M_model.potentialSpace()->element();
        auto p_mean = M_model.potentialSpace()->element();

        M_model.solve();

        M_u = M_model.fluxField();
        M_p = M_model.potentialField();

        double errU = normL2(_range=elements(M_mesh), _expr=idv(*M_u)-M_u_exact);
        double mean_p_exact = mean( elements(M_mesh), M_p_exact)(0,0);
        double mean_p = mean( elements(M_mesh), idv(*M_p))(0,0);
        double errP = normL2(_range=elements(M_mesh), _expr=idv(*M_p)-cst(mean_p)-M_p_exact+cst(mean_p_exact));

        cout << "h" << "\t" << "nDofU" << "\t" << "errU"
             << "\t" << "DofP" << "\t" << "errP" << std::endl;
        cout << h << "\t" << nDofU << "\t" << errU << "\t"
             << nDofP << "\t" << errP << std::endl;
        cvg_u << h << "\t" << nDofU << "\t" << errU << std::endl;
        cvg_p << h << "\t" << nDofP << "\t" << errP << std::endl;

        u_ex.on( elements(M_mesh), M_u_exact);
        p_ex.on( elements(M_mesh), M_p_exact-cst(mean_p_exact));
        p_mean.on( elements(M_mesh), idv(*M_p)-cst(mean_p));

        e->step(i)->setMesh(M_mesh);
        e->step(i)->add("u", *M_u);
        e->step(i)->add("p", p_mean);
        e->step(i)->add("u_ex", u_ex);
        e->step(i)->add("p_ex", p_ex);
        e->save();

        h/=2.;
    }

    cvg_u.close();
    cvg_p.close();
}

}
