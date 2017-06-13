#include "nlthermoelectric.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="nl-thermoelectric"),
                     _desc=makeOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    using rb_model_type = Thermoelectric;
    using rb_model_ptrtype = boost::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = boost::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = boost::shared_ptr<crb_type>;
    using wn_type = typename crb_type::wn_type;
    using vectorN_type = Eigen::VectorXd;
    using export_vector_wn_type = typename crb_type::export_vector_wn_type;

    // init
    rb_model_ptrtype model = boost::make_shared<rb_model_type>();
    crb_model_ptrtype crbModel = boost::make_shared<crb_model_type>(model);
    crb_ptrtype crb = boost::make_shared<crb_type>("thermoelectric_crb", crbModel);

    // offline
    crb->offline();

    // parameter from option
    double sigma = doption("thermoelectric.sigma");
    double alpha = doption("thermoelectric.alpha");
    double L = doption("thermoelectric.L");
    double current = doption("thermoelectric.current");
    double h = doption("thermoelectric.h");
    double Tw = doption("thermoelectric.Tw");
    auto paramSpace = crbModel->parameterSpace();
    auto mu = paramSpace->element();
    mu.setParameterNamed("sigma", sigma);
    mu.setParameterNamed("current", current);
    mu.setParameterNamed("alpha", alpha);
    mu.setParameterNamed("k", L);
    mu.setParameterNamed("h", h);
    mu.setParameterNamed("Tw", Tw);
    Feel::cout << "using parameter:" << std::endl << mu << std::endl;

    // online
    int N = ioption("thermoelectric.N");
    if( N > crb->dimension() || N < 1 )
        N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);
    crb->fixedPointPrimal(N, mu, uNs, uNolds, outputs);
    vectorN_type uN = uNs[0];
    Feel::cout << "uN:" << std::endl << uN << std::endl;

    // export basis
    wn_type WN = crb->wn();
    std::vector<wn_type> WN_vec = std::vector<wn_type>();
    WN_vec.push_back(WN);
    std::vector<std::string> name_vec = std::vector<std::string>(1, "primal");
    export_vector_wn_type exportWn = boost::make_tuple(WN_vec, name_vec);
    crb->exportBasisFunctions(exportWn);

    // export
    auto u = crb->expansion( uN, uN.size(), WN );
    auto V = u.template element<0>();
    auto T = u.template element<1>();
    auto e = Exporter<rb_model_type::mesh_type>::New( "thermoelectric" );
    e->setMesh( model->mesh() );
    e->add( "V", V );
    e->add( "T", T );

    // export EIM basis
    auto eimG = model->scalarDiscontinuousEim()[0];
    int M = eimG->mMax();
    auto betaEim = eimG->beta(mu);
    auto jEim = eimG->functionSpace()->element();
    for( int m = 0; m < M; ++m )
    {
        e->add( (boost::format("q_%1%" ) % m).str(), eimG->q(m) );
        jEim += vf::project( _range=elements(model->mesh()), _space=eimG->functionSpace(),
                             _expr=betaEim(m)*idv(eimG->q(m)) );
    }
    e->add( "j2Eim", jEim );

    auto j = vf::project( _range=elements(model->mesh()), _space=eimG->functionSpace(),
                          _expr=mu.parameterNamed("sigma")*inner(gradv(V)) );
    e->add( "j2", j );

    // export FE
    auto VTFE = model->solve(mu);
    auto VFE = VTFE.template element<0>();
    auto TFE = VTFE.template element<1>();
    e->add( "VFE", VFE );
    e->add( "TFE", TFE );

    e->save();

    return 0;
}
