#include "poissonCRB.hpp"

using namespace Feel;

int main( int argc, char** argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions()
                     .add(crbOptions())
                     .add(crbSEROptions())
                     .add(eimOptions())
                     .add(podOptions())
                     .add(backend_options("backend-primal"))
                     .add(backend_options("backend-dual"))
                     .add(backend_options("backend-l2"))
                     .add(bdf_options("PoissonCRB")) );

    using rb_model_type = Poisson;
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
    crb_ptrtype crb = boost::make_shared<crb_type>("poissonmodel_crb", crbModel);

    // offline
    crb->offline();

    // parameter from option
    double kappa = doption("poisson.kappa");
    double flux = doption("poisson.flux");
    auto paramSpace = crbModel->parameterSpace();
    auto muMin = paramSpace->min();
    auto muMax = paramSpace->max();
    auto mu = paramSpace->element();
    mu << kappa, flux;
    Feel::cout << "using parameter:" << std::endl << mu << std::endl;

    // online
    int N = ioption("poisson.N");
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
    auto e = Exporter<rb_model_type::mesh_type>::New( "poisson" );
    e->setMesh( model->mesh() );
    e->add( "u", u );
    e->save();

    return 0;
}
