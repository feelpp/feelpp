#include <fmt/core.h>
#include <fmt/compile.h>
#include <feel/feelmor/toolboxmor.hpp>
#include <feel/feelmodels/heat/heat.hpp>

using namespace Feel;

void writeErrors(std::ostream& out, std::vector<std::vector<double> > const& err)
{
    if( out && Environment::isMasterRank() )
    {
        int N = err.size();
        int size = err[0].size();
        out << std::setw(5) << "N";
        for(int i = 0; i < size; ++i)
            out << std::setw(24) << "mu_" << i;
        out << std::endl;
        for(int n = 0; n < N; ++n)
        {
            out << std::setw(5) << n+1;
            for(int i = 0; i < size; ++i)
                out << std::setw(25) << err[n][i];
            out << std::endl;
        }
    }
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>
computeStats(std::vector<std::vector<double>> errs)
{
    int N = errs.size();
    std::vector<double> min(N), max(N), mean(N), stdev(N);
    if( N == 0 )
        return std::make_tuple(min, max, mean, stdev);
    int size = errs[0].size();
    for(int n = 0; n < N; ++n)
    {
        min[n] = *std::min_element(errs[n].begin(), errs[n].end());
        max[n] = *std::max_element(errs[n].begin(), errs[n].end());
        double s = std::accumulate(errs[n].begin(), errs[n].end(), 0.0);
        mean[n] = s/size;
        double accum = std::accumulate(errs[n].begin(), errs[n].end(), 0.0,
                                       [s,size](double a, double b) {
                                           return a + (b-s/size)*(b-s/size);
                                       });
        stdev[n] = accum/size;
    }
    return std::make_tuple(min, max, mean, stdev);
}

void writeStats(std::ostream& out, std::vector<double> min, std::vector<double> mean, std::vector<double> max, std::vector<double> stdev)
{
    if( out && Environment::isMasterRank() )
    {
        int N = min.size();
        out << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
            << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
        for(int n = 0; n < N; ++n)
            out << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                << std::setw(25) << mean[n] << std::setw(25) << stdev[n] << std::endl;
    }
}

std::vector<double> computeOutputs(std::vector<std::vector<std::vector<Eigen::VectorXd>>> const& Lqm_pr, 
                                   std::vector<std::vector<std::vector<double>>> const& beta,
                                   Eigen::VectorXd const& uN)
{
    std::vector<double> outputs;
    for(int i = 1; i < beta.size(); ++i )
    {
        Eigen::VectorXd F_pr = Eigen::VectorXd::Zero(uN.size());
        for(int q = 0; q < beta[i].size(); ++q)
            for(int m = 0; m < beta[i][q].size(); ++m)
                F_pr += beta[i][q][m]*Lqm_pr[i-1][q][m].head(uN.size());
        outputs.push_back(F_pr.dot(uN));
    }
    return outputs;
}

template<int Dim, int Order>
int runSimulation()
{
    using convex_type = Simplex<Dim>;
    using base_type = Lagrange<Order, Scalar, Continuous, PointSetFekete>;
    using heat_tb_type = FeelModels::Heat<convex_type, base_type>;
    using heat_tb_ptrtype = std::shared_ptr<heat_tb_type>;
    using space_type = typename heat_tb_type::space_temperature_type;

    using rb_model_type = ToolboxMor<space_type>;
    using rb_model_ptrtype = std::shared_ptr<rb_model_type>;
    using crb_model_type = CRBModel<rb_model_type>;
    using crb_model_ptrtype = std::shared_ptr<crb_model_type>;
    using crb_type = CRB<crb_model_type>;
    using crb_ptrtype = std::shared_ptr<crb_type>;
    using wn_type = typename crb_type::wn_type;
    using vectorN_type = Eigen::VectorXd;
    using export_vector_wn_type = typename crb_type::export_vector_wn_type;
    using mesh_type = typename rb_model_type::mesh_type;
    using mesh_ptrtype = typename rb_model_type::mesh_ptrtype;
    using parameter_type = typename rb_model_type::parameter_type;
    using sampling_type = typename crb_type::sampling_type;
    using sampling_ptrtype = std::shared_ptr<sampling_type>;

    using deim_function_type = typename rb_model_type::deim_function_type;
    using mdeim_function_type = typename rb_model_type::mdeim_function_type;

    auto heatBox = heat_tb_type::New(_prefix="heat");
    heatBox->init();
    heatBox->printAndSaveInfo();

    rb_model_ptrtype model = std::make_shared<rb_model_type>(soption("toolboxmor.name"));
    model->setFunctionSpaces(heatBox->spaceTemperature());
    auto heatBoxModel = DeimMorModelToolbox<heat_tb_type>::New(heatBox);
    model->initToolbox(heatBoxModel);

    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model);
    crb_ptrtype crb = crb_type::New(soption("toolboxmor.name"), crbModel, crb::stage::offline);

    tic();
    crb->offline();
    toc("offline");

    if( !boption("toolboxmor.do-cvg") )
        return 0;

    // convergence study
    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outs(timeSteps, 0);

    auto allOutputs = model->modelProperties()->outputs();
    auto outputs = allOutputs.ofTypes({"integrate","mean","sensor","point"});
    auto Fqm = model->getFqm();
    std::vector<std::vector<std::vector<Eigen::VectorXd>>> Lqm_pr(Fqm.size()-1);
    for( int i = 1; i < Fqm.size(); ++i)
    {
        Lqm_pr[i-1].resize(Fqm[i].size());
        for( int q = 0; q < Fqm[i].size(); ++q )
        {
            Lqm_pr[i-1][q].resize(Fqm[i][q].size());
            for( int m = 0; m < Fqm[i][q].size(); ++m )
            {
                Lqm_pr[i-1][q][m] = Eigen::VectorXd(N);
                for( int n = 0; n < N; ++n )
                    Lqm_pr[i-1][q][m](n) = crbModel->Fqm( i, q, m, crbModel->rBFunctionSpace()->primalBasisElement(n) );
            }
        }
    }

    sampling_ptrtype sampling( new sampling_type( crbModel->parameterSpace() ) );
    int size = ioption("toolboxmor.sampling-size");
    sampling->clear();
    sampling->randomize( size, true );

    std::vector<std::vector<double> > errs(N, std::vector<double>(size)), errsRel(N, std::vector<double>(size));
    std::vector<std::vector<std::vector<double>>> errsOutput(Lqm_pr.size(), std::vector<std::vector<double>>(N, std::vector<double>(size)));
    std::vector<std::vector<std::vector<double>>> errsOutputRel(Lqm_pr.size(), std::vector<std::vector<double>>(N, std::vector<double>(size)));
    std::vector<std::vector<double>> errsOutputRef(N, std::vector<double>(size)), errsOutputRefRel(N, std::vector<double>(size));
    std::vector<double> output;
    auto Xh = model->functionSpace();
    auto mesh = Xh->mesh();
    auto rangeT = elements(support(Xh));
    auto TFE = Xh->element();
    auto TRB = Xh->element();

    Feel::cout << "starting convergence study with " << size << " random parameters" << std::endl;
    int j = 0;
    for( auto const& mu : *sampling )
    {
        Feel::cout << "cvg for parameter mu=" << mu.toString() << std::endl;
        for( int i = 0; i < mu.size(); ++i )
            heatBox->addParameterInModelProperties(mu.parameterName(i), mu(i));
        heatBox->updateParameterValues();
        heatBox->solve();
        TFE = heatBox->fieldTemperature();
        auto normT = normL2( _range=rangeT, _expr=idv(TFE) );
        int k = 1;
        std::vector<double> outputsFE;
        model->computeBetaQm(mu);
        for( auto& [name, output] : outputs )
            outputsFE.push_back(model->output(k++, mu, TFE));

        auto betaFqm = model->computeBetaQm(mu).template get<1>();
        for(int n = 0; n < N; ++n)
        {
            crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outs);
            vectorN_type uN = uNs[0];
            TRB = crb->expansion( uN, n+1 );
            errs[n][j] = normL2( _range=rangeT, _expr=idv(TRB)-idv(TFE) );
            errsRel[n][j] = errs[n][j]/normT;

            auto outputsRB = computeOutputs(Lqm_pr, betaFqm, uN);
            for( int k = 0; k < outputsRB.size(); ++k )
            {
                errsOutput[k][n][j] = std::abs(outputsRB[k]-outputsFE[k]);
                errsOutputRel[k][n][j] = errsOutput[k][n][j]/std::abs(outputsFE[k]);
            }
        }
        ++j;
    }

    fs::ofstream cvgErr( "err.dat" );
    writeErrors(cvgErr, errs);
    cvgErr.close();
    fs::ofstream cvgErrR( "errR.dat" );
    writeErrors(cvgErrR, errsRel);
    cvgErrR.close();

    auto [min,max,mean,stdev] = computeStats(errsRel);
    fs::ofstream cvgStat( "stat.dat" );
    writeStats(cvgStat, min, mean, max, stdev);
    cvgStat.close();
    Feel::cout << "stats on field over size of RB" << std::endl;
    writeStats(std::cout, min, mean, max, stdev);

    double tol = doption("toolboxmor.tolerance");
    bool status = mean[N-1] < tol;

    int k = 0;
    for( auto const& [name,output] : outputs)
    {
        fs::ofstream cvgErrO( "err_"+name+".dat" );
        writeErrors(cvgErrO, errsOutput[k]);
        cvgErrO.close();
        fs::ofstream cvgErrOR( "errR_"+name+".dat" );
        writeErrors(cvgErrOR, errsOutputRel[k]);
        cvgErrOR.close();
        auto [minO,maxO,meanO,stdevO] = computeStats(errsOutputRel[k]);
        fs::ofstream cvgStatO( "stat_"+name+".dat" );
        writeStats(cvgStatO, minO, meanO, maxO, stdevO);
        cvgStatO.close();
        Feel::cout << "stats on output " << name << " over size of RB" << std::endl;
        writeStats(std::cout, minO, meanO, maxO, stdevO);
        status = status && meanO[N-1] < tol;
        k++;
    }

    auto e = exporter(_mesh=mesh);
    e->add("TFE", TFE);
    e->add("TRB", TRB);
    e->save();

    return !status;
}


int main( int argc, char** argv)
{
    using namespace Feel;
    try
    {
        po::options_description opt("options");
        opt.add_options()
            ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
            ("case.discretization", Feel::po::value<std::string>()->default_value( "P1" ), "discretization : P1,P2,P3 ")
            ( "toolboxmor.name", po::value<std::string>()->default_value( "toolboxmor" ), "Name of the db directory" )
            ( "toolboxmor.do-cvg", po::value<bool>()->default_value( true ), "do convergence test" )
            ( "toolboxmor.sampling-size", po::value<int>()->default_value(10), "size of the sampling" )
            ( "toolboxmor.tolerance", po::value<double>()->default_value(5e-2), "tolerance" )
            ;

        Environment env( _argc=argc, _argv=argv,
                         _desc=opt.add(makeToolboxMorOptions())
                         .add(toolboxes_options("heat")) );

        int dimension = ioption(_name="case.dimension");
        std::string discretization = soption(_name="case.discretization");
        

        int status = 0;
        hana::for_each( Pc_t<>,
                        [&discretization, &dimension, &status]( auto const& d ) {
                            constexpr int _dim = std::decay_t<decltype( hana::at_c<0>( d ) )>::value;
                            constexpr int _torder = std::decay_t<decltype( hana::at_c<1>( d ) )>::value;
                            std::string const& _discretization = hana::at_c<2>( d );
                            if ( dimension == _dim && discretization == _discretization )
                                status = runSimulation<_dim, _torder>();
                        } );
        return status;
    }
    catch( ... )
    {
        handleExceptions();
    }
    return EXIT_FAILURE;
}
