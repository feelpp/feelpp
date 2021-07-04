#include <feel/feelcrb/toolboxmor.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes.hpp>
#include <feel/feelmodels/coefficientformpdes/coefficientformpdes_registered_type.hpp>

using namespace Feel;

void writeErrors(fs::ofstream& out, std::vector<std::vector<double> > const& err)
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
        out.close();
    }
}

template<typename ConvexType, typename BasisType>
int runSimulation(std::shared_ptr<FeelModels::coefficient_form_PDEs_t<ConvexType>>& cfpdes)
{
    using model_type = FeelModels::coefficient_form_PDEs_t<ConvexType>;
    using model_ptrtype = std::shared_ptr<model_type>;
    using space_type = FunctionSpace<Mesh<ConvexType>, bases<BasisType>>;
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

    rb_model_ptrtype model = std::make_shared<rb_model_type>();
    auto cfpde = cfpdes->coefficientFormPDE( std::get<1>(cfpdes->pdes()[0])->equationName(), hana::type_c<BasisType> );
    model->setFunctionSpaces(cfpde->spaceUnknown());
    auto rhs = cfpdes->algebraicFactory()->rhs()->clone();
    auto mat = cfpdes->algebraicFactory()->matrix();
    deim_function_type assembleDEIM =
        [&cfpdes,&rhs,&mat](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    cfpdes->addParameterInModelProperties(mu.parameterName(i), mu(i));
                cfpdes->updateParameterValues();
                rhs->zero();
                cfpdes->algebraicFactory()->applyAssemblyLinear( cfpdes->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.lhs"} );
                return rhs;
            };
    model->setAssembleDEIM(assembleDEIM);
    mdeim_function_type assembleMDEIM =
        [&cfpdes,&rhs,&mat](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    cfpdes->addParameterInModelProperties(mu.parameterName(i), mu(i));
                cfpdes->updateParameterValues();
                mat->zero();
                cfpdes->algebraicFactory()->applyAssemblyLinear( cfpdes->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhs, {"ignore-assembly.rhs"} );
                return mat;
            };
    model->setAssembleMDEIM(assembleMDEIM);
    model->initModel();

    auto deimToolbox = std::make_shared<model_type>("cfpdes");
    deimToolbox->setMesh(model->getDEIMReducedMesh());
    deimToolbox->init();
    deimToolbox->printAndSaveInfo();
    deim_function_type assembleOnlineDEIM =
        [deimToolbox](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    deimToolbox->addParameterInModelProperties(mu.parameterName(i), mu(i));
                deimToolbox->updateParameterValues();
                auto rhs = deimToolbox->algebraicFactory()->rhs()->clone();
                rhs->zero();
                auto matTMP = deimToolbox->algebraicFactory()->matrix()->clone();
                deimToolbox->algebraicFactory()->applyAssemblyLinear( deimToolbox->algebraicBlockVectorSolution()->vectorMonolithic(), matTMP, rhs, {"ignore-assembly.lhs"} );
                return rhs;
            };
    model->setOnlineAssembleDEIM(assembleOnlineDEIM);
    // model->setOnlineAssembleDEIM(assembleDEIM);

    auto mdeimToolbox = std::make_shared<model_type>("cfpdes");
    mdeimToolbox->setMesh(model->getMDEIMReducedMesh());
    mdeimToolbox->init();
    mdeimToolbox->printAndSaveInfo();
    mdeim_function_type assembleOnlineMDEIM =
        [mdeimToolbox](parameter_type const& mu)
            {
                for( int i = 0; i < mu.size(); ++i )
                    mdeimToolbox->addParameterInModelProperties(mu.parameterName(i), mu(i));
                mdeimToolbox->updateParameterValues();
                auto mat = mdeimToolbox->algebraicFactory()->matrix()->clone();
                mat->zero();
                auto rhsTMP = mdeimToolbox->algebraicFactory()->rhs()->clone();
                mdeimToolbox->algebraicFactory()->applyAssemblyLinear( mdeimToolbox->algebraicBlockVectorSolution()->vectorMonolithic(), mat, rhsTMP, {"ignore-assembly.rhs"} );
                return mat;
            };
    model->setOnlineAssembleMDEIM(assembleOnlineMDEIM);
    // model->setOnlineAssembleMDEIM(assembleMDEIM);

    model->postInitModel();
    model->setInitialized(true);

    crb_model_ptrtype crbModel = std::make_shared<crb_model_type>(model);
    crb_ptrtype crb = crb_type::New("toolboxmor", crbModel, crb::stage::offline);

    tic();
    crb->offline();
    toc("offline");

    int N = crb->dimension();
    int timeSteps = 1;
    std::vector<vectorN_type> uNs(timeSteps, vectorN_type(N)), uNolds(timeSteps, vectorN_type(N));
    std::vector<double> outputs(timeSteps, 0);

    sampling_ptrtype sampling( new sampling_type( crbModel->parameterSpace() ) );
    int size = ioption("toolboxmor.sampling-size");
    sampling->clear();
    sampling->randomize( size, true );

    std::vector<std::vector<double> > errs(N, std::vector<double>(size));
    std::vector<std::vector<double> > errsRel(N, std::vector<double>(size));
    auto Xh = model->functionSpace();
    auto mesh = Xh->mesh();
    auto rangeU = elements(mesh);//Xh->dof()->meshSupport()->rangeElements();
    auto UFE = Xh->element();
    auto URB = Xh->element();

    int j = 0;
    for( auto const& mu : *sampling )
    {
        for( int i = 0; i < mu.size(); ++i )
            cfpdes->addParameterInModelProperties(mu.parameterName(i), mu(i));
        cfpdes->updateParameterValues();
        cfpdes->solve();
        UFE = cfpde->fieldUnknown();
        auto normU = normL2( rangeU, idv(UFE) );
        for(int n = 0; n < N; ++n)
        {
            crb->fixedPointPrimal(n+1, mu, uNs, uNolds, outputs);
            vectorN_type uN = uNs[0];
            URB = crb->expansion( uN, n+1 );
            errs[n][j] = normL2( rangeU, idv(URB)-idv(UFE) );
            errsRel[n][j] = errs[n][j]/normU;
        }
        ++j;
    }

    std::vector<double> min(N), max(N), mean(N), stdev(N);
    for(int n = 0; n < N; ++n)
    {
        min[n] = *std::min_element(errsRel[n].begin(), errsRel[n].end());
        max[n] = *std::max_element(errsRel[n].begin(), errsRel[n].end());
        double s = std::accumulate(errsRel[n].begin(), errsRel[n].end(), 0.0);
        mean[n] = s/size;
        double accum = std::accumulate(errsRel[n].begin(), errsRel[n].end(), 0.0,
                                       [s,size](double a, double b) {
                                           return a + (b-s/size)*(b-s/size);
                                       });
        stdev[n] = accum/size;
    }

    fs::ofstream cvgErr( "err.dat" );
    fs::ofstream cvgErrR( "errR.dat" );
    fs::ofstream cvgStat( "stat.dat" );
    writeErrors(cvgErr, errs);
    writeErrors(cvgErrR, errsRel);
    if( cvgStat && Environment::isMasterRank() )
    {
        cvgStat << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
                << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
        for(int n = 0; n < N; ++n)
            cvgStat << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                    << std::setw(25) << mean[n] << std::setw(25) << stdev[n] << std::endl;
        cvgStat.close();
    }
    Feel::cout << std::setw(5) << "N" << std::setw(25) << "min" << std::setw(25) << "max"
               << std::setw(25) << "mean" << std::setw(25) << "stdev" << std::endl;
    for(int n = 0; n < N; ++n)
        Feel::cout << std::setw(5) << n+1 << std::setw(25) << min[n] << std::setw(25) << max[n]
                   << std::setw(25) << mean[n] << std::setw(25) << stdev[n] << std::endl;

    auto e = exporter(mesh);
    e->add("UFE", UFE);
    e->add("URB", URB);
    e->save();

    return 0;
}

int main( int argc, char** argv)
{
    po::options_description opt("options");
    opt.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ( "toolboxmor.sampling-size", po::value<int>()->default_value(10), "size of the sampling" )
        ;
    Environment env( _argc=argc, _argv=argv,
                     _desc=opt.add(makeToolboxMorOptions())
                     .add(toolboxes_options( "coefficient-form-pdes", "cfpdes" ) ) );

    int dimension = ioption(_name="case.dimension");
    auto dimt = hana::make_tuple(hana::int_c<2>/*,hana::int_c<3>*/);
    int status = 0;
    hana::for_each( dimt,
                    [&dimension,&status]( auto const& d ) {
                        constexpr int _dim = std::decay_t<decltype(d)>::value;
                        if ( dimension == _dim  ) {
                            using model_type = FeelModels::coefficient_form_PDEs_t< Simplex<_dim> >;
                            using model_ptrtype = std::shared_ptr<model_type>;

                            std::shared_ptr<model_type> cfpdes( new model_type("cfpdes") );
                            cfpdes->init();
                            cfpdes->printAndSaveInfo();
                            auto cfpde = std::get<1>(cfpdes->pdes()[0]);
                            auto unknownBasis = cfpde->unknownBasis();

                            using cfpdes_type = unwrap_ptr_t<std::decay_t<decltype(cfpdes)>>;
                            hana::for_each( cfpdes_type::tuple_type_unknown_basis,
                                            [&cfpdes,&_dim,&status,&unknownBasis]( auto & e ) {
                                                if ( cfpdes_type::unknowBasisTag( e ) == unknownBasis )
                                                {
                                                    status = runSimulation<Simplex<_dim>, typename std::decay_t<decltype(e)>::type>(cfpdes);
                                                }
                                            });

                        }
                    } );
    return status;
}
