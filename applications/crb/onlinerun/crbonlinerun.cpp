
#include <heat3d/heat3d.hpp>
#include <linearelasticity3d/linearelasticity3d.hpp>
//#include <BenchmarkGrepl/benchmarkgrepl-linear-elliptic.hpp>
//#include <BenchmarkGrepl/benchmarkgrepl-nonlinear-elliptic.hpp>

std::string
loadModelName( std::string const& filename )
{
    using namespace Feel;
    if ( !fs::exists( filename ) )
    {
        LOG(INFO) << "Could not find " << filename << std::endl;
        return std::string("");
    }

    auto json_str_wo_comments = removeComments(readFromFile(filename));
    //LOG(INFO) << "json file without comment:" << json_str_wo_comments;

    boost::property_tree::ptree ptree;
    std::istringstream istr( json_str_wo_comments );
    boost::property_tree::read_json( istr, ptree );

    auto const& ptreeCrbModel = ptree.get_child( "crbmodel" );
    std::string modelName = ptreeCrbModel.template get<std::string>( "model-name" );
    return modelName;
}

template <typename ModelType>
boost::shared_ptr<Feel::CRB<Feel::CRBModel< ModelType> > >
loadCrbOnline( std::string const& filename )
{
    using namespace Feel;
    typedef Feel::CRBModel< ModelType > crbmodel_type;
    typedef Feel::CRB<crbmodel_type> crb_type;
    boost::shared_ptr<crb_type> crb;

    if ( !fs::exists( filename ) )
        return crb;

    boost::shared_ptr<ModelType> model( new ModelType );
    model->loadJson( filename, "crbmodel" );
    boost::shared_ptr<crbmodel_type> crbmodel( new crbmodel_type( model, Feel::CRBModelMode::CRB, false ) );
    crbmodel->loadJson( filename, "crbmodel" );

    crb.reset( new crb_type );
    crb->setTruthModel( crbmodel );
    crb->loadJson( filename );
//crb->setOfflineStep( false );

    return crb;
}
template <typename ModelType>
bool
runCrbOnline( std::string const& filename )
{
    using namespace Feel;

    auto crb = loadCrbOnline<ModelType>( Environment::expand( soption(_name="db.filename") ) );
    if ( !crb )
    {
        std::cout << "failure in crb loading -> exit program\n";
        return false;
    }

    Eigen::VectorXd/*typename crb_type::vectorN_type*/ time_crb;
    double online_tol = 1e-2;//Feel::doption(Feel::_name="crb.online-tolerance");
    bool print_rb_matrix = false;//boption(_name="crb.print-rb-matrix");
    auto muspace = crb->Dmu();

    auto mysampling = muspace->sampling();
    int nSample = ioption(_name="sampling.size");
    std::string sampler = soption("sampling.type");
    mysampling->sample( nSample, sampler );

    bool loadFiniteElementDatabase = boption(_name="crb.load-elements-database");
    boost::shared_ptr<Exporter<typename ModelType::mesh_type> > fieldExporter;
    if ( loadFiniteElementDatabase )
        fieldExporter = exporter( _mesh=crb->model()->rBFunctionSpace()->mesh() );


    int nSamples = mysampling->size();
    for ( int k=0;k<nSamples;++k )
    {
        auto const& mu = (*mysampling)[k];
        std::ostringstream ostrmu;
        for ( uint16_type d=0;d<muspace->dimension();++d)
            ostrmu << mu(d) << " ";
        std::cout << "--------------------------------------\n";
        std::cout << "mu["<<k<<"] : " << ostrmu.str() << "\n";
        // std::cout << "parameterSpace->min: "<< muspace->min() << "\n";
        // std::cout << "parameterSpace->max: "<< muspace->max() << "\n";
        //auto mu = crb->Dmu()->element();
        //std::cout << "input mu\n" << mu << "\n";
        auto crbResult = crb->run( mu, time_crb, online_tol, -1, print_rb_matrix);
        auto resOuptut = boost::get<0>( crbResult );
        auto resError = boost::get<0>( boost::get<6>( crbResult ) );
        std::cout << "output " << resOuptut.back() << "\n";
        std::cout << "err " << resError.back() << "\n";

        if ( loadFiniteElementDatabase )
        {
            auto const& solutions = crbResult.template get<2>();
            auto uN = crb->model()->rBFunctionSpace()->element();
            //auto const& uN = solutions.template get<0>();
            uN.container() = solutions.template get<0>().back();
            auto sol = crb->model()->rBFunctionSpace()->expansion( uN );
            fieldExporter->add( (boost::format("sol-%1%")%k).str(), sol );
        }
    }
    if ( loadFiniteElementDatabase )
        fieldExporter->save();

    return true;

}


int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description crbonlinerunoptions( "crb online run options" );
	crbonlinerunoptions.add_options()
        ( "db.filename", po::value<std::string>()->default_value( "" ), "database filename" )
        ( "sampling.size", po::value<int>()->default_value( 10 ), "size of sampling" )
        ( "sampling.type", po::value<std::string>()->default_value( "random" ), "type of sampling" )
	 	;
	po::options_description crbonlinerunliboptions( "crb online run lib options" );
    crbonlinerunliboptions.add(crbOptions())
        .add(crbSEROptions())
        .add(eimOptions())
        .add(podOptions())
        .add(backend_options("backend-primal"))
        .add(backend_options("backend-dual"))
        .add(backend_options("backend-l2"))
        ;

	Environment env( _argc=argc, _argv=argv,
                     _desc=crbonlinerunoptions,
                     _desc_lib=crbonlinerunliboptions.add( feel_options() ),
                     _about=about(_name="crbonlinerun",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    std::string jsonfilename = Environment::expand( soption(_name="db.filename") );
    std::string modelName = loadModelName( jsonfilename );
    bool successrun = false;
    if ( modelName == "Heat3d" )
        successrun = runCrbOnline<Heat3d>( jsonfilename );
    else if ( modelName == "LinearElasticity" )
        successrun = runCrbOnline<LinearElasticity3d>( jsonfilename );
#if 0
    else if ( modelName == "BenchMarkGreplLinearElliptic1" )
        successrun = runCrbOnline<BenchmarkGreplLinearElliptic<1>>( jsonfilename );
    else if ( modelName == "BenchMarkGreplNonlinearElliptic1" )
        successrun = runCrbOnline<BenchmarkGreplNonlinearElliptic<1,2>>( jsonfilename );
#endif
    else
    {
        std::cout << "invalid modelName : " <<  modelName << "\n";
        return 0;
    }

    return 0;
}
