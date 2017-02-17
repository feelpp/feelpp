
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feelvf/vf.hpp>


template<typename DatabaseType,typename SpaceType>
void
runPOD( DatabaseType & myDb, boost::shared_ptr<SpaceType> const& space, std::string const& fieldPod )
{
    using namespace Feel;
    bool useScalarProductL2 = soption(_name="scalar-product") == "L2";
    auto ui = space->element();
    auto uj = space->element();
    auto const& fieldsInfo = myDb.fieldsInfo();
    auto const& timeSetDb = myDb.timeSet();

    double ti = timeSetDb.front();
    double tf = timeSetDb.back();
    if ( Environment::vm().count("time-initial") )
        ti = doption(_name="time-initial");
    if ( Environment::vm().count("time-final") )
        tf = doption(_name="time-final");
    std::vector<double> timeSetIndex;
    for (int k=0; k<timeSetDb.size();++k)
    {
        if ( ( timeSetDb[k] >= (ti - 1e-9) ) &&
             ( timeSetDb[k] <= (tf + 1e-9 ) ) )
            timeSetIndex.push_back( timeSetDb[k] );
    }
    int nTimeStep = timeSetIndex.size();
    if ( nTimeStep == 0 )
    {
        std::cout << "exit pod because timeset is empty\n";
        return;
    }

    if ( myDb.worldComm().isMasterRank() )
        std::cout << "start build matrix pod of size : " << nTimeStep << "," << nTimeStep << "\n";

    //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pod( timeSet.size(), timeSet.size() );
    MatrixEigenDense<double> matrixPodFeel( nTimeStep,nTimeStep );
    auto & pod = matrixPodFeel.mat();

    for ( int i=0;i<nTimeStep;++i )
    {
        myDb.load( timeSetIndex[i],fieldPod,ui );
        for ( int j=0;j<i;++j )
        {
            myDb.load( timeSetIndex[j],fieldPod,uj );
            if ( useScalarProductL2 )
                pod( i,i ) = integrate( _range=elements(space->mesh()),
                                        _expr=inner(idv(ui),idv(uj) ) ).evaluate(false)(0,0);
            else
                pod( i,j ) = inner_product( ui, uj );
            pod( j,i ) = pod( i,j );
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "("<<i<<","<<j<<")"<<std::flush;
        }
        if ( useScalarProductL2 )
            pod( i,i ) = integrate( _range=elements(space->mesh()),
                                    _expr=inner(idv(ui),idv(ui), mpl::int_<InnerProperties::IS_SAME/*|InnerProperties::SQRT*/>() ) ).evaluate(false)(0,0);
        else
            pod( i,i ) = inner_product( ui,ui );
        if ( myDb.worldComm().isMasterRank() )
            std::cout << "("<<i<<","<<i<<")\n"<<std::flush;
    }

    if ( myDb.worldComm().isMasterRank() )
    {
        std::cout << "\npod=\n" << pod << "\n";
        std::ofstream file("pod.txt");
        if (file.is_open())
        {
            file << pod;
            file.close();
        }
        matrixPodFeel.printMatlab( "pod.m" );
    }

    Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
    eigenSolver.compute( pod );
    auto eigenValues = eigenSolver.eigenvalues();
    //double myeigen = real(eigenValues[k]);
    auto eigenVectors = eigenSolver.eigenvectors();
    if ( myDb.worldComm().isMasterRank() )
        std::cout << "eigenValues=\n" << eigenValues << "\n";
}

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description exportdboptions( "export database options" );
	exportdboptions.add_options()
        ( "ifile", po::value<std::string>(), "input database file" )
        ( "field", po::value<std::string>(), "field loaded in database " )
        ( "scalar-product", po::value<std::string>()->default_value( "L2" ), "scalar product used in pod : L2,euclidian " )
        ( "time-initial", po::value<double>(), "initial time used for pod" )
        ( "time-final", po::value<double>(), "final time used for pod" )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=exportdboptions,
                   _about=about(_name="export_database",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string ifile;
    if ( Environment::vm().count("ifile") )
        ifile = Environment::vm()["ifile"].as<std::string>();
    std::string fieldPod;
    if ( Environment::vm().count("field") )
        fieldPod = Environment::vm()["field"].as<std::string>();

    if ( ifile.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile is missing\n";
        return 0;
    }
    if ( fieldPod.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --field is missing\n";
        return 0;
    }

    typedef Mesh<Simplex<3>> mesh_type;
    FeelppDatabase<mesh_type> myDb;
    myDb.setFilename( ifile );
    myDb.loadInfo();

    if ( !myDb.hasField( fieldPod ) )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because fieldPod does not found in databse\n";
        return 0;

   }

    CHECK( myDb.fieldInfo( fieldPod ).basisName() == "lagrange" ) << "only lagrange instantiated : " << myDb.fieldInfo( fieldPod ).basisName();
    CHECK( myDb.fieldInfo( fieldPod ).basisOrder() == 1 ) << "only order 1 instantiated : " << myDb.fieldInfo( fieldPod ).basisOrder();

    auto mesh = myDb.loadMesh( MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES );


    if ( myDb.isScalarField( fieldPod ) )
    {
        auto Vh = Pch<1>( mesh );
        runPOD( myDb, Vh, fieldPod );
    }
    else if ( myDb.isScalarField( fieldPod ) )
    {
        auto Vh = Pchv<1>( mesh );
        runPOD( myDb, Vh, fieldPod );
    }

    return 0;
}
