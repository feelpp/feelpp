
#include <feel/feelcore/environment.hpp>
#include <feel/feelalg/matrixeigendense.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feelvf/vf.hpp>


template<typename DatabaseType,typename SpaceType>
void
runPOD( DatabaseType & myDb, std::shared_ptr<SpaceType> const& space, std::string const& fieldPod )
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
            timeSetIndex.push_back( k );
    }
    int nTimeStep = timeSetIndex.size();
    if ( nTimeStep == 0 )
    {
        std::cout << "exit pod because timeset is empty\n";
        return;
    }

    if ( myDb.worldComm().isMasterRank() )
        std::cout << "start build matrix pod of size : " << nTimeStep << "," << nTimeStep << "\n";

    std::shared_ptr<MatrixSparse<double>> scalarProductOperator;
    bool useScalarProductOperator = boption(_name="scalar-product.use-operator");
    if ( useScalarProductL2 && useScalarProductOperator )
    {
        scalarProductOperator = backend()->newMatrix(_test=space,_trial=space);
        form2(_test=space,_trial=space,_matrix=scalarProductOperator ) =
            integrate(_range=elements(space->mesh()),_expr=inner(idt(ui),id(ui)) );
    }

    //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pod( timeSet.size(), timeSet.size() );
    MatrixEigenDense<double> matrixPodFeel( nTimeStep,nTimeStep );
    auto & pod = matrixPodFeel.mat();

    for ( int i=0;i<nTimeStep;++i )
    {
        if ( myDb.worldComm().isMasterRank() )
            std::cout << "compute values of pod matrix at row "<<i+1<<"/"<<nTimeStep << std::flush;
        myDb.load( timeSetIndex[i],fieldPod,ui );
        for ( int j=0;j<i;++j )
        {
            myDb.load( timeSetIndex[j],fieldPod,uj );
            if ( useScalarProductL2 )
            {
                if ( useScalarProductOperator )
                    pod( i,j ) = scalarProductOperator->energy( ui,uj );
                else
                    pod( i,j ) = integrate( _range=elements(space->mesh()),
                                            _expr=inner(idv(ui),idv(uj) ) ).evaluate()(0,0);
            }
            else
            {
                pod( i,j ) = inner_product( ui, uj );
            }
            pod( j,i ) = pod( i,j );
        }
        if ( useScalarProductL2 )
        {
            if ( useScalarProductOperator )
                pod( i,i ) = scalarProductOperator->energy( ui,ui );
            else
                pod( i,i ) = integrate( _range=elements(space->mesh()),
                                        _expr=inner(idv(ui),idv(ui), mpl::int_<InnerProperties::IS_SAME/*|InnerProperties::SQRT*/>() ) ).evaluate()(0,0);
        }
        else
        {
            pod( i,i ) = inner_product( ui,ui );
        }
        if ( myDb.worldComm().isMasterRank() )
            std::cout << " : done\n"<<std::flush;
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

    std::string eigenSolverType = soption(_name="eigen-solver");
    if ( eigenSolverType == "general" )
    {
        Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
        eigenSolver.compute( pod );
        auto eigenValues = eigenSolver.eigenvalues();
        //double myeigen = real(eigenValues[k]);
        auto eigenVectors = eigenSolver.eigenvectors();
        /*if ( myDb.worldComm().isMasterRank() )
            std::cout << "eigenValues=\n" << eigenValues << "\n";
        if ( myDb.worldComm().isMasterRank() )
         std::cout << "eigenVectors=\n" << eigenVectors << "\n";*/
        //std::cout << "eigenValues.size() " << eigenValues.size() << " and eigenVectors.size() " << eigenVectors.size() << " " << eigenVectors.rows() << " " << eigenVectors.cols()  << "\n";
        int nEigenValues = eigenValues.size();
        std::vector<std::pair<double,int> > eigenValuesSortedWithInitialId;
        for (int k=0;k<nEigenValues;++k)
            eigenValuesSortedWithInitialId.push_back( std::make_pair(real(eigenValues[k]),k) );
        std::sort(eigenValuesSortedWithInitialId.begin(), eigenValuesSortedWithInitialId.end(),
                  [](std::pair<double,int> const& a, std::pair<double,int>const& b) { return b.first > a.first; });
        Eigen::VectorXd eigenValuesSorted(nEigenValues);
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigenVectorsSorted(nEigenValues,nEigenValues);
        for(int k1=0;k1<nEigenValues;++k1)
        {
            int k1Initial = eigenValuesSortedWithInitialId[k1].second;
            eigenValuesSorted(k1) = real(eigenValues[k1Initial]);
            for(int k2=0;k2<nEigenValues;++k2)
                eigenVectorsSorted(k2,k1) = real(eigenVectors( k2,k1Initial ));
        }
        if ( false )
        {
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "eigenValues (sorted)=\n" << eigenValuesSorted << "\n";
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "eigenVectors (sorted)=\n" << eigenVectorsSorted << "\n";
        }
    }

    if ( eigenSolverType == "self-adjoint" )
    {
        Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
        eigenSolver.compute( pod );
        auto eigenValues = eigenSolver.eigenvalues();
        auto eigenVectors = eigenSolver.eigenvectors();
        if ( false )
        {
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "eigenValues=\n" << eigenValues << "\n";
            if ( myDb.worldComm().isMasterRank() )
                std::cout << "eigenVectors=\n" << eigenVectors << "\n";
        }
    }
}

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description exportdboptions( "export database options" );
	exportdboptions.add_options()
        ( "ifile", po::value<std::string>(), "input database file" )
        ( "field", po::value<std::string>(), "field loaded in database " )
        ( "scalar-product", po::value<std::string>()->default_value( "L2" ), "scalar product used in pod : L2,euclidian " )
        ( "scalar-product.use-operator", po::value<bool>()->default_value( true ), "build  scalar product operator  " )
        ( "time-initial", po::value<double>(), "initial time used for pod" )
        ( "time-final", po::value<double>(), "final time used for pod" )
        ( "eigen-solver", po::value<std::string>()->default_value( "self-adjoint" ), "self-adjoint, general" )
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
    else if ( myDb.isVectorialField( fieldPod ) )
    {
        auto Vh = Pchv<1>( mesh );
        runPOD( myDb, Vh, fieldPod );
    }

    return 0;
}
