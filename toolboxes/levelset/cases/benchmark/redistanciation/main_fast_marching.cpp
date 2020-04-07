#include <feel/feelmodels/levelset/levelset.hpp>

namespace Feel
{

template <uint16_type OrderLevelset, uint16_type OrderLevelsetPN = LEVELSET_PN_ORDER >
void
runLevelsetApplication()
{
  tic();
    using namespace Feel;

    typedef Simplex<FEELPP_DIM,1> convex_type;
    typedef Lagrange<OrderLevelset, Scalar, Continuous, PointSetFekete> basis_type;
    // This line was changed, see main_ls.cpp
    //typedef FunctionSpace<Mesh<convex_type>, Feel::detail::bases<Lagrange<OrderLevelset, Vectorial, Continuous, PointSetFekete>>, double, Periodicity<NoPeriodicity>, mortars<NoMortar>> space_advection_velocity_type;
    typedef FunctionSpace<Mesh<convex_type>, bases<typename FeelModels::detail::ChangeBasisPolySet<Vectorial, basis_type>::type> > space_advection_velocity_type;
    typedef Lagrange<OrderLevelsetPN, Scalar, Continuous, PointSetFekete> basis_PN_type;

    typedef FeelModels::LevelSet<
        convex_type,
        basis_type,
        NoPeriodicity,
        space_advection_velocity_type,
        basis_PN_type
        > model_type;
    
    auto LS = model_type::New("levelset");
    
    Feel::cout << "============================================================\n";
    Feel::cout << "Level set benchmark : redistanciation\n";
    Feel::cout << "============================================================\n";

    toc("Before init: ");
    tic();
    LS->init();
    toc("init(): ");
    tic();
    LS->printAndSaveInfo();
    toc("print(AndSaveInfo(): ");
    if( !LS->doRestart() )
      {
        tic();
	LS->exportResults(0.);
	toc("exportResults(0.): ");
      }

    std::shared_ptr<Exporter<typename model_type::mesh_type>> myExporter;

    double hSize = doption( _name="levelset.gmsh.hsize" );
    int nDof = LS->functionSpace()->nDof();
    int nElements = LS->mesh()->numGlobalElements();
    double redistTime;
    // TODO : create meaningful base names : 
    std::string timeFileBaseName = "timeFileBaseName";
    std::string geoBaseName = "geoBaseName";
    std::string jsonBaseName = "jsonBaseName";
    std::string timeFileExtension = "csv";
    std::string timeFileName = timeFileBaseName + "_" + geoBaseName + "_" + jsonBaseName + "_" + std::to_string(hSize) + "." + timeFileExtension;
    
    Feel::cout << "Redistanciation started. Will write results to " << timeFileName << std::endl;
    tic();
    LS->redistanciate();
    redistTime = toc();
    Feel::cout << "Redistanciation time : " << redistTime << std::endl;
    Feel::cout << "h size is : " << hSize << std::endl;

    tic();
    if ( Environment::isMasterRank() )
      {
	std::ofstream timerFile;
	timerFile.open( timeFileName, std::ofstream::app );
	// write the number of MPI processes, the hsize, the number of dof, the number of elements, and the redistanciation time
	timerFile << Environment::worldComm().globalSize() << ", " << hSize << ", " << nDof << ", " << nElements << ", " << redistTime << std::endl;
	Feel::cout << "Timing results written to " << timeFileName <<  std::endl;
	timerFile.close();
	// TODO handle errors :
	//Feel::cout << "Could not write to file : " << timeFileName <<  std::endl;
      }
    toc("Write to file block :");
    tic();
    LS->exportResults();
    toc("Final LS->exportResults() :");
}

}  // namespace Feel

int main( int argc, char** argv )
{
    using namespace Feel;

    po::options_description levelsetoptions( "application levelset options" );
    levelsetoptions.add( toolboxes_options( "levelset" ) );
    levelsetoptions.add_options()
        ("fe-approximation", Feel::po::value<std::string>()->default_value( "P1" ), "fe-approximation : P2, P1" )
        ("export-dist-to-boundary", Feel::po::value<bool>()->default_value( false ), "compute and export the distance to the boundary" )
        ("levelset.redistanciation-timer-file", Feel::po::value<std::string>()->default_value( "redistanciation_FEELPP_DIMd_timer.csv" ), "Path to a file to write the redistanciation time" )
        ;

    Environment env( 
            _argc=argc, _argv=argv,
            _desc=levelsetoptions,
            _about=about(_name="application_levelset_benchmark_fast_marching",
                _author="Feel++ Consortium",
                _email="feelpp-devel@feelpp.org")
            );
    std::string feapprox = soption(_name="fe-approximation");
    if ( feapprox == "P2" )
        runLevelsetApplication<2>();
    else if ( feapprox == "P1" )
        runLevelsetApplication<1>();
    
    else CHECK( false ) << "invalid feapprox " << feapprox;

    return 0;
}
