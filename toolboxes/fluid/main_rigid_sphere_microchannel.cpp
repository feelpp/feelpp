#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feells/disttoentityrange.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feelvf/print.hpp>
#include <feel/feelmodels/fluid/fluidmechanics.hpp>
#include <feel/options.hpp>

namespace Feel
{

// OK
template <typename type_element>
auto computeElementFacesAreasRatio( type_element elt )
{
    auto elementFacesAreas = elt.faceMeasures();
    auto minFaceArea = elementFacesAreas.minCoeff();
    auto maxFaceArea = elementFacesAreas.maxCoeff();
    auto areasRatio = maxFaceArea/maxFaceArea;
    return areasRatio;
}

// OK
template <typename type_mesh, typename type_threshold>
bool mustRemesh( type_mesh const& mesh, type_threshold ratioThreshold )
{
    bool willRemesh = false;

    for ( auto const& eltWrap : elements( mesh ) )
    {
        auto const& elt = unwrap_ref( eltWrap );
        auto areasRatio = computeElementFacesAreasRatio( elt );
        if ( areasRatio > ratioThreshold )
        {
            Feel::cout   << "mustRemesh(): Degenerate tetrahedron at " 
                        <<  elt.barycenter() 
                        << " faces areas ratio is" 
                        << areasRatio 
                        << std::endl;
            willRemesh = true;
            break;
        }
    }
    return willRemesh;
}


// auto initialRemeshedMesh = myRemesh( referenceFM->meshALE()->movingMesh(), referenceFM->meshALE() );
// Seems OK, check metric.on(...)
template <typename typemesh, typename typealemesh>
typemesh const myRemesh( typemesh const& mesh, typealemesh const& alemesh, typemesh parent = {} )
{
    // Mesh before remeshing
    auto originalMesh = mesh;
    // Get P1 function space defined on the mesh
    auto Xh = Pch<1>( originalMesh );
    // Create metric field 
    auto metric = Xh->element();
    // Original mesh statistics, we could use the average element size
    // which should be close to the desired target size if element size is constant on the mesh
    // or use the minimun element size if the average is already too large 
    // (but shouldn't the degenerate tetrahedra/triangle criterion avoid this issue?)
    auto [havg, hmin, hmax] = hMeasures( originalMesh );
    Feel::cout   << "[havg, hmin, hmax] = " << "[" << havg << ", " << hmin << ", " << hmax << "]" << std::endl;    
    // For now, let's try with the min since the test case's initial h-size is constant
    // Update, I switch to havg because mmg segfaults, so I try with havg (update : same result).
    // We set this target element size on the whole mesh 
    // TODO maybe set this only on the "solvent" marked elements? Can MMG handle this?
    Feel::cout   << "Before metric" << std::endl;
    metric.on( _range = elements( mesh ),
               _expr = cst( havg ) );
    // Create remesher.
    // Not sure what the vector, the flag and the parent are for ? Maybe for partial (submesh) remeshing ?
    // TODO require non moving parts to remain stable (inlet, oulet, walls etc)
    Feel::cout   << "Before remesher() with membrane" << std::endl;
    auto myRemesher = remesher( originalMesh, std::vector<int>{}, alemesh->aleFactory()->flagSet( "moving" ), parent );
    // This segfaults, maybe the rigid body's boundary is not set to "moving" before the first time step?
    // I try to specify the "Membrane" marker for the elements :
    //auto myRemesher = remesher( originalMesh, std::vector<int>{}, markedfaces( alemesh, "Membrane" ), parent );
    // With the following line -> Check failed: it != vm.end() Invalid option remesh.verbose
    //auto myRemesher = remesher( originalMesh, std::vector<int>{}, std::vector<int>{}, parent );
    // Set remeshing target element size
    Feel::cout   << "Before setMetric()" << std::endl;
    myRemesher.setMetric( metric );
    // Create new remeshed mesh
    Feel::cout   << "Before myRemesher.execute()" << std::endl;
    auto newMesh = myRemesher.execute();
    // ?
    Feel::cout   << "Before updateForUse" << std::endl;
    newMesh->updateForUse();
    // Copy original mesh markers on new mesh 
    Feel::cout   << "Before setMarkerNames" << std::endl;
    newMesh->setMarkerNames( originalMesh->markerNames() );
    // Return new, remeshed mesh
    return newMesh;
}


// To add all the measures to a single, global CSV file
// because fluid.measures.csv is overwritten
void updateMeasureFile( std::string fileName, std::string rootRepository )
{
    std::string filename = rootRepository + "/" + fileName;
    std::ofstream of_c( filename, std::ios_base::binary | std::ios_base::app );
    std::string measureName = rootRepository + "/fluid.measures.csv";
    std::ifstream if_b( measureName, std::ios_base::binary );
    of_c.seekp( 0, std::ios_base::end );
    of_c << if_b.rdbuf();
}

template <int nDim,uint16_type OrderVelocity,uint16_type OrderPressure, uint16_type OrderGeo = 1>
int
runApplicationFluid()
{
    using namespace Feel;

    // Option : max/min areas threshold
    double maxminfacesareasratio = doption( _name="remesh.maxminareasratiothreshold" );
    Feel::cout << "maxminfacesareasratio " << maxminfacesareasratio << std::endl;

    typedef FeelModels::FluidMechanics< Simplex<nDim,OrderGeo>,
                                        Lagrange<OrderVelocity, Vectorial,Continuous,PointSetFekete>,
                                        Lagrange<OrderPressure, Scalar,Continuous,PointSetFekete> > model_type;
    
    // Used for displacement field
    using trace_space_component_type = typename model_type::space_trace_velocity_type::component_functionspace_type;

    auto referenceFM = model_type::New("fluid");

    // referenceFM was FM_ref
    referenceFM->init();
    
    auto initialGlobalMesh = referenceFM->mesh();
    
    auto globalBDF = bdf( _space = referenceFM->functionSpaceVelocity(), _name = "mybdf", _prefix = "mybdf" );

    // Exporters for reference and moving mesh
    // Define a moving mesh exporter
    auto movExporter = exporter( _mesh = initialGlobalMesh, _name = "move", _geo = "change" );
    // Exporter save function
    auto movExporterSave = [&movExporter]( double t, auto alemesh, auto vel, auto press, auto w ) {
        movExporter->step( t )->setMesh( alemesh->movingMesh() );
        movExporter->step( t )->add( "disp",alemesh->displacement() );
        movExporter->step( t )->add( "vel", vel );
        movExporter->step( t )->add( "press",press);
        movExporter->step( t )->add( "test_interp", w );
        movExporter->save(); 
    };
    // Define a reference mesh exporter
    auto referenceExporter = exporter( _mesh = initialGlobalMesh, _name = "reference", _geo = "change" );
    // Exporter save function
    auto referenceExporterSave = [&referenceExporter]( double t, auto alemesh ) {
        referenceExporter->step( t )->setMesh( alemesh->referenceMesh() );
        referenceExporter->step( t )->add( "disp", alemesh->displacement() );
        referenceExporter->save();
    };

    //
    bool initialConfigurationStep = true;
    // New FM
    std::shared_ptr<model_type> currentFM( new model_type( "fluid" ) );


    // delta_CM and delta_CM_init are used for center of mass management
    double delta_CM, delta_CM_init;
    int number_of_remesh = 0;


    if ( initialConfigurationStep == true )
    {
        Feel::cout << "Forced first remeshing" << std::endl;

        // Initial remeshing (why not, this way I can see if it crashes or not)
        auto initialRemeshedMesh = myRemesh( referenceFM->meshALE()->movingMesh(), referenceFM->meshALE() );
        Feel::cout << "done." << std::endl; 
        // Saving remeshed mesh as json and msh
        initialRemeshedMesh->saveHDF5( "inital_remesh.json" );
        saveGMSHMesh( _mesh = initialRemeshedMesh, _filename = "inital_remesh.msh" );
        auto remeshing_time = referenceFM->time();
        
        // Creating new fluid model (toolbox) 
        // currentFM was FM
        auto currentFM = model_type::New( "fluid" );
        currentFM->setMesh( initialRemeshedMesh );
        currentFM->setTimeInitial( remeshing_time );
        // initial "t=0" time (this will change later)
        currentFM->addParameterInModelProperties( "T0", 0 );
        currentFM->init( true );

        // Quality field - actually useless since we rely on an element tetrahedron
        // TODO This should be re-enabled if we switch to a mesh quality field instead
        // Xh_Quality = space_quality_type::New( _mesh = currentFM->mesh() );
        // quality_field = Xh_Quality->element();

        // Time step management
        currentFM->timeStepBase()->setTimeInitial( globalBDF->time() - globalBDF->timeStep() );
        currentFM->printAndSaveInfo();
        currentFM->startTimeStep();

        // Center of mass management
        auto center_of_mass_init = integrate( _range = markedelements( currentFM->mesh(), "Cell" ), _expr = P() ).evaluate();
        auto mass = integrate( _range = markedelements( currentFM->mesh(), "Cell" ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
        center_of_mass_init /= mass;
        delta_CM_init = center_of_mass_init( 0, 0 );
        initialConfigurationStep = false;

    }// End of the initial configuration step

    Feel::cout << "Entering global BDF loop " << std::endl;
    for ( globalBDF->start(); !globalBDF->isFinished(); globalBDF->next() )
    {

        if ( referenceFM->worldComm().isMasterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "globalBDF: " << globalBDF->time() << "s \n";
            std::cout << "referenceFM time" << referenceFM->time() << "s \n";
            std::cout << "============================================================\n";
        }

        if ( currentFM->worldComm().isMasterRank() )
        {
            std::cout << "============================================================\n";
            std::cout << "time simulation (currentFM): " << currentFM->time() << "s \n";
            std::cout << "============================================================\n";
            std::cout << "CHECK IF THE MESH QUALITY IS GOOD\n";

            // Analyse each element's quality for the whole mesh.
            // If any degenerate element is found (according to a certain criterion), trigger the remeshing.
            if( mustRemesh(currentFM->meshALE()->movingMesh(), maxminfacesareasratio) )
            {
                number_of_remesh += 1;

                delta_CM = delta_CM_init;
                std::cout << "Remeshing..." << std::endl;
                auto new_mesh = myRemesh( currentFM->meshALE()->movingMesh(), currentFM->meshALE(), referenceFM->mesh() );
                new_mesh->saveHDF5( "new_remesh.json" );
                saveGMSHMesh( _mesh = new_mesh, _filename = "new_mesh.msh" );
                auto remeshing_time = currentFM->time();
                
                // Creating new fluid model (toolbox) after remeshing
                // remeshedFM was FM_after_remeshing
                std::shared_ptr<model_type> remeshedFM( new model_type( "fluid" ) );
                remeshedFM->setMesh( new_mesh );
                remeshedFM->setTimeInitial( remeshing_time - currentFM->timeStep() );
                remeshedFM->addParameterInModelProperties( "T0", remeshing_time - currentFM->timeStep() );
                auto center_of_mass_curent = integrate( _range = markedelements( currentFM->mesh(), "Cell" ), _expr = P() ).evaluate();
                // This should remain constant if cell domain is not remeshed
                auto mass = integrate( _range = markedelements( currentFM->mesh(), "Cell" ), _expr = cst( 1.0 ) ).evaluate()( 0, 0 );
                center_of_mass_curent /= mass;
                delta_CM -= center_of_mass_curent( 0, 0);
                std::cout << "delta_CM " << delta_CM << std::endl;
                
                remeshedFM->init( true );
                remeshedFM->updateParameterValues();
                auto Yh = Pch<1>( new_mesh );
                auto w = Yh->element();
                w.on( _range = elements( support( Yh ) ), _expr = Px() );
                std::cout << "OK after init" << std::endl;
                
                remeshedFM->printAndSaveInfo();

                //
                std::vector<std::string> markersInterpolate;
                markersInterpolate = { "Solvent" };
                remeshedFM->init( currentFM, markersInterpolate );
                remeshedFM->startTimeStep();
                currentFM = remeshedFM;

            }// End if mustRemesh

            currentFM->updateParameterValues();

            // I think I need this because of the rigid body, 
            // which imposes a moving boundary of the fluid domain
            // see feel/feelmodels//fluid/fluidmechanicsothers.cpp:879
            currentFM->setApplyMovingMeshBeforeSolve( false );

            // // SWIMMING not needed here for now
            // // ========
            // auto expr_swimming_f = [&FM](){
            // [...]
            // currentFM->updateALEmesh( se );
            // ========

            // SOLVE
            currentFM->solve();
            currentFM->exportResults();
            updateMeasureFile( "fluid_global.measures.csv", currentFM->rootRepository() );
            {
                auto Yh = Pch<1>( currentFM->mesh() );
                auto w = Yh->element();
                w.on( _range = elements( support( Yh ) ), _expr = Px() );
                movExporterSave( globalBDF->time(), currentFM->meshALE(), currentFM->fieldVelocity(), currentFM->fieldPressure(), w );
                referenceExporterSave( globalBDF->time(), currentFM->meshALE() );
            }
            currentFM->updateTimeStep();
        }
    }
    return 0;
}



}// namespace Feel

int
main( int argc, char** argv )
{
    using namespace Feel;
	po::options_description fluidmecoptions( "application fluid-mechanics options" );
    fluidmecoptions.add( toolboxes_options("fluid") );
    fluidmecoptions.add( bdf_options( "mybdf" ) ).add( ts_options( "mybdf" ) );
    fluidmecoptions.add_options()
        ("case.dimension", Feel::po::value<int>()->default_value( 3 ), "dimension")
        ("case.discretization", Feel::po::value<std::string>()->default_value( "P2P1G1" ), "discretization : P2P1G1,P2P1G2")
        ("remesh.maxminareasratiothreshold", Feel::po::value<double>()->default_value( 1.5 ), "Trigger remeshing whenever, for any element of the mesh, the ratio of the largest to the smallest face areas is larger than this threshold")
        ;
    fluidmecoptions.add( backend_options( "Ip" ) ); // TODO check if I need these 4 lines
    fluidmecoptions.add( backend_options( "Iv" ) );
    fluidmecoptions.add( backend_options( "Idisp" ) );
    fluidmecoptions.add( alemesh_options( "alemesh_init" ) );

	Environment env( _argc=argc, _argv=argv,
                     _desc=fluidmecoptions,
                   _about=about(_name="application_fluid",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    int dimension = ioption(_name="case.dimension");
    std::string discretization = soption(_name="case.discretization");
    if ( discretization == "P2P1" )
        discretization = "P2P1G1";

    auto dimt = hana::make_tuple(hana::int_c<2>,hana::int_c<3>);

    auto discretizationt = hana::make_tuple( hana::make_tuple("P2P1G1", hana::make_tuple( hana::int_c<2>,hana::int_c<1>,hana::int_c<1>) ) );

    int status = 0;
    hana::for_each( hana::cartesian_product(hana::make_tuple(dimt,discretizationt)), [&discretization,&dimension,&status]( auto const& d )
                    {
                        constexpr int _dim = std::decay_t<decltype(hana::at_c<0>(d))>::value;
                        std::string const& _discretization = hana::at_c<0>( hana::at_c<1>(d) );
                        constexpr int _uorder = std::decay_t<decltype(hana::at_c<0>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _porder = std::decay_t<decltype(hana::at_c<1>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        constexpr int _gorder = std::decay_t<decltype(hana::at_c<2>(hana::at_c<1>( hana::at_c<1>(d)) ))>::value;
                        if ( dimension == _dim && discretization == _discretization )
                            status = runApplicationFluid<_dim,_uorder,_porder,_gorder>();
                    } );
    return status;
}
