
#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/feelppdatabase.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/exporter.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	po::options_description exportdboptions( "export database options" );
	exportdboptions.add_options()
        ( "ifile", po::value<std::string>(), "input database file" )
        ( "add-partitioning", po::value<bool>()->default_value( false ), "add partitioning " )
		;

	Environment env( _argc=argc, _argv=argv,
                     _desc=exportdboptions,
                   _about=about(_name="export_database",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    std::string ifile;
    if ( Environment::vm().count("ifile") )
        ifile = Environment::vm()["ifile"].as<std::string>();
    if ( ifile.empty() )
    {
        if ( Environment::isMasterRank() )
            std::cout << "do nothing because --ifile is missing\n";
        return 0;
    }
    bool addPartitioning = boption(_name="add-partitioning" );

    FeelppDatabase<Mesh<Simplex<3>>> myDb;
    myDb.setFilename( ifile );
    myDb.loadInfo();

    auto mesh = myDb.loadMesh( MESH_UPDATE_FACES_MINIMAL|MESH_NO_UPDATE_MEASURES );

    auto VhScal = Pch<1>( mesh );
    auto VhVec = Pchv<1>( mesh );
    auto uScal = VhScal->element();
    auto uVec = VhVec->element();
    std::string geoExportType = "static";
    auto e = exporter( _mesh=mesh,_geo=geoExportType );
    auto const& fieldsInfo = myDb.fieldsInfo();
    auto const& timeSet = myDb.timeSet();

    for ( int i=0;i<timeSet.size();++i )
    {
        double time = timeSet[i];
        if ( myDb.worldComm().isMasterRank() )
            std::cout << "reload db at time " << time << "\n";
        for (auto const& fieldsInfo : fieldsInfo )
        {
            std::string fieldName = fieldsInfo.first;
            int nComp = std::get<2>( fieldsInfo.second );
            if ( nComp == 1 )
            {
                myDb.load( i,fieldName,uScal );
                e->step(time)->add( fieldName, uScal );
            }
            else if ( nComp == mesh->nRealDim )
            {
                myDb.load( i,fieldName,uVec );
                e->step(time)->add( fieldName, uVec );
            }
        }
        if ( addPartitioning )
            e->step(time)->addRegions();
        e->save();
    }

    return 0;
}
