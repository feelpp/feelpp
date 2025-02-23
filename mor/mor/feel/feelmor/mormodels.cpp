
#include <feel/feelmor/mormodels.hpp>
#include <feel/feelmor/zip.hpp>

namespace Feel
{

void
MORModel::loadPlugin()
{
    dbroot_path = fs::path( Environment::expand( this->dbroot ) );
    std::string crbmodelName = Environment::expand( this->name );
    CRBModelDB crbmodelDB{ this->output, uuids::nil_uuid(), dbroot_path.string() };

    std::string attribute = this->attribute;
    std::string attribute_data;
    if ( attribute == "id" || attribute == "name" )
    {
        attribute_data = Environment::expand( soption( _name = fmt::format( "crbmodel.db.{}", attribute ) ) );
    }
    else if ( attribute == "last_created" || attribute == "last_modified" )
    {
        std::vector<std::string> split_;
        boost::split( split_, attribute, boost::is_any_of( "_" ) );
        attribute_data = split_[1];
    }
    else
    {
        throw std::runtime_error( "no crbmodel selection, crbmodel.db.id or crbmodel.db.last should be defined" );
    }
    auto meta = crbmodelDB.loadDBMetaData( attribute, attribute_data );
    VLOG(2) << "-- crbmodelDB::dbRepository() = " << crbmodelDB.dbRepository() << std::endl;
#if 0
    if ( boption( _name = "export-solution" ) )
        p_ = crbmodelDB.loadDBPlugin( meta, "all" );
    else
#endif
    p_ = crbmodelDB.loadDBPlugin( meta, this->load );
}
MORModels::MORModels( fs::path const& fpp )
{
    auto _load = [this]( fs::path const& fpp ) {
            std::ifstream f( fpp / "mormodels.json" );
            nl::json j = nl::json::parse( f );
            f.close();
            *this = j["mormodels"].get<MORModels>();
            this->load();
    };
    if ( !fs::exists( fpp ) )
    {
        throw std::invalid_argument( fmt::format( "[MORModels::MORModels] resource {} does not exist", fpp.string() ) );
    }
    if ( fs::is_directory(fpp) && fs::exists(fpp/"mormodels.json") )
    {
        // if there is a cfg file in the directory, use it
        fs::path cfg = fs::path(fpp.filename())/fs::path(fpp.stem().string() + ".cfg");
        Environment::setConfigFile( cfg.string() );
        VLOG(1) << fmt::format("-- MORModels::MORModels: setConfigFile : {} done", cfg.string() ) << std::endl;
        _load( fpp );
    }
    else if ( fpp.extension() != ".fpp")
    {
        throw std::invalid_argument( fmt::format( "[MORModels::MORModels] the zip file should have a .fpp extension", fpp.string() ) );
    }
    else
    {
        extractZipFile(fpp.string(), fpp.parent_path().string() );
        VLOG(2) << fmt::format("-- MORModels::MORModels: extractZipFile : {} done", fpp.string()) << std::endl;
        fs::path cfg = fpp.parent_path() / fpp.stem() / fs::path(fpp.stem().string() + ".cfg");
        Environment::setConfigFile( cfg.string() );
        fs::path desc = fpp.parent_path() / fpp.stem() / "mormodels.json";
        VLOG(2) << fmt::format("-- MORModels::MORModels: Looking for  {} ", desc.string() ) << std::endl;
        if ( !fs::exists(desc) )
        {
            throw std::invalid_argument( fmt::format( "[MORModels::MORModels] FPP file {} does not contain a mormodels.json file : {}", fpp.string(), desc.string() ) );
        }
        _load( fpp.parent_path() / fpp.stem() );
        cleanupTemporaryDirectory( ( fpp.parent_path() / fpp.stem() ).string() );
    }
    if ( Environment::isMasterRank() )
        std::cout << fmt::format( "-- model {} loaded from {} with {} outputs", this->operator[](0).name, fpp.string(), this->size() ) << std::endl;
    //throw std::invalid_argument( fmt::format( "[MORModels::MORModels] FPP file {} does not exist", fppfilename ) );
}
std::vector<std::vector<CRBResults>>
MORModels::run( std::shared_ptr<ParameterSpaceX::Sampling> const& sampling, nl::json const& data ) const
{
    using namespace Feel;

    auto N = data.value("N",-1);
    Eigen::VectorXd /*typename crb_type::vectorN_type*/ time_crb;
    auto online_tol = data.value( "tolerance", 1e-2 );
    auto print_rb_matrix = data.value( "print_rb_matrix", false );

    std::vector<std::vector<CRBResults>> results( sampling->size() );
    for ( auto const& [k, mu] : enumerate( *sampling ) )
    {
        results[k].reserve( this->size() );
        for ( auto const& p : *this )
        {
            tic();
            results[k].emplace_back( p.run( mu, time_crb, online_tol, N, print_rb_matrix ) );
            double t = toc( fmt::format( "rb-online-{}-{}", p.name, p.output ), FLAGS_v > 0 );
        }
        for ( auto const& o : observers_ )
        {
            o->update( std::pair{ k, mu }, results[k] );
        }
    } // end for sample k

    return results;
}

} // namespace Feel
