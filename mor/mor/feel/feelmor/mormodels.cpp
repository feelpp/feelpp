
#include <feel/feelmor/mormodels.hpp>

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
    std::cout << "-- crbmodelDB::dbRepository() = " << crbmodelDB.dbRepository() << std::endl;
#if 0
    if ( boption( _name = "export-solution" ) )
        p_ = crbmodelDB.loadDBPlugin( meta, "all" );
    else
#endif    
    p_ = crbmodelDB.loadDBPlugin( meta, this->load );
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
