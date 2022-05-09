#include <feel/feelcore/environment.hpp>
#include <feel/feelmor/crbmodelproperties.hpp>
#include <feel/feelmor/crbmodelparameters.hpp>
#include <feel/feelmor/crbmodeloutputs.hpp>

namespace Feel {


nl::json
CRBModelProperties::read_json( std::string const& _filename, worldcomm_ptr_t const& world )
{
    std::string filename = Environment::expand( _filename );
    if ( !fs::exists( filename ) )
    {
        if ( world->isMasterRank() )
        {
            std::cout << "[CRBmodelProperties]  Could not find :\"" << filename << "\"" <<std::endl;
        }
        LOG(INFO) << "Could not find " << filename << std::endl;
        return {};
    }
    else
    {
        if ( world->isMasterRank() )
        {
            std::cout << "[CRBmodelProperties] Loading Model Properties : \"" << filename << "\"" << std::endl;
        }
        LOG(INFO) << "Loading " << filename << std::endl;
    }

    std::ifstream ifs(filename);
    return nl::json::parse(ifs,nullptr,true,true);
}

CRBModelProperties::CRBModelProperties( std::string const& directoryLibExpr, worldcomm_ptr_t const& world, std::string const& prefix, po::variables_map const& vm )
    :
    super( world ),
    M_prefix( prefix ),
    M_directoryLibExpr( directoryLibExpr ),
    M_params( world ),
    M_outputs( world )
{
    if ( countoption( _name="json.merge_patch",_prefix=M_prefix,_vm=vm ) > 0 )
    {
        for ( std::string const& patch : vsoption( _name="json.merge_patch",_prefix=M_prefix,_vm=vm ) )
            M_json_merge_patch.merge_patch( json::parse( patch ) );
    }
    if ( countoption( _name="json.patch",_prefix=M_prefix,_vm=vm ) > 0 )
    {
        for ( std::string const& patch : vsoption( _name="json.patch", _prefix=M_prefix,_vm=vm ) )
            M_json_patch.push_back( json::parse( patch ) );
    }
}

void
CRBModelProperties::setup( std::vector<std::string> const& filenames )
{
    nl::json j;
    for ( std::string const& filename : filenames )
        j.merge_patch( CRBModelProperties::read_json( filename, this->worldCommPtr() ) );
    this->setup( std::move( j ) );
}

void
CRBModelProperties::setupFromFilenameOption( po::variables_map const& vm )
{
    if ( countoption( _name="json.filename",_prefix=M_prefix,_vm=vm ) > 0 )
        this->setup( vsoption( _name="json.filename",_prefix=M_prefix,_vm=vm ) );
}

void
CRBModelProperties::setupImpl()
{
    if ( !M_json_merge_patch.is_null() )
        M_jsonData.merge_patch( M_json_merge_patch );
    if ( !M_json_patch.empty() )
        M_jsonData = M_jsonData.patch( M_json_patch );
    nl::json & jarg = M_jsonData;
    //std::cout << jarg.dump(1) << std::endl;
    if ( jarg.contains("Name") )
    {
        auto const& j_name = jarg.at("Name");
        if ( j_name.is_string() )
            M_name = j_name.get<std::string>();
    }
    if ( jarg.contains("ShortName") )
    {
        auto const& j_shortname = jarg.at("ShortName");
        if ( j_shortname.is_string() )
            M_shortname = j_shortname.get<std::string>();
    }

    if ( jarg.contains("CRBParameters") )
    {
        LOG(INFO) << "CRBModel with parameters\n";
        M_params.setPTree( jarg.at("CRBParameters") );
    }

    if ( jarg.contains("CRBOutputs") )
    {
        LOG(INFO) << "Model with outputs\n";
        if ( !M_directoryLibExpr.empty() )
            M_outputs.setDirectoryLibExpr( M_directoryLibExpr );
        M_outputs.setPTree( jarg.at("CRBOutputs") );
    }
}

}

