#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelmor/crbmodelparameters.hpp>

namespace Feel {

CRBModelParameters::CRBModelParameters( worldcomm_ptr_t const& world )
    :
    super( world )
{}

CRBModelParameters::~CRBModelParameters()
{}

void
CRBModelParameters::setPTree( nl::json const& jarg )
{
    M_p = jarg;
    setup();
}


void
CRBModelParameters::setup()
{
    for ( auto const& [jargkey,jargval] : M_p.items() )
    {
        std::string const& t = jargkey; // parameter name
        LOG(INFO) << "reading parameter " << t;

        if ( jargval.is_object() )
        {
            std::string name = t;
            if ( jargval.contains("name") )
            {
                auto const& j_name = jargval.at("name");
                if ( j_name.is_string() )
                    name = j_name.get<std::string>();
            }
            std::string desc;
            if ( jargval.contains("description") )
            {
                auto const& j_description = jargval.at("description");
                if ( j_description.is_string() )
                    desc = j_description.get<std::string>();
            }
            std::string sampling;
            if ( jargval.contains("sampling") )
            {
                auto const& j_sampling = jargval.at("sampling");
                if ( j_sampling.is_string() )
                    sampling = j_sampling.get<std::string>();
            }

            double min = 0.;
            double max = 0.;
            if ( jargval.contains("min") )
            {
                auto const& j_min = jargval.at("min");
                if ( j_min.is_number() )
                    min = j_min.get<double>();
                else if ( j_min.is_string() )
                    min = std::stod( j_min.get<std::string>() );
                LOG( INFO ) << fmt::format( "parameter {} has min value {}", name, min );
            }
            else
            {
                LOG(INFO) << "parameter " << name << " does not have min, skipping it";
                continue;
            }
            if ( jargval.contains("max") )
            {
                auto const& j_max = jargval.at("max");
                if ( j_max.is_number() )
                    max = j_max.get<double>();
                else if ( j_max.is_string() )
                    max = std::stod( j_max.get<std::string>() );
                LOG( INFO ) << fmt::format( "parameter {} has max value {}", name, max );
            }
            else
            {
                LOG(INFO) << "parameter " << name << " does not have max, skipping it";
                continue;
            }
 
            double value = 0.;
            if ( jargval.contains("value") )
            {
                auto const& j_val = jargval.at("value");
                if ( j_val.is_number() )
                    value = j_val.get<double>();
                else if ( j_val.is_string() )
                    value = std::stod( j_val.get<std::string>() );
                LOG(INFO) << fmt::format("parameter {} has value {}", name, value);
            }
            int sample = 0;
            if ( jargval.contains("sampling-size") )
            {
                auto const& j_sample = jargval.at("sampling-size");
                if ( j_sample.is_number() )
                    sample = j_sample.get<int>();
                else if ( j_sample.is_string() )
                    sample = std::stoi( j_sample.get<std::string>() );
            }
            this->operator[](name) = CRBModelParameter( name, min, max, value, sample, sampling, desc );
        }
    }
}


} // namespace Feel
