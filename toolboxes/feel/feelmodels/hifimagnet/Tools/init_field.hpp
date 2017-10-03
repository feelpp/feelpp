#ifndef __INIT_FIELD_HPP
#define __INIT_FIELD_HPP

#include <feel/feelcore/environment.hpp>

namespace Feel
{
    template <typename element_type> bool init_field(element_type & M_V,  const std::string Field, const std::string Init)
    {
        tic();
        bool is_double = false;
        try
            {
                double val = stod(Init);
                is_double = true;
                Feel::cout << "Initialize " << Field << " to " << val << std::endl;
                M_V = vf::project( M_V.functionSpace(), elements(M_V.mesh()), cst(val) );
            }
        catch (const std::invalid_argument& ia)
            {
                // check if Init is an expression (find x:y:z in string) or a filename
                if ( !Init.empty() )
                    {
                        std::string datafile = ( boost::format( "%1%" ) % Environment::expand(Init) ).str();
                        Feel::cout << "Loading " << Field << " init from " << datafile << std::endl;
                        fs::path VInitPath(datafile);
                        if ( fs::exists( VInitPath ) )
                            M_V.loadHDF5( datafile );
                        else
                            throw std::logic_error( "No such file : " + datafile + " for " + Field + " initialization" );
                    }
            }
        return is_double;
    }
}
#endif /* __INIT_FIELD_HPP */
