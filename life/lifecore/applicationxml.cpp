/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-24

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file applicationxml.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-03-24
 */

#include <applicationxml.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

namespace Life
{
namespace fs = boost::filesystem;
namespace detail
{
po::options_description
makeOptions( po::options_description const & opt )
{
    po::options_description xml( "XML application options" );
    xml.add_options()
        ("capabilities", "generate xml file describing the capabilities");

    return xml.add( opt );
}
}
ApplicationXML::ApplicationXML( int argc,
                                char** argv,
                                AboutData const& ad,
                                po::options_description const& od )
    :
    super( argc, argv, ad, detail::makeOptions( od ) ),
    M_params(),
    M_outputs(),
    M_parameter_values(),
    M_output_values()
{
}
ApplicationXML::ApplicationXML( ApplicationXML const& app )
    :
    super( app ),
    M_params( app.M_params ),
    M_outputs( app.M_outputs ),
    M_parameter_values( app.M_parameter_values ),
    M_output_values( app.M_output_values )
{}
ApplicationXML::~ApplicationXML()
{}
ApplicationXML&
ApplicationXML::operator=( ApplicationXML const& app )
{
    if (this != &app )
        {
            M_params = app.M_params;
            M_outputs = app.M_outputs;
            M_parameter_values = app.M_parameter_values;
            M_output_values = app.M_output_values;
        }
    return *this;
}

ApplicationXML::RunStatus
ApplicationXML::preProcessing()
{
    Debug( 1000 ) << "start preprocessing\n";
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return RUN_EXIT;
        }

    if ( this->vm().count( "capabilities" ) )
		{
			Debug( 1000 ) << "Writing capabilities..." << "\n";
            fs::path rep_path;
            std::string fmtstr = (boost::format( "%1%/" ) % "xml").str();
            rep_path = rootRepository();
            rep_path = rep_path / fmtstr;
            if ( !fs::exists( rep_path ) )
                fs::create_directory( rep_path );
            rep_path = rep_path / "xml_response.xml";
			xmlParser::writeResponse( rep_path.string(),
                                      this->about().appName(),
                                      M_params,
                                      M_outputs);
            std::string rep="";
            for (unsigned int i=0; i<M_params.size(); i++) {
                Debug( 1000 ) << "rep = " << rep << "\n";
                rep+=M_params[i].getName();
                rep+="_";
                rep+=M_parameter_values[i];
                rep+="/";
            }
            this->changeRepository( boost::format( "%1%/%2%" )
                                    % this->about().appName()
                                    % rep
                                    );
            Debug( 1000 ) << "Capabilities writed..." << "\n";
			return RUN_EXIT;
		}

    std::string rep="";
    for (unsigned int i=0; i<M_params.size(); i++) {
        Debug( 1000 ) << "rep = " << rep << "\n";
        rep+=M_params[i].getName();
        rep+="_";
        rep+=M_parameter_values[i];
        rep+="/";
    }
	this->changeRepository( boost::format( "%1%/%2%" )
							% this->about().appName()
							% rep
                            );

    return RUN_CONTINUE;
}
void
ApplicationXML::postProcessing()
{
    fs::path rep_path;
    std::string fmtstr = (boost::format( "%1%/" ) % "xml").str();
    rep_path = rootRepository();
    rep_path = rep_path / fmtstr;
    if ( !fs::exists( rep_path ) )
        fs::create_directory( rep_path );
    rep_path = rep_path / "xml_result.xml";
	xmlParser::writeResult( rep_path.string(),
                            this->about().appName(),
                            M_params,
                            M_outputs,
                            M_parameter_values,
                            M_output_values );
    /*string rep="";
    for (unsigned int i=0; i<M_params.size(); i++) {
        Debug( 1000 ) << "rep = " << rep << "\n";
        rep+=M_params[i].getName();
        rep+="_";
        rep+=M_parameter_values[i];
        rep+="/";
    }
	this->changeRepository( boost::format( "%1%/%2%" )
							% this->about().appName()
							% rep
                            );
    */
}



}

