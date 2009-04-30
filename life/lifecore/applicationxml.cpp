/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-03-24

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

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


namespace Life
{
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
    std::cout << "start preprocessing\n";
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return RUN_EXIT;
        }

    if ( this->vm().count( "capabilities" ) )
		{
			std::cout << "Writing capabilities..." << "\n";
			this->changeRepository( boost::format( "%1%/" )
									% "xml" );
			xmlParser::writeResponse("xml_response.xml",
                                     this->about().appName(),
                                     M_params,
                                     M_outputs);
            string rep="";
            for (unsigned int i=0; i<M_params.size(); i++) {
                std::cout << "rep = " << rep << "\n";
                rep+=M_params[i].getName();
                rep+="_";
                rep+=M_parameter_values[i];
                rep+="/";
            }
            this->changeRepository( boost::format( "%1%/%2%" )
                                    % this->about().appName()
                                    % rep
                                    );
			return RUN_EXIT;
		}

    string rep="";
    for (unsigned int i=0; i<M_params.size(); i++) {
        std::cout << "rep = " << rep << "\n";
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
    this->changeRepository( boost::format( "%1%/")
                            % "xml"
                            );


	xmlParser::writeResult( "xml_result.xml",
                            this->about().appName(),
                            M_params,
                            M_outputs,
                            M_parameter_values,
                            M_output_values );
    string rep="";
    for (unsigned int i=0; i<M_params.size(); i++) {
        std::cout << "rep = " << rep << "\n";
        rep+=M_params[i].getName();
        rep+="_";
        rep+=M_parameter_values[i];
        rep+="/";
    }
	this->changeRepository( boost::format( "%1%/%2%" )
							% this->about().appName()
							% rep
                            );
}



}

