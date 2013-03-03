/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-24

  Copyright (C) 2009-2011 Université de Grenoble 1

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-24
 */

#include <feel/feelcore/applicationxml.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

namespace Feel
{
namespace fs = boost::filesystem;
namespace detail
{
po::options_description
makeOptions()
{
    po::options_description xml( "XML application options" );
    xml.add_options()
    ( "capabilities", "generate xml file describing the capabilities" );

    return xml;
}
}
ApplicationXML::ApplicationXML( int argc,
                                char** argv,
                                AboutData const& ad,
                                po::options_description const& od )
    :
    super( argc, argv, ad, detail::makeOptions().add( od ) ),
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
    if ( this != &app )
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
    DVLOG(2) << "start preprocessing\n";

    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return RUN_EXIT;
    }

    if ( this->vm().count( "capabilities" ) )
    {
        DVLOG(2) << "Writing capabilities..." << "\n";
        std::cerr << "Writing capabilities..." << "\n";
        fs::path rep_path;
        std::string fmtstr = ( boost::format( "%1%/" ) % "xml" ).str();
        rep_path = rootRepository();
        rep_path = rep_path / fmtstr;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );

        rep_path = rep_path / "xml_response.xml";
        std::cerr << "Writing xml..." << "\n";
        xmlParser::writeResponse( rep_path.string(),
                                  this->about().appName(),
                                  M_params,
                                  M_outputs );
        std::cerr << "preparing repository..." << "\n";
        std::string rep="";

        for ( unsigned int i=0; i<M_params.size(); i++ )
        {
            DVLOG(2) << "rep = " << rep << "\n";
            rep+=M_params[i].getName();
            rep+="_";
            rep+=M_parameter_values[i];
            rep+="/";
        }

        std::cerr << "changing to repository" << rep << "\n";
        this->changeRepository( boost::format( "%1%/%2%" )
                                % this->about().appName()
                                % rep
                              );
        DVLOG(2) << "Capabilities written..." << "\n";
        std::cerr << "Capabilities written..." << "\n";
        return RUN_EXIT;
    }

    std::string rep="";

    for ( unsigned int i=0; i<M_params.size(); i++ )
    {
        DVLOG(2) << "rep = " << rep << "\n";
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
    std::string fmtstr = ( boost::format( "%1%/" ) % "xml" ).str();
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
        DVLOG(2) << "rep = " << rep << "\n";
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

