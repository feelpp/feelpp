/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-07-21

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exporter.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-07-21
 */
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefilters/exportergmsh.hpp>
#include <life/lifefilters/exporterensight.hpp>

namespace Life
{
template<typename MeshType>
Exporter<MeshType>::Exporter( std::string const& __type, std::string const& __prefix, int __freq  )
    :
    super1(),
    super2(),
    M_do_export( true ),
    M_type( __type ),
    M_prefix( __prefix ),
    M_freq( __freq ),
    M_ft( ASCII ),
    M_path( "." )
{

}

template<typename MeshType>
Exporter<MeshType>::Exporter( po::variables_map const& vm, std::string const& exp_prefix )
    :
    super1(),
    super2(),
    M_do_export( true ),
    M_type(),
    M_prefix( exp_prefix ),
    M_freq(1),
    M_ft( ASCII ),
    M_path( "." )
{
    std::cout << "[exporter::exporter] do export = " << doExport() << std::endl;
}

template<typename MeshType>
Exporter<MeshType>::Exporter( Exporter const & __ex )
    :
    super1(),
    super2(),
    M_do_export( __ex.M_do_export ),
    M_type( __ex.M_type ),
    M_prefix( __ex.M_prefix ),
    M_freq( __ex.M_freq ),
    M_ft( __ex.M_ft ),
    M_path( __ex.M_path )
{

}

template<typename MeshType>
Exporter<MeshType>::~Exporter()
{}

template<typename MeshType>
Exporter<MeshType>*
Exporter<MeshType>::New( std::string const& exportername, std::string prefix )
{
    Exporter<MeshType>* exporter =  Factory::type::instance().createObject( exportername  );

    exporter->addTimeSet( timeset_ptrtype( new timeset_type( prefix ) ) );
    exporter->setPrefix( prefix );
    return exporter;
}

template<typename MeshType>
Exporter<MeshType>*
Exporter<MeshType>::New( po::variables_map const& vm, std::string prefix )
{
    std::string estr = vm["exporter"].template as<std::string>();
    Exporter<MeshType>* exporter =  Factory::type::instance().createObject( estr  );
    exporter->setOptions( vm );
    std::cout << "[exporter::New] do export = " << exporter->doExport() << std::endl;
    exporter->addTimeSet( timeset_ptrtype( new timeset_type( prefix ) ) );
    exporter->setPrefix( prefix );
    return exporter;
}

template<typename MeshType>
Exporter<MeshType>*
Exporter<MeshType>::setOptions( po::variables_map const& vm, std::string const& exp_prefix )
{
    std::string _prefix = exp_prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    M_do_export = vm[_prefix+"export"].template as<bool>();
    M_type =  vm[_prefix+"exporter"].template as<std::string>();
    if ( vm.count ( _prefix+"exporter-prefix" ) )
        M_prefix = vm[_prefix+"exporter-prefix"].template as<std::string>();
    M_freq = vm[_prefix+"exporter-freq"].template as<int>();
    M_ft = file_type( vm[_prefix+"exporter-file-type"].template as<int>() );

    Debug() << "[Exporter] type:  " << M_type << "\n";
    Debug() << "[Exporter] prefix:  " << M_prefix << "\n";
    Debug() << "[Exporter] freq:  " << M_freq << "\n";
    Debug() << "[Exporter] ft:  " << M_ft << "\n";
    return this;
}

template<typename MeshType>
Exporter<MeshType>*
Exporter<MeshType>::addPath( boost::format fmt )
{
    fs::path rep_path = ".";
    typedef std::vector< std::string > split_vector_type;

    split_vector_type dirs; // #2: Search for tokens
    std::string fmtstr = fmt.str();
    boost::split( dirs, fmtstr, boost::is_any_of("/") );

    BOOST_FOREACH( std::string const& dir, dirs )
        {
            //Debug( 1000 ) << "[Application::Application] option: " << s << "\n";
            rep_path = rep_path / dir;
            if (!fs::exists( rep_path ) )
                fs::create_directory( rep_path );
        }

    M_path = rep_path.string();

    return this;
}

//
// Exporter Options
//
/**
 * \return the command lines options for the exporter
 */
po::options_description exporter_options( std::string const& prefix )
{
    std::string _prefix = prefix;
    if ( !_prefix.empty() )
        _prefix += "-";

    po::options_description _options( "Exporter " + prefix + " options");
    _options.add_options()
        // do export
        ((_prefix+"export").c_str(), Life::po::value<bool>()->default_value( true ), "true if export, false otherwise")

        // exporter type
        ((_prefix+"exporter").c_str(), Life::po::value<std::string>()->default_value( "ensight" ), "type of exporter")

        // prefix options
        ((_prefix+"exporter-prefix").c_str(), Life::po::value<std::string>()->default_value( prefix ), "prefix for exported files")

        // frequency options
        ((_prefix+"exporter-freq").c_str(), Life::po::value<int>()->default_value( 1 ), "frequency at which results are exported")

        // file type options
        ((_prefix+"exporter-file-type").c_str(), Life::po::value<int>()->default_value( ASCII ), "file type in which the results are exported ('ascii' = 0 or 'binary' = 1)")
        ;
    return _options;
}

}
