/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-21

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file exporterimpl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-21
 */
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#include <feel/feelfilters/exportergmsh.hpp>
#include <feel/feelfilters/exporterensight.hpp>

namespace Feel
{
template<typename MeshType, int N> class ExporterEnsight;
template<typename MeshType, int N> class ExporterGmsh;

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( std::string const& __type, std::string const& __prefix, int __freq, WorldComm const& worldComm  )
    :
    super1(),
    super2(),
    M_worldComm( worldComm ),
    M_do_export( true ),
    M_type( __type ),
    M_prefix( __prefix ),
    M_freq( __freq ),
    M_cptOfSave( 0 ),
    M_ft( ASCII ),
    M_path( "." )
{

}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( po::variables_map const& vm, std::string const& exp_prefix, WorldComm const& worldComm )
    :
    super1(),
    super2(),
    M_worldComm( worldComm ),
    M_do_export( true ),
    M_type(),
    M_prefix( exp_prefix ),
    M_freq( 1 ),
    M_cptOfSave( 0 ),
    M_ft( ASCII ),
    M_path( "." )
{
    VLOG(1) << "[exporter::exporter] do export = " << doExport() << "\n";
}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( Exporter const & __ex )
    :
    super1(),
    super2(),
    M_worldComm( __ex.M_worldComm ),
    M_do_export( __ex.M_do_export ),
    M_type( __ex.M_type ),
    M_prefix( __ex.M_prefix ),
    M_freq( __ex.M_freq ),
    M_cptOfSave( __ex.M_cptOfSave ),
    M_ft( __ex.M_ft ),
    M_path( __ex.M_path )
{

}

template<typename MeshType, int N>
Exporter<MeshType, N>::~Exporter()
{}

template<typename MeshType, int N>
Exporter<MeshType, N>*
Exporter<MeshType, N>::New( std::string const& exportername, std::string prefix, WorldComm const& worldComm )
{
    Exporter<MeshType, N>* exporter =  0;//Factory::type::instance().createObject( exportername  );

    if ( N == 1 && ( exportername == "ensight" || Environment::numberOfProcessors() > 1 ) )
        exporter = new ExporterEnsight<MeshType, N>;
    else if ( N > 1 || ( exportername == "gmsh" ) )
        exporter = new ExporterGmsh<MeshType,N>;
    else // fallback
        exporter = new ExporterEnsight<MeshType, N>;

    exporter->addTimeSet( timeset_ptrtype( new timeset_type( prefix ) ) );
    exporter->setPrefix( prefix );
    return exporter;
}

template<typename MeshType, int N>
Exporter<MeshType, N>*
Exporter<MeshType, N>::New( po::variables_map const& vm, std::string prefix, WorldComm const& worldComm )
{
    std::string estr = vm["exporter.format"].template as<std::string>();
    Exporter<MeshType, N>* exporter =  0;//Factory::type::instance().createObject( estr  );

    if ( N == 1 && ( estr == "ensight"  || Environment::numberOfProcessors() > 1 ) )
        exporter = new ExporterEnsight<MeshType, N>( vm,prefix,worldComm );
    else if ( N > 1 || estr == "gmsh" )
        exporter = new ExporterGmsh<MeshType,N>;
    else // fallback
        exporter = new ExporterEnsight<MeshType, N>( vm,prefix,worldComm );

    exporter->setOptions( vm );
    //std::cout << "[exporter::New] do export = " << exporter->doExport() << std::endl;
    exporter->addTimeSet( timeset_ptrtype( new timeset_type( prefix ) ) );
    exporter->setPrefix( prefix );
    return exporter;
}

template<typename MeshType, int N>
Exporter<MeshType, N>*
Exporter<MeshType, N>::setOptions( po::variables_map const& vm, std::string const& exp_prefix )
{
    std::string _prefix = exp_prefix;

    if ( !_prefix.empty() )
        _prefix += "-";

    //M_do_export = vm[_prefix+"export"].template as<bool>();
    M_do_export = vm[_prefix+"exporter.export"].template as<bool>();
    M_type =  vm[_prefix+"exporter.format"].template as<std::string>();

    if ( vm.count ( _prefix+"exporter.prefix" ) )
        M_prefix = vm[_prefix+"exporter.prefix"].template as<std::string>();

    M_freq = vm[_prefix+"exporter.freq"].template as<int>();
    M_ft = file_type( vm[_prefix+"exporter.file-type"].template as<int>() );

    VLOG(1) << "[Exporter] type:  " << M_type << "\n";
    VLOG(1) << "[Exporter] prefix:  " << M_prefix << "\n";
    VLOG(1) << "[Exporter] freq:  " << M_freq << "\n";
    VLOG(1) << "[Exporter] ft:  " << M_ft << "\n";
    return this;
}

template<typename MeshType, int N>
Exporter<MeshType, N>*
Exporter<MeshType, N>::addPath( boost::format fmt )
{
    fs::path rep_path = ".";
    typedef std::vector< std::string > split_vector_type;

    split_vector_type dirs; // #2: Search for tokens
    std::string fmtstr = fmt.str();
    boost::split( dirs, fmtstr, boost::is_any_of( "/" ) );

    BOOST_FOREACH( std::string const& dir, dirs )
    {
        //Debug( 1000 ) << "[Application::Application] option: " << s << "\n";
        rep_path = rep_path / dir;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
    }

    M_path = rep_path.string();

    return this;
}


}

