/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-21

  Copyright (C) 2010-2014 Feel++ Consortium

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
   \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
   \date 2010-04-21
 */
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/factory.hpp>
#include <feel/feelcore/singleton.hpp>

#ifdef FEELPP_HAS_GMSH
#include <feel/feelfilters/exportergmsh.hpp>
#endif
#include <feel/feelfilters/exporterensight.hpp>

#ifdef FEELPP_HAS_MPIIO
#include <feel/feelfilters/exporterensightgold.hpp>
#endif

#if defined(FEELPP_HAS_VTK)
#include <feel/feelfilters/exportervtk.hpp>
#endif

#if defined(FEELPP_HAS_HDF5)
#include <feel/feelfilters/exporterxdmf.hpp>
#endif

#include <feel/feelfilters/exporterexodus.hpp>



namespace Feel
{
template<typename MeshType, int N> class ExporterEnsight;

#if defined(FEELPP_HAS_MPIIO)
template<typename MeshType, int N> class ExporterEnsightGold;
#endif

#ifdef FEELPP_HAS_GMSH
template<typename MeshType, int N> class ExporterGmsh;
#endif
#if defined(FEELPP_HAS_HDF5)
//template<typename MeshType, int N> class ExporterXDMF;
#endif

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    super1(),
    super2(),
    M_do_export( true ),
    M_use_single_transient_file( false ),
    M_type(),
    M_prefix( Environment::about().appName() ),
    M_freq( 1 ),
    M_ft( ASCII ),
    M_path( "." ),
    M_ex_geometry( EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
{
    VLOG(1) << "[exporter::exporter] do export = " << doExport() << "\n";
}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( std::string const& __type, std::string const& __prefix, int __freq, worldcomm_ptr_t const& worldComm  )
    :
    super( worldComm ),
    super1(),
    super2(),
    M_do_export( true ),
    M_use_single_transient_file( false ),
    M_type( __type ),
    M_prefix( __prefix ),
    M_freq( __freq ),
    M_ft( ASCII ),
    M_path( "." ),
    M_ex_geometry( EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
{

}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( po::variables_map const& vm, std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    super1(),
    super2(),
    M_do_export( true ),
    M_use_single_transient_file( false ),
    M_type(),
    M_prefix( exp_prefix ),
    M_freq( 1 ),
    M_ft( ASCII ),
    M_path( "." ),
    M_ex_geometry( EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
{
    VLOG(1) << "[exporter::exporter] do export = " << doExport() << "\n";
}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( std::string const& exp_prefix, worldcomm_ptr_t const& worldComm )
    :
    super( worldComm ),
    super1(),
    super2(),
    M_do_export( true ),
    M_use_single_transient_file( false ),
    M_type(),
    M_prefix( exp_prefix ),
    M_freq( 1 ),
    M_ft( ASCII ),
    M_path( "." ),
    M_ex_geometry( EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
{
    VLOG(1) << "[exporter::exporter] do export = " << doExport() << "\n";
}

template<typename MeshType, int N>
Exporter<MeshType, N>::Exporter( Exporter const & __ex )
    :
    super( __ex ),
    super1(),
    super2(),
    M_do_export( __ex.M_do_export ),
    M_use_single_transient_file( __ex.M_use_single_transient_file ),
    M_type( __ex.M_type ),
    M_prefix( __ex.M_prefix ),
    M_freq( __ex.M_freq ),
    M_ft( __ex.M_ft ),
    M_path( __ex.M_path ),
    M_ex_geometry( EXPORTER_GEOMETRY_CHANGE_COORDS_ONLY )
{

}

template<typename MeshType, int N>
Exporter<MeshType, N>::~Exporter()
{}

template<typename MeshType, int N>
std::shared_ptr<Exporter<MeshType, N> >
Exporter<MeshType, N>::New( std::string const& exportername, std::string prefix, worldcomm_ptr_t const& worldComm )
{
    Exporter<MeshType, N>* exporter =  0;//Factory::type::instance().createObject( exportername  );

    if ( N == 1 && ( exportername == "ensight" ) )
        exporter = new ExporterEnsight<MeshType, N>( worldComm );
#if defined(FEELPP_HAS_MPIIO)
    else if ( N <= 2 && ( exportername == "ensightgold"  ) )
        exporter = new ExporterEnsightGold<MeshType, N>( worldComm );
#endif
    else if ( N == 1 && ( exportername == "exodus"  ) )
        exporter = new ExporterExodus<MeshType, N>( worldComm );
#if defined(FEELPP_HAS_HDF5)
    else if ( N == 1 && ( exportername == "xdmf" ))
        exporter = new ExporterXDMF<MeshType, N> ( worldComm ) ;
#endif
#if defined(FEELPP_HAS_VTK)
    else if ( N == 1 && ( exportername == "vtk"  ) )
        exporter = new ExporterVTK<MeshType, N>( worldComm );
#endif
#ifdef FEELPP_HAS_GMSH
    else if ( N > 1 || ( exportername == "gmsh" ) )
        exporter = new ExporterGmsh<MeshType,N>( worldComm );
#endif
    else // fallback
    {
        LOG(INFO) << "[Exporter] The exporter format " << exportername << " Cannot be found. Falling back to Ensight exporter." << std::endl;
        exporter = new ExporterEnsight<MeshType, N>( worldComm );
    }

    exporter->addTimeSet( timeset_ptrtype( new timeset_type( prefix ) ) );
    exporter->setPrefix( prefix );
    return std::shared_ptr<Exporter<MeshType, N> >(exporter);
}

template<typename MeshType, int N>
std::shared_ptr<Exporter<MeshType, N> >
Exporter<MeshType, N>::New( po::variables_map const& vm, std::string prefix, worldcomm_ptr_t const& worldComm )
{
    return New( prefix, worldComm );
}
template<typename MeshType, int N>
std::shared_ptr<Exporter<MeshType, N> >
Exporter<MeshType, N>::New( std::string prefix, worldcomm_ptr_t const& worldComm )
{
    std::string estr = soption("exporter.format");
    std::shared_ptr<Exporter<MeshType, N> > exporter;

    LOG(INFO) << "[Exporter] format :  " << estr << "\n";
    LOG(INFO) << "[Exporter] N      :  " << N << "\n";
    if( N > 1 && estr != "gmsh" )
        LOG(WARNING) << "[Exporter] format " << estr << " is not available for mesh order > 1 - using gmsh exporter instead\n";

    if ( N == 1 && ( estr == "ensight"   ) )
        exporter = std::make_shared<ExporterEnsight<MeshType, N>>( prefix, worldComm );
#if defined(FEELPP_HAS_MPIIO)
    else if ( N <= 2 && ( estr == "ensightgold"   ) )
        exporter = std::make_shared<ExporterEnsightGold<MeshType, N>>( prefix, worldComm );
#endif
    else if ( N == 1 && ( estr == "exodus"   ) )
        exporter = std::make_shared<ExporterExodus<MeshType, N>>( prefix, worldComm );
#if defined(FEELPP_HAS_HDF5)
    else if ( N == 1 && ( estr == "xdmf" ) )
        exporter = std::make_shared<ExporterXDMF<MeshType, N>> ( prefix, worldComm ) ;
#endif
#if defined(FEELPP_HAS_VTK)
    else if ( N == 1 && ( estr == "vtk"  ) )
        exporter = std::make_shared<ExporterVTK<MeshType, N>>( prefix, worldComm );
#endif
#ifdef FEELPP_HAS_GMSH
    else if ( N > 1 || estr == "gmsh" )
        exporter = std::make_shared<ExporterGmsh<MeshType,N>>( prefix, worldComm );
#endif
    else // fallback
    {
        LOG(INFO) << "[Exporter] The exporter format " << estr << " Cannot be found. Falling back to Ensight exporter." << std::endl;
        exporter = std::make_shared<ExporterEnsight<MeshType, N>>( prefix, worldComm );
    }


    exporter->setOptions();
    //std::cout << "[exporter::New] do export = " << exporter->doExport() << std::endl;
    exporter->addTimeSet( std::make_shared<timeset_type>( prefix ) );
    exporter->setPrefix( prefix );
    return exporter;
}

template<typename MeshType, int N>
Exporter<MeshType, N>*
Exporter<MeshType, N>::setOptions( std::string const& exp_prefix )
{
    //M_do_export = Environment::vm(_prefix+"export"].template as<bool>();
    M_do_export = Environment::vm(_name="exporter.export",_prefix=exp_prefix).template as<bool>();
    M_type =  Environment::vm(_name="exporter.format",_prefix=exp_prefix).template as<std::string>();

    std::string p = exp_prefix;
    if ( !p.empty() )
        p += ".";
    if ( Environment::vm().count( p+"exporter.prefix" ) )
        M_prefix = Environment::vm(_name="exporter.prefix",_prefix=exp_prefix).template as<std::string>();

    M_freq = Environment::vm(_name="exporter.freq",_prefix=exp_prefix).template as<int>();
    std::string ftstr = soption(_name="exporter.file-type",_prefix=exp_prefix);
    if ( ftstr == "binary" )
        M_ft = BINARY;
    else
        M_ft = ASCII;

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

    for( auto const& dir: dirs )
    {
        //DVLOG(2) << "[Application::Application] option: " << s << "\n";
        rep_path = rep_path / dir;

        if ( !fs::exists( rep_path ) )
            fs::create_directory( rep_path );
    }

    M_path = rep_path.string();

    return this;
}

template<typename MeshType, int N>
void
Exporter<MeshType, N>::serve() const
{

#if 0
    ./bin/startRenderingServer.js --pvpython /Applications/ParaView-5.5.0.app/Contents/bin/pvpython --data-directory-path data --data-load-signature-decoder b\
        > in/decodeDataLoadSignature.js
#endif
}


}
