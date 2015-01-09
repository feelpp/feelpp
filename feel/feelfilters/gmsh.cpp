/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4 tw=0

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-02-10

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)
  Copyright (C) 2008-2014 Feel++ Consortium

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
   \file gmsh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-02-10
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include  <boost/preprocessor/punctuation/paren.hpp>
#include  <boost/preprocessor/punctuation/comma.hpp>
#include <boost/regex.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/concept_check.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelmesh/hypercube.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshsimplexdomain.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelfilters/gmshellipsoiddomain.hpp>

#include <feel/feelcore/fmemopen.h>

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */


#if defined( FEELPP_HAS_GMSH_H )

#include <GmshConfig.h>
#include <Gmsh.h>
#include <GModel.h>
#include <OpenFile.h>
#include <GmshDefines.h>
#include <Context.h>


int PartitionMesh( GModel *const model, meshPartitionOptions &options );

int gmsh_yyparse();
int gmsh_yylex();
void gmsh_yyflush();
int gmsh_yy_scan_string ( const char *str );
extern std::string gmsh_yyname;
extern int gmsh_yylineno;
extern FILE* gmsh_yyin;
extern int gmsh_yyerrorstate;
extern int gmsh_yyviewindex;
extern char *gmsh_yytext;
void PrintParserSymbols(bool help, std::vector<std::string> &vec);

#endif // FEELPP_HAS_GMSH_H


namespace Feel
{
namespace fs = boost::filesystem;

const char* FEELPP_GMSH_FORMAT_VERSION = "2.2";
#if defined(HAVE_METIS)
const GMSH_PARTITIONER GMSH_PARTITIONER_DEFAULT = GMSH_PARTITIONER_METIS;
#else
const GMSH_PARTITIONER GMSH_PARTITIONER_DEFAULT = GMSH_PARTITIONER_CHACO;
#endif


int
ParseGeoFromMemory( GModel* model, std::string const& name, std::string const& geo  )
{
#if defined(FEELPP_HAS_GMSH_H)
    gmsh_yyname = name;
    gmsh_yylineno = 1;
    gmsh_yyerrorstate = 0;
    gmsh_yyviewindex = 0;
    gmsh_yyin = fmemopen( (void*)geo.c_str(), geo.size()*sizeof(char), "r");

    while(!feof(gmsh_yyin))
    {
        gmsh_yyparse();
        if(gmsh_yyerrorstate > 20)
        {
          if(gmsh_yyerrorstate != 999) // 999 is a volontary exit
            std::cerr << "Too many errors: aborting parser...\n";
          gmsh_yyflush();
          break;
        }
    }
    fclose(gmsh_yyin);

    int imported = model->importGEOInternals();
#else
    LOG(ERROR) << "Invalid call: Gmsh is not available on this system.";
    int imported = 0;
#endif /* FEELPP_HAS_GMSH_H */

    return imported;
}
Gmsh::Gmsh( int nDim, int nOrder, WorldComm const& worldComm )
    :
    M_worldComm( worldComm ),
    M_dimension( nDim ),
    M_order( nOrder ),
    M_version( FEELPP_GMSH_FORMAT_VERSION ),
    M_format( GMSH_FORMAT_ASCII ),
    M_in_memory( false ),
    M_I( nDim ),
    M_h( 0.1 ),
    M_addmidpoint( true ),
    M_usePhysicalNames( false ),
    M_partitioner( (GMSH_PARTITIONER)GMSH_PARTITIONER_DEFAULT ),
    M_partitions( worldComm.size() ),
    M_partition_file( 0 ),
    M_shear( 0 ),
    M_recombine( 0 ),
    //M_structured( false ),
    M_structured( 2 ),
    M_refine_levels( 0 ),
    M_substructuring( false ),
    M_periodic(),
    M_gmodel( new GModel() )
{
    this->setReferenceDomain();
    setX( std::make_pair( 0., 1.) );
    if ( nDim > 1 )
        setY( std::make_pair( 0., 1.) );
    if ( nDim > 2 )
        setZ( std::make_pair( 0., 1.) );
}
Gmsh::Gmsh( Gmsh const & __g )
    :
    M_worldComm( __g.M_worldComm ),
    M_dimension( __g.M_dimension ),
    M_order( __g.M_order ),
    M_version( __g.M_version ),
    M_format( __g.M_format ),
    M_geoParamMap( __g.M_geoParamMap ),
    M_in_memory( __g.M_in_memory ),
    M_I( __g.M_I ),
    M_h( __g.M_h ),
    M_addmidpoint( __g.M_addmidpoint ),
    M_usePhysicalNames( __g.M_usePhysicalNames ),
    M_partitioner( __g.M_partitioner ),
    M_partitions( __g.M_partitions ),
    M_partition_file( __g.M_partition_file ),
    M_shear( __g.M_shear ),
    M_recombine( __g.M_recombine ),
    M_structured( __g.M_structured ),
    M_refine_levels( __g.M_refine_levels ),
    M_substructuring( __g.M_substructuring ),
    M_periodic( __g.M_periodic ),
    M_gmodel( __g.M_gmodel )
{}
Gmsh::~Gmsh()
{
    M_gmodel->deleteMesh();
    M_gmodel->destroy();
    delete M_gmodel;
}

Gmsh&
Gmsh::operator=( Gmsh const& __g )
{
    if (  this != &__g )
    {
        M_dimension = __g.M_dimension;
        M_order = __g.M_order;
        M_version = __g.M_version;
        M_format = __g.M_format;
        M_geoParamMap = __g.M_geoParamMap;
        M_in_memory = __g.M_in_memory;
        M_addmidpoint = __g.M_addmidpoint;
        M_usePhysicalNames = __g.M_usePhysicalNames;
        M_shear = __g.M_shear;
        M_refine_levels = __g.M_refine_levels;
        M_periodic = __g.M_periodic;
        M_worldComm = __g.M_worldComm;
        M_gmodel = __g.M_gmodel;
    }

    return *this;
}

boost::shared_ptr<Gmsh>
Gmsh::New( po::variables_map const& vm )
{
    std::ostringstream ostr;
    ostr << vm["convex"].as<std::string>() << "(" << vm["dim"].as<int>() << "," << vm["order"].as<int>() << ")";
    boost::shared_ptr<Gmsh> gmsh_ptr( Gmsh::Factory::type::instance().createObject( ostr.str() ) );
    return gmsh_ptr;
}

boost::shared_ptr<Gmsh>
Gmsh::New( std::string const& shape, uint16_type d, uint16_type o, std::string const& ct, WorldComm const& wc )
{
    std::ostringstream ostr;

    if ( shape != "hypercube" && ct != "hypercube" )
        ostr << shape << "(" << d << "," << o << ")";

    else
        ostr << shape << "(" << d << "," << o << "," << boost::to_lower_copy( ct ) << ")";

    boost::shared_ptr<Gmsh> gmsh_ptr( Gmsh::Factory::type::instance().createObject( ostr.str() ) );
    gmsh_ptr->setWorldComm( wc );
    return gmsh_ptr;
}

std::string
Gmsh::prefix( std::string const& __name, uint16_type __N ) const
{
    std::ostringstream __p;
    __p << __name << "-" << __N << "-" << this->order();

    return __p.str();
}

std::string
Gmsh::getDescriptionFromFile( std::string const& file ) const
{
    if ( !fs::exists( file ) )
    {
        std::ostringstream ostr;
        ostr << "File " << file << " does not exist";
        throw std::invalid_argument( ostr.str() );
    }

    std::ifstream __geoin( file.c_str() );

    std::ostringstream __geostream;
    std::istreambuf_iterator<char> src( __geoin.rdbuf() );
    std::istreambuf_iterator<char> end;
    std::ostream_iterator<char> dest( __geostream );

    std::copy( src,end,dest );

    __geoin.close();
    return  __geostream.str();
}

//
// We split the geo parameters list in two times. First at the `:`, then
// for each character `=`. We add each key and value to the map.
// Example : "key1=val1:key2=val2:key3=val3"
//
std::map<std::string, std::string>
Gmsh::gpstr2map( std::string const& _geopars )
{
    std::map<std::string, std::string> geopm;
    std::string geopars = _geopars;
    if(!geopars.empty())
    {
        boost::algorithm::trim( geopars );
        boost::char_separator<char> sep(":");
        boost::char_separator<char> sep2("=");
        boost::tokenizer< boost::char_separator<char> > kvlist(geopars, sep);
        for( const auto& ikvl : kvlist )
        {
            boost::tokenizer< boost::char_separator<char> > kv( ikvl, sep2);

            assert( distance( kv.begin(), kv.end() ) == 2 ); //! TODO modify message !

            try{
                //geopm[ *(kv.begin()) ] = boost::lexical_cast<double>( *(++(kv.begin())) );
                geopm[ *(kv.begin()) ] = *(++(kv.begin()));
            }
            catch( boost::bad_lexical_cast& e )
            {
                std::cerr<< "Error : " << e.what() << std::endl;
            }
        }
    }
    return geopm;
}

std::map<std::string, std::string>
Gmsh::retrieveGeoParameters( std::string const& __geo ) const
{
    std::map<std::string, std::string> __geopm;
    // (TODO should strip C/CPP comments)
    // Regex for a `keyword=value;` expression. We capture only [keyword]
    // and the [value] (ex : `h_2=-1,3e+4`).
    boost::regex kvreg("([[:word:]]+)[[:blank:]]*=[[:blank:]]*([+-]?(?:(?:(?:[[:digit:]]*\\.)?[[:digit:]]*(?:[eE][+-]?[[:digit:]]+)?)));" );
    boost::sregex_token_iterator iRex( __geo.begin(), __geo.end(), kvreg, 0 );
    boost::sregex_token_iterator end;

    for( ; iRex != end; ++iRex )
    {
        boost::smatch what;
        boost::regex_search( (*iRex).str(), what, kvreg );
        try
        {
            auto par = std::string( what[1].first, what[1].second );
            //auto val = boost::lexical_cast<double>( std::string( what[2].first, what[2].second ) );
            std::string val = std::string( what[2].first, what[2].second );


            __geopm[ par ] = val;
            LOG(INFO) << "[Gmsh::retrieveGeoParameter] New geometry parameter : "<< par << " = " << __geopm[par] << std::endl;
        }
        catch( boost::bad_lexical_cast& e )
        {
            std::cerr<< "Error : " << e.what() << std::endl;
        }
    }

    return __geopm;
}

bool
Gmsh::generateGeo( std::string const& __name, std::string const& __geo, bool const modifGeo ) const
{
    std::string _geo = __geo;

    if ( modifGeo )
    {
        // Create a new geo description from mapped parameters.
        for( const auto& iGpm : M_geoParamMap )
        {
            // Check any regular expression `mykey=myvalue;` in the description (see
            // retrieveGeoParameters()).
            boost::regex regex1( "(?:(" + iGpm.first  + "))[[:blank:]]*=[[:blank:]]*[+-]?(?:(?:(?:[[:digit:]]*\\.)?[[:digit:]]*(?:[eE][+-]?[[:digit:]]+)?));" );
            std::ostringstream _ostr;
            try{
                _ostr << "(?1$1) = " << iGpm.second << ";";
                LOG(INFO) << "[Gmsh::generateGeo] Geo geometry parameter "
                          << ( ( regex_search( __geo, regex1, boost::match_default) )?
                             ( iGpm.first + "=" + iGpm.second  + " now !" )
                             : iGpm.first + " not found ! " )
                          << std::endl;
            }
            catch( boost::bad_lexical_cast& e )
            {
                std::cerr<< "Error : " << e.what() << std::endl;
            }

            _geo = boost::regex_replace( _geo, regex1, _ostr.str(), boost::match_default | boost::format_all );
        }

        // Get the 'h' for hsize and modify its value in the geo file. (TODO could be included in the
        // geo-variables-list option (previous loop).
        // -----------
        boost::regex regex2( "(?:(lc|h))[[:blank:]]*=[[:blank:]]*[+-]?(?:(?:(?:[[:digit:]]*\\.)?[[:digit:]]*(?:[eE][+-]?[[:digit:]]+)?));" );
        std::ostringstream hstr;
        hstr << "(?1$1) = " << M_h << ";";

        LOG(INFO) << "[Gmsh::generateGeo] Parameter "
                  << ( ( regex_search( __geo, regex2, boost::match_default) )?
                     ( "hsize  = " + hstr.str() + " now ! (overwrite geo param h !)" )
                     : "h parameter (hsize) not found ! " )
                  << std::endl;

        _geo = boost::regex_replace( _geo, regex2, hstr.str(), boost::match_default | boost::format_all );
        // -----------
    }

    bool geochanged = true;
    // save geo file on disk if in-memory is false
    if ( M_in_memory == false )
    {
        std::ostringstream __geoname;
        __geoname << __name << ".geo";
        fs::path __path( __geoname.str() );
        geochanged = false;

        // Create a new .geo.
        if ( !fs::exists( __path ) )
        {
            LOG(INFO) << "[Gmsh::generateGeo] file :" << __path << "/" << __geoname.str() << std::endl;
            std::ofstream __geofile( __geoname.str().c_str() );
            __geofile << _geo;
            __geofile.close();
            geochanged = true;
        }
        // Overwrite .geo if it exists and differs.
        else
        {
            std::string s = this->getDescriptionFromFile( __geoname.str() );

            if ( s != _geo )
            {
                LOG(INFO) << __path << " exists but is different from the expected geometry, we overwrite it now";
                std::ofstream __geofile( __geoname.str().c_str() );
                __geofile << _geo;
                __geofile.close();
                geochanged = true;
            }
            else
            {
                LOG(INFO) << __path << " exists and its content correspond to the expected geometry";
            }
        }
    }
    M_geo = std::make_pair( __name, _geo );
    return geochanged;
}

std::string
Gmsh::generate( std::string const& name ) const
{
    std::string descr = this->getDescription();
    return this->generate( name, descr ).get<0>();
}

boost::tuple<std::string,bool>
Gmsh::generate( std::string const& __name, std::string const& __geo, bool const __forceRebuild, bool const parametric,bool const modifGeo ) const
{
    std::string fname;
    bool generated = false;
    if ( !mpi::environment::initialized() || ( mpi::environment::initialized()  && this->worldComm().globalRank() == this->worldComm().masterRank() ) )
    {


        LOG(INFO) << "Generate mesh on processor " <<  this->worldComm().globalRank() << "/" << this->worldComm().globalSize() << "\n";
        bool geochanged = generateGeo( __name,__geo,modifGeo );
        std::ostringstream __geoname;
        __geoname << __name << ".geo";


        // generate mesh
        std::ostringstream __meshname;
        __meshname << __name << ".msh";
        LOG( INFO ) << "Mesh file " << __meshname.str() << " generated from geometry " << __geoname.str();
        LOG( INFO ) << " - does mesh file name exists : " << (fs::exists( __meshname.str() )?"true":"false") << "\n";
        fs::path __meshpath( __meshname.str() );

        if ( geochanged || __forceRebuild || !fs::exists( __meshpath ) )
        {
            LOG( INFO ) << "Generating " << __meshname.str() << "...\n";

            generate( __geoname.str(), this->dimension(), parametric );

            LOG( INFO ) << "Generating " << __meshname.str() << " done.\n";
            generated = true;
        }
        fname=__meshname.str();

    }
    google::FlushLogFiles(INFO);

    auto ret = boost::make_tuple(fname,generated);
	if ( mpi::environment::initialized() )
    {
		this->worldComm().barrier();
        LOG(INFO) << "Broadcast mesh file " << fname << " to all other mpi processes, (created: " << generated << ")";
        mpi::broadcast( this->worldComm().globalComm(), ret, 0 );
    }
    return ret;
}
std::string
Gmsh::refine( std::string const& name, int level, bool parametric  ) const
{
#if FEELPP_HAS_GMSH
    std::ostringstream filename;
	std::string _name;

    if ( !mpi::environment::initialized() || ( mpi::environment::initialized()  && this->worldComm().globalRank() == this->worldComm().masterRank() ) )
    {
#if BOOST_FILESYSTEM_VERSION == 3
		filename << fs::path( name ).stem().string() << "-refine-" << level << ".msh";
		boost::system::error_code ec;
        auto na = fs::path( name );
        auto fi= fs::path( filename.str() );
		fs::copy_file( na, fi, fs::copy_option::overwrite_if_exists, ec );
#elif BOOST_FILESYSTEM_VERSION == 2
		filename << fs::path( name ).stem() << "-refine-" << level << ".msh";
		fs::copy_file( fs::path( name ), fs::path( filename.str() ), fs::copy_option::overwrite_if_exists );
#endif

#if !defined(FEELPP_HAS_GMSH_LIBRARY)
        // generate mesh
        std::ostringstream __str;

        if ( parametric )
            __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                  << " -parametric -refine " << filename.str();

        else
            __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                  << " -refine " << filename.str();

		for ( int l = 0; l < level; ++l )
		{
		    auto err = ::system( __str.str().c_str() );
		}

#else
		M_gmodel->readMSH(filename.str());

		CTX::instance()->mesh.order = M_order;
		CTX::instance()->mesh.secondOrderIncomplete = 0;
		CTX::instance()->mesh.secondOrderLinear = 1; // has to 1 to work


        int partitions = M_gmodel->getMeshPartitions().size();
		LOG(INFO) << "[Gmsh::refine] Original mesh : " << filename.str() << "\n";
		//LOG(INFO) << "[Gmsh::refine] vertices : " << M_gmodel->getNumMeshVertices() << "\n";
		LOG(INFO) << "[Gmsh::refine] elements : " << M_gmodel->getNumMeshElements() << "\n";
		LOG(INFO) << "[Gmsh::refine] partitions : " << partitions << "\n";
		//std::cout << "secondOrderLinear=" << CTX::instance()->mesh.secondOrderLinear << std::endl << std::flush;
		CTX::instance()->partitionOptions.num_partitions =  M_partitions;
		CTX::instance()->partitionOptions.partitioner =  M_partitioner;
		CTX::instance()->partitionOptions.mesh_dims[0] = M_partitions;
		CTX::instance()->mesh.mshFilePartitioned = M_partition_file;
		CTX::instance()->mesh.mshFileVersion = std::atof( this->version().c_str() );

		for ( int l = 0; l < level; ++l )
		{
		    M_gmodel->refineMesh( CTX::instance()->mesh.secondOrderLinear );
		}

		if ( partitions )
            {
                LOG(INFO) << "[Gmsh::refine] Repartioning mesh : " << filename.str() << "\n";
                PartitionMesh( M_gmodel, CTX::instance()->partitionOptions );
            }
        M_gmodel->writeMSH( filename.str() );
		LOG(INFO) << "[Gmsh::refine] Refined mesh : " << filename.str() << "\n";
		//LOG(INFO) << "[Gmsh::refine] vertices : " << M_gmodel->getNumMeshVertices() << "\n";
		LOG(INFO) << "[Gmsh::refine] elements : " << M_gmodel->getNumMeshElements() << "\n";
		LOG(INFO) << "[Gmsh::refine] partitions : " << M_gmodel->getMeshPartitions().size() << "\n";

		_name = filename.str();

		M_gmodel->destroy();
#endif

	}

	if ( mpi::environment::initialized() )
    {
		this->worldComm().barrier();

        _name = filename.str();
        mpi::broadcast( this->worldComm().globalComm(), _name, this->worldComm().masterRank() );
        LOG(INFO) << "[Gmsh::refine] broadcast mesh filename : " << _name << " to all other processes\n";
    }

    return _name;

#else
    throw std::invalid_argument( "Gmsh is not available on this system" );
#endif
}
void
Gmsh::generate( std::string const& __geoname, uint16_type dim, bool parametric  ) const
{
#if FEELPP_HAS_GMSH
#if !defined(FEELPP_HAS_GMSH_LIBRARY)
    // generate mesh
    std::ostringstream __str;

    //__str << "gmsh -algo tri -" << dim << " " << "-order " << this->order() << " " << __geoname;
    if ( parametric )
        __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
              << " -parametric -" << dim << " " << __geoname;

    else
        __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
              << " -" << dim << " -part " << M_partitions  << " " << __geoname;


    LOG(INFO) << "[Gmsh::generate] execute '" <<  __str.str() << "\n";
    LOG(INFO) << "[Gmsh::generate] partitions: " <<  M_partitions << "\n";
    LOG(INFO) << "[Gmsh::generate] partitioner: " <<  M_partitioner << "\n";

    auto err = ::system( __str.str().c_str() );
#else
    fs::path gp = __geoname;
    std::string _name = gp.stem().string();

    static bool gmshIsInit =false;

    if ( ! gmshIsInit )
    {
        gmshIsInit=true;
        GmshInitialize();
    }
    LOG(INFO) << "[Gmsh::generate] env.part: " <<  Environment::numberOfProcessors() << "\n";
    LOG(INFO) << "[Gmsh::generate] env.part: " <<  Environment::worldComm().size() << "\n";
    LOG(INFO) << "[Gmsh::generate] partitions: " <<  M_partitions << "\n";
    LOG(INFO) << "[Gmsh::generate] partitioner: " <<  M_partitioner << "\n";
    CTX::instance()->partitionOptions.num_partitions =  M_partitions;
    CTX::instance()->partitionOptions.partitioner =  M_partitioner;
    CTX::instance()->partitionOptions.mesh_dims[0] = M_partitions;

    CTX::instance()->mesh.mshFileVersion = std::atof( this->version().c_str() );
    CTX::instance()->mesh.lcExtendFromBoundary = 1;
    CTX::instance()->mesh.lcFromPoints = 1;
    CTX::instance()->mesh.order = M_order;
    CTX::instance()->mesh.secondOrderIncomplete = 0;

		if(doption("gmsh.randFactor") > 0.)
			CTX::instance()->mesh.randFactor = doption("gmsh.randFactor");
    //if ( M_recombine )
    if ( 0 )
    {
        CTX::instance()->mesh.algo2d = 5;
        CTX::instance()->mesh.algoRecombine = 1;
        CTX::instance()->mesh.recombineAll = 1;
    }

    else
    {
        CTX::instance()->mesh.algo2d = ( Environment::vm().count("gmsh.algo2d") )? ioption("gmsh.algo2d") : ALGO_2D_FRONTAL;
#if defined(HAVE_TETGEN)
        CTX::instance()->mesh.algo3d = ( Environment::vm().count("gmsh.algo3d") )? ioption("gmsh.algo3d") : ALGO_3D_DELAUNAY;
#else
        CTX::instance()->mesh.algo3d = ( Environment::vm().count("gmsh.algo3d") )? ioption("gmsh.algo3d") : ALGO_3D_FRONTAL;
#endif
    }
    // disable heap checking if enabled
    {
#if defined(FEELPP_HAS_GPERFTOOLS)
        HeapLeakChecker::Disabler disabler;
#endif /* FEELPP_HAS_GPERFTOOLS */


        CTX::instance()->mesh.mshFilePartitioned = M_partition_file;

        M_gmodel = new GModel;
        M_gmodel->setName( _name );
#if !defined( __APPLE__ )
        M_gmodel->setFileName( _name );
#endif
        if ( M_in_memory )
            ParseGeoFromMemory( M_gmodel, _name, geo().second );
        else
            M_gmodel->readGEO( __geoname );

        M_gmodel->mesh( dim );
        CHECK(M_gmodel->getMeshStatus() == dim)  << "Invalid Gmsh Mesh. Something went wrong with Gmsh.  Gmsh status : " << M_gmodel->getMeshStatus()
                                                 << " should be == " << dim;
        LOG(INFO) << "Mesh refinement levels : " << M_refine_levels << "\n";
        for( int l = 0; l < M_refine_levels; ++l )
        {
            M_gmodel->refineMesh( CTX::instance()->mesh.secondOrderLinear );
        }
        PartitionMesh( M_gmodel, CTX::instance()->partitionOptions );
        LOG(INFO) << "Mesh partitions : " << M_gmodel->getMeshPartitions().size() << "\n";

        CHECK(M_gmodel->getMeshStatus() == dim)  << "Invalid Gmsh Mesh. Something went wrong with Gmsh.  Gmsh status : " << M_gmodel->getMeshStatus()
                                                 << " should be == " << dim;
        if ( M_in_memory == false )
        {
            CTX::instance()->mesh.binary = M_format;
            LOG(INFO) << "Writing GMSH file " << _name+".msh" << " in " << (M_format?"binary":"ascii") << " format\n";
            M_gmodel->writeMSH( _name+".msh", 2.2, CTX::instance()->mesh.binary );
        }

    }
#endif
#else
    throw std::invalid_argument( "Gmsh is not available on this system" );
#endif
}

void
Gmsh::rebuildPartitionMsh( std::string const& nameMshInput,std::string const& nameMshOutput ) const
{
#if defined(FEELPP_HAS_GMSH_LIBRARY)

    std::string _name;

    if ( !mpi::environment::initialized() || ( mpi::environment::initialized()  && this->worldComm().globalRank() == this->worldComm().masterRank() ) )
    {
#if BOOST_FILESYSTEM_VERSION == 3
		std::string _name = fs::path( nameMshInput ).stem().string();
#elif BOOST_FILESYSTEM_VERSION == 2
        std::string _name = fs::path( nameMshInput ).stem();
#endif

        GModel* M_gmodel=new GModel();
        M_gmodel->readMSH( nameMshInput );

        meshPartitionOptions newPartionOption;
        newPartionOption.num_partitions = M_partitions;
        newPartionOption.mesh_dims[0] = M_partitions;
        if (M_partitions==1)
            newPartionOption.partitioner=GMSH_PARTITIONER_METIS;
        else
            newPartionOption.partitioner =  M_partitioner;

        CTX::instance()->mesh.mshFilePartitioned = M_partition_file;
        CTX::instance()->mesh.mshFileVersion = std::atof( this->version().c_str() );
        PartitionMesh( M_gmodel, newPartionOption );

        CTX::instance()->mesh.binary = M_format;
        //M_gmodel.mesh.binary = M_format;

        //M_gmodel->writeMSH( nameMshOutput );
        M_gmodel->writeMSH( nameMshOutput, 2.2, CTX::instance()->mesh.binary );

        M_gmodel->destroy();
    }

    if ( mpi::environment::initialized() )
    {
        this->worldComm().barrier();
        mpi::broadcast( this->worldComm().globalComm(), _name, this->worldComm().masterRank() );
        LOG(INFO) << "[Gmsh::rebuildPartitionMsh] broadcast mesh filename : " << _name << " to all other processes\n";
    }


#else
    throw std::invalid_argument( "Gmsh library is not available on this system" );
#endif

}

/* if Gmsh API is not detected, some variables need to be define
   see gmsh/Common/GmshDefines.h
*/

#ifndef FEELPP_HAS_GMSH_H

// 2D meshing algorithms (numbers should not be changed)
#define ALGO_2D_MESHADAPT      1
#define ALGO_2D_AUTO           2
#define ALGO_2D_MESHADAPT_OLD  4
#define ALGO_2D_DELAUNAY       5
#define ALGO_2D_FRONTAL        6
#define ALGO_2D_BAMG           7
#define ALGO_2D_FRONTAL_QUAD   8

// 3D meshing algorithms (numbers should not be changed)
#define ALGO_3D_DELAUNAY       1
#define ALGO_3D_FRONTAL        4
#define ALGO_3D_FRONTAL_DEL    5
#define ALGO_3D_FRONTAL_HEX    6
#define ALGO_3D_MMG3D          7
#define ALGO_3D_RTREE          9

#endif

    std::string
    Gmsh::preamble() const
    {
        std::ostringstream ostr;

        ostr << "Mesh.MshFileVersion = " << this->version() << ";\n"
             << "Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
             << "Mesh.CharacteristicLengthFromPoints=1;\n"
             << "Mesh.ElementOrder=" << M_order << ";\n"
             << "Mesh.SecondOrderIncomplete = 0;\n";

        if ( M_recombine )
            ostr << "Mesh.Algorithm = 5;\n";
        else
        {
            ostr << "Mesh.Algorithm = " << ALGO_2D_FRONTAL << ";\n";
#if defined(HAVE_TETGEN)
            ostr << "Mesh.Algorithm3D = " << ALGO_3D_DELAUNAY << ";\n";
#else
            ostr << "Mesh.Algorithm3D = " << ALGO_3D_FRONTAL << ";\n";
#endif
        }

        ostr << "//Mesh.OptimizeNetgen=1;\n";

        if ( this->worldComm().globalSize() != 1 )
        {
            ostr << "// partitioning data\n"
                 << "Mesh.Partitioner=" << M_partitioner << ";\n"
                 << "Mesh.NbPartitions=" << M_partitions << ";\n"
                 << "Mesh.MshFilePartitioned=" << M_partition_file << ";\n";
        }
        //ostr << "Mesh.Optimize=1;\n"
        //<< "Mesh.CharacteristicLengthFromCurvature=1;\n"

        // if (this->structuredMesh() == 3)
        // {
        //     ostr << "nx=" << M_nx << ";\n"
        //          << "ny=" << M_ny << ";\n";
        // }
        // else
        // {
        ostr << "h=" << M_h << ";\n";
        //}

        if ( M_recombine )
        {
            ostr << "Mesh.RecombinationAlgorithm=1;//blossom\n"
                 << "Mesh.RecombineAll=1; //all\n";
        }
        else
        {
            ostr << "Mesh.RecombinationAlgorithm=0;\n";
        }

        return ostr.str();
    }

/// \cond detail
namespace detail
{

struct HypercubeDomain
{
    HypercubeDomain( int _Dim, int _Order )
        :
        Dim( _Dim ), Order( _Order ), RDim( _Dim ), Hyp( false )
        {}
    HypercubeDomain( int _Dim, int _Order, int _RDim, std::string const& hyp )
        :
        Dim( _Dim ), Order( _Order ), RDim( _RDim ), Hyp( hyp == "hypercube" )
        {
        }
    Gmsh* operator()()
        {
            return new GmshHypercubeDomain( Dim, Order, RDim, Hyp );
        }
    int Dim, Order,RDim;
    bool Hyp;
};

struct SimplexDomain
{
    SimplexDomain( int _Dim, int _Order )
        :
        Dim( _Dim ), Order( _Order )
        {}
    Gmsh* operator()()
        {
            return new GmshSimplexDomain( Dim, Order );
        }
    int Dim, Order;
};

struct EllipsoidDomain
{
    EllipsoidDomain( int _Dim, int _Order )
        :
        Dim( _Dim ), Order( _Order )
        {}
    Gmsh* operator()()
        {
            return new GmshEllipsoidDomain( Dim, Order );
        }
    int Dim, Order;
};

} // detail



# define DIMS BOOST_PP_TUPLE_TO_LIST(3,(1,2,3))
# define ORDERS BOOST_PP_TUPLE_TO_LIST(5,(1,2,3,4,5))
# define SHAPES1 BOOST_PP_TUPLE_TO_LIST(3, ((2,(simplex, Simplex)) ,    \
                                            (2,(ellipsoid, Ellipsoid)) , \
                                            (2,(hypercube, Hypercube)) ))
# define SHAPES2 BOOST_PP_TUPLE_TO_LIST(2, ((3,(hypercube, Hypercube, Simplex)), \
                                            (3,(hypercube, Hypercube, Hypercube)) ) )


#define FACTORY1NAME( LDIM, LORDER, LSHAPE )                            \
    BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0,LSHAPE) BOOST_PP_LPAREN() LDIM BOOST_PP_COMMA() LORDER BOOST_PP_RPAREN())

# define FACTORY1(LDIM,LORDER,LSHAPE )                                  \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( mesh, LDIM ), LORDER), BOOST_PP_ARRAY_ELEM(1,LSHAPE))  = \
                Gmsh::Factory::type::instance().registerProduct( boost::to_lower_copy(boost::algorithm::erase_all_copy( std::string( FACTORY1NAME(LDIM, LORDER, LSHAPE ) ), " " ) ), \
                                                                 *new Feel::detail::BOOST_PP_CAT(BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER) );

# define FACTORY1_OP(_, GDO) FACTORY1 GDO


#define FACTORY2NAME( LDIM, LORDER, LSHAPE )                            \
    BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0,LSHAPE) BOOST_PP_LPAREN() LDIM BOOST_PP_COMMA() LORDER BOOST_PP_COMMA() BOOST_PP_ARRAY_ELEM(2,LSHAPE) BOOST_PP_RPAREN())

# define FACTORY2(LDIM,LORDER,LSHAPE )                                  \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( mesh, LDIM ), LORDER), BOOST_PP_ARRAY_ELEM(1,LSHAPE)), BOOST_PP_ARRAY_ELEM(2,LSHAPE))   = \
                Gmsh::Factory::type::instance().registerProduct( boost::to_lower_copy( boost::algorithm::erase_all_copy( std::string( FACTORY2NAME(LDIM, LORDER, LSHAPE ) ), " " ) ), \
                                                                 *new Feel::detail::BOOST_PP_CAT(BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER,LDIM,boost::to_lower_copy(std::string(BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(2,LSHAPE)))) ));

# define FACTORY2_OP(_, GDO) FACTORY2 GDO

// only up to 4 for mesh data structure not supported for higher order in Gmsh
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY1_OP, 3, ( DIMS, ORDERS, SHAPES1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY2_OP, 3, ( DIMS, ORDERS, SHAPES2 ) )

const bool meshs213s = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,3,simplex)", *new Feel::detail::HypercubeDomain( 2, 1, 3, "simplex" ) );
const bool meshs213ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,3,hypercube)", *new Feel::detail::HypercubeDomain( 2, 1, 3, "hypercube" ) );
const bool meshs112s = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1,2,simplex)", *new Feel::detail::HypercubeDomain( 1, 1, 2, "simplex" ) );
const bool meshs112ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1,2,yypercube)", *new Feel::detail::HypercubeDomain( 1, 1, 2, "hypercube" ) );

/// \endcond detail







} // Feel
