/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2005-02-10

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-02-10
 */
#include <cstdlib>

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

#if defined( FEELPP_HAS_GMSH_H )
#include <GmshConfig.h>
#include <Gmsh.h>
#include <GModel.h>
#include <GmshDefines.h>
#include <Context.h>
//#include <meshPartition.h>
int PartitionMesh( GModel *const model, meshPartitionOptions &options );
#endif

namespace Feel
{
namespace fs = boost::filesystem;

const char* FEELPP_GMSH_FORMAT_VERSION = "2.2";

Gmsh::Gmsh( int nDim, int nOrder, WorldComm const& worldComm )
    :
    M_worldComm( worldComm ),
    M_dimension( nDim ),
    M_order( nOrder ),
    M_version( FEELPP_GMSH_FORMAT_VERSION ),
    M_I( nDim ),
    M_h( 0.1 ),
    M_addmidpoint( true ),
    M_usePhysicalNames( false ),
    M_partitioner( GMSH_PARTITIONER_CHACO ),
    M_partitions( 1 ),
    M_partition_file( 0 ),
    M_shear( 0 ),
    M_refine_levels( 0 )
{
    this->setReferenceDomain();
}
Gmsh::Gmsh( Gmsh const & __g )
    :
    M_worldComm( __g.M_worldComm ),
    M_dimension( __g.M_dimension ),
    M_order( __g.M_order ),
    M_version( __g.M_version ),
    M_I( __g.M_I ),
    M_h( __g.M_h ),
    M_addmidpoint( __g.M_addmidpoint ),
    M_usePhysicalNames( __g.M_usePhysicalNames ),
    M_partitioner( __g.M_partitioner ),
    M_partitions( __g.M_partitions ),
    M_partition_file( __g.M_partition_file ),
    M_shear( __g.M_shear ),
    M_refine_levels( __g.M_refine_levels )
{}
Gmsh::~Gmsh()
{}

boost::shared_ptr<Gmsh>
Gmsh::New( po::variables_map const& vm )
{
    std::ostringstream ostr;
    ostr << vm["convex"].as<std::string>() << "(" << vm["dim"].as<int>() << "," << vm["order"].as<int>() << ")";
    boost::shared_ptr<Gmsh> gmsh_ptr( Gmsh::Factory::type::instance().createObject( ostr.str() ) );
    return gmsh_ptr;
}

boost::shared_ptr<Gmsh>
Gmsh::New( std::string const& shape, uint16_type d, uint16_type o, std::string const& ct )
{
    std::ostringstream ostr;

    if ( shape != "hypercube" && ct != "hypercube" )
        ostr << shape << "(" << d << "," << o << ")";

    else
        ostr << shape << "(" << d << "," << o << "," << boost::to_lower_copy( ct ) << ")";

    boost::shared_ptr<Gmsh> gmsh_ptr( Gmsh::Factory::type::instance().createObject( ostr.str() ) );
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
bool
Gmsh::generateGeo( std::string const& __name, std::string const& __geo,bool const modifGeo ) const
{
    std::string _geo;

    if ( modifGeo )
    {
        boost::regex regex( "(?:(lc|h))[[:blank:]]*=[[:blank:]]*[+-]?(?:(?:(?:[[:digit:]]*\\.)?[[:digit:]]*(?:[eE][+-]?[[:digit:]]+)?));" );
        std::ostringstream hstr;
        hstr << "(?1$1) = " << M_h << ";";
        Debug( 10000 ) << "found hsize: " << regex_search(__geo, regex, boost::match_default) << "\n";
        Debug( 10000 ) << "hstr: " << hstr.str() << "\n";

        _geo = boost::regex_replace( __geo, regex, hstr.str(), boost::match_default | boost::format_all );
    }

    else
    {
        _geo = __geo;
    }

    // generate geo
    std::ostringstream __geoname;
    __geoname << __name << ".geo";
    fs::path __path( __geoname.str() );
    bool geochanged = false;

    if ( !fs::exists( __path ) )
    {
        Debug( 10000 ) << "generating: " << __geoname.str() << "\n";
        std::ofstream __geofile( __geoname.str().c_str() );
        __geofile << _geo;
        __geofile.close();
        geochanged = true;
    }

    else
    {
        std::string s = this->getDescriptionFromFile( __geoname.str() );

        if ( s != _geo )
        {
            std::ofstream __geofile( __geoname.str().c_str() );
            __geofile << _geo;
            __geofile.close();
            geochanged = true;
        }
    }

    return geochanged;
}
std::string
Gmsh::generate( std::string const& name ) const
{
    std::string descr = this->getDescription();
    return this->generate( name, descr );
}

std::string
Gmsh::generate( std::string const& __name, std::string const& __geo, bool const __forceRebuild, bool const parametric,bool const modifGeo ) const
{
    std::string fname;

    if ( !mpi::environment::initialized() || ( mpi::environment::initialized()  && this->worldComm().globalRank() == this->worldComm().masterRank() ) )
    {
        Log() << "[Gmsh::generate] generate on processor " <<  this->worldComm().globalRank() << "/" << this->worldComm().globalSize() << "\n";
        bool geochanged ( generateGeo( __name,__geo,modifGeo ) );
        std::ostringstream __geoname;
        __geoname << __name << ".geo";

        // generate mesh
        std::ostringstream __meshname;
        __meshname << __name << ".msh";
        Debug( 10000 ) << "mesh file name: " << __meshname.str() << "\n";
        Debug( 10000 ) << "does mesh file name exists ?: " << fs::exists( __meshname.str() ) << "\n";
        fs::path __meshpath( __meshname.str() );

        if ( geochanged || __forceRebuild || !fs::exists( __meshpath ) )
        {
            Debug( 10000 ) << "generating: " << __meshname.str() << "\n";
#if 0

            if ( __geo.find( "Volume" ) != std::string::npos )
                generate( __geoname.str(), 3, parametric );

            else if ( __geo.find( "Surface" ) != std::string::npos )
                generate( __geoname.str(), 2, parametric );

            else //if ( __geo.find( "Line" )  != std::string::npos )
                generate( __geoname.str(), 1, parametric );

            //else
            //generate( __geoname.str(), 3, parametric );
#else
            generate( __geoname.str(), this->dimension(), parametric );
#endif
        }

        Log() << "[Gmsh::generate] meshname = " << __meshname.str() << "\n";
        fname=__meshname.str();
    }

    if ( mpi::environment::initialized() )
    {
        mpi::broadcast( this->worldComm().globalComm(), fname, 0 );
        Log() << "[Gmsh::generate] broadcast mesh filename : " << fname << " to all other processes\n";

    }

    return fname;
}
std::string
Gmsh::refine( std::string const& name, int level, bool parametric  ) const
{
#if FEELPP_HAS_GMSH
    std::ostringstream filename;
    filename << fs::path( name ).stem() << "-refine-" << level << ".msh";

#if BOOST_FILESYSTEM_VERSION == 3
    boost::system::error_code ec;
    fs::copy_file( fs::path( name ), fs::path( filename.str() ), fs::copy_option::overwrite_if_exists, ec );
#elif BOOST_FILESYSTEM_VERSION == 2
    fs::copy_file( fs::path( name ), fs::path( filename.str() ), fs::copy_option::overwrite_if_exists );
#endif

    for ( int l = 0; l < level; ++l )
    {
        // generate mesh
        std::ostringstream __str;

        if ( parametric )
            __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                  << " -parametric -refine " << filename.str();

        else
            __str << BOOST_PP_STRINGIZE( GMSH_EXECUTABLE )
                  << " -refine " << filename.str();

        ::system( __str.str().c_str() );
    }

    return filename.str();
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

    Log() << "[Gmsh::generate] execute '" <<  __str.str() << "\n";

    auto err = ::system( __str.str().c_str() );
#else
    std::string _name = fs::path( __geoname ).stem().string();

    static bool gmshIsInit =false;

    if ( ! gmshIsInit )
    {
        gmshIsInit=true;
        GmshInitialize();
    }

    CTX::instance()->partitionOptions.num_partitions =  M_partitions;
    CTX::instance()->partitionOptions.partitioner =  M_partitioner;


    CTX::instance()->mesh.mshFileVersion = std::atof( this->version().c_str() );
    CTX::instance()->mesh.lcExtendFromBoundary = 1;
    CTX::instance()->mesh.lcFromPoints = 1;
    CTX::instance()->mesh.order = M_order;
    CTX::instance()->mesh.secondOrderIncomplete = 0;

    //if ( M_recombine )
    if ( 0 )
    {
        CTX::instance()->mesh.algo2d = 5;
        CTX::instance()->mesh.algoRecombine = 1;
        CTX::instance()->mesh.recombineAll = 1;
    }

    else
    {
        CTX::instance()->mesh.algo2d = ALGO_2D_FRONTAL;
#if defined(HAVE_TETGEN)
        CTX::instance()->mesh.algo3d = ALGO_3D_DELAUNAY;
#else
        CTX::instance()->mesh.algo3d = ALGO_3D_FRONTAL;
#endif
    }

    CTX::instance()->mesh.mshFilePartitioned = M_partition_file;

    new GModel();
    GModel::current()->setName( _name );
    GModel::current()->setFileName( _name );
    GModel::current()->readGEO( _name+".geo" );
    GModel::current()->mesh( dim );
    for( int l = 0; l < M_refine_levels-1; ++l )
        GModel::current()->refineMesh( M_order==1 );
    PartitionMesh( GModel::current(), CTX::instance()->partitionOptions );
    //std::cout << "size : " << GModel::current()->getMeshPartitions().size() << "\n";
    GModel::current()->writeMSH( _name+".msh" );
    //GModel::current()->destroy();
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

        std::string _name = fs::path( nameMshInput ).stem().string();

        GModel* newGmshModel=new GModel();
        newGmshModel->readMSH( nameMshInput );

        meshPartitionOptions newPartionOption;
        newPartionOption.num_partitions = M_partitions;
        newPartionOption.mesh_dims[0] = M_partitions;
        if (M_partitions==1)
            newPartionOption.partitioner=GMSH_PARTITIONER_METIS;
        else
            newPartionOption.partitioner =  M_partitioner;

        CTX::instance()->mesh.mshFilePartitioned = M_partition_file;
        CTX::instance()->mesh.mshFileVersion = std::atof( this->version().c_str() );
        PartitionMesh( newGmshModel, newPartionOption );

        newGmshModel->writeMSH( nameMshOutput );

        newGmshModel->destroy();
        delete newGmshModel;

    }

    if ( mpi::environment::initialized() )
    {
        mpi::broadcast( this->worldComm().globalComm(), _name, this->worldComm().masterRank() );
        Log() << "[Gmsh::rebuildPartitionMsh] broadcast mesh filename : " << _name << " to all other processes\n";
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

    ostr << "//Mesh.OptimizeNetgen=1;\n"
         << "// partitioning data\n"
         << "Mesh.Partitioner=" << M_partitioner << ";\n"
         << "Mesh.NbPartitions=" << M_partitions << ";\n"
         << "Mesh.MshFilePartitioned=" << M_partition_file << ";\n"
        //<< "Mesh.Optimize=1;\n"
        //<< "Mesh.CharacteristicLengthFromCurvature=1;\n"
         << "h=" << M_h << ";\n";

    if ( M_recombine )
    {
        ostr << "Mesh.RecombinationAlgorithm=1;//blossom\n"
             << "Mesh.RecombineAll=1; //all\n";
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
                                                                 *new detail::BOOST_PP_CAT(BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER) );

# define FACTORY1_OP(_, GDO) FACTORY1 GDO


#define FACTORY2NAME( LDIM, LORDER, LSHAPE )                            \
    BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(0,LSHAPE) BOOST_PP_LPAREN() LDIM BOOST_PP_COMMA() LORDER BOOST_PP_COMMA() BOOST_PP_ARRAY_ELEM(2,LSHAPE) BOOST_PP_RPAREN())

# define FACTORY2(LDIM,LORDER,LSHAPE )                                  \
    const bool BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( BOOST_PP_CAT( mesh, LDIM ), LORDER), BOOST_PP_ARRAY_ELEM(1,LSHAPE)), BOOST_PP_ARRAY_ELEM(2,LSHAPE))   = \
                Gmsh::Factory::type::instance().registerProduct( boost::to_lower_copy( boost::algorithm::erase_all_copy( std::string( FACTORY2NAME(LDIM, LORDER, LSHAPE ) ), " " ) ), \
                                                                 *new detail::BOOST_PP_CAT(BOOST_PP_ARRAY_ELEM(1,LSHAPE),Domain)(LDIM,LORDER,LDIM,boost::to_lower_copy(std::string(BOOST_PP_STRINGIZE(BOOST_PP_ARRAY_ELEM(2,LSHAPE)))) ));

# define FACTORY2_OP(_, GDO) FACTORY2 GDO

// only up to 4 for mesh data structure not supported for higher order in Gmsh
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY1_OP, 3, ( DIMS, ORDERS, SHAPES1 ) )
BOOST_PP_LIST_FOR_EACH_PRODUCT( FACTORY2_OP, 3, ( DIMS, ORDERS, SHAPES2 ) )

const bool meshs213s = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,3,simplex)", *new detail::HypercubeDomain( 2, 1, 3, "simplex" ) );
const bool meshs213ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,3,hypercube)", *new detail::HypercubeDomain( 2, 1, 3, "hypercube" ) );
const bool meshs112s = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1,2,simplex)", *new detail::HypercubeDomain( 1, 1, 2, "simplex" ) );
const bool meshs112ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1,2,yypercube)", *new detail::HypercubeDomain( 1, 1, 2, "yypercube" ) );

/// \endcond detail
} // Feel
