/* -*- mode: c++ -*-

  This file is part of the Life library

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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/concept_check.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifemesh/simplexproduct.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/gmshsimplexdomain.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>

namespace Life
{
namespace fs = boost::filesystem;

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
    ostr << shape << "(" << d << "," << o << ")";
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
Gmsh::generateCircle( std::string const& __name, double __h )
{
    uint16_type N=( uint16_type )(1./__h);

    Debug(10000) << "prefix file name: " << __name << "\n";
    Debug(10000) << "h: " << __h << "\n";
    Debug(10000) << "N: " << N << "\n";

    // generate geo
    std::ostringstream __geoname;
    __geoname << this->prefix( __name, N ) << ".geo";
    fs::path __path( __geoname.str() );
    if ( !fs::exists( __path ) )
        {
            Debug(10000) << "generating: " << __geoname.str() << "\n";
            std::ofstream __geo( __geoname.str().c_str() );
            __geo << "Mesh.MshFileVersion = 2;\n"
                  << "h=" << __h << ";\n"
                  << "Point(1) = {0.0,0.0,0.0,h};\n"
                  << "Point(2) = {1,0.0,0.0,h};\n"
                  << "Point(3) = {0.2,0.0,0.0,h};\n"
                  << "Point(4) = {0,1,0,h};\n"
                  << "Point(5) = {0,0.2,0,h};\n"
                  << "Point(6) = {0,-0.2,0,h};\n"
                  << "Point(7) = {0,-1,0,h};\n"
                  << "Point(8) = {-0.2,0,0,h};\n"
                  << "Point(9) = {-1,0,0,h};\n"
                  << "Circle(1) = {2,1,4};\n"
                  << "Circle(2) = {4,1,9};\n"
                  << "Circle(3) = {9,1,7};\n"
                  << "Circle(4) = {7,1,2};\n"
                  << "Circle(5) = {3,1,5};\n"
                  << "Circle(6) = {5,1,8};\n"
                  << "Circle(7) = {8,1,6};\n"
                  << "Circle(8) = {6,1,3};\n"
                  << "Line Loop(9) = {1,2,3,4};\n"
                  << "Line Loop(10) = {6,7,8,5};\n"
                  << "Plane Surface(11) = {9,10};\n";

            __geo.close();
        }
    // generate mesh
    std::ostringstream __meshname;
    __meshname << this->prefix( __name, N ) << ".msh";
    Debug(10000) << "mesh file name: " << __meshname.str() << "\n";
    Debug(10000) << "does mesh file name exists ?: " << fs::exists(__meshname.str() ) << "\n";
    fs::path __meshpath( __meshname.str() );
    if ( !fs::exists( __meshpath ) )
        {
            Debug(10000) << "generating: " << __meshname.str() << "\n";
            generate( __geoname.str(), 2 );
        }
    return __meshname.str();
}
std::string
Gmsh::generateLine( std::string const& __name, double __h )
{
    uint16_type N=( uint16_type )(1./__h);

    Debug(10000) << "prefix file name: " << __name << "\n";
    Debug(10000) << "h: " << __h << "\n";
    Debug(10000) << "N: " << N << "\n";

    // generate geo
    std::ostringstream __geoname;
    __geoname << this->prefix( __name, N ) << ".geo";
    fs::path __path( __geoname.str() );
    if ( !fs::exists( __path ) )
        {
            Debug(10000) << "generating: " << __geoname.str() << "\n";
            std::ofstream __geo( __geoname.str().c_str() );
            __geo << "Mesh.MshFileVersion = 2;\n"
                  << "h=" << __h << ";\n"
                  << "Point(1) = {0.0,0.0,0.0,h};\n"
                  << "Point(2) = {1.0,0,0.0,h};\n"
                  << "Line(1) = {1,2};\n";
            __geo.close();
        }
    // generate mesh
    std::ostringstream __meshname;
    __meshname << this->prefix( __name, N ) << ".msh";
    Debug(10000) << "mesh file name: " << __meshname.str() << "\n";
    Debug(10000) << "does mesh file name exists ?: " << fs::exists(__meshname.str() ) << "\n";
    fs::path __meshpath( __meshname.str() );
    if ( !fs::exists( __meshpath ) )
        {
            Debug(10000) << "generating: " << __meshname.str() << "\n";
            generate( __geoname.str(), 1 );
        }
    return __meshname.str();
}
std::string
Gmsh::generateSquare( std::string const& __name, double __h )
{
    uint16_type N=( uint16_type )(1./__h);

    Debug(10000) << "prefix file name: " << __name << "\n";
    Debug(10000) << "h: " << __h << "\n";
    Debug(10000) << "N: " << N << "\n";

    // generate geo
    std::ostringstream __geoname;
    __geoname << this->prefix( __name, N ) << ".geo";
    fs::path __path( __geoname.str() );
    if ( !fs::exists( __path ) )
        {
            Debug(10000) << "generating: " << __geoname.str() << "\n";
            std::ofstream __geo( __geoname.str().c_str() );
            __geo << "Mesh.MshFileVersion = 2;\n"
                  << "h=" << __h << ";\n"
                  << "Point(1) = {0.0,0.0,0.0,h};\n"
                  << "Point(2) = {0.0,1,0.0,h};\n"
                  << "Point(3) = {1,1,0.0,h};\n"
                  << "Point(4) = {1,0,0.0,h};\n"
                  << "Line(1) = {1,4};\n"
                  << "Line(2) = {4,3};\n"
                  << "Line(3) = {3,2};\n"
                  << "Line(4) = {2,1};\n"
                  << "Line Loop(5) = {3,4,1,2};\n"
                  << "Plane Surface(6) = {5};\n"
                  << "Physical Line(10) = {1,2,4};\n"
                  << "Physical Line(20) = {3};\n"
                  << "Physical Surface(30) = {6};\n";
            __geo.close();
        }
    // generate mesh
    std::ostringstream __meshname;
    __meshname << this->prefix( __name, N ) << ".msh";
    Debug(10000) << "mesh file name: " << __meshname.str() << "\n";
    Debug(10000) << "does mesh file name exists ?: " << fs::exists(__meshname.str() ) << "\n";
    fs::path __meshpath( __meshname.str() );
    if ( !fs::exists( __meshpath ) )
        {
            Debug(10000) << "generating: " << __meshname.str() << "\n";
            generate( __geoname.str(), 2 );
        }
    return __meshname.str();
}

std::string
Gmsh::generateCube( std::string const& __name, double __h )
{
    uint16_type N=( uint16_type )(1./__h);

    Debug(10000) << "prefix file name: " << __name << "\n";
    Debug(10000) << "h: " << __h << "\n";
    Debug(10000) << "N: " << N << "\n";

    // generate geo
    std::ostringstream __geoname;
    __geoname << this->prefix( __name, N ) << ".geo";
    fs::path __path( __geoname.str() );
    if ( !fs::exists( __path ) )
        {
            Debug(10000) << "generating: " << __geoname.str() << "\n";
            std::ofstream __geo( __geoname.str().c_str() );
            __geo << "Mesh.MshFileVersion = 2;\n"
                  << "h=" << __h << ";\n"
                  << "Point(1) = {0.0,0.0,0.0,h};\n"
                  << "Point(2) = {0.0,1,0.0,h};\n"
                  << "Point(3) = {1,1,0.0,h};\n"
                  << "Point(4) = {1,0,0.0,h};\n"
                  << "Line(1) = {1,4};\n"
                  << "Line(2) = {4,3};\n"
                  << "Line(3) = {3,2};\n"
                  << "Line(4) = {2,1};\n"
                  << "Line Loop(5) = {3,4,1,2};\n"
                  << "Plane Surface(6) = {5};\n"
                  << "Extrude Surface {6, {0.0,0.0,1.0}};\n"
                  << "Physical Surface(10) = {19,27,15,23,6};\n"
                  << "Physical Surface(20) = {28};\n"
                  << "Surface Loop(31) = {28,15,-6,19,23,27};\n"
                  << "\n"
                  << "// volume \n"
                  << "Volume(1) = {31};\n"
                  << "Physical Volume(2) = {1};\n";
            __geo.close();
        }
    // generate mesh
    std::ostringstream __meshname;
    __meshname << this->prefix( __name, N ) << ".msh";
    Debug(10000) << "mesh file name: " << __meshname.str() << "\n";
    Debug(10000) << "does mesh file name exists ?: " << fs::exists(__meshname.str() ) << "\n";
    fs::path __meshpath( __meshname.str() );
    if ( !fs::exists( __meshpath ) )
        {
            Debug(10000) << "generating: " << __meshname.str() << "\n";
            generate( __geoname.str(), 3 );
        }
    return __meshname.str();
}

bool
Gmsh::generateGeo( std::string const& __name, std::string const& __geo_ ) const
{
    std::ostringstream ofsgeo;
    ofsgeo << "Mesh.ElementOrder=" << _M_order << ";\n"
        //<< "Mesh.SecondOrderExperimental = 1;\n"//Ne semble plus indispensable
           << "Mesh.SecondOrderIncomplete = 0;\n"
           << "Mesh.Algorithm = 6;\n"
           << __geo_;
    std::string __geo = ofsgeo.str();

    // generate geo
    std::ostringstream __geoname;
    __geoname << __name << ".geo";
    fs::path __path( __geoname.str() );
    bool geochanged = false;
    if ( !fs::exists( __path ) )
        {
            Debug(10000) << "generating: " << __geoname.str() << "\n";
            std::ofstream __geofile( __geoname.str().c_str() );
            __geofile << __geo.c_str();
            __geofile.close();
            geochanged = true;
        }
    else
        {
            // has the file changed ?
            std::ifstream __geoin( __geoname.str().c_str() );

            std::ostringstream __geostream;
            std::istreambuf_iterator<char> src(__geoin.rdbuf());
            std::istreambuf_iterator<char> end;
            std::ostream_iterator<char> dest(__geostream);

            std::copy(src,end,dest);

            std::string s = __geostream.str();
            __geoin.close();
            if ( s != __geo )
                {
                    std::ofstream __geofile( __geoname.str().c_str() );
                    __geofile << __geo.c_str();
                    __geofile.close();
                    geochanged = true;
                }
        }

    return geochanged;
}

std::string
Gmsh::generate( std::string const& __name, std::string const& __geo, bool const __forceRebuild ) const
{
    std::string fname;
    if ( !mpi::environment::initialized() || (mpi::environment::initialized()  && M_comm.rank() == 0 ) )
        {
            bool geochanged (generateGeo(__name,__geo));
            std::ostringstream __geoname;
            __geoname << __name << ".geo";

            // generate mesh
            std::ostringstream __meshname;
            __meshname << __name << ".msh";
            Debug(10000) << "mesh file name: " << __meshname.str() << "\n";
            Debug(10000) << "does mesh file name exists ?: " << fs::exists(__meshname.str() ) << "\n";
            fs::path __meshpath( __meshname.str() );
            if ( geochanged || __forceRebuild || !fs::exists( __meshpath ) )
                {
                    Debug(10000) << "generating: " << __meshname.str() << "\n";
                    if ( __geo.find( "Volume" ) != std::string::npos )
                        generate( __geoname.str(), 3 );
                    else if ( __geo.find( "Surface" ) != std::string::npos )
                        generate( __geoname.str(), 2 );
                    else if ( __geo.find( "Line" )  != std::string::npos )
                        generate( __geoname.str(), 1 );
                    else
                        generate( __geoname.str(), 3 );
                }
            Debug( 10000 ) << "[Gmsh::generate] meshname = " << __meshname.str() << "\n";
            fname=__meshname.str();
        }
    if ( mpi::environment::initialized() )
        mpi::broadcast( M_comm, fname, 0 );
    return fname;
}
void
Gmsh::generate( std::string const& __geoname, uint16_type dim  ) const
{
#if HAVE_GMSH
    // generate mesh
    std::ostringstream __str;
    //__str << "gmsh -algo tri -" << dim << " " << "-order " << this->order() << " " << __geoname;
    __str << "gmsh -" << dim << " " << __geoname;
    ::system( __str.str().c_str() );
#else
    throw std::invalid_argument("Gmsh is not available on this system");
#endif
}

namespace detail {
template<int Dim, int Order>
Gmsh* createSimplexDomain() { return new GmshSimplexDomain<Dim, Order>; }
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
Gmsh* createTensorizedDomain() { return new GmshTensorizedDomain<Dim, Order, RDim, Entity>; }

} // detail

//
// Simplex 1,1
//
const bool meshs11s = Gmsh::Factory::type::instance().registerProduct( "simplex(1,1)", &detail::createSimplexDomain<1,1> );
const bool meshs12s = Gmsh::Factory::type::instance().registerProduct( "simplex(1,2)", &detail::createSimplexDomain<1,2> );

const bool meshs21s = Gmsh::Factory::type::instance().registerProduct( "simplex(2,1)", &detail::createSimplexDomain<2,1> );
const bool meshs22s = Gmsh::Factory::type::instance().registerProduct( "simplex(2,2)", &detail::createSimplexDomain<2,2> );

const bool meshs31s = Gmsh::Factory::type::instance().registerProduct( "simplex(3,1)", &detail::createSimplexDomain<3,1> );
const bool meshs32s = Gmsh::Factory::type::instance().registerProduct( "simplex(3,2)", &detail::createSimplexDomain<3,2> );

const bool meshs112ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1,2,Simplex)",
                                                                         &detail::createTensorizedDomain<1,1,2,Simplex> );

const bool meshs11tsxx = Gmsh::Factory::type::instance().registerProduct( "hypercube(1,1)",
                                                                          &detail::createTensorizedDomain<1,1,1,Simplex> );
const bool meshs21tsxx = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1)",
                                                                          &detail::createTensorizedDomain<2,1,2,Simplex> );

const bool meshs21ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,Simplex)",
                                                                        &detail::createTensorizedDomain<2,1,2,Simplex> );
const bool meshs213ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,3,Simplex)",
                                                                         &detail::createTensorizedDomain<2,1,3,Simplex> );
const bool meshs22ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,2,Simplex)",
                                                                        &detail::createTensorizedDomain<2,2,2,Simplex> );

const bool meshs21tsp = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,1,SimplexProduct)",
                                                                         &detail::createTensorizedDomain<2,1,2,SimplexProduct> );
const bool meshs22tsp = Gmsh::Factory::type::instance().registerProduct( "hypercube(2,2,SimplexProduct)",
                                                                         &detail::createTensorizedDomain<2,2,2,SimplexProduct> );


const bool meshs31ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(3,1,Simplex)",
                                                                        &detail::createTensorizedDomain<3,1,3,Simplex> );
const bool meshs32ts = Gmsh::Factory::type::instance().registerProduct( "hypercube(3,2,Simplex)",
                                                                        &detail::createTensorizedDomain<3,2,3,Simplex> );

const bool meshs31tsp = Gmsh::Factory::type::instance().registerProduct( "hypercube(3,1,SimplexProduct)",
                                                                         &detail::createTensorizedDomain<3,1,3,SimplexProduct> );
const bool meshs32tsp = Gmsh::Factory::type::instance().registerProduct( "hypercube(3,2,SimplexProduct)",
                                                                         &detail::createTensorizedDomain<3,2,3,SimplexProduct> );


} // Life
