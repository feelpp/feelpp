/*
  This file is part of the Feel library.

  Author: Christophe Prud'homme (christophe.prudhomme@feelpp.org)

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2004 EPFL

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
#include <cstring>
#include <cstdlib>
#include <errno.h>

#include <list>
#include <map>
#include <vector>

#include <iterator>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>

#include <feel/feelconfig.h>

#ifdef FEELPP_HAS_BACKTRACE
# include <execinfo.h>
#endif

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/debug.hpp>

namespace Feel
{
namespace posix_time = boost::posix_time;

/*!

  *\ingroup Core
  *\brief Area debugging tool

  \c VLOG(1) provides a debug stream to which you can pass a number, say 100, associated
  to an area of the code, say a class \c A.  In the implementation of the class \c A, you use
  debug statement like

  void A::f()
  {
    DVLOG(2) << "A::f() is called.\n";
    // do something here
  }

  Now the debug message "A::f() is called." will be seen only if the area 100
  is defined in the environment(shell) variable \c DEBUG while executing a program
  \c A::f() is called \c runWithA that makes use of our class \c A.

  > runwithA
    --> no debug message related to A
  > export DEBUG="100"
  > runwithA
    A::f() is called.

   With this tool you can select the area you want to debug explicitly while keeping the
   others hidden.

  @author Christophe Prud'homme (christophe.prudhomme@feelpp.org)
*/
enum DebugLevels
{
    DEBUG_INFO  = 0,
    DEBUG_WARN  = 1,
    DEBUG_ERROR = 2,
    DEBUG_FATAL = 3
};
struct DebugStream::Private
{
    Private()
        :
        debug( false ),
        __flush_function( 0 )
    {}
    bool debug;
    std::ostringstream M_output;

    stprintf __flush_function;

    static bool _S_attached;
    static std::ofstream _S_logfile;
    //static std::map<int, std::ostream*> _S_logfile_per_area;

};
//
// getDescription
//
static std::map<unsigned int, std::string>* DebugAreas = 0;
static std::string* StringNull = 0;
static std::list<int>* AREAS;
static std::string* DEBUG_AREA = 0;

namespace
{
#define DEBUG_ADD_AREA( area, areastring ) \
 DebugAreas->insert( std::make_pair ( area, areastring ) )

// this Function makes sure that the static variables are initialized
// properly before being used
void
initDebugAreas ()
{
    static bool alloc = false;

    if ( alloc == false )
    {
        DEBUG_AREA = new std::string ( "" );
        AREAS = new std::list<int>;
        StringNull = new std::string ( "" );
        DebugAreas = new std::map<unsigned int, std::string>;
        alloc = true;

        DEBUG_ADD_AREA( 1000, "Feel::Application" );
        DEBUG_ADD_AREA( 1010, "Feel::Application" );
        DEBUG_ADD_AREA( 1020, "Feel::Application" );
        DEBUG_ADD_AREA( 1030, "Feel::Application" );
        DEBUG_ADD_AREA( 1100, "notyetassigned" );
        DEBUG_ADD_AREA( 1200, "notyetassigned" );
        DEBUG_ADD_AREA( 2200, "Feel::factory" );
        DEBUG_ADD_AREA( 3000, "Feel::array" );
        DEBUG_ADD_AREA( 3100, "Feel::exporters" );
        DEBUG_ADD_AREA( 4000, "Feel::mesh" );
        DEBUG_ADD_AREA( 4005, "Feel::meshEntity" );
        DEBUG_ADD_AREA( 4010, "Feel::RegionTree" );
        DEBUG_ADD_AREA( 4011, "Feel::KDTree" );
        DEBUG_ADD_AREA( 4012, "Feel::SubFaceOf" );
        DEBUG_ADD_AREA( 4015, "Feel::Mesh" );
        DEBUG_ADD_AREA( 4020, "Feel::PartitionerMetis" );
        DEBUG_ADD_AREA( 4021, "Feel::PartitionerMetis" );
        DEBUG_ADD_AREA( 4100, "Feel::mesh_util_base" );
        DEBUG_ADD_AREA( 5000, "Feel::fem" );
        DEBUG_ADD_AREA( 5005, "Feel::Dof" );
        DEBUG_ADD_AREA( 5010, "Feel::FunctionSpace" );
        DEBUG_ADD_AREA( 5015, "Feel::FunctionSpace::Element" );
        DEBUG_ADD_AREA( 5017, "Feel::BDF" );
        DEBUG_ADD_AREA( 5020, "Feel::Polynomial" );
        DEBUG_ADD_AREA( 5025, "Feel::Polynomial" );
        DEBUG_ADD_AREA( 5030, "Feel::FiniteElement" );
        DEBUG_ADD_AREA( 5032, "Operator" );
        DEBUG_ADD_AREA( 5033, "OperatorLinear" );
        DEBUG_ADD_AREA( 5034, "OperatorInterpolation" );
        DEBUG_ADD_AREA( 5035, "OperatorLagrangeP1" );
        DEBUG_ADD_AREA( 5040, "Feel::PolynomialSet" );
        DEBUG_ADD_AREA( 5042, "Feel::PolynomialSet" );
        DEBUG_ADD_AREA( 5045, "Feel::Lagrange" );
        DEBUG_ADD_AREA( 5046, "Feel::GeoMap" );
        DEBUG_ADD_AREA( 5047, "Feel::GeoMap" );
        DEBUG_ADD_AREA( 5048, "Feel::QuadRule" );
        DEBUG_ADD_AREA( 5050, "Feel::VF::Type" );
        DEBUG_ADD_AREA( 5051, "Feel::VF::Expr" );
        DEBUG_ADD_AREA( 5055, "Feel::VF::BilinearForm" );
        DEBUG_ADD_AREA( 5060, "Feel::VF::LinearForm" );
        DEBUG_ADD_AREA( 5065, "Feel::VF::Integrator" );
        DEBUG_ADD_AREA( 5066, "Feel::VF::IntegratorOn" );
        DEBUG_ADD_AREA( 5067, "Feel::VF::IntegratorOn" );
        DEBUG_ADD_AREA( 5080, "Feel::gauss1d" );
        DEBUG_ADD_AREA( 5100, "Feel::SolverUMFPACK" );
        DEBUG_ADD_AREA( 5600, "Feel::VectorUblas" );
        DEBUG_ADD_AREA( 5620, "Feel::LU" );
        DEBUG_ADD_AREA( 5800, "Feel::ImplicitFunction" );
        DEBUG_ADD_AREA( 6000, "Feel::solver" );
        DEBUG_ADD_AREA( 6010, "Feel::NavierStokesWithFlux" );
        DEBUG_ADD_AREA( 6020, "Feel::NavierStokesSolverPC" );
        DEBUG_ADD_AREA( 6100, "Feel::DarcySolver" );
        DEBUG_ADD_AREA( 6200, "Feel::operFS" );
        DEBUG_ADD_AREA( 6201, "Feel::operFS" );
        DEBUG_ADD_AREA( 6205, "Feel::exactJacobian" );
        DEBUG_ADD_AREA( 6210, "Feel::fixedPoint" );
        DEBUG_ADD_AREA( 6215, "Feel::steklovPoincare" );
        DEBUG_ADD_AREA( 6220, "Feel::FSISolver" );
        DEBUG_ADD_AREA( 7000, "Feel::alg" );
        DEBUG_ADD_AREA( 7002, "Feel::cg" );
        DEBUG_ADD_AREA( 7004, "Feel::gmres" );
        DEBUG_ADD_AREA( 7005, "Feel::BackendPetsc" );
        DEBUG_ADD_AREA( 7006, "Feel::BackendGmm" );
        DEBUG_ADD_AREA( 7010, "Feel::SolverLinearPetsc" );
        DEBUG_ADD_AREA( 7011, "Feel::VectorPetsc" );
        DEBUG_ADD_AREA( 7013, "Feel::MatrixPetsc" );
        DEBUG_ADD_AREA( 7015, "Feel::MatrixGmm" );
        DEBUG_ADD_AREA( 7020, "Feel::SolverNonLinearPetsc" );

        DEBUG_ADD_AREA( 7050, "Feel::SolverPardiso" );
        DEBUG_ADD_AREA( 8000, "Feel::TimeSet" );
        DEBUG_ADD_AREA( 8005, "Feel::TimeSet::Step" );
        DEBUG_ADD_AREA( 8005, "Feel::Exporter" );
        DEBUG_ADD_AREA( 8006, "Feel::ExporterEnsight" );
        DEBUG_ADD_AREA( 8007, "Feel::ExporterGmsh" );
        DEBUG_ADD_AREA( 8010, "Feel::ReadMesh1D" );
        DEBUG_ADD_AREA( 8011, "Feel::ReadMesh2D" );
        DEBUG_ADD_AREA( 8012, "Feel::ReadMesh3D" );
        DEBUG_ADD_AREA( 8013, "Feel::ImporterGambit" );
        DEBUG_ADD_AREA( 8098, "Feel::PointSetToMesh" );
        DEBUG_ADD_AREA( 8099, "Feel::FilterFromVtk" );
        DEBUG_ADD_AREA( 10000, "testsuite" );


        char * __env = getenv( "DEBUG" );

        if ( __env )
        {
            *DEBUG_AREA = __env;
        }

        std::istringstream __is ( *DEBUG_AREA );


        std::copy ( std::istream_iterator<int,char> ( __is ),
                    std::istream_iterator<int,char> (),
                    std::back_inserter ( *AREAS ) );

    }
}
std::string
getDescription ( unsigned int __area )
{
    if ( DebugAreas->empty() )
        return std::string( "Area " ) + boost::lexical_cast<std::string>( __area );

    std::map<unsigned int, std::string>::iterator entry_it = DebugAreas->find ( __area );

    if ( entry_it != DebugAreas->end() )
        return entry_it->second;

    else
        return std::string( "Area " ) + boost::lexical_cast<std::string>( __area );


}
}

//
// DebugStream
//
DebugStream::DebugStream( int area, int level, bool print )
    :
    __p( new Private )
{
    initDebugAreas ();

    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug =  ( ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() && print ) ||
                        !area );
    }

    else
    {
        __p->debug =  ( print && !area );
    }

    if ( __p->debug && level == DEBUG_INFO )
    {
        posix_time::ptime __time( posix_time::second_clock::local_time() );

        if ( area )
            __p->M_output << "[" << getDescription ( area ) << "] ";

        //<< posix_time::to_simple_string( __time )<< ") : ";
    }

}
DebugStream::DebugStream( const char* initialString, int area, int level, bool print )
    :
    __p( new Private )
{
    initDebugAreas ();

    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug =  ( ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() &&
                          print ) || !area );
    }

    else
    {
        __p->debug =  ( print && !area );
    }

    if ( __p->debug && level == DEBUG_INFO )
    {
        posix_time::ptime __time( posix_time::second_clock::local_time() );

        if ( area )
            __p->M_output << "[" << getDescription ( area ) << "] "
                           //<< posix_time::to_simple_string( __time )<< ") : "
                           << initialString;
    }
}
DebugStream::DebugStream( const DebugStream& sd )
    :
    __p( new Private )
{
    __p->debug = sd.__p->debug;
    __p->__flush_function = sd.__p->__flush_function;
}
DebugStream::~DebugStream()
{
    delete __p;
}

bool
DebugStream::doPrint() const
{
    return __p->debug;
}

DebugStream&
DebugStream::operator<<( double s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
DebugStream&
DebugStream::operator<<( std::complex<double> s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
#if defined(FEELPP_HAS_QD_H)
DebugStream&
DebugStream::operator<<( dd_real s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
DebugStream&
DebugStream::operator<<( qd_real s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
#endif /* FEELPP_HAS_QD_H */
DebugStream&
DebugStream::operator<<( bool s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}

DebugStream&
DebugStream::operator<<( uint16_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}

DebugStream&
DebugStream::operator<<( uint32_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
#if defined (__s390x__) || defined( __s390__ ) || defined( __APPLE__ )
DebugStream&
DebugStream::operator<<( size_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
#endif
#if defined( __APPLE__ )
DebugStream&
DebugStream::operator<<( ptrdiff_t s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
#endif
DebugStream&
DebugStream::operator<<( uint64_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}

DebugStream&
DebugStream::operator<<( int16_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}

DebugStream&
DebugStream::operator<<( int32_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}
DebugStream&
DebugStream::operator<<( int64_type s )
{
    if ( __p->debug )
    {
        __p->M_output  << s;
        flush();
    }

    return *this;
}

DebugStream&
DebugStream::operator<<( const char* s )
{
    if ( __p->debug )
        __p->M_output  << s;

    flush();
    return *this;
}
DebugStream&
DebugStream::operator<<( std::string const& s )
{
    if ( __p->debug )
        __p->M_output  << s;

    size_t found = s.find( '\n' );

    if ( found != std::string::npos )
        flush();

    return *this;
}


DebugStream&
DebugStream::operator<<( Feel::LManipFunction __f )
{
    if ( __p->debug )
    {
        ( *__f )( *this );
    }

    return *this;
}

void
DebugStream::setFlush( stprintf func )
{
    __p->__flush_function = func;
}
void
DebugStream::flush(  )
{
    if ( !__p->M_output.str().empty() )
    {
        if ( Private::_S_attached )
        {
            Private::_S_logfile << __p->M_output.str();
            Private::_S_logfile.flush();
        }

        else if ( __p->__flush_function == 0 )
        {
            std::cerr << __p->M_output.str();
        }

        else
        {
            __p->__flush_function( "%s", __p->M_output.str().c_str() );
        }

        __p->M_output.str( "" );
    }

}
void
DebugStream::addDebugArea( uint16_type area, std::string const& description )
{
    DEBUG_ADD_AREA( area, description );
}
void
DebugStream::showDebugAreas( std::string const& areas )
{
    // make sure that we get a space between each area when we concatenate
    *DEBUG_AREA += " ";
    *DEBUG_AREA += areas;
    std::cout << "DEBUG_AREA = " << *DEBUG_AREA << std::endl;

    std::istringstream __is ( *DEBUG_AREA );

    // clear first
    AREAS->clear();
    // inset in AREAS
    std::copy ( std::istream_iterator<int,char> ( __is ),
                std::istream_iterator<int,char> (),
                std::back_inserter ( *AREAS ) );
}

bool DebugStream::Private::_S_attached = false;
std::ofstream DebugStream::Private::_S_logfile;
//std::map<int, std::ostream*> DebugStream::_S_logfile_per_area;

void DebugStream::attach( std::string const& __logfile )
{
    std::ostringstream __filename;
    __filename <<  __logfile;

    if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile.close();
    }

    Private::_S_logfile.open( __filename.str().c_str(), std::ios::out );

    if ( Private::_S_logfile.fail() )
    {
        Warning() << "DebugStream::attach( " << __logfile.c_str() << " ) failed to open "  << __filename.str() << "\n";
        Warning() << "Redirecting to default output\n";
        Private::_S_attached = false;
    }

    else if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile << __filename.str() << " is opened for debug" << std::endl;
        Private::_S_attached = true;
    }
}
void
DebugStream::attach( std::string const& /*__logfile*/, int /*__area*/ )
{

}
void
DebugStream::detach( std::string const& /*__logfile*/, int /*__area*/ )
{}

void
DebugStream::detachAll()
{
    if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile.close();
        Private::_S_attached = false;
    }
}

DebugStream
Log( int area, DebugStream::stprintf func )
{
    DebugStream s( area, DEBUG_INFO );
    s.setFlush( func );
    return s;
}

DebugStream
Log( bool cond, int area, DebugStream::stprintf /*func*/ )
{
    if ( cond )
        return DebugStream( area, DEBUG_INFO );

    else
        return DebugStream( 0, 0, false );
}

#ifndef NDEBUG
DebugStream
Debug( int area, DebugStream::stprintf func )
{
    DebugStream s( area, DEBUG_INFO );
    s.setFlush( func );
    return s;
}

DebugStream
Debug( bool cond, int area, DebugStream::stprintf /*func*/ )
{
    if ( cond )
        return DebugStream( area, DEBUG_INFO );

    else
        return DebugStream( 0, 0, false );
}
#endif

DebugStream
Warning( int area )
{
    return DebugStream( "WARNING: ", area, DEBUG_WARN );
}

DebugStream
Warning( bool cond, int area )
{
    if ( cond )
        return DebugStream( "WARNING: ", area, DEBUG_WARN );

    else
        return DebugStream( 0, 0, false );

}

DebugStream
Error( int area )
{
    //DVLOG(2) << LBacktrace() << "\n";
    return DebugStream( "ERROR: ", area, DEBUG_ERROR );
}

DebugStream
Error( bool cond, int area )
{
    //DVLOG(2) << LBacktrace() << "\n";
    if ( cond )
        return DebugStream( "ERROR: ", area, DEBUG_ERROR );

    else
        return DebugStream( 0, 0, false );

}

DebugStream
Fatal( int area )
{
    //LBacktrace();
    return DebugStream( "FATAL: ", area, DEBUG_FATAL );
}

DebugStream
Fatal( bool cond, int area )
{
    //LBacktrace();
    if ( cond )
        return DebugStream( "FATAL: ", area, DEBUG_FATAL );

    else
        return DebugStream( 0, 0, false );
}

std::string
backtrace ()
{
    // show all backtrace
    return backtrace( -1 );
}
std::string
backtrace ( int /*__levels*/ )
{
    std::ostringstream os;
#ifdef FEELPP_HAS_BACKTRACE

    void* trace[256];
    int n = backtrace ( trace, 256 );
    char** strings = backtrace_symbols ( trace, n );

    if ( __levels != -1 )
        n = ( std::min ) ( n, __levels );

    os << "[\n";

    for ( int i = 0; i < n; ++i )
        os << i << ": " << strings[i] << "\n";

    os << "]\n";
    free ( strings );
#endif
    return os.str();
}
}

Feel::DebugStream&
perror( Feel::DebugStream& s )
{
    s << " " << strerror( errno );
    return s;
}
Feel::DebugStream&
endl( Feel::DebugStream& s )
{
    s << "\n";
    return s;
}
Feel::DebugStream&
flush( Feel::DebugStream& s )
{
    s.flush();
    return s;
}
