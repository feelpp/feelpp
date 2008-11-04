/*
  This file is part of the Life library.

  Author: Christophe Prud'homme (christophe.prudhomme@ujf-grenoble.fr)

  Copyright (C) 2004 EPFL

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

#include <lifeconfig.h>

#ifdef HAVE_BACKTRACE
# include <execinfo.h>
#endif

#include <life/lifecore/life.hpp>

namespace Life
{
namespace posix_time = boost::posix_time;

/*!
  \class Debug
  *\ingroup Core
  *\brief Area debugging tool

  \c Debug() provides a debug stream to which you can pass a number, say 100, associated
  to an area of the code, say a class \c A.  In the implementation of the class \c A, you use
  debug statement like

  void A::f()
  {
    Debug(100) << "A::f() is called.\n";
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

  @author Christophe Prud'homme (christophe.prudhomme@ujf-grenoble.fr)
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
    std::ostringstream _M_output;

    stprintf __flush_function;

    static bool _S_attached;
    static std::ofstream _S_logfile;
    //static std::map<int, std::ostream*> _S_logfile_per_area;

};
//
// getDescription
//
static std::map<uint, std::string>* DebugAreas = 0;
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
            DebugAreas = new std::map<uint, std::string>;
            alloc = true;

            DEBUG_ADD_AREA( 1000, "Life::Application" );
            DEBUG_ADD_AREA( 1010, "Life::Application" );
            DEBUG_ADD_AREA( 1020, "Life::Application" );
            DEBUG_ADD_AREA( 1030, "Life::Application" );
            DEBUG_ADD_AREA( 1100, "notyetassigned" );
            DEBUG_ADD_AREA( 1200, "notyetassigned" );
            DEBUG_ADD_AREA( 2200, "Life::factory" );
            DEBUG_ADD_AREA( 3000, "Life::array" );
            DEBUG_ADD_AREA( 3100, "Life::exporters" );
            DEBUG_ADD_AREA( 4000, "Life::mesh" );
            DEBUG_ADD_AREA( 4005, "Life::meshEntity" );
            DEBUG_ADD_AREA( 4010, "Life::RegionTree" );
            DEBUG_ADD_AREA( 4011, "Life::KDTree" );
            DEBUG_ADD_AREA( 4012, "Life::SubFaceOf" );
            DEBUG_ADD_AREA( 4015, "Life::Mesh" );
            DEBUG_ADD_AREA( 4020, "Life::PartitionerMetis" );
            DEBUG_ADD_AREA( 4021, "Life::PartitionerMetis" );
            DEBUG_ADD_AREA( 4100, "Life::mesh_util_base" );
            DEBUG_ADD_AREA( 5000, "Life::fem" );
            DEBUG_ADD_AREA( 5005, "Life::Dof" );
            DEBUG_ADD_AREA( 5010, "Life::FunctionSpace" );
            DEBUG_ADD_AREA( 5015, "Life::FunctionSpace::Element" );
            DEBUG_ADD_AREA( 5017, "Life::BDF" );
            DEBUG_ADD_AREA( 5020, "Life::Polynomial" );
            DEBUG_ADD_AREA( 5025, "Life::Polynomial" );
            DEBUG_ADD_AREA( 5030, "Life::FiniteElement" );
            DEBUG_ADD_AREA( 5032, "Operator" );
            DEBUG_ADD_AREA( 5033, "OperatorLinear" );
            DEBUG_ADD_AREA( 5034, "OperatorInterpolation" );
            DEBUG_ADD_AREA( 5035, "OperatorLagrangeP1" );
            DEBUG_ADD_AREA( 5040, "Life::PolynomialSet" );
            DEBUG_ADD_AREA( 5042, "Life::PolynomialSet" );
            DEBUG_ADD_AREA( 5045, "Life::Lagrange" );
            DEBUG_ADD_AREA( 5046, "Life::GeoMap" );
            DEBUG_ADD_AREA( 5047, "Life::GeoMap" );
            DEBUG_ADD_AREA( 5048, "Life::QuadRule" );
            DEBUG_ADD_AREA( 5050, "Life::VF::Type" );
            DEBUG_ADD_AREA( 5051, "Life::VF::Expr" );
            DEBUG_ADD_AREA( 5055, "Life::VF::BilinearForm" );
            DEBUG_ADD_AREA( 5060, "Life::VF::LinearForm" );
            DEBUG_ADD_AREA( 5065, "Life::VF::Integrator" );
            DEBUG_ADD_AREA( 5066, "Life::VF::IntegratorOn" );
            DEBUG_ADD_AREA( 5067, "Life::VF::IntegratorOn" );
            DEBUG_ADD_AREA( 5080, "Life::gauss1d" );
            DEBUG_ADD_AREA( 5100, "Life::SolverUMFPACK" );
            DEBUG_ADD_AREA( 5600, "Life::VectorUblas" );
            DEBUG_ADD_AREA( 5620, "Life::LU" );
            DEBUG_ADD_AREA( 5800, "Life::ImplicitFunction" );
            DEBUG_ADD_AREA( 6000, "Life::solver" );
            DEBUG_ADD_AREA( 6010, "Life::NavierStokesWithFlux" );
            DEBUG_ADD_AREA( 6020, "Life::NavierStokesSolverPC" );
            DEBUG_ADD_AREA( 6100, "Life::DarcySolver" );
            DEBUG_ADD_AREA( 6200, "Life::operFS" );
            DEBUG_ADD_AREA( 6201, "Life::operFS" );
            DEBUG_ADD_AREA( 6205, "Life::exactJacobian" );
            DEBUG_ADD_AREA( 6210, "Life::fixedPoint" );
            DEBUG_ADD_AREA( 6215, "Life::steklovPoincare" );
            DEBUG_ADD_AREA( 6220, "Life::FSISolver" );
            DEBUG_ADD_AREA( 7000, "Life::alg" );
            DEBUG_ADD_AREA( 7002, "Life::cg" );
            DEBUG_ADD_AREA( 7004, "Life::gmres" );
            DEBUG_ADD_AREA( 7005, "Life::BackendPetsc" );
            DEBUG_ADD_AREA( 7006, "Life::BackendGmm" );
            DEBUG_ADD_AREA( 7010, "Life::SolverLinearPetsc" );
            DEBUG_ADD_AREA( 7011, "Life::VectorPetsc" );
            DEBUG_ADD_AREA( 7013, "Life::MatrixPetsc" );
            DEBUG_ADD_AREA( 7015, "Life::MatrixGmm" );
            DEBUG_ADD_AREA( 7020, "Life::SolverNonLinearPetsc" );

            DEBUG_ADD_AREA( 7050, "Life::SolverPardiso" );
            DEBUG_ADD_AREA( 8000, "Life::TimeSet" );
            DEBUG_ADD_AREA( 8005, "Life::TimeSet::Step" );
            DEBUG_ADD_AREA( 8005, "Life::Exporter" );
            DEBUG_ADD_AREA( 8006, "Life::ExporterEnsight" );
            DEBUG_ADD_AREA( 8007, "Life::ExporterGmsh" );
            DEBUG_ADD_AREA( 8010, "Life::ReadMesh1D" );
            DEBUG_ADD_AREA( 8011, "Life::ReadMesh2D" );
            DEBUG_ADD_AREA( 8012, "Life::ReadMesh3D" );
            DEBUG_ADD_AREA( 8013, "Life::ImporterGambit" );
            DEBUG_ADD_AREA( 8098, "Life::PointSetToMesh" );
            DEBUG_ADD_AREA( 8099, "Life::FilterFromVtk" );
            DEBUG_ADD_AREA( 10000, "testsuite" );


            char * __env = getenv("DEBUG");
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
getDescription ( uint __area )
{
    if ( DebugAreas->empty() )
        return std::string( "Area " ) + boost::lexical_cast<std::string>(__area);

    std::map<uint, std::string>::iterator entry_it = DebugAreas->find ( __area );

    if ( entry_it != DebugAreas->end() )
        return entry_it->second;
    else
        return std::string( "Area " ) + boost::lexical_cast<std::string>(__area);


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
        __p->debug =  ( (std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() && print) ||
                        !area );
    }
    else
    {
        __p->debug =  ( print && !area );
    }
    if ( __p->debug && level == DEBUG_INFO )
    {
        posix_time::ptime __time( posix_time::second_clock::local_time() );
        __p->_M_output << "[" << getDescription ( area ) << "] ";
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
        __p->_M_output << "[" << getDescription ( area ) << "] "
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
DebugStream::operator<<( double s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
DebugStream&
DebugStream::operator<<( std::complex<double> s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
#if defined(HAVE_QD_H)
DebugStream&
DebugStream::operator<<( dd_real s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
DebugStream&
DebugStream::operator<<( qd_real s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
#endif /* HAVE_QD_H */
DebugStream&
DebugStream::operator<<( bool s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( uint16_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( uint32_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
DebugStream&
DebugStream::operator<<( uint64_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( int16_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( int32_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
DebugStream&
DebugStream::operator<<( int64_type s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( const char* s)
{
    if ( __p->debug )
        __p->_M_output  << s;
    flush();
    return *this;
}
DebugStream&
DebugStream::operator<<( std::string const& s)
{
    if ( __p->debug )
        __p->_M_output  << s;
    if ( s[s.size() -1] == '\n')
        flush();
    return *this;
}


DebugStream&
DebugStream::operator<<( Life::LManipFunction __f )
{
    if ( __p->debug )
    {
        (*__f)( *this );
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
    if ( !__p->_M_output.str().empty() )
    {
        if ( Private::_S_attached )
        {
            Private::_S_logfile << __p->_M_output.str();
            Private::_S_logfile.flush();
        }
        else if ( __p->__flush_function == 0 )
        {
            std::cerr << __p->_M_output.str();
        }
        else
        {
            __p->__flush_function( "%s", __p->_M_output.str().c_str() );
        }
        __p->_M_output.str( "" );
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

#ifndef NDEBUG_OLD
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
    //Debug () << LBacktrace() << "\n";
    return DebugStream( "ERROR: ", area, DEBUG_ERROR );
}

DebugStream
Error( bool cond, int area )
{
    //Debug () << LBacktrace() << "\n";
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
#ifdef HAVE_BACKTRACE

    void* trace[256];
    int n = backtrace ( trace, 256 );
    char** strings = backtrace_symbols ( trace, n );

    if ( __levels != -1 )
        n = ( std::min ) ( n, __levels );
    os << "[\n";

    for (int i = 0; i < n; ++i)
        os << i << ": " << strings[i] << "\n";
    os << "]\n";
    free (strings);
#endif
    return os.str();
}
}

Life::DebugStream&
perror( Life::DebugStream& s )
{
    s << " " << strerror( errno ); return s;
}
Life::DebugStream&
endl( Life::DebugStream& s )
{
    s << "\n"; return s;
}
Life::DebugStream&
flush( Life::DebugStream& s )
{
    s.flush(); return s;
}

