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
#ifndef __Debug_H
#define __Debug_H 1

#include <cstdio>
#include <iosfwd>

#include <string>
#include <sstream>

#include <feel/feelcore/feelmacros.hpp>


namespace Feel
{
class DebugStream;
class NdebugStream;

typedef DebugStream & ( *LManipFunction )( DebugStream & ); // manipulator function
typedef NdebugStream & ( *LNManipFunction )( NdebugStream& ); // manipulator function

#ifdef __GNUC__
# define FEELPP_FUNCINFO "[" << __PRETTY_FUNCTION__ << "] "
#else
# define FEELPP_FUNCINFO "[" << __FILE__ << ":" << __LINE__ << "] "
#endif

#define FEELPP_LINEINFO "[" << __FILE__ << ":" << __LINE__ << "] "

namespace detail
{
#if 0
class print
{
public:
    explicit print ( std::ostream& __os ) : M_os ( __os ) {}

    template<typename T>
    std::ostream& operator<< ( T const& __t )
    {
        __print ( __t, St::SInt2Type<St::STypeTraits<T>::isFundamental>() );
        return M_os;
    }
    template<typename T>
    std::ostream& operator<< ( T const* __t )
    {
        //__print ( __t, St::SInt2Type<St::STypeTraits<T>::isFundamental>() );
        M_os << __t;
        return M_os;
    }
private:
    template<typename T>
    void __print ( T const& __t, St::SInt2Type<true> )
    {
        M_os << __t;
    }
    template<typename T>
    void __print ( T const& __t, St::SInt2Type<false> )
    {

    }
private:
    std::ostream& M_os;
};
#endif
} // end namespace detail

class DebugStream
{
public:


    /** @name Internal Structures
     */
    //@{

    struct Private;

    typedef int ( *stprintf )( const char* format, ... );
    //@}

    /** @name Constructors, destructor
     */
    //@{

    DebugStream( int area = 0, int level = 1, bool print = true );
    DebugStream( const char* initialString, int area = 0, int level = 1, bool print = true );
    DebugStream( DebugStream const& );
    ~DebugStream();

    //@}

    /** @name  Methods
     */
    //@{

    /** \return true if the stream prints the debug output, false otherwise */
    bool doPrint() const;

    //@}
    void setFlush( stprintf = 0 );

    /** @name  Methods
     */
    //@{

    static void addDebugArea( uint16_type area, std::string const& description );
    static void showDebugAreas( std::string const& areas );
    static void attach( std::string const& __logfile );
    static void attach( std::string const& __logfile, int area );
    static void detach( std::string const& __logfile, int area );
    static void detachAll();
    void flush();


    DebugStream& operator<<( bool );
    DebugStream& operator<<( int16_type );
    DebugStream& operator<<( int32_type );
    DebugStream& operator<<( uint16_type );
    DebugStream& operator<<( uint32_type );
#if defined (__s390x__) || defined( __s390__ )
    DebugStream& operator<<( size_type );
#endif
#if defined( __APPLE__ )
    DebugStream& operator<<( size_type );
    DebugStream& operator<<( ptrdiff_t );
#endif

#if !defined( BOOST_NO_INT64_T )
    DebugStream& operator<<( int64_type );
    DebugStream& operator<<( uint64_type );
#endif

    DebugStream& operator<<( double );
    DebugStream& operator<<( std::complex<double> );
#if defined(FEELPP_HAS_QD_H)
    DebugStream& operator<<( dd_real );
    DebugStream& operator<<( qd_real );
#endif /* FEELPP_HAS_QD_H */

    DebugStream& operator<<( const char* );
    DebugStream& operator<<( std::string const& );
    DebugStream& operator<<( LManipFunction f );
    //@}

protected:

private:
    Private* __p;

};

template<typename T>
DebugStream& operator<< ( DebugStream& __s, T const* __t )
{
    std::ostringstream __os;
    __os << __t;
    __s << __os.str();
    return __s;
}
std::string backtrace ();
std::string backtrace ( int );

class NdebugStream
{
public:
    /** @name Constructors, destructor
     */
    //@{
    typedef int ( *stprintf )( const char* format, ... );

    NdebugStream() {}
    ~NdebugStream() {}

    //@}

    /** @name  Methods
     */
    //@{
    static void attach( std::string const& ) {}
    static void attach( std::string const&, int ) {}
    static void detach( std::string const&, int ) {}
    static void detachAll() {}
    void flush( stprintf = 0 ) {}
    NdebugStream& operator<<( char const* )
    {
        return *this;
    }

    NdebugStream& operator<<( bool )
    {
        return *this;
    }
    NdebugStream& operator<<( int16_type )
    {
        return *this;
    }
    NdebugStream& operator<<( int32_type )
    {
        return *this;
    }
    NdebugStream& operator<<( uint16_type )
    {
        return *this;
    }
    NdebugStream& operator<<( uint32_type )
    {
        return *this;
    }
#if defined (__s390x__) || defined( __s390__ )
    NdebugStream& operator<<( size_type )
    {
        return *this;
    }
#endif
#if defined( __APPLE__ )
    NdebugStream& operator<<( size_type )
    {
        return *this;
    }
    NdebugStream& operator<<( ptrdiff_t )
    {
        return *this;
    }
#endif

#if !defined( BOOST_NO_INT64_T )
    NdebugStream& operator<<( uint64_type )
    {
        return *this;
    }
    NdebugStream& operator<<( int64_type )
    {
        return *this;
    }
#endif

    NdebugStream& operator<<( double )
    {
        return *this;
    }
    NdebugStream& operator<<( std::complex<double> )
    {
        return *this;
    }
#if defined(FEELPP_HAS_QD_H)
    NdebugStream& operator<<( dd_real )
    {
        return *this;
    }
    NdebugStream& operator<<( qd_real )
    {
        return *this;
    }
#endif /* FEELPP_HAS_QD_H */

    NdebugStream& operator<<( std::string const& )
    {
        return *this;
    }
    NdebugStream& operator<<( LManipFunction  )
    {
        return *this;
    }

    //@}
};

inline NdebugStream& perror( NdebugStream& s )
{
    return s;
}
inline NdebugStream& endl( NdebugStream& s )
{
    return s;
}
inline NdebugStream& flush( NdebugStream& s )
{
    return s;
}

DebugStream Log( int area = 0, DebugStream::stprintf = 0 ) FEELPP_DEPRECATED;
DebugStream Log( bool cond, int area = 0, DebugStream::stprintf = 0 ) FEELPP_DEPRECATED;


#ifndef NDEBUG
DebugStream Debug( int area = 0, DebugStream::stprintf = 0 ) FEELPP_DEPRECATED;
DebugStream Debug( bool cond, int area = 0, DebugStream::stprintf = 0 ) FEELPP_DEPRECATED;
#else
#define Debug Ndebug
inline NdebugStream Ndebug( int = 0, NdebugStream::stprintf = &printf )
{
    return NdebugStream();
}
inline NdebugStream Ndebug( bool /*cond*/, int = 0, NdebugStream::stprintf = &printf )
{
    return NdebugStream();
}
#endif


DebugStream Warning( int area = 0 );
DebugStream Warning( bool cond, int area = 0 );

DebugStream Error( int area = 0 );
DebugStream Error( bool cond, int area = 0 );

DebugStream Fatal( int area = 0 );
DebugStream Fatal( bool cond, int area = 0 );

/// Standard function announcer
#define FEELPP_DEBUG_FUNC_INFO(area) Debug(area) << FEELPP_FUNCINFO << "\n";

/// Use these to introduce and extroduce functions
#define FEELPP_DEBUG_BEGIN(area) Debug(area) << "BEGIN: " << __PRETTY_FUNCTION__ << "\n";
#define FEELPP_DEBUG_END(area)   Debug(area) << "END: " << __PRETTY_FUNCTION__ << "\n";

}




Feel::DebugStream& perror( Feel::DebugStream& s );
Feel::DebugStream& endl( Feel::DebugStream& s );
Feel::DebugStream& flush( Feel::DebugStream& );


#endif /* __Debug_H */
