/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-16

  Copyright (C) 2009 Université de Grenoble 1
  Copyright (C) 2005,2006 EPFL

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
   \file SmartAssert.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-16
 */
#if !defined(SMART_ASSERT_H)
#define SMART_ASSERT_H

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <map>
#if defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning(disable:780)
#endif
#include <glog/logging.h>
#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

namespace Feel
{
enum
{
    lvl_info = google::GLOG_INFO,

    // default behavior - just loggs this assert
    // (a message is shown to the user to the console)
    lvl_warn = google::GLOG_WARNING,

    // default behavior - asks the user what to do:
    // Ignore/ Retry/ etc.
    lvl_debug = google::GLOG_INFO,

    // default behavior - throws a SmartAssert_error
    lvl_error = google::GLOG_ERROR,

    // default behavior - dumps all assert context to console,
    // and aborts
    lvl_fatal = google::GLOG_FATAL
};



/**
   \class AssertContext
  *\ingroup Core
  *\brief contains details about a failed assertion
*/
class AssertContext
{
    typedef std::string string;
public:
    AssertContext() : M_level( lvl_info )
    {}

    // where the assertion failed: file & line
    void setFileLine( const char * file, int line )
    {
        M_file = file;
        M_line = line;
    }
    const string & getContextFile() const
    {
        return M_file;
    }
    int getContextLine() const
    {
        return M_line;
    }

    // get/ set expression
    void setExpression( const string & str )
    {
        M_expression = str;
    }
    const string & expression() const
    {
        return M_expression;
    }

    typedef std::pair< string, string> val_and_str;
    typedef std::vector< val_and_str> vals_array;
    // return values array as a vector of pairs:
    // [Value, corresponding string]
    const vals_array & get_vals_array() const
    {
        return M_vals;
    }
    // adds one value and its corresponding string
    void add_val( const string & val, const string & str )
    {
        M_vals.push_back( val_and_str( val, str ) );
    }

    // get/set level of assertion
    void setLevel( int nLevel )
    {
        M_level = nLevel;
    }
    int get_level() const
    {
        return M_level;
    }

    // get/set (user-friendly) message
    void setLevelMsg( const char * strMsg )
    {
        if ( strMsg )
            M_msg = strMsg;

        else
            M_msg.erase();
    }
    const string & get_level_msg() const
    {
        return M_msg;
    }

private:
    // where the assertion occured
    string M_file;
    int M_line;

    // expression and values
    string M_expression;
    vals_array M_vals;

    // level and message
    int M_level;
    string M_msg;
};


namespace SmartAssert
{

typedef void ( *assert_function_type )( const AssertContext & context );

// helpers
std::string getTypeofLevel( int nLevel );
void dumpContextSummary( const AssertContext & context, std::ostream & out );
void dumpContextDetail( const AssertContext & context, std::ostream & out );

// defaults
void defaultWarnHandler( const AssertContext & context );
void defaultDebugHandler( const AssertContext & context );
void defaultErrorHandler( const AssertContext & context );
void defaultFatalHandler( const AssertContext & context );
void defaultLogger( const AssertContext & context );

} // namespace SmartAssert

namespace Private
{
void initAssert();
void setDefaultLogStream( std::ostream & out );
void setDefaultLogName( const char * str );

// allows finding if a value is of type 'const char *'
// and is null; if so, we cannot print it to an ostream
// directly!!!
template< class T>
struct isNullFinder
{
    bool is( const T & ) const
    {
        return false;
    }
};

template<>
struct isNullFinder< char*>
{
    bool is( char * const & val )
    {
        return val == 0;
    }
};

template<>
struct isNullFinder< const char*>
{
    bool is( const char * const & val )
    {
        return val == 0;
    }
};


} // namespace Private


struct Assert
{
    typedef SmartAssert::assert_function_type assert_function_type;

    // helpers, in order to be able to compile the code
    Assert & SMART_ASSERT_A;
    Assert & SMART_ASSERT_B;

    Assert( const char * expr )
        : SMART_ASSERT_A( *this ),
          SMART_ASSERT_B( *this ),
          M_needs_handling( true )
    {
        M_context.setExpression( expr );

        if ( ( logger() == 0 ) || handlers().size() < 4 )
        {
            // used before main!
            Private::initAssert();
        }
    }

    Assert( const Assert & other )
        : SMART_ASSERT_A( *this ),
          SMART_ASSERT_B( *this ),
          M_context( other.M_context ),
          M_needs_handling( true )
    {
        other.M_needs_handling = false;
    }

    ~Assert()
    {
        if ( M_needs_handling )
            handleAssert();
    }

    template< class type>
    Assert & printCurrentValue( const type & val, const char * msg );

    Assert & printContext( const char * file, int line )
    {
        M_context.setFileLine( file, line );
        return *this;
    }

    Assert & msg( const char * strMsg )
    {
        M_context.setLevelMsg( strMsg );
        return *this;
    }

    Assert & level( int nLevel, const char * strMsg = 0 )
    {
        M_context.setLevel( nLevel );
        M_context.setLevelMsg( strMsg );
        return *this;
    }

    Assert & info( const char * strMsg = 0 )
    {
        return level( lvl_info, strMsg );
    }
    Assert & warn( const char * strMsg = 0 )
    {
        return level( lvl_warn, strMsg );
    }

    Assert & debug( const char * strMsg = 0 )
    {
        return level( lvl_debug, strMsg );
    }

    Assert & error( const char * strMsg = 0 )
    {
        return level( lvl_error, strMsg );
    }

    Assert & fatal( const char * strMsg = 0 )
    {
        return  level( lvl_fatal, strMsg );
    }

    // in this case, we set the default logger, and make it
    // write everything to this file
    static void setLog( const char * strFileName )
    {
        Private::setDefaultLogName( strFileName );
        logger() = &SmartAssert::defaultLogger;
    }

    // in this case, we set the default logger, and make it
    // write everything to this log
    static void setLog( std::ostream & out )
    {
        Private::setDefaultLogStream( out );
        logger() = &SmartAssert::defaultLogger;
    }

    static void setLog( assert_function_type log )
    {
        logger() = log;
    }

    static void setHandler( int nLevel, assert_function_type handler )
    {
        handlers()[ nLevel] = handler;
    }

private:
    // handles the current assertion.
    void handleAssert()
    {
        logger()( M_context );
        get_handler( M_context.get_level() )( M_context );
    }

    /*
      IMPORTANT NOTE:
      The only reason logger & handlers are functions, are
      because you might use SMART_ASSERT before main().

      In this case, since they're statics, they might not
      be initialized. However, making them functions
      will make it work.
    */

    // the log
    static assert_function_type & logger()
    {
        static assert_function_type inst;
        return inst;
    }

    // the handler
    typedef std::map< int, assert_function_type> handlers_collection;
    static handlers_collection & handlers()
    {
        static handlers_collection inst;
        return inst;
    }

    static assert_function_type get_handler( int nLevel )
    {
        handlers_collection::const_iterator found = handlers().find( nLevel );

        if ( found != handlers().end() )
            return found->second;

        else
            // we always assume the debug handler has been set
            return handlers().find( lvl_debug )->second;
    }

private:
    AssertContext M_context;
    mutable bool M_needs_handling;

};

template<class type>
Assert &
Assert::printCurrentValue( const type & _val, const char * _msg )
{
    std::ostringstream out;

    Private::isNullFinder< type> f;
    bool bIsNull = f.is( _val );

    if ( !bIsNull )
        out << _val;

    else
        // null string
        out << "null";

    M_context.add_val( out.str(), _msg );
    return *this;
}


namespace SmartAssert
{
inline ::Feel::Assert make_assert( const char * expr )
{
    return ::Feel::Assert( expr );
}
} // namespace SmartAssert

} // Feel Namespace

#ifdef FEELPP_SMART_ASSERT_DEBUG_MODE

#if FEELPP_SMART_ASSERT_DEBUG_MODE == 1
#define FEELPP_SMART_ASSERT_DEBUG
#else
#undef FEELPP_SMART_ASSERT_DEBUG
#endif

#else

// defaults
#ifndef NDEBUG
#define FEELPP_SMART_ASSERT_DEBUG
#else
#undef FEELPP_SMART_ASSERT_DEBUG
#endif

#endif


#ifdef FEELPP_SMART_ASSERT_DEBUG
// "debug" mode
#define FEELPP_SMART_ASSERT( expr) \
    if ( (expr) ) ; \
    else ::Feel::SmartAssert::make_assert( #expr).printContext( __FILE__, __LINE__).SMART_ASSERT_A \
    /**/

#else
// "release" mode
#define FEELPP_SMART_ASSERT( expr) \
    if ( true ) ; \
    else ::Feel::SmartAssert::make_assert( "").SMART_ASSERT_A \
    /**/

#endif // ifdef FEELPP_SMART_ASSERT_DEBUG

// FEELPP_ASSERT is a equivalent to FEELPP_SMART_ASSERT
#define FEELPP_ASSERT( expr) FEELPP_SMART_ASSERT(expr)


#define FEELPP_SMART_VERIFY( expr) \
    if ( (expr) ) ; \
    else ::Feel::SmartAssert::make_assert( #expr).error().printContext( __FILE__, __LINE__).SMART_ASSERT_A \
    /**/
#define FEELPP_VERIFY( expr) FEELPP_SMART_VERIFY(expr)


#define SMART_ASSERT_A(x) FEELPP_SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) FEELPP_SMART_ASSERT_OP(x, A)

#define FEELPP_SMART_ASSERT_OP(x, next)         \
    SMART_ASSERT_A.printCurrentValue((x), #x).SMART_ASSERT_ ## next \
    /**/



#endif
