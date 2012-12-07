/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the FeelV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-06-22

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file logs.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-06-22
 */
#include <boost/log/functions.hpp>
#include <boost/limits.hpp>
#include <iostream>
#include <fstream>
#include <time.h>

#include "logs.hpp"

BOOST_DEFINE_LOG( app, "app" )
BOOST_DEFINE_LOG( dbg, "app.dbg" )
BOOST_DEFINE_LOG( err, "app.err" )
BOOST_DEFINE_LOG( warn, "app.warn" )
BOOST_DEFINE_LOG( info, "info" )


// *** Appenders
void write_to_cout( const std::string &, const std::string &msg )
{
    std::cout << msg;
}


// *** Modifiers
void prefix_time( const std::string &, std::string & msg )
{
    char time_buff[ 20];
    time_t t = time( 0 );
    tm details = *localtime( &t );
    sprintf( time_buff, "%02d:%02d:%02d  ", details.tm_hour, details.tm_min, details.tm_sec );
    msg = time_buff + msg;
}

void init_logs()
{
    using namespace boost::logging;

    // Modifiers for all:
    // [type_of_message] original_message append_enter_if_needed
    //add_modifier("*", &prefix_time, INT_MAX );
    add_modifier( "*", &append_enter );
    // Modifiers for app and its ascendants
    // <time> [type_of_message] original_message append_enter_if_needed
    add_modifier( "*", &prepend_prefix );
    // Modifiers for "app" only
    // <time> [Thread ID] [type_of_message] original_message append_enter_if_needed
    //add_modifier("app", &prepend_thread_id, 0);

    // Log Functions
    // all messages are written to cout
    //add_appender("*", write_to_cout);
    // "app*" messages are written to file as well
    add_appender( "app*", write_to_file( "report.txt" ) );
    // 'app' only and dbg messages are written to Output Debug Window as well

    flush_log_cache();
}



