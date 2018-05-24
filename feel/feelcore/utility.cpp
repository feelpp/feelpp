/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 27 sept. 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <assert.h>

#include <feel/feelcore/utility.hpp>

#ifdef FEELPP_OS_LINUX
#include <termios.h>
#include <unistd.h>
#endif

namespace Feel {

std::string 
readFromFile(std::string const& infile)
{
    std::ifstream instream(infile.c_str());
    if (!instream.is_open()) {
        std::cerr << "Couldn't open file: " << infile << std::endl;
        exit(-1);
    }
    instream.unsetf(std::ios::skipws);      // No white space skipping!
    return std::string(std::istreambuf_iterator<char>(instream.rdbuf()),
                       std::istreambuf_iterator<char>());
}

// Get a unique char from standard input (no return cariage required).
int getOneChar() {
      int c=10;
#ifdef FEELPP_OS_LINUX
      struct termios org_opts, new_opts;
      int res = 0;

      // Store settings.
//      res = tcgetattr(STDIN_FILENO, &org_opts);
      res = tcgetattr(ORTE_IOF_STDIN, &org_opts);
      if( res!=0 )
      {
          std::cerr << "tcgetattr: store failed " << std::strerror(errno) << std::endl;
          std::exit(1);
      }
      // New params.
      memcpy(&new_opts, &org_opts, sizeof(new_opts));
      new_opts.c_lflag &= ~(ICANON | ECHO | ECHOE | ECHOK | ECHONL | ECHOPRT | ECHOKE | ICRNL);
      tcsetattr(STDIN_FILENO, TCSANOW, &new_opts);
      c=getchar();
      // Restore settings.
      res = tcsetattr(STDIN_FILENO, TCSANOW, &org_opts);
      if( res!=0 )
      {
          std::cerr << "tcgetattr: restore failed " << std::strerror(errno) << std::endl;
          std::exit(1);
      }
#endif
      return(c);
}


// Ask to the user to fill a password.
// This function support password corrections.
// Note: Sequential only!
std::string askPassword( const std::string& msg )
{
    const char RETURN=10;
    const char ERASE=8;
    std::string pw;
    pw.reserve(100);
    std::cout << msg;
    char ch = 0;
#ifdef FEELPP_OS_LINUX
    while( (ch = getOneChar()) != RETURN ) // return key
    {
        switch(ch)
        {
            case ERASE:
                if( pw.size() )
                {
                    std::cout << "\b";
                    std::cout << " ";
                    std::cout << "\b" << std::flush;
                    pw.pop_back();
                }
                break;

            case RETURN:
                break;

            default:
                std::cout << '*';
                pw.push_back(ch);
        }
    }
    std::cout << std::endl;
#endif
    return pw;
}

}
