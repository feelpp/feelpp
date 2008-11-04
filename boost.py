# -*- Mode: python coding: latin-1 -*-
###  boost.py --- scons boost checks

#  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
#       Date: 2006-05-25
#
#  Copyright (C) 2006 EPFL
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
import copy

def CheckBoost(context,libs):
    custom_boost_tests = {
        'filesystem' : CheckBoostFilesystem,
        'serialization' : CheckBoostSerialization,
        'program_options' : CheckBoostProgramOptions,
        'date_time' : CheckBoostDateTime,
        'signals' : CheckBoostSignals,
        'test' : CheckBoostTest,
        }
    lastCPPPATH= context.env['CPPPATH']
    lastLIBPATH = context.env['LIBPATH']

    context.env.Append(LIBPATH = context.env['boost_lib_dir'],
                       CPPPATH = context.env['boost_inc_dir'] )

    ret = 1
    for lib in libs:
        ret &= custom_boost_tests[lib](context)

    if not ret :
        context.env.Replace(LIBPATH=lastLIBPATH, CPPPATH=lastCPPPATH)
    context.Result( ret )
    return ret

def CheckBoostDateTime(context):
    context.Message( 'Checking for boost/date_time... ' )
    lastLIBS = context.env['LIBS']

    context.env.Append(LIBS = 'boost_date_time' )
    ret = context.TryLink("""
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <string>
#include "boost/date_time/posix_time/posix_time.hpp"

int main(int argc, char **argv) {
    using namespace boost::posix_time;
    using namespace boost::gregorian;

    ptime t(second_clock::local_time());
    return 0;
}
""",extension='.cpp')
    if not ret:
        context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret

def CheckBoostSignals(context):
    context.Message( 'Checking for boost/signals... ' )
    lastLIBS = context.env['LIBS']

    context.env.Append(LIBS = 'boost_signals' )
    ret = context.TryLink("""
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/signal.hpp>
struct HelloWorld
{
  void operator()() const
  {
    std::cout << "Hello" << std::endl;
  }
};
int main(int argc, char **argv) {
    boost::signal<void()> sig;
    HelloWorld hello;
    sig.connect(hello);
    sig();
}
""",extension='.cpp')
    if not ret:
        context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret

def CheckBoostProgramOptions(context):
    context.Message( 'Checking for boost/program_options... ' )
    lastLIBS = context.env['LIBS']

    context.env.Append(LIBS = 'boost_program_options' )
    ret = context.TryLink("""
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <string>
#include "boost/program_options.hpp"
int main(int argc, char **argv) {
    int Argc = 1;
    char* Argv[] = { "toto" };
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("prefix,p", po::value<std::string>(), "filename prefix")
        ("start,s", po::value<int>(), "start iteration")
        ("finish,f", po::value<int>(), "finish iteration")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(Argc, Argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\\n";
        return 1;
    }
    return 0;
}
""",extension='.cpp')
    if not ret:
        context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret

def CheckBoostFilesystem(context):
    context.Message( 'Checking for boost/filesystem... ' )
    lastLIBS = context.env['LIBS']

    context.env.Append(LIBS = 'boost_filesystem' )
    ret = context.TryLink("""
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include <iostream>
int main(int argc, char **argv) {
boost::filesystem::remove_all( "foobar" );
    boost::filesystem::create_directory( "foobar" );
    boost::filesystem::ofstream file( "foobar/cheeze" );
    file << "tastes good!\\n";
    file.close();
    if ( !boost::filesystem::exists( "foobar/cheeze" ) )
    {
      std::cout << "Something is rotten in foobar\\n";
      return 1;
    }
    boost::filesystem::remove_all( "foobar" );
    return 0;
}
""",extension='.cpp')
    if not ret:
        context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret


def CheckBoostSerialization(context):
    context.Message( 'Checking for boost/serialization... ' )
    lastLIBS = context.env['LIBS']

    context.env.Append(LIBS = 'boost_serialization' )
    ret = context.TryLink("""

#include <iostream>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
class gps_position
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & degrees;
        ar & minutes;
        ar & seconds;
    }
    int degrees;
    int minutes;
    float seconds;
public:
    gps_position(){};
    gps_position(int d, int m, float s) :
        degrees(d), minutes(m), seconds(s)
    {}
};
int main(int argc, char **argv) {
    std::ofstream ofs("filename");
    boost::archive::text_oarchive oa(ofs);
    const gps_position g(35, 59, 24.567f);
    oa << g;
    ofs.close();
    std::ifstream ifs("filename", std::ios::binary);
    boost::archive::text_iarchive ia(ifs);
    gps_position newg;
    ia >> newg;
    ifs.close();
    return 0;
}
""",extension='.cpp')
    if not ret:
        context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret

def CheckBoostTest(context):
    context.Message( 'Checking for boost/test... ' )
    lastLIBS = copy.copy(context.env['LIBS'])
    context.env.Append(LIBS = 'boost_unit_test_framework' )
    ret = context.TryLink("""
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/test/test_tools.hpp>
int init_unit_test_suite(int argc, char **argv) {
   BOOST_CHECK( 1==1);
}
""",extension='.cpp')
    context.env.Replace(LIBS = lastLIBS)
    context.Result( ret )
    return ret
