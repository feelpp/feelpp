// This file is part of the Feel library
//
// Author(s): Feel++ Contortium
//      Date: 2017-07-10
//
// @copyright (C) 2017 University of Strasbourg
// @copyright (C) 2012-2017 Feel++ Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#define BOOST_TEST_MODULE test_observers
#include <testsuite.hpp>

#include <iostream>
#include <cassert>
#include <feel/feelconfig.h>

#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/observer.hpp>

using namespace Feel;


// Object used as simulation info manager.
// This class will owns the signal where each watched object will send
// notifications.
class ManagerTest
: virtual public Observer::SimInfoManager
{
};

// Object observed. A simInfoNotify has to be defined to send notifications
// to the simulation info manager.
// Note that a "typename" and a "name" key has to be set in the property tree
class ProbeTest1
:   virtual public Observer::SimInfoWatcher // observe the simulation
{
    public:
        ProbeTest1( std::string name = "default" ): M_name( name ) {}

        // Notification for SimInfo observer.
        const pt::ptree simInfoNotify() const
        {
            pt::ptree p;
            p.put("typename","ProbeTest1"); // required
            p.put("name", M_name); // required
            p.put("a","1");
            p.put("b","2");
            p.put("c.d","3");
            return p;
        }

    private:
        int M_val;
        std::string M_name;
};


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( observers )

BOOST_AUTO_TEST_CASE( siminfo_basic )
{

    ManagerTest m;
    ProbeTest1 p1( "p1" );
    ProbeTest1 p2( "p2" );

    // Example to show list of signals and list of slots.
    if( Environment::isMasterRank() )
    {
        m.signalStaticShow();
        p1.slotShow(); // non static.
    }
    // Example how to retrieve a signal using signalhandler.
    // (The signals template arguments are required to cast into the proper signal.)
    const auto& sigptr = ManagerTest::signalStatic< pt::ptree(), Observer::SimInfoMerge >( "simInfoManager" );
    std::cout << "number of connected slot: " << sigptr->num_slots() << std::endl;

    p1.simInfoConnect();
    p2.simInfoConnect();

    std::cout << "number of connected slot: " << sigptr->num_slots() << std::endl;

    // Merge p1, p2 simulation info property tree into one using a call from
    // the manager.

    const auto& res = ManagerTest::simInfoPull();
    // This also works "m.simInfoPull();"

    // Save into a json file.
    ManagerTest::simInfoSave();

    auto t = res.get_child( "ProbeTest1.p1" );
    auto a = t.get<int>("a");
    auto b = t.get<int>("b");
    auto c = t.get_child("c");
    auto d = c.get<int>("d");
    if( Environment::isMasterRank() )
        std::cout << a << " " << b <<  " " << d << std::endl;

    CHECK( a == 1 );
    CHECK( b == 2 );
    CHECK( d == 3 );
}

BOOST_AUTO_TEST_SUITE_END()

// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
