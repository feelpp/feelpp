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
#include <feel/feelcore/testsuite.hpp>

#include <iostream>
#include <cassert>
#include <feel/feelconfig.h>

#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/observer.hpp>

using namespace Feel;


// This test shows how to feed the benchmark system.
// One work here with the simulation info observer "Journal".
// An observer is defined by its managers and its watchers.
//
//  JournalManager  JournalWatcher
//   +-------+        +-------+
//   |object1|        |object2|
//   +-------+        +-------+
//   | Sig1  +-------->slot1  |
//   +-------+   |    +-------+
//               |
//               |    +-------+
//               |    |object3|
//               |    +-------+
//               +---->slot1  |
//                    +-------+
//
// One have to define one or several manager which can ask
// to their watchers (to be defined also) to send their information.
// A Journal observer


// Object1 is defined as the simulation info manager.
class Object1
: virtual public Observer::JournalManager
{
};

// Object2 is an object to be watched.
// Note: A journalNotify has to be defined with the notifications to be send
// to the simulation info manager.
// Note that a "typename" and a "name" key are mandatory
class Object2
: virtual public Observer::JournalWatcher // observe the simulation
{
    public:
        Object2( std::string name = "default" ): M_name( name ) {}

        // Notification for Journal observer.
        const pt::ptree journalNotify() const
        {
            pt::ptree p;
            p.put( M_name + ".a","1" );
            p.put( M_name + ".b","2" );
            p.put( M_name + ".c.d","3" );
            return p;
        }

    private:
        int M_val;
        std::string M_name;
};


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( observers )

BOOST_AUTO_TEST_CASE( journal_basic )
{
    Object1::journalAutoMode(false);
    Object1 m;
    Object2 p1( "p1" );
    Object2 p2( "p2" );

    // Example to show list of signals and list of slots.
    if( Environment::isMasterRank() )
    {
        m.signalStaticShow();
        p1.slotShow(); // non static.
    }
    // Example how to retrieve a signal using signalhandler.
    // (The signals template arguments are required to cast into the proper signal.)
    const auto& sigptr = Object1::signalStatic< pt::ptree(), Observer::JournalMerge >( "journalManager" );
    BOOST_TEST_MESSAGE( "number of connected slot: " << sigptr->num_slots() );

    p1.journalConnect();
    p2.journalConnect();

    BOOST_TEST_MESSAGE( "number of connected slot: " << sigptr->num_slots() );

    // Merge p1, p2 simulation info property tree into one using a call from
    // the manager.

    // Retrieve the merged property tree (no MPI).
    const auto& res = Object1::journalPull(false);
    // This also works
    // m.journalPull();

    // Save into a json file.
    Object1::journalSave();

    auto t = res.get_child( "p1" );
    auto a = t.get<int>("a");
    auto b = t.get<int>("b");
    auto c = t.get_child("c");
    auto d = c.get<int>("d");
    if( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( a << " " << b <<  " " << d );

    CHECK( a == 1 );
    CHECK( b == 2 );
    CHECK( d == 3 );
}

BOOST_AUTO_TEST_CASE( journal_auto )
{
    // Object1 class is a journal manager with environment class (aside).
    Object1::journalAutoMode(true);
    
    // We retrieve a signal of Object1 called JournalManager (used for the feel++ journal
    // system). The signals template arguments are required to cast into the proper signal type.
    const auto& sigptr = Object1::journalSignal();
    int na = sigptr->num_slots();
    BOOST_TEST_MESSAGE( "Number of slots: " << na );

    // The three objects are connected automatically via their constructor.
    std::shared_ptr<Object2> p1 = std::make_shared<Object2>( "p1" );
    std::shared_ptr<Object2> p2 = std::make_shared<Object2>( "p2" );
    std::shared_ptr<Object2> p3 = std::make_shared<Object2>( "p3" );
    if(1)
    {
      std::shared_ptr<Object2> p4 = std::make_shared<Object2>( "p4" );
      BOOST_TEST_MESSAGE( "Number of slots: " << sigptr->num_slots() );
    }
    // p4 is destructed here, thus disconnected from the journal by the destructor.

    p3->journalDisconnect();

    int nb = sigptr->num_slots();

    BOOST_TEST_MESSAGE( "Number of slots: " << sigptr->num_slots() );

    // Only p1 and p2 remains.
    CHECK( nb-na == 2 );

    Object1::journalFinalize();
}

BOOST_AUTO_TEST_SUITE_END()

// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
