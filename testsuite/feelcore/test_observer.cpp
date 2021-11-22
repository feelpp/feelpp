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


// Object2 is an object to be watched.
// Note: A journalNotify has to be defined with the notifications to be send
// to the simulation info manager.
// Note that a "typename" and a "name" key are mandatory
class Object2
: virtual public JournalWatcher // observe the simulation
{
    public:
        Object2( std::string name = "" ) :
            JournalWatcher( "Object2", name ) {}

    void updateInformationObject( nl::json & p ) const override
        {
            p.emplace( "a",1 );
            p.emplace( "b",2 );
            p["/c/d"_json_pointer] =3;
        }
};


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( observers )

BOOST_AUTO_TEST_CASE( journal_basic )
{
    Environment::setJournalAutoMode(false);
    Object2 p1( "p1" );
    Object2 p2( "p2" );

    // Example to show list of signals and list of slots.
    if( Environment::isMasterRank() )
    {
        Environment::signalStaticShow();
        p1.slotShow(); // non static.
    }
    // Example how to retrieve a signal using signalhandler.
    // (The signals template arguments are required to cast into the proper signal.)
    const auto& sigptr = Environment::signalStatic< void ( nl::json & ) >( "journalManager" );
    BOOST_TEST_MESSAGE( "number of connected slot: " << sigptr->num_slots() );

    p1.journalConnect();
    p2.journalConnect();

    BOOST_TEST_MESSAGE( "number of connected slot: " << sigptr->num_slots() );

    // Merge p1, p2 simulation info property tree into one using a call from
    // the manager.

    // Retrieve the merged property tree (no MPI).
    const auto& res = Environment::journalPull();

    // Save into a json file.
    Environment::journalSave();

    auto t = res.at( "/Object2/p1"_json_pointer );
    // std::cout << "t="<<t.dump(1) << std::endl;
    auto a = t.at("a").get<int>();
    auto b = t.at("b").get<int>();
    auto c = t.at("c");
    auto d = c.at("d").get<int>();
    if( Environment::isMasterRank() )
        BOOST_TEST_MESSAGE( a << " " << b <<  " " << d );

    BOOST_CHECK( a == 1 );
    BOOST_CHECK( b == 2 );
    BOOST_CHECK( d == 3 );
}

BOOST_AUTO_TEST_CASE( journal_auto )
{
    Environment::setJournalAutoMode(true);

    // We retrieve a signal of Object1 called JournalManager (used for the feel++ journal
    // system). The signals template arguments are required to cast into the proper signal type.
    const auto& sigptr = Environment::journalSignal();
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
    BOOST_CHECK( nb-na == 2 );

    //Environment::journalFinalize();
}

BOOST_AUTO_TEST_SUITE_END()

// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
