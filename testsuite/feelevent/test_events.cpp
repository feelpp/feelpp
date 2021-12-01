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

#define BOOST_TEST_MODULE test_events
#include <feel/feelcore/testsuite.hpp>

#include <iostream>
#include <cassert>
#include <feel/feelconfig.h>
#include <feel/feelevent/events.hpp>

using namespace Feel;

// Object which is signalable and slotable.
class MyObject
: virtual public Event::SignalHandler,
  virtual public Event::SlotHandler
{
    MyObject& operator=( MyObject && ) = delete;
};

class Slot1
{
public:
    using result_type = int;
    int operator()( int i ){return i;}
};

class Slot2
{
public:
    using result_type = int;
    int operator()(){return 0;}
};

// Signal event.
struct AddEvent{
    using result_type=int;
    template<typename InputIterator>
        result_type operator()( InputIterator first, InputIterator last ) const
        {
            int sum=0;
            while(first != last) {
                sum += *first;
                ++first;
            }
            return sum;
        }
};

int slot_test1() { return 1; }
int slot_test2() { return 2; }
int slot_test3() { return 3; }


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( events )

BOOST_AUTO_TEST_CASE( create_signal )
{
    MyObject obj1, obj2;
    // MyObject instance own their signals.
    const auto& sig1 = obj1.signalNew< int() >( "sig1" );
    const auto& sig2 = obj1.signalNew< int() >( "sig1" ); // is sig1
    const auto& sig3 = obj2.signalNew< int() >( "sig3" );
    const auto& sig4 = obj2.signalNew< int() >( "sig4" );

    // Static signal is shared between MyObject instances.
    const auto& sig5 = obj2.signalStaticNew< int() >( "sig5" );
    const auto& sig6 = MyObject::signalStaticNew< double() >( "sig6" );
    const auto& sig7 = MyObject::signalStaticNew< int() >( "sig5" ); // is sig5

    std::cout << "* OBJECT1: " << std::endl;
    obj1.signalShow();
    obj1.signalStaticShow();

    std::cout << "* OBJECT2: " << std::endl;
    obj2.signalShow();
    obj2.signalStaticShow();

    CHECK( obj1.signals().size() != obj2.signals().size() );
    CHECK( obj1.signalsStatic().size() == obj2.signalsStatic().size() );
    CHECK( obj1.signalsStatic().size() == MyObject::signalsStatic().size() );
}

BOOST_AUTO_TEST_CASE( create_slot )
{
    MyObject obj1, obj2;
    const auto& slot1 = obj1.slotNew< int() >( "slot1", &slot_test1 );
    const auto& slot2 = obj2.slotNew< int() >( "slot1", Slot2() ); // is slot1.
    const auto& slot3 = obj2.slotNew< int() >( "slot3", [](){ return 3;} );
    const auto& slot4 = obj2.slotNew< int() >( "slot4", [](){ return 3;} );

    const auto& slot5 = obj2.slotStaticNew< int() >( "slot5", [](){ return 3;} );
    const auto& slot6 = obj2.slotStaticNew< int() >( "slot5", [](){ return 3;} ); // is slot5
    const auto& slot7 = MyObject::slotStaticNew< int() >( "slot7", [](){ return 3;} );

    std::cout << "* OBJECT1: " << std::endl;
    obj1.slotShow();
    obj1.slotStaticShow();

    std::cout << "* OBJECT2: " << std::endl;
    obj2.slotShow();
    obj2.slotStaticShow();

    CHECK( obj1.slots().size() != obj2.slots().size() );
    CHECK( obj1.slotsStatic().size() == obj2.slotsStatic().size() );
    CHECK( obj1.slotsStatic().size() == MyObject::slotsStatic().size() );
}

BOOST_AUTO_TEST_CASE( access_sig_slot )
{
    MyObject obj;
    // Cast to the correct signal/slot.
    const auto& sig1 = obj.signalNew< int() >( "sig1" );
    const auto& slot1 = obj.slotNew< int() >( "slot1", &slot_test1 );

    decltype(auto) g = obj.signal< int() >( "sig1" );
    const auto& t = obj.slot< int() >( "slot1" );
}

BOOST_AUTO_TEST_CASE( connect_sig_to_slot )
{
    // Connect obj1 to obj2.
    MyObject obj1, obj2;
    // Signals
    obj1.signalNew< int(), AddEvent >( "sig1" );
    obj1.signalStaticNew< int(), AddEvent >( "sig2" );
    // Slots
    obj2.slotNew< int() >( "slot1", &slot_test1 );
    obj2.slotStaticNew< int() >( "slot2", &slot_test2 );
    
    // TODO Fix this case !!
//    obj1.signalConnect< int(), AddEvent >( "sig1", obj2, "slot1" );
    MyObject::signalStaticConnect< int(), AddEvent >( "sig2", obj2, "slot2", Event::SLOT_STATIC );

    MyObject::signalStaticConnect< int(), AddEvent >( "sig2", obj2, "slot1" );
}

BOOST_AUTO_TEST_SUITE_END()

// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
