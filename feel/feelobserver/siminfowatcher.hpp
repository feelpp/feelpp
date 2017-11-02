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

#ifndef FEELPP_SIMINFOWATCHER_HPP
#define FEELPP_SIMINFOWATCHER_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/siminfomanager.hpp>
#include <feel/feelobserver/functors/siminfomerge.hpp>

// All Feel++ observers are defined here.

namespace Feel
{
namespace Observer
{

namespace pt =  boost::property_tree;

// Watcher
class SimInfoWatcher
    : public Event::SlotHandler
{
public:
    // Type alias.
    using notify_type = pt::ptree;

    //! Watch child properties and notify the manager.
    virtual const pt::ptree simInfoNotify() const = 0;

    //! Constructors
    //! @{

    //! Default constructor.
    //!
    //! When this constructor is called by a child class, the new child observer
    //! is connected a new simulation info slot is created and can be connected
    //! to the simulation info manager
    //!
    //! \see SimInfoManager
    SimInfoWatcher()
    {
        slotNew< notify_type () >( "simInfoWatcher", std::bind( &SimInfoWatcher::simInfoNotify, this) );
        if( M_siminfo_autoconnect )
            simInfoConnect();
    }

    //! @}

    //! Setters
    //! @{
    void simInfoAutoconnect( bool b = true )
    {
        M_siminfo_autoconnect = b;
    }

    //! @}

    //! Misc
    //! @{

    //! Connect the derived object to the simulation info manager.
    void simInfoConnect()
    {
        SimInfoManager::signalStaticConnect< notify_type (), SimInfoMerge >( "simInfoManager", *this, "simInfoWatcher" );
    }

    //! Disconnect the derived object from the simulation info manager.
    void simInfoDisConnect()
    {
        SimInfoManager::signalStaticDisconnect< notify_type (), SimInfoMerge >( "simInfoManager", *this, "simInfoWatcher" );
    }
private:
    static bool M_siminfo_autoconnect;
    //! @}
};


} // Observer namespace
} // Feel namespace.

#endif // FEELPP_SIMINFOWATCHER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
