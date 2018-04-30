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

#ifndef FEELPP_JOURNALWATCHER_HPP
#define FEELPP_JOURNALWATCHER_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/journalmanager.hpp>
#include <feel/feelobserver/functors/journalmerge.hpp>

// All Feel++ observers are defined here.

namespace Feel
{
namespace Observer
{

namespace pt =  boost::property_tree;

// Watcher
class JournalWatcher
    : public Event::SlotHandler
{
public:
    // Type alias.
    using notify_type = pt::ptree;

    //! Watch child properties and notify the manager.
    virtual const pt::ptree journalNotify() const = 0;

    //! Constructors
    //! @{

    //! Default constructor.
    //!
    //! When this constructor is called by a child class, the new child observer
    //! is connected a new simulation info slot is created and can be connected
    //! to the simulation info manager
    //!
    //! \see JournalManager
    JournalWatcher()
    {
        slotNew< notify_type () >( "journalWatcher", std::bind( &JournalWatcher::journalNotify, this) );
        if( M_journal_autoconnect )
            journalConnect();
    }

    //! Default destructor.
    virtual ~JournalWatcher()
    {
        if( M_journal_autoconnect )
            journalDisconnect();
    }

    //! @}

    //! Setters
    //! @{
    void journalAutoconnect( bool b = true )
    {
        M_journal_autoconnect = b;
    }

    //! @}

    //! Misc
    //! @{

    //! Connect the derived object to the simulation info manager.
    void journalConnect()
    {
        JournalManager::signalStaticConnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
    }

    //! Disconnect the derived object from the simulation info manager.
    void journalDisconnect()
    {
        JournalManager::signalStaticDisconnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
    }
private:
    static bool M_journal_autoconnect;
    //! @}
};


} // Observer namespace.
} // Feel namespace.

#include <feel/feelobserver/detail/journalfeed.hpp>

#endif // FEELPP_JOURNALWATCHER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
