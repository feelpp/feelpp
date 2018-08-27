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

    //! Constructors
    //! @{

    //! Default constructor.
    //! 
    //! \param force Force the connection to the journal.
    //!
    //! When this constructor is called by a child class, the new child observer
    //! is automatically connected (in "auto" mode) to the journal system.
    //! A slot is created and and is connected to a JournalManager.
    //!
    //! \see JournalManager
    JournalWatcher( bool force = false )
        : M_journal_is_connected( false )
    {
        slotNew< notify_type () >( "journalWatcher", std::bind( &JournalWatcher::journalNotify, this ) );
        if( JournalManager::journalEnabled()
            and JournalManager::journalAutoMode() )
        {
           journalConnect( force );
        }
    }

    //! Default destructor.
    //! The (inherited) object is always disconnected from the journal during the
    //! destruction.
    virtual ~JournalWatcher()
    {
        // store info in the global ptree (No MPI comm!).
        if( JournalManager::journalAutoPullAtDelete() )
            JournalManager::journalLocalPull();
        journalDisconnect();
    }

    //! @}
    
    //! Getters
    //! @{
    //! Check if the object is connected to the journal.
    //! \return true i
    const bool journalIsConnected()
    {
        return M_journal_is_connected;
    }

    //! @}

    //! Misc
    //! @{

    //! Connect the derived object to the simulation info manager.
    void journalConnect( bool force = false )
    {
        if( ( JournalManager::journalEnabled() and not journalIsConnected() )
            or force )
        {
            JournalManager::signalStaticConnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
            M_journal_is_connected=true;
        }
    }

    //! Disconnect the derived object from the simulation info manager.
    //! The disconnection is safe.
    void journalDisconnect( bool force = false )
    {
        if( journalIsConnected()
            or force )
        {
            JournalManager::signalStaticDisconnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
            M_journal_is_connected=false;
        }
    }

private:
    //! Private Methods
    //! @{

    //! Watch child properties and notify the manager.
    //! Note: Only this class can call journalNotify!
    virtual const pt::ptree journalNotify() const
    {
        pt::ptree p;
        // return empty tree by default;
        return p;
    }

    //! @}
private:
    bool M_journal_is_connected;
};


} // Observer namespace.
} // Feel namespace.

#include <feel/feelobserver/detail/journalfeed.hpp>

#endif // FEELPP_JOURNALWATCHER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
