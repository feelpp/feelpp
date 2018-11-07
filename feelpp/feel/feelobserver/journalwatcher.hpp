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

    using function_update_information_type = std::function<void ( pt::ptree & )>;
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
    explicit JournalWatcher( std::string const& category = "", std::string const& name = "", bool connect = JournalManager::journalAutoMode(), bool force = false )
        :
        JournalWatcher( std::bind( &JournalWatcher::updateInformationObject, this, std::placeholders::_1 ), category, name, connect, force )
        {}

    explicit JournalWatcher( function_update_information_type const& func, std::string const& category = "", std::string const& name = "", bool connect = JournalManager::journalAutoMode(), bool force = false )
        :
        M_journal_is_connected( false ),
        M_category( category ),
        M_name( name ),
        M_number( S_call_counter ),
        M_function_updateInformationObject( func )
    {
        if ( M_name.empty() )
            M_name += "object-" + std::to_string( M_number );

        slotNew< notify_type () >( "journalWatcher", std::bind( &JournalWatcher::journalNotify, this ) );

        if( connect )
        {
            this->journalConnect( force );
        }

        S_call_counter++;
    }

    //! Default destructor.
    //! The (inherited) object is always disconnected from the journal during the
    //! destruction.
    virtual ~JournalWatcher()
    {
        this->journalFinalize();
    }
    //! @}

    //! Getters
    //! @{
    //! Check if the object is connected to the journal.
    //! \return true
    bool journalIsConnected() const
    {
        return M_journal_is_connected;
    }

    //! Return the full instance name for this watched object.
    //! The name is composed with the base name and a suffix index corresponding
    //! to the nth call.
    //! \param isauto If true, the instance name is suffixed by the call number
    std::string const& journalWatcherInstanceName() const
    {
        return M_name;
    }

    pt::ptree const& informationObject() const { return M_informationObject; }
    //! @}

    //! Setters
    //! @{

    //! Set journal watcher a name for the instance. If the journal is set
    //! in autommatic mode, then a number will be added to this nam
    //! \param automode If true, a index suffix is added to the name
    void setJournalWatcherName( std::string const& name, bool automode = false )
    {
        M_name = name;
        if ( automode )
            M_name += "-" + std::to_string( M_number );
    }

    //! @}

    //! Misc
    //! @{

    //! Connect the derived object to the simulation info manager
    void journalConnect( bool force = false )
    {
        DVLOG(2) << "[Journal] Connect: " << journalWatcherInstanceName();
        if( ( JournalManager::journalEnabled() and not journalIsConnected() )
            or force )
        {
            JournalManager::signalStaticConnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
            M_journal_is_connected=true;
            DVLOG(2) << "[Journal] " << journalWatcherInstanceName() << " Connected!";
        }
    }

    //! Disconnect the derived object from the simulation info manager.
    //! The disconnection is safe.
    void journalDisconnect( bool force = false )
    {
        DVLOG(2) << "[Journal] Disconnect: " << journalWatcherInstanceName();
        if( journalIsConnected()
            or force )
        {
            JournalManager::signalStaticDisconnect< notify_type (), JournalMerge >( "journalManager", *this, "journalWatcher" );
            M_journal_is_connected=false;
            DVLOG(2) << "[Journal] " << journalWatcherInstanceName() << " Disconnected!";
        }
    }

    //! set information object
    void setInformationObject( pt::ptree const& informationObject )
        {
            M_informationObject = informationObject;
        }
    //! put information object
    void putInformationObject( pt::ptree const& informationObject )
        {
            for( const auto& p : informationObject )
                M_informationObject.put_child( p.first, p.second );
        }
    //! add information object
    void addInformationObject( pt::ptree const& informationObject )
        {
            for( const auto& p : informationObject )
                M_informationObject.add_child( p.first, p.second );
        }

    //! finalize journal publication
    void journalFinalize()
        {
            if ( !JournalManager::journalEnabled() || !this->journalIsConnected() )
                return;
            // store info in the global ptree (No MPI comm!).
            if ( JournalManager::journalAutoPullAtDelete() )
            {
                DVLOG(2) << "[JournalManager] Destructor call. Nofification send (signal exec)!";
                JournalManager::journalLocalPull( this->journalNotify() );
            }
            this->journalDisconnect();
        }

protected:
    //! Protected Methods
    //! @{

    void applyUpdateInformationObject()
        {
            M_function_updateInformationObject( M_informationObject );
        }
    virtual void updateInformationObject( pt::ptree & p ) {}

    //! Watch child properties and notify the manager.
    //! Note: Only this class can call journalNotify!
    pt::ptree const journalNotify()
    {
        this->applyUpdateInformationObject();

        std::string prefix;
        if ( !M_category.empty() && !M_name.empty() )
            prefix = M_category + "." + M_name;
        else if ( !M_category.empty() )
            prefix = M_category;
        else
            prefix = M_name;

        if ( prefix.empty() || M_informationObject.empty() )
            return M_informationObject;
        else
        {
            pt::ptree pt;
            pt.put_child( prefix, M_informationObject );
            return pt;
        }
    }

    //! @}

private:
    //! Private attributes
    //! @{

    bool M_journal_is_connected;

    //! category of object
    std::string M_category;

    //! Unique instance name for the watched object.
    std::string M_name;

    //! Counter for instance call of this object.
    static uint16_t S_call_counter;
    //! Unique instance number for the watched object.
    uint16_t M_number;

    function_update_information_type M_function_updateInformationObject;
    pt::ptree M_informationObject;
    //! @}
};


} // Observer namespace.
} // Feel namespace.

#include <feel/feelobserver/detail/journalfeed.hpp>

#endif // FEELPP_JOURNALWATCHER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
