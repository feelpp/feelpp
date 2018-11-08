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

#include <feel/feelevent/events.hpp>
#include <feel/feelcore/journalmanager.hpp>

namespace Feel
{

// Watcher
class JournalWatcher : public Event::SlotHandler
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
    explicit JournalWatcher( std::string const& category = "", std::string const& name = "", bool connect = JournalManager::journalAutoMode() );

    explicit JournalWatcher( std::string const& category, bool useDefaultName, bool connect = JournalManager::journalAutoMode() );

    explicit JournalWatcher( function_update_information_type const& func, std::string const& category = "", std::string const& name = "",
                             bool useDefaultNameIfEmpty = true, bool connect = JournalManager::journalAutoMode() );

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

    //! @}

    //! Misc
    //! @{

    //! Connect the derived object to the simulation info manager
    void journalConnect();

    //! Disconnect the derived object from the simulation info manager.
    //! The disconnection is safe.
    void journalDisconnect();

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
    void journalFinalize();

protected:
    //! Protected Methods
    //! @{

    void applyUpdateInformationObject()
        {
            M_function_updateInformationObject( M_informationObject );
        }
    virtual void updateInformationObject( pt::ptree & p ) {}

private :
    //! Watch child properties and notify the manager.
    //! Note: Only this class can call journalNotify!
    void journalNotify( notify_type & journalData );

    //! @}

private:
    //! Private attributes
    //! @{

    bool M_journal_is_connected;

    //! category of object
    std::string M_category;

    //! Unique instance name for the watched object.
    std::string M_name;

    //! Counter of object by category.
    static std::map<std::string,uint16_type> S_counterObjectByCategory;

    //! function which update the information object
    function_update_information_type M_function_updateInformationObject;

    //! information object data structure
    pt::ptree M_informationObject;
    //! @}
};


} // Feel namespace.

#endif // FEELPP_JOURNALWATCHER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
