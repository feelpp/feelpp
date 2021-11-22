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

#include <feel/feelcore/journalwatcher.hpp>

namespace Feel
{

JournalWatcher::JournalWatcher( std::string const& category, std::string const& name, bool connect )
    :
    JournalWatcher( std::bind( &JournalWatcher::updateInformationObject, this, std::placeholders::_1 ), category, name, true, connect )
{}

JournalWatcher::JournalWatcher( std::string const& category, bool useDefaultName, bool connect )
    :
    JournalWatcher( std::bind( &JournalWatcher::updateInformationObject, this, std::placeholders::_1 ), category, "", useDefaultName, connect )
{}

JournalWatcher::JournalWatcher( function_update_information_type const& func, std::string const& category, std::string const& name,
                                bool useDefaultNameIfEmpty, bool connect )
    :
    M_journal_is_connected( false ),
    M_category( category ),
    M_name( name ),
    M_function_updateInformationObject( func )
{
    if ( M_name.empty() && useDefaultNameIfEmpty )
        M_name += "object-" + std::to_string( S_counterObjectByCategory[category] );

    slotNew< void ( notify_type & ) >( "journalWatcher", std::bind( &JournalWatcher::journalNotify, this, std::placeholders::_1 ) );

    if ( connect )
        this->journalConnect();

    S_counterObjectByCategory[category]++;
}

JournalWatcher::JournalWatcher( JournalWatcher const& jw )
    :
    M_journal_is_connected( false ),
    M_category( jw.M_category ),
    M_name( jw.M_name ),
    M_function_updateInformationObject( std::bind( &JournalWatcher::updateInformationObject, this, std::placeholders::_1 ) ),
    M_informationObject( jw.M_informationObject )
{
    slotNew< void ( notify_type & ) >( "journalWatcher", std::bind( &JournalWatcher::journalNotify, this, std::placeholders::_1 ) );

    if ( jw.journalIsConnected() )
        this->journalConnect();

    S_counterObjectByCategory[M_category]++;
}

typename nl::json::json_pointer
JournalWatcher::journalSection() const
{
    std::string res;
    if ( !M_category.empty() )
        res += "/" + M_category;
    if ( !M_name.empty() )
        res += "/" + M_name;
    return nl::json::json_pointer( res );
}

void
JournalWatcher::journalConnect()
{
    DVLOG(2) << "[Journal] Connect: " << journalWatcherInstanceName();
    if( JournalManager::journalEnabled() && !this->journalIsConnected())
    {
        JournalManager::signalStaticConnect< void ( notify_type & ) >( "journalManager", *this, "journalWatcher" );
        M_journal_is_connected=true;
        DVLOG(2) << "[Journal] " << journalWatcherInstanceName() << " Connected!";
    }
}

void
JournalWatcher::journalDisconnect()
{
    DVLOG(2) << "[Journal] Disconnect: " << journalWatcherInstanceName();
    if( this->journalIsConnected() )
    {
        JournalManager::signalStaticDisconnect< void ( notify_type & ) >( "journalManager", *this, "journalWatcher" );
        M_journal_is_connected=false;
        DVLOG(2) << "[Journal] " << journalWatcherInstanceName() << " Disconnected!";
    }
}

void
JournalWatcher::journalFinalize()
{
    if ( !JournalManager::journalEnabled() || !this->journalIsConnected() )
        return;

    // store info in the global ptree
    DVLOG(2) << "[JournalManager] Destructor : send info to the journal";
    this->journalNotify( JournalManager::journalData() );

    this->journalDisconnect();
}

void
JournalWatcher::journalNotify( notify_type & journalData )
{
    this->applyUpdateInformationObject();

    if ( M_informationObject.empty() )
        return;

    journalData[this->journalSection()].update( M_informationObject.begin(), M_informationObject.end() );
}


// Init static variable.
std::map<std::string,uint16_type> JournalWatcher::S_counterObjectByCategory = {};

} // Feel namespace.


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
