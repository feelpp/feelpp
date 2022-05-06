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

#ifndef FEELPP_JOURNALMANAGER_HPP
#define FEELPP_JOURNALMANAGER_HPP 1

#include <feel/feelcore/feel.hpp>
#include <feel/feelevent/events.hpp>
#include <feel/feelcore/mongocxx.hpp>
#include <feel/feelcore/json.hpp>

#define FEELPP_DB_JOURNAL_VERSION "0.1.0"

namespace Feel
{

//! JournalManager that manage the  journal system handles all journalWatchers.
//! \see JournalManager JournalWatcher Environment
class JournalManager : public Event::SignalHandler
{
public:
    // Types alias.
    using notify_type = nl::json;

    //! Constructor
    //! @{

    //! Default constructor
    //! This constructor create a new signal waiting for called
    //! 'journalManager'. journalWatcher object connect a specific slot to this
    //! signal.
    //! \see JournalWatcher
    JournalManager();

    //! @}

    //! @{
    //! Destructor
    virtual ~JournalManager() = default;

    //! @}

    //! Setters
    //! @{

    //! Set JSON file name.
    //! \param name the file name.
    static void setJournalFilename( const std::string& name )
        {
            S_journal_filename = name;
        }

    //! Set the journal mode.
    //! \param b automatic mode true or false.
    static void setJournalAutoMode( bool b )
        {
            S_automode = b;
        }

    //! Activate or deactivate the journal.
    static void setJournalEnable( bool m )
        {
            S_enable = m;
        }

    //! @}

    //! Getters
    //! @{

    //! Get the journal filename (no ext).
    static const std::string& journalFilename()
        {
            return S_journal_filename;
        }

    static nl::json & journalData()
        {
            return S_journal_ptree;
        }

    //! Get the current journal mode status.
    //! \return true in automatic mode.
    static const bool journalAutoMode()
        {
            return S_automode;
        }

    //! Get the journal status.
    //! \return true if journal is enabled.
    static const bool journalEnabled()
        {
            return S_enable;
        }

    //! Get the current checkpoint counter.
    //! \return the counter value.
    static const uint32_type
    journalCurrentCheckpoint()
        {
            return S_journal_checkpoint;
        }

    //! Get the journal underlying static signal.
    //! \return
    static decltype(auto)
        journalSignal()
        {
            return signalStatic< void ( notify_type & ) >( "journalManager" );
        }

    //! @}

    //! Journal gather/save.
    //! @{

    //! Fetch and merge notifications coming from all observed objects.
    //! \return The global journal property tree.
    //! \see JournalWatcher
    static notify_type const&
    journalPull()
        {
            auto thesig = signalStatic< void ( notify_type & ) >( "journalManager" );
            (*thesig)( S_journal_ptree );
            return S_journal_ptree;
        }

    //! Save the global property tree into a json file.
    static void
    journalSave( std::string const& filename = "" );

    //! Generate a checkpoint
    //! \param save enable intermediate save.
    //! \param filename name of the intermediate file to save in.
    //! The filename "journal" is overwritten by default. Use
    //! journalCurrentCheckpoint to add an index to the files.
    //!
    //! \see journalCurrentCheckpoint
    static void
    journalCheckpoint( bool save = true, const std::string& filename="" );

    //! Create the journal.
    static void
    journalFinalize()
        {
            journalCheckpoint();
            setJournalEnable( false );
        }

    //! @}

    //! Setters
    //! @{

    //! Set the mongodb database name.
    static void journalSetDBName( const std::string& dbname = "feelpp")
        {
            S_journal_db_config.name = dbname;
        }

    //! Set the mongodb host.
    static void journalSetDBHost( const std::string&  host = "localhost")
        {
            S_journal_db_config.host = host;
        }

    //! Set the mongodb user name.
    static void journalSetDBUsername( const std::string& user = "")
        {
            S_journal_db_config.user = user;
        }

    //! Set the mongodb user password.
    static void journalSetDBPassword( const std::string& password = "")
        {
            S_journal_db_config.password = password;
        }

    //! Set the mongodb port.
    static void journalSetDBPort( const std::string& port = "27017")
        {
            S_journal_db_config.port = port;
        }

    //! Set the collection used to authenticate.
    static void journalSetDBAuthsrc( const std::string& authsrc = "admin")
        {
            S_journal_db_config.authsrc = authsrc;
        }

    //! Set mongodb collection to use for the journal.
    static void journalSetDBCollection( const std::string& dbname = "journal")
        {
            S_journal_db_config.name = dbname;
        }

    //! Set mongodb collection to use for the journal.
    static void journalDBConfig( const MongoConfig& mc )
        {
            S_journal_db_config = mc;
        }
    //! @}

private:
    //! Private methods
    //! @{

    //! Save the global property tree into a json file.
    static void
    journalJSONSave( const std::string& filename, nl::json const& pt );

    //! Save the json in a mongodb database.
    //! This function read a json file in a bson format. Then this bson entry
    //! is send in the mongodb database. The database has to be configured
    //! beforehand.
    static void
    journalDBSave( const std::string& filename );

    //! @}

private:
    //! JSON filename.
    static std::string S_journal_filename;
    //! Global property tree.
    static nl::json S_journal_ptree;

    //! MongoDB specific attribute. These attributes are set via options.
    //! @{
    static MongoConfig S_journal_db_config;
    //! @}

    //! checkpoint number
    static uint32_type S_journal_checkpoint;

    //! journal enable or disable
    static bool S_enable;

    //! connect automaticaly an observer to the journal if true
    static bool S_automode;
};

} // Feel namespace.

#endif // FEELPP_JOURNALMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
