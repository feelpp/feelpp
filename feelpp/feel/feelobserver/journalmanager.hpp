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
#include <feel/feelobserver/functors/journalmerge.hpp>
#include <feel/feelcore/mongocxx.hpp>

#define FEELPP_DB_JOURNAL_VERSION "0.1.0"

namespace Feel
{

namespace Observer
{

namespace pt =  boost::property_tree;

//! JournalManager handles all journalWatchers.
//! The purpose of this class is to be inherited by class that manage the
//! journal system
//!
//! \note Journal manager class should inherit from this class
//!
//! \remark Environment is the favored journal manager (child class). However,
//! environment (child) is initialized after JournalManager (mother).
//! Thus a  template parameter is defined in order to use environment static
//! members.
//!
//! \see JournalManager JournalWatcher Environment
class JournalManager : public Event::SignalHandler
{

public:
        // Types alias.
        using notify_type = pt::ptree;

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
            Options::automode = b;
        }

        //! Set the journal delete mode.
        //! The journal manager allows its watchers to execute the signal and
        //! send their information before they are deleted.
        //! \param m pull at delete mode true or false
        //! \return true if journal if delete pull is allowed.
        static void setJournalAutoPullAtDelete( bool m )
        {
            Options::allow_destructor_call = m;
        }

        //! Activate or deactivate the journal.
        static void setJournalEnable( bool m )
        {
            Options::enable = m;
        }

        //! @}

        //! Getters
        //! @{

        //! Get the journal filename (no ext).
        static const std::string& journalFilename()
        {
            return S_journal_filename;
        }

        //! Get the current journal mode status.
        //! \return true in automatic mode.
        static const bool journalAutoMode()
        {
            return Options::automode;
        }

        //! Get the journal status.
        //! \return true if journal is enabled.
        static const bool journalEnabled()
        {
            return Options::enable;
        }

        //! Get the journal pull at delete status.
        //! The journal manager allows its watchers to execute the signal.
        //! \return true if journal if delete pull is allowed.
        static const bool journalAutoPullAtDelete()
        {
            return Options::allow_destructor_call;
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
            return signalStatic< notify_type (), JournalMerge >( "journalManager" );
        }

        //! @}

        //! Journal gather/save.
        //! @{

        //! Merge global property_tree form a ptree
        //! \return The global journal property tree.
        static
        notify_type const&
        journalPull( pt::ptree const& pt )
        {
            ptMerge( S_journal_ptree, pt );
            return S_journal_ptree;
        }

        //! Fetch and merge notifications coming from all observed objects.
        //! \return The global journal property tree.
        //! \see JournalWatcher
        static notify_type const&
        journalPull( std::string const& signame = "journalManager" )
        {
            const pt::ptree& pt_merged = signalStaticPull< notify_type (), JournalMerge >( signame );
            ptMerge( S_journal_ptree, pt_merged );
            return S_journal_ptree;
        }

        //! Save the global property tree into a json file.
        static void
        journalSave( std::string filename = "" );

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
        journalJSONSave( const std::string& filename, pt::ptree const& pt );

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
        static pt::ptree S_journal_ptree;

        //! MongoDB specific attribute. These attributes are set via options.
        //! @{
        static MongoConfig S_journal_db_config;
        //! @}

        //! checkpoint number
        static uint32_type S_journal_checkpoint;

public:
        // Options for the journal.
        // These default options are setup via the Feel++ commandline options.
        // See environment.
        class Options
        {
            public:
                //! Journal automatic mode enable or disable.
                static bool enable;
                //! Journal automatic mode enable or disable.
                static bool automode;
                //! Journal automatic mode enable or disable.
                static bool allow_destructor_call;
        };
};

} // Observer namespace.
} // Feel namespace.

#endif // FEELPP_JOURNALMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
