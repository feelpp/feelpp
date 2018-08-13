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

#include <chrono>
#include <ctime>
#include <iomanip>
#include <boost/mpi.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/uuid/uuid.hpp>    
#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/functors/journalmerge.hpp>
#include <feel/feelcore/mongocxx.hpp>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#endif

#define FEELPP_DB_JOURNAL_VERSION "0.1.0"

namespace Feel
{
// Forward declaration.
class Environment;

namespace Observer
{

namespace pt =  boost::property_tree;
namespace mpi = boost::mpi;
namespace uuids =  boost::uuids;

//! JournalManagerBase handles all journalWatchers.
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
template< typename Env = Feel::Environment >
class JournalManagerBase
   : public Event::SignalHandler
{

public:
        // Types alias.
        using notify_type = pt::ptree;
        using env_type = Env;

        //! Constructor
        //! @{

        //! Default constructor
        //! This constructor create a new signal waiting for called
        //! 'journalManager'. journalWatcher object connect a specific slot to this
        //! signal.
        //! \see JournalWatcher
        JournalManagerBase()
        {
#if defined(FEELPP_HAS_MONGOCXX )
            // Create mongo unique instance!
            Feel::MongoCxx::instance();
#endif 
            //M_journal_uuid = Env::randomUUID();
            std::time_t t = std::time(nullptr);
            const std::string tag = "journal.time";
            S_journal_ptree.put( "journal.version", FEELPP_DB_JOURNAL_VERSION );
            S_journal_ptree.put( tag + ".time_t", t );
            S_journal_ptree.put( tag + ".gm", std::put_time(std::gmtime(&t), "%c %Z") );
            S_journal_ptree.put( tag + ".local", std::put_time(std::localtime(&t), "%c %Z") );
            // Create a signal for simulation info.
            signalStaticNew< notify_type (), JournalMerge >( "journalManager" );
        }

        //! @}

        //! @{
        //! Destructor
        virtual ~JournalManagerBase() = default;

        //! @}

        //! Setters
        //! @{

        //! Set JSON file name.
        //! \param name the file name.
        static void journalFilename( const std::string& name )
        {
            S_journal_filename = name;
        }
        
        //! Set the journal mode.
        //! \param m automatic mode true or false.
        static void journalAutoMode( bool m )
        {
            S_journal_auto = m;
        }

        //! Activate or deactivate the journal.
        static void journalEnable( bool m )
        {
            S_journal_enabled = m;
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
            return S_journal_auto;
        }

        //! Get the journal status.
        //! \return true if journal is enabled.
        static const bool journalEnabled()
        {
            return S_journal_enabled;
        }

        //! Get the current checkpoint counter.
        //! \return the counter value.
        static const uint32_t
        journalCurrentCheckpoint()
        {
            return S_journal_checkpoint;        
        }

        //! @}

        //! Journal gather/save.
        //! @{

        //! Fetch and merge notifications coming from all observed objects into
        //! a global property tree.
        //! \see JournalWatcher
        static const notify_type
        journalPull( bool parallel = true )
        {
            const pt::ptree& pt_merged = signalStaticPull< notify_type (), JournalMerge >( "journalManager" );
            ptMerge( S_journal_ptree, pt_merged );

            if( Env::isParallel() and parallel )
            {
                pt::ptree p;
                mpi::reduce( Env::worldComm(), S_journal_ptree, p, ptMerge, Env::masterRank() );
                S_journal_ptree = p;
            }

            return S_journal_ptree;
        }

        //! Save the global property tree into a json file.
        static void
        journalSave( const std::string& filename = "" )
        {
            if( not filename.empty() )
                journalFilename( filename );

            if( Env::isMasterRank() )
            {
                journalJSONSave( S_journal_filename );
                journalDBSave( S_journal_filename );
            }
        }
        
        //! Generate a checkpoint
        //! \param save enable intermediate save.
        //! \param filename name of the intermediate file to save in.
        //! \param parallel gather information from every proc.
        //! The filename "journal" is overwritten by default. Use
        //! journalCurrentCheckpoint to add an index to the files. By default,
        //! a property tree is generated by gathering information from every process
        //! in parallel. This option might have a non neglectible cost.
        //!
        //! Remark: Performing parallel operation during intermediate checkpoints is
        //! not recommended.
        //! 
        //! \see journalCurrentCheckpoint
        static const void
        journalCheckpoint( bool parallel = true,
                           bool save = false,
                           const std::string& filename="" )
        {
            if( S_journal_enabled )
            {
                S_journal_checkpoint++;

                journalPull( parallel );

                std::time_t t = std::time(nullptr);
                const std::string tag = "journal.checkpoints.checkpoint-" + std::to_string( S_journal_checkpoint ) + ".time";
                S_journal_ptree.put( "journal.checkpoints.number", S_journal_checkpoint);
                S_journal_ptree.put( tag + ".time_t", t );
                S_journal_ptree.put( tag + ".gm", std::put_time(std::gmtime(&t), "%c %Z") );
                S_journal_ptree.put( tag + ".local", std::put_time(std::localtime(&t), "%c %Z") );

                if( save ) 
                    journalSave( filename );
            }
        }

        //! Create the journal.
        static const void
        journalFinalize()
        {
            journalCheckpoint( true, true );
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
        journalJSONSave( const std::string& filename = "" )
        {
            if( S_journal_enabled )
            {
                if( not filename.empty() )
                    S_journal_filename = filename;

                //std::cout << "[Observer: Journal] generate report (JSON).";
                if( not S_journal_ptree.empty() )
                    write_json( S_journal_filename + ".json", S_journal_ptree );
            }
        }

        //! Save the json in a mongodb database. 
        //! This function read a json file in a bson format. Then this bson entry 
        //! is send in the mongodb database. The database has to be configured
        //! beforehand.
        static void
        journalDBSave( const std::string& filename = "" )
        {
            if( S_journal_enabled )
            {
                if( not filename.empty() )
                    S_journal_filename = filename;
                auto vm =  Env::vm();
                bool enable = vm["journal.database"].template as<bool>();
                if( enable )
                {
#if defined(FEELPP_HAS_MONGOCXX )
                    using bsoncxx::builder::stream::close_array;
                    using bsoncxx::builder::stream::close_document;
                    using bsoncxx::builder::stream::document;
                    using bsoncxx::builder::stream::finalize;
                    using bsoncxx::builder::stream::open_array;
                    using bsoncxx::builder::stream::open_document;
                    auto uri_str = S_journal_db_config();
                    mongocxx::uri uri( uri_str );
                    mongocxx::client client(uri);
                    mongocxx::database journaldb = client[S_journal_db_config.name];
                    auto journal = journaldb[S_journal_db_config.collection];
                    //auto builder = bsoncxx::builder::stream::document{};
                    std::ifstream json( S_journal_filename + ".json");
                    std::stringstream jsonbuff;
                    jsonbuff << json.rdbuf();
                    // TODO json is read from file. An improvement would be to extract to add
                    // from_ptree method to avoid disk access.
                    bsoncxx::document::value document = bsoncxx::from_json(jsonbuff.str());
                    bsoncxx::stdx::optional<mongocxx::result::insert_one> result = journal.insert_one(document.view());
#endif
                }
            }
        }

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
 
        //! Journal automatic mode enable or disable.
        static bool S_journal_enabled;

        //! Journal automatic mode enable or disable.
        static bool S_journal_auto;

        //! checkpoint number
        static uint32_t S_journal_checkpoint;
};

// Extern explicit instanciation.
template<> std::string JournalManagerBase<>::S_journal_filename;
template<> pt::ptree JournalManagerBase<>::S_journal_ptree;
template<> MongoConfig JournalManagerBase<>::S_journal_db_config;
template<> bool JournalManagerBase<>::S_journal_enabled;
template<> bool JournalManagerBase<>::S_journal_auto;
template<> uint32_t JournalManagerBase<>::S_journal_checkpoint;

// Manager class should be derived from this alias class.
using JournalManager = JournalManagerBase<>;

} // Observer namespace.
} // Feel namespace.

#endif // FEELPP_JOURNALMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
