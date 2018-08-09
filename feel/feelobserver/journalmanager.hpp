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

//! JournalManagerBase handles all journalWatchers.
//! The purpose of this class is to be inherited by class that manage the 
//! journal system
//!
//! \note Managers should inherit from type JournalManager directly.
//!
//! \remark Environment is the favored journal manager (child class). However,
//! as mpi, options, etc.. are initialized in Environment constructor, an Env
//! template parameter has to be define to use static members.
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
            std::time_t t = std::time(nullptr);
            const std::string tag = "database.time";
            M_journal_ptree.put( "database.version", FEELPP_DB_JOURNAL_VERSION );
            M_journal_ptree.put( tag + ".time_t", t );
            M_journal_ptree.put( tag + ".gm", std::put_time(std::gmtime(&t), "%c %Z") );
            M_journal_ptree.put( tag + ".local", std::put_time(std::localtime(&t), "%c %Z") );
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
        void JournalSetFilename( std::string name )
        {
            M_journal_filename = name;
        }
        //! @}

        //! Getters
        //! @{
        
        static const bool journalAutoMode()
        {
            return M_journal_auto;
        }

        static const uint32_t
        journalCurrentCheckpoint()
        {
            return M_journal_checkpoint;        
        }

        //! @}

        //! Mutators
        //! @{
        
        static void journalAutoMode( bool m )
        {
            M_journal_auto = m;
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
            ptMerge( M_journal_ptree, pt_merged );

            pt::ptree p;
            if( Env::isParallel() and parallel )
                mpi::reduce( Env::worldComm(), M_journal_ptree, p, ptMerge, Env::masterRank() );
            M_journal_ptree = p;

            return M_journal_ptree;
        }

        //! Save the global property tree into a json file.
        static void
        journalSave( const std::string& filename = "journal" )
        {
            if( Env::isMasterRank() )
            {
                journalJSONSave( filename );
                journalDBSave( filename );
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
                           const std::string& filename="journal" )
        {
            journalPull( parallel );

            std::time_t t = std::time(nullptr);
            const std::string tag = "checkpoints.checkpoint-" + std::to_string( M_journal_checkpoint ) + ".time";
            M_journal_ptree.put( "checkpoints.number", M_journal_checkpoint);
            M_journal_ptree.put( tag + ".time_t", t );
            M_journal_ptree.put( tag + ".gm", std::put_time(std::gmtime(&t), "%c %Z") );
            M_journal_ptree.put( tag + ".local", std::put_time(std::localtime(&t), "%c %Z") );
            
            M_journal_checkpoint++;
            if( save ) 
                journalSave( filename );
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
            M_journal_db_config.name = dbname;
        }

        //! Set the mongodb host.
        static void journalSetDBHost( const std::string&  host = "localhost")
        {
            M_journal_db_config.host = host;
        }

        //! Set the mongodb user name.
        static void journalSetDBUsername( const std::string& user = "")
        {
            M_journal_db_config.user = user;
        }

        //! Set the mongodb user password.
        static void journalSetDBPassword( const std::string& password = "")
        {
            M_journal_db_config.password = password;
        }

        //! Set the mongodb port.
        static void journalSetDBPort( const std::string& port = "27017")
        {
            M_journal_db_config.port = port;
        }

        //! Set the collection used to authenticate.
        static void journalSetDBAuthsrc( const std::string& authsrc = "admin")
        {
            M_journal_db_config.authsrc = authsrc;
        }
        
        //! Set mongodb collection to use for the journal.
        static void journalSetDBCollection( const std::string& dbname = "journal")
        {
            M_journal_db_config.name = dbname;
        }
        
        //! Set mongodb collection to use for the journal.
        static void journalDBConfig( const MongoConfig& mc )
        {
            M_journal_db_config = mc;
        }
        //! @}

private:
        //! Private methods
        //! @{

        //! Save the global property tree into a json file.
        static void
        journalJSONSave( const std::string& filename = "journal" )
        {
            std::string fname = M_journal_filename;
            if( filename != "journal" )
                fname = filename;

            //std::cout << "[Observer: Journal] generate report (JSON).";
            if( not M_journal_ptree.empty() )
                write_json( fname + ".json", M_journal_ptree );
        }

        //! Save the json in a mongodb database. 
        //! This function read a json file in a bson format. Then this bson entry 
        //! is send in the mongodb database. The database has to be configured
        //! beforehand.
        static void
        journalDBSave( const std::string& filename = "journal" )
        {
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
                auto uri_str = M_journal_db_config();
                mongocxx::uri uri( uri_str );
                mongocxx::client client(uri);
                mongocxx::database journaldb = client[M_journal_db_config.name];
                auto journal = journaldb[M_journal_db_config.collection];
                //auto builder = bsoncxx::builder::stream::document{};
                std::ifstream json( filename + ".json");
                std::stringstream jsonbuff;
                jsonbuff << json.rdbuf();
                // TODO json is read from file. An improvement would be to extract to add
                // from_ptree method to avoid disk access.
                bsoncxx::document::value document = bsoncxx::from_json(jsonbuff.str());
                bsoncxx::stdx::optional<mongocxx::result::insert_one> result = journal.insert_one(document.view());
#endif
            }
        }

        //! @}

private:
        //! JSON filename.
        static std::string M_journal_filename;
        //! Global property tree.
        static pt::ptree M_journal_ptree;

        //! MongoDB specific attribute. These attributes are set via options.
        //! @{
        static MongoConfig M_journal_db_config;
        //! @}
 
        //! Journal automatic mode enable or disable.
        static bool M_journal_auto;

        //! checkpoint number
        static uint32_t M_journal_checkpoint;
};

// Extern explicit instanciation.
template<> std::string JournalManagerBase<>::M_journal_filename;
template<> pt::ptree JournalManagerBase<>::M_journal_ptree;
template<> MongoConfig JournalManagerBase<>::M_journal_db_config;
template<> bool JournalManagerBase<>::M_journal_auto;
template<> uint32_t JournalManagerBase<>::M_journal_checkpoint;


// Manager class should be derived from this alias class.
using JournalManager = JournalManagerBase<>;

} // Observer namespace.
} // Feel namespace.

#endif // FEELPP_JOURNALMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
