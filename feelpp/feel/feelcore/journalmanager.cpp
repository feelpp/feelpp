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

#include <feel/feelcore/journalmanager.hpp>

//#include <boost/property_tree/json_parser.hpp>
//#include <boost/property_tree/ptree_serialization.hpp>
#include <feel/feelcore/environment.hpp>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#endif

namespace Feel
{

JournalManager::JournalManager()
{
#if defined(FEELPP_HAS_MONGOCXX )
    // Create mongo unique instance!
    Feel::MongoCxx::instance();
#endif
    //M_journal_uuid = Env::randomUUID();
    std::time_t t = std::time(nullptr);
    const std::string tag = "journal.time";
    S_journal_ptree["journal"]["version"] = FEELPP_DB_JOURNAL_VERSION;
    S_journal_ptree["journal"]["time"]["time_t"] = t;
    std::ostringstream ostr; ostr << std::put_time(std::gmtime(&t), "%c %Z");
    S_journal_ptree["journal"]["time"]["gm"] = ostr.str();
    ostr.str(""); ostr << std::put_time(std::localtime(&t), "%c %Z");
    S_journal_ptree["journal"]["time"]["local"] = ostr.str();

    // Create a signal for simulation info.
    signalStaticNew< void ( notify_type& ) >( "journalManager" );
}


//! Save the global property tree into a json file.
void
JournalManager::journalSave( std::string const& filename )
{
    std::string const& filenameUsed = filename.empty()? journalFilename() : filename;

    if( Environment::isMasterRank() )
    {
        journalJSONSave( filenameUsed, S_journal_ptree );
        journalDBSave( filenameUsed );
    }
}

void
JournalManager::journalCheckpoint( bool save,
                                   const std::string& filename )
{
    if ( journalEnabled() )
    {
        S_journal_checkpoint++;

        journalPull();

        std::time_t t = std::time(nullptr);
        auto & jJournalCheckpoints = S_journal_ptree["journal"]["checkpoints"];
        auto & jJournalCheckpointsCurrentTime = jJournalCheckpoints["checkpoint-" + std::to_string( S_journal_checkpoint )]["time"];
        jJournalCheckpoints["number"] = S_journal_checkpoint;
        jJournalCheckpointsCurrentTime["time_t"] = t;
        std::ostringstream ostr;  ostr << std::put_time(std::gmtime(&t), "%c %Z");
        jJournalCheckpointsCurrentTime["gm"] = ostr.str();
        ostr.str(""); ostr << std::put_time(std::localtime(&t), "%c %Z");
        jJournalCheckpointsCurrentTime["local"] = ostr.str();

        if ( save )
            journalSave( filename );
    }
}

void
JournalManager::journalJSONSave( const std::string& filename, nl::json const& pt )
{
    if( !journalEnabled() )
        return;

    //std::cout << "[Observer: Journal] generate report (JSON).";
    //write_json( filename, pt );
    std::ofstream o( filename );
    o << pt.dump(1);
}

void
JournalManager::journalDBSave( const std::string& filename )
{
    if( journalEnabled() )
    {
        bool enable = boption(_name="journal.database");
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
            std::ifstream json( filename );
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


// Init static variables.
std::string JournalManager::S_journal_filename = "journal.json";
nl::json JournalManager::S_journal_ptree = {};
MongoConfig JournalManager::S_journal_db_config = {};
uint32_type JournalManager::S_journal_checkpoint = 0;

bool JournalManager::S_enable = false;
bool JournalManager::S_automode = false;

} // Feel namespace.


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
