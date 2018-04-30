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
#include <boost/property_tree/ptree.hpp>
#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/functors/journalmerge.hpp>

#if defined(FEELPP_HAS_MONGOCXX )
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#endif

#define FEELPP_DB_JOURNAL_VERSION "0.1.0"

namespace Feel
{
namespace Observer
{

namespace pt =  boost::property_tree;

// Manager
class JournalManager
   : public Event::SignalHandler
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
        JournalManager()
        {
            std::time_t t = std::time(nullptr);
            M_journal_ptree.put( "database.version", FEELPP_DB_JOURNAL_VERSION );
            M_journal_ptree.put( "database.time.gm", std::put_time(std::gmtime(&t), "%c %Z") );
            M_journal_ptree.put( "database.time.local", std::put_time(std::localtime(&t), "%c %Z") );
            // Create a signal for simulation info.
            signalStaticNew< notify_type (), JournalMerge >( "journalManager" );
        }

        //! @}
        //
        //! Setters
        //! @{

        //! Set JSON file name.
        void JournalSetFilename( std::string name )
        {
            M_journal_filename = name;
        }
        //! @}

        //! Misc
        //! @{

        //! Fetch and merge notifications coming from all observed objects into
        //! a global property tree.
        //! \see JournalWatcher
        static const notify_type
        journalPull()
        {
            //std::cout << "[Observer: Journal] Pull watcher information.";
            const pt::ptree& pt_merged = signalStaticPull< notify_type (), JournalMerge >( "journalManager" );
//            for( const auto& it : pt_merged )
//                M_journal_ptree.put_child( it.first, it.second );
            ptMerge( M_journal_ptree, pt_merged );
            return M_journal_ptree;
        }

        //! Save the simulation info global property tree into a json file.
        static void
        journalSave( const std::string& filename = "journal" )
        {
            std::string fname = M_journal_filename;
            if( filename != "journal" )
                fname = filename;

            //std::cout << "[Observer: Journal] generate report (JSON).";
            if( not M_journal_ptree.empty() )
                write_json( fname + ".json", M_journal_ptree );
        }

        //! @}

    private:
        static void
        journalDBsave()
        {
#if defined(FEELPP_HAS_MONGOCXX )
            using bsoncxx::builder::stream::close_array;
            using bsoncxx::builder::stream::close_document;
            using bsoncxx::builder::stream::document;
            using bsoncxx::builder::stream::finalize;
            using bsoncxx::builder::stream::open_array;
            using bsoncxx::builder::stream::open_document;
#endif
        }

    private:
        //! JSON filename.
        static std::string M_journal_filename;
        //! Global property tree.
        static pt::ptree M_journal_ptree;
};

} // Observer namespace
} // Feel namespace.

#endif // FEELPP_JOURNALMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
