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

#ifndef FEELPP_SIMINFOMANAGER_HPP
#define FEELPP_SIMINFOMANAGER_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <feel/feelevent/events.hpp>
#include <feel/feelobserver/functors/siminfomerge.hpp>

namespace Feel
{
namespace Observer
{

namespace pt =  boost::property_tree;

// Manager
class SimInfoManager
   : public Event::SignalHandler
{
    public:
        // Types alias.
        using notify_type = pt::ptree;

        //! Constructor
        //! @{

        //! Default constructor
        //! This constructor create a new signal waiting for simulation info
        //! watcher object to connect.
        //! \see SimInfoWatcher
        SimInfoManager()
        {
            // Create a signal for simulation info.
            signalStaticNew< notify_type (), SimInfoMerge >( "simInfoManager" );
        }

        //! @}
        //
        //! Setters
        //! @{

        //! Set JSON file name.
        void SimInfoSetFilename( std::string name )
        {
            M_siminfo_filename = name;
        }
        //! @}

        //! Misc
        //! @{

        //! Fetch and merge notifications coming from all observed objects into
        //! a global property tree.
        //! \see SimInfoWatcher
        static const notify_type
        simInfoPull()
        {
            M_siminfo_ptree = signalStaticPull< notify_type (), SimInfoMerge >( "simInfoManager" );
            return M_siminfo_ptree;
        }

        //! Save the simulation info global property tree into a json file.
        static void
        simInfoSave()
        {
            if( not M_siminfo_ptree.empty() )
                write_json( M_siminfo_filename , M_siminfo_ptree );
        }

        //! @}

    private:
        //! JSON filename.
        static std::string M_siminfo_filename;
        //! Global property tree.
        static pt::ptree M_siminfo_ptree;
};

} // Observer namespace
} // Feel namespace.

#endif // FEELPP_SIMINFOMANAGER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
