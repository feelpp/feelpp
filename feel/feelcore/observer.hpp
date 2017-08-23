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

#ifndef FEELPP_OBSERVER_HPP
#define FEELPP_OBSERVER_HPP 1

//#include <feel/feelcore/environment.hpp>
//#include <feel/feelcore/signal.hpp>
#include <boost/signals2.hpp>
#include <boost/any.hpp>
#include <boost/property_tree/ptree.hpp>

namespace Feel
{
namespace Observer
{

namespace pt =  boost::property_tree;

//! Enable signals & slots to inherited objects.
//! Add new functionnalities to connect objects between each others.
//! Note: This class is used to generate benchmark database.
class SignalManager
{
public:
    template< typename Notify, typename Event >
        using sig_type = boost::signals2::signal< Notify, Event >;


    // Create a new Feel++ signal.
    template< typename Notify, typename Event >
    sig_type<Notify, Event>
    newSignal( std::string s )
    {
        const auto v = std::make_shared< sig_type< Notify, Event > >;
        M_sigs.insert(
            std::pair<std::string,
                      sig_type< Notify, Event > > ( s, v )
            );
    };

    // Get an existing Feel++ signal. Use typeid for information, if you
    // desire to dynamic casting to the signal real type.
    // Otherwise, prefer using createSignal.
    boost::any signal( std::string s )
    {
        return M_sigs[s];
    };

private:
    //! Map containing signal name and object.
    std::map< std::string, boost::any > M_sigs;
};


class SlotManager
{
};


//! Retrieve information from child class (observer pattern)
//!
//! Each child must implement its inform() method. It returns a property tree
//! containing objects information, aka simulations input/output parameters
//! important for reproducibility.
//class SimulationInfo
//    :
//        public ObserverBase< pt::ptree ( std::string ), pt::ptree >
//{
//public:
//    //! Watch child properties.
//    virtual pt::ptree& notify() const = 0;
//};


} // Observer namespace
} // Feel namespace.

#endif // FEELPP_OBSERVER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
