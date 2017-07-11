// This file is part of the Feel library
//
// Author(s): Guillaume Dolle <gdolle@unistra.fr>
//      Date: 2012-02-14
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

#include <feel/feelcore/environment.hpp>

namespace Feel
{
namespace pt =  boost::property_tree;

//! Retrieve information from child class (observer pattern)
//!
//! Each child must implement its inform() method. It returns a property tree
//! containing objects information, aka simulations input/output parameters
//! important for reproducibility.
class Observer
{
protected:
    //! Default constructor.
    //!
    //! When this constructor is called by a child class, the new child observer
    //! is connected to the global signal contained by Environment, to inform
    //! Environment to watch him.
    Observer()
    {
        Environment::simDataSignals().connect(std::bind(&Observer::watch,this));
    }
public:
    //! Watch child properties.
    virtual pt::ptree& watch() const = 0;
};

} // Feel namespace.

#endif // FEELPP_OBSERVER_HPP


// MODELINES
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t
// -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:
