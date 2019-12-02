/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Jean-Baptiste Wahl
            Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Dec 2019

 Copyright (C) 2018-2019 Feel++ Consortium

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#ifndef FEELPP_SOLVERBASE_HPP
#define FEELPP_SOLVERBASE_HPP 1

#include <feel/feelfmi/fmumodelbase.hpp>

namespace Feel
{

/**
 * This virtual class provide a generic interface for all FMI solvers
 */
class SolverBase
{
public :
    typedef FmuModelBase fmumodel_type;
    typedef std::shared_ptr<fmumodel_type> fmumodel_ptrtype;

    SolverBase()
    {}

    SolverBase( fmumodel_ptrtype model ) :
        M_model( model ),
        M_tcur( 0 )
    {}

    virtual ~SolverBase()
    {}

    void setTimeStep( double const& step )
    {
        M_step = step;
    }

    void setRelativeTol( double const& tol )
    {
        M_tol = tol;
    }

    double currentTime() const
    {
        return M_tcur;
    }

    virtual void initialize( double const& t_init, double const& t_final, double const& tol )=0;
    virtual void simulate()=0;
    virtual void doSteps( double t_stop )=0;

protected :
    fmumodel_ptrtype M_model;
    double M_step, M_tol, M_tcur, M_tfinal;

}; //class SovlerBase


} // namespace Feel

#endif
