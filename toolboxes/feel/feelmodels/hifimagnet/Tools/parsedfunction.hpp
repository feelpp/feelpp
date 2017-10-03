/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   \file parsedfunction.h
   \author Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2012-06-15
*/

#ifndef __PARSEDFUNCTION_HPP
#define __PARSEDFUNCTION_HPP 1

#include <feel/options.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>
using namespace Feel;

#include <matheval.h>
#include <string>

//namespace ublas = boost::numeric::ublas;


class parsedfunc
{
public :
    /* introduce matheval function */
    void *f, *f_prim;	/* Evaluators for function and function derivative.  */
    char **names;		/* Function variables names. */
    int count;			/* Number of function variables. */

    parsedfunc() : f((void *) NULL), names((char **) NULL), count(0)
    {
        LOG(INFO) << "create empty function\n";
    };

    parsedfunc(const std::string & function)
    {
        f = evaluator_create (const_cast<char *> (function.c_str()) );
        assert(f);
        evaluator_get_variables (f, &names, &count);
        LOG(INFO) << "regular constructor ";
        LOG(INFO) << "f()=" <<evaluator_get_string(f) << " ---- " << count << " var\n";
    };

    parsedfunc(const parsedfunc & func)
    {
        if ( !func.f)
            LOG(INFO) << "copy constructor with empty matheval\n";

        f = evaluator_create (evaluator_get_string(func.f));
        names = func.names;
        count = func.count;

        LOG(INFO) << "copy constructor ";
        LOG(INFO) << "f()=" <<evaluator_get_string(f) << " ---- " << count << " var\n";
    };

    ~parsedfunc()
    {
        if (f)
            {
                LOG(INFO) << "destructor : f() " << evaluator_get_string(f) << "\n";
                evaluator_destroy (f);
                f = (void *) NULL;
                names = (char **) NULL;
            }
    };

    parsedfunc & operator= (const parsedfunc & func)
    {
        if ( this != &func)
            {
                f = evaluator_create (evaluator_get_string(func.f));
                evaluator_get_variables (f, &names, &count);
            };

        LOG(INFO) << "operator= constructor ";
        LOG(INFO) << "f()=" << evaluator_get_string(f) << " ---- " << count << " var\n";
        return *this;
    };


    static const size_type context = vm::JACOBIAN|vm::POINT;
    typedef double value_type;
    typedef Feel::uint16_type uint16_type;
    static const uint16_type rank = 0;
    static const uint16_type imorder = 2;
    static const bool imIsPoly = true;

    double operator()( uint16_type, uint16_type, ublas::vector<double> const& x, ublas::vector<double> const& /*n*/ ) const
    {

        int nparams = count > x.size() ? count : x.size();
        double *params = new double [nparams];
        for (int i=0; i<x.size(); i++)
             params[i] = x[i];

        double val = evaluator_evaluate(f, count, names, params);
        delete[] params;

        return val;

    };
};
#endif /* __PARSEDFUNCTION_HPP */
