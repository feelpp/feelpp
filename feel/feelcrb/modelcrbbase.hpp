/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Stephane Veys <stephane.veys@imag.fr>
 Date: 2013-02-22

 Copyright (C) 2008-2012 Universite Joseph Fourier (Grenoble I)

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 3.0 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file model.hpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-02-22
*/
#ifndef ModelCrbBase_H
#define ModelCrbBase_H 1

#include <feel/feel.hpp>
#include <feel/feelcrb/eim.hpp>

namespace Feel
{


class ParameterDefinitionBase
{
public :
    typedef ParameterSpace<1> parameterspace_type ;
};

template <typename ParameterDefinition>
class EimDefinitionBase
{

public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    //warning
    //mesh, basis, space and eim objects are here only to have fun_type and fund_type

    /*mesh*/
    typedef Simplex<1> entity_type ;
    typedef Mesh<entity_type > mesh_type ;

    /*basis*/
    typedef Lagrange<1,Scalar> basis_type ;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;

    /* EIM */
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fun_type ;
    typedef EIMFunctionBase<space_type , space_type  , parameterspace_type > fund_type ;

};

template <typename ParameterDefinition = ParameterDefinitionBase, typename EimDefinition = EimDefinitionBase<ParameterDefinition> >
class ModelCrbBase
{

public :
    typedef typename EimDefinition::fun_type fun_type;
    typedef typename EimDefinition::fund_type fund_type;

    typedef boost::shared_ptr<fun_type> fun_ptrtype;
    typedef std::vector<fun_ptrtype> funs_type;

    typedef boost::shared_ptr<fund_type> fund_ptrtype;
    typedef std::vector<fund_ptrtype> funsd_type;

    virtual funs_type scalarContinuousEim () const
    {
        return M_funs;
    }

    virtual funsd_type scalarDiscontinuousEim () const
    {
        return M_funs_d;
    }

protected :

    funs_type M_funs;
    funsd_type M_funs_d;

};

}//Feel
#endif /* __Model_H */

