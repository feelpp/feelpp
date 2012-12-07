/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-03-24

  Copyright (C) 2009-2011 Université de Grenoble 1

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
   \file applicationxml.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-03-24
 */
#ifndef __ApplicationXML_H
#define __ApplicationXML_H 1

#include <boost/lexical_cast.hpp>

#include <feel/feelcore/application.hpp>

#include <feel/feelcore/xmlparser.hpp>

namespace Feel
{
/**
 * \class ApplicationXML
 * \brief XML application
 *
 * @author Christophe Prud'homme
 * @see
 */
class ApplicationXML : public Application
{
    //! super class
    typedef Application super;
public:


    /** @name Constants
     */
    //@{

    enum RunStatus
    {
        RUN_CONTINUE,
        RUN_EXIT
    };


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    ApplicationXML( int argc, char** argv, AboutData const& ad, po::options_description const& od );
    //! copy constructor
    ApplicationXML( ApplicationXML const & );
    //! destructor
    ~ApplicationXML();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    ApplicationXML& operator=( ApplicationXML const & o );
    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    ApplicationXML& addOutput( Output const& o )
    {
        M_outputs.push_back( o );
        return *this;
    }
    ApplicationXML& addParameter( Parameter const& o )
    {
        M_params.push_back( o );
        return *this;
    }

    template<typename T>
    ApplicationXML& addOutputValue( T const& val )
    {
        std::ostringstream oss;
        oss << val;
        M_output_values.push_back( oss.str() );
        return *this;
    }

    template<typename T>
    ApplicationXML& addParameterValue( T const& val )
    {
        char sci_val[11];
        sprintf( sci_val,"%.5e",( double )val );
        M_parameter_values.push_back( std::string( sci_val ) );
        return *this;
    }

    virtual RunStatus preProcessing();

    virtual void postProcessing();

    //@}



protected:

    //! parameters
    std::vector<Parameter> M_params;

    //! outputs
    std::vector<Output> M_outputs;

    //! parameter values
    std::vector<std::string> M_parameter_values;

    //! output values
    std::vector<std::string> M_output_values;


private:

};
} // Feel
#endif /* __ApplicationXML_H */
