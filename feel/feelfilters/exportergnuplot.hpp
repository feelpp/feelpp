/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-25

  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file exportergnuplot.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-25
 */
#ifndef __ExporterGnuplot_H
#define __ExporterGnuplot_H 1

#include <iostream>
#include <fstream>


#include <boost/lambda/lambda.hpp>

#include <feel/feelcore/debug.hpp>

#include <feel/feelfilters/exporter.hpp>

namespace Feel
{
/**
 * \class ExporterGnuplot
 * \brief Exporter to GNUPLOT format
 *
 * This class implements exporting meshes and functions using GNUplot
 * in 1D
 *
 * \ingroup Exporter
 * @author Christophe Prud'homme
 */
template<typename MeshType>
class ExporterGnuplot
    :
public Exporter<MeshType>
{
public:

    /**
     * Define enumerations to set plotting properties on construction
     */
    enum PlottingProperties
    {
        GRID_ON    = 1,
        PNG_OUTPUT = 2
    };

    /** @name Typedefs
     */
    //@{

    typedef MeshType mesh_type;

    typedef Exporter<MeshType> super;
    typedef typename mesh_type::value_type value_type;
    typedef typename super::timeset_type timeset_type;
    typedef typename super::timeset_ptrtype timeset_ptrtype;
    typedef typename super::timeset_iterator timeset_iterator;
    typedef typename super::timeset_const_iterator timeset_const_iterator;

    typedef typename matrix_node<value_type>::type matrix_node_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Constructor.
     * To set the properties, we input a bitwise OR of the
     * ExporterGnuplot::PlottingProperties enumerations,
     * e.g. ExporterGnuplot::GRID_ON | ExporterGnuplot::PNG_OUTPUT
     */
    ExporterGnuplot( std::string const& __p = "default", int freq = 1, int properties = 0 );

    ExporterGnuplot( po::variables_map const& vm, std::string const& exp_prefix = "", int properties = 0 );

    ExporterGnuplot( ExporterGnuplot const & __ex );

    ~ExporterGnuplot();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * Set title of plot
     */
    void setTitle( std::string const& title )
    {
        M_title = title;
    }

    /**
     * Turn grid on or off.
     */
    void useGrid( bool grid )
    {
        M_grid = grid;
    }


    /**
     * Write output to a .png file useing gnuplot
     */
    void setPngOutput( bool png_output )
    {
        M_png_output = png_output;
    }

    Exporter<MeshType>* setOptions( po::variables_map const& vm, std::string const& exp_prefix = "" )
    {
        super::setOptions( vm, exp_prefix );

        return this;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
       save the timeset
    */
    void save() const;

    /**
     * export mesh
     */
    void visit( mesh_type* mesh );

    //@}



protected:


private:

    std::string M_title;
    bool M_grid;
    bool M_png_output;
};

} // Feel

#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feelfilters/exportergnuplot.cpp>
#endif // FEELPP_INSTANTIATION_MODE

#endif /* __ExporterGnuplot_H */


