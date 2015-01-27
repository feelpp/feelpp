/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-21

  Copyright (C) 2008 Universite Joseph Fourier (Grenoble I)

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
   \file exporterquick.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-21
 */
#ifndef __FEELPP_EXPORTERQUICK_HPP__
#define __FEELPP_EXPORTERQUICK_HPP__ 1
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/exporterensight.hpp>
#include <feel/feelfilters/exportergmsh.hpp>

namespace Feel
{
/**
 * \class ExporterQuick
 * \brief simple interface to exporter
 *
 */
template<typename MeshType>
class ExporterQuick
{
public:
    typedef MeshType mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef Exporter<mesh_type,1> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;
    typedef typename export_type::timeset_ptrtype timeset_ptrtype;

    ExporterQuick( std::string const& name, po::variables_map& vm )
        :
        exporter( new ExporterEnsight<mesh_type,1> ( vm )),//Exporter<mesh_type>::New( vm["exporter"].template as<std::string>() )->setOptions( vm ) ),
        timeSet( new timeset_type( name ) )
    {
        exporter->setOptions( vm );
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( name );

    }
    ExporterQuick( std::string const& name, std::string const& exp )
        :
        exporter( new ExporterEnsight<mesh_type,1> ( exp,1 ) ),//Exporter<mesh_type>::New( exp ) ),
        timeSet( new timeset_type( name ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( name );

    }
    void save( mesh_ptrtype const& mesh )
    {
        if ( !mesh )
        {
            FEELPP_ASSERT( mesh ).error( "[ExporterQuick] invalid mesh (=0)" );
            return;
        }

        typename timeset_type::step_ptrtype timeStep = timeSet->step( 0 );
        timeStep->setMesh( mesh );

        exporter->save();
    }

    template<typename F1>
    void save( double time, F1 const& f1 )
    {
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( f1.functionSpace()->mesh() );

        timeStep->add( f1.name(), f1 );

        exporter->save();
    }

    template<typename F1, typename F2>
    void save( double time, F1 const& f1, F2 const& f2 )
    {
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( f1.functionSpace()->mesh() );

        timeStep->add( f1.name(), f1 );
        timeStep->add( f2.name(), f2 );

        exporter->save();
    }


    template<typename F1, typename F2, typename F3>
    void save( double time, F1 const& f1, F2 const& f2, F3 const& f3 )
    {
        typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
        timeStep->setMesh( f1.functionSpace()->mesh() );

        timeStep->add( f1.name(), f1 );
        timeStep->add( f2.name(), f2 );
        timeStep->add( f3.name(), f3 );

        exporter->save();
    }

private:
    export_ptrtype exporter;
    timeset_ptrtype timeSet;
};

}
#endif // __FEELPP_EXPORTERQUICK_HPP__
