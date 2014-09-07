/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-25

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file stvenant_kirchhoff_main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-25
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/importergmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>


#include <structurebase.hpp>




namespace Feel
{
/**
 * StVenantKirchhoffApp Application
 *
 */
class StructureApp   : public Application
{
    typedef Application super;
public:

    typedef StructureBase structure_type;
    typedef boost::shared_ptr<structure_type> structure_ptrtype;

    StructureApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        this->changeRepository( boost::format( "benchmarks/convergence/%1%/%2%/Simplex_%3%_1/P%4%/h_%5%" )
                                % this->about().appName()
                                % this->vm()["material"].as<std::string>()
                                % 2 //% this->vm()["dim"].as<int>()
                                % this->vm()["sorder"].as<int>()
                                % this->vm()["hsize"].as<double>()
                              );


        M_structure = structure_type::New( this->vm() );

        // print data
        M_structure->print();

    }

    void run()
    {
        M_structure->run();
    }

private:

    structure_ptrtype M_structure;

}; // StructureApp

} // Feel




int
main( int argc, char** argv )
{
    Feel::StructureApp app( argc, argv, Feel::StructureBase::makeAbout(), Feel::StructureBase::makeOptions() );

    app.run();
}







