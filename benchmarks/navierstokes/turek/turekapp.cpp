/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-10

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
   \file turekapp.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-10
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <data.hpp>


namespace Feel
{

/**
 * Turek Model
 *
 */
class TurekApp   : public Application
{
    typedef Application super;
public:

    typedef Data turek_type;
    typedef boost::shared_ptr<turek_type> turek_ptrtype;

    TurekApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        this->changeRepository( boost::format( "%1%/Re_%8%/Simplex_%2%_%3%/P%4%P%5%/h_%6%_%7%/to_%9%_dt_%10%" )
                                % this->about().appName()
                                % this->vm()["d"].as<int>()
                                % this->vm()["order-geo"].as<int>()
                                % this->vm()["order-u"].as<int>() % ( this->vm()["order-u"].as<int>() - 1 )
                                % this->vm()["hsize"].as<double>()
                                % this->vm()["h-cyl-scale"].as<double>()
                                % this->vm()["Re"].as<double>()
                                % this->vm()["bdf.time-order"].as<int>()
                                % this->vm()["bdf.time-step"].as<double>()
                              );


        M_turek = Data::New( this->vm() );

        // print data
        M_turek->print();

    }

    void run()
    {
        M_turek->run();
    }

private:

    turek_ptrtype M_turek;

}; // Turek

} // Feel




int
main( int argc, char** argv )
{
#if 0

    try
    {
        Feel::TurekApp app( argc, argv, Data::makeAbout(), Data::makeOptions() );

        app.run();
    }

    catch ( std::logic_error const& e )
    {
        std::cout << "std::logic_error caught: " << std::endl;
    }

    catch ( ... )
    {
        std::cout << "exception caught: " << std::endl;
    }

#else
    Feel::TurekApp app( argc, argv, Data::makeAbout(), Data::makeOptions() );

    app.run();
#endif
}







