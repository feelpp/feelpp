/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
/**
   \file heat1d.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelmodels/Rbheat/Cabineheat.hpp>


namespace Feel
{
/**
 * \fn makeAbout()
 * \brief Create the About data of the OpusApp
 *
 */
AboutData
makeAbout()
{
    Feel::AboutData about( "cabineheatfem" ,
                           "cabineheatfem" ,
                           "0.1",
                           "CabineHeat FEM model",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 * \class Heat1dcrbApp
 * \brief Opus application
 *
 * This class implements the Opus application, getting the command
 * line arguments and running the actual code.
 *
 * @author Christophe Prud'homme
 */
class CabineHeatFemApp   : public Application
{
    typedef Application super;
public:

    typedef CRBModel<CabineHeat> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;

    CabineHeatFemApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                              );

        M_crbmodel = crbmodel_ptrtype( new crbmodel_type( this->vm() ) );

    }

    void run()
    {
        //M_crbmodel->solve();
        Feel::ParameterSpace<3>::Element mu( M_crbmodel->parameterSpace() );
        mu <<
           this->vm()["mu1"].as<double>(),
                this->vm()["mu2"].as<double>(),
                this->vm()["mu3"].as<double>(),

                std::cout << "mu = " << mu << "\n";
        M_crbmodel->solve( mu );

    }

private:

    crbmodel_ptrtype M_crbmodel;
};

} // Feel


int
main( int argc, char** argv )
{
    Feel::CabineHeatFemApp app( argc, argv, Feel::makeAbout(), Feel::makeOptions() );
    app.run();
}

