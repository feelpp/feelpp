/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-10

  Copyright (C) 2009-2011 Université Joseph Fourier (Grenoble I)

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
   \file opusscm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-10
 */
#include <eads.hpp>
#include <eadsscmapp.hpp>



namespace Feel
{
template<int OrderU, int OrderP, int OrderT> class OpusModelRB;

po::options_description
makeEadsSCMOptions()
{
    return Feel::makeEadsOptions();

}

EadsSCMApp::EadsSCMApp( AboutData const& ad, po::options_description const& od )
    :
    super( ad, od )
{
    this->init();
}

EadsSCMApp::EadsSCMApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od )
{
    this->init();
}
void EadsSCMApp::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

#if 0

    if ( this->vm()["steady"].as<bool>() )
        this->changeRepository( boost::format( "%1%/P%2%P%3%P%4%/D_%5%/h_%6%/stab_%7%/steady" )
                                % this->about().appName()
                                % this->vm()["order-u"].as<int>() % this->vm()["order-p"].as<int>() % this->vm()["order-temp"].as<int>()
                                % this->vm()["fluid.flow-rate"].as<double>()
                                % this->vm()["hsize"].as<double>()
                                % this->vm()["stab"].as<bool>()
                              );

    else
        this->changeRepository( boost::format( "%1%/P%2%P%3%P%4%/D_%5%/h_%6%/stab_%7%/to_%8%_dt_%9%" )
                                % this->about().appName()
                                % this->vm()["order-u"].as<int>() % this->vm()["order-p"].as<int>()  % this->vm()["order-temp"].as<int>()
                                % this->vm()["fluid.flow-rate"].as<double>()
                                % this->vm()["hsize"].as<double>()
                                % this->vm()["stab"].as<bool>()
                                % this->vm()["bdf.order"].as<int>()
                                % this->vm()["bdf.time-step"].as<double>()
                              );

#endif
    M_opusmodel = opusmodel_ptrtype( new opusmodel_type( this->vm() ) );
    M_scm = scm_ptrtype( new scm_type( this->about().appName(), this->vm() ) );
    M_scm->setTruthModel( M_opusmodel );



}

void
EadsSCMApp::run( std::ofstream& os, scm_type::parameter_type const& mu, int K )
{
    std::cout << "------------------------------------------------------------\n";
    double lb,lbti;
    boost::tie( lb, lbti ) = M_scm->lb( mu, K );
    double ub,ubti;
    boost::tie( ub, ubti ) = M_scm->ub( mu, K );
    double ex,exti;
    boost::tie( ex, exti ) = M_scm->ex( mu );
    std::cout << "lb=" << lb << " ub=" << ub << " ex=" << ex << "\n";
    std::cout << ( ex-lb )/( ub-lb ) << "\n";
    os << K << " "
       << std::setprecision( 16 ) << lb << " "
       << std::setprecision( 3 ) << lbti << " "
       << std::setprecision( 16 ) << ub << " "
       << std::setprecision( 3 ) << ubti << " "
       << std::setprecision( 16 ) << ex << " "
       << std::setprecision( 16 ) << ( ub-lb )/( ub ) << " "
       << std::setprecision( 16 ) << ( ex-lb )/( ex ) << " "
       << std::setprecision( 16 ) << ( ub-ex )/( ex ) << " "
       << "\n";
    std::cout << "------------------------------------------------------------\n";
}
void
EadsSCMApp::run()
{
    std::vector<boost::tuple<double,double,double> > ckconv = M_scm->offline();

    std::ofstream osck( ( boost::format( "ckconv_K_%1%_Mp_%2%_Ma_%3%_Xi_%4%_L_%5%.dat" )
                          % M_scm->KMax()
                          % M_scm->Mplus()
                          % M_scm->Malpha()
                          % this->vm()["crb-scm-sampling-size"].as<int>()
                          % this->vm()["crb-scm-level"].as<int>() ).str().c_str() );

    for ( size_type k = 0; k < ckconv.size(); ++k )
    {
        osck << k << "  "  << std::setprecision( 16 ) << ckconv[k] << "\n";
    }
}

}

