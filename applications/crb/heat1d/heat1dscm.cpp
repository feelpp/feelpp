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
#include <boost/tuple/tuple_io.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

//#include <opusdata.hpp>
//#include <opusmodelrb.hpp>
#include <heat1d.hpp>

#include <feel/feelcrb/crbscm.hpp>
#include <feel/feelcrb/crbmodel.hpp>


namespace Feel
{
template<int OrderU, int OrderP, int OrderT> class OpusModelRB;

/**
 * \fn makeAbout()
 * \brief Create the About data of the OpusApp
 *
 */
AboutData
makeAbout()
{
    Feel::AboutData about( "heat1dscm" ,
                           "heat1dscm" ,
                           "0.1",
                           "1D/2D OPUS/EADS Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010-2012 Universite de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 * \class OpusSCMApp
 * \brief Opus application
 *
 * This class implements the Opus application, getting the command
 * line arguments and running the actual code.
 *
 * @author Christophe Prud'homme
 */
class OpusSCMApp   : public Application
{
    typedef Application super;
public:

    //typedef CRBModel<OpusModelRB<2,1,2> > opusmodel_type;
    typedef CRBModel<Heat1D> opusmodel_type;
    typedef boost::shared_ptr<opusmodel_type> opusmodel_ptrtype;
    typedef CRBSCM<opusmodel_type> scm_type;
    typedef boost::shared_ptr<scm_type> scm_ptrtype;

    OpusSCMApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
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
                                    % this->vm()["fluid-flow-rate"].as<double>()
                                    % this->vm()["hsize"].as<double>()
                                    % this->vm()["stab"].as<bool>()
                                  );

        else
            this->changeRepository( boost::format( "%1%/P%2%P%3%P%4%/D_%5%/h_%6%/stab_%7%/to_%8%_dt_%9%" )
                                    % this->about().appName()
                                    % this->vm()["order-u"].as<int>() % this->vm()["order-p"].as<int>()  % this->vm()["order-temp"].as<int>()
                                    % this->vm()["fluid-flow-rate"].as<double>()
                                    % this->vm()["hsize"].as<double>()
                                    % this->vm()["stab"].as<bool>()
                                    % this->vm()["bdf.time-order"].as<int>()
                                    % this->vm()["bdf.time-step"].as<double>()
                                  );

#endif
        M_opusmodel = opusmodel_ptrtype( new opusmodel_type( this->vm() ) );
        M_scm = scm_ptrtype( new scm_type( this->about().appName(), this->vm() ) );
        M_scm->setTruthModel( M_opusmodel );



    }

    void run( std::ofstream& os, scm_type::parameter_type const& mu, int K )
    {
        std::cout << "------------------------------------------------------------\n";
        double lb,lbti;
        boost::tie( lb, lbti ) = M_scm->lb( mu, K );
        double ub,ubti;
        boost::tie( ub, ubti ) = M_scm->ub( mu, K );
        double ex, exti;
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
    void run()
    {

#if 0

        if ( this->vm()["scm-generate"].as<bool>() )
            M_scm->offline();

        scm_type::bounds_type bounds;
        bounds =  M_scm->online();
#else
        M_opusmodel->init();
        std::vector<boost::tuple<double,double,double> > ckconv = M_scm->offline();

        std::ofstream osck( ( boost::format( "ckconv_K_%1%_Mp_%2%_Ma_%3%_Xi_%4%_L_%5%.dat" )
                              % M_scm->KMax()
                              % M_scm->Mplus()
                              % M_scm->Malpha()
                              % this->vm()["crb-scm-sampling-size"].as<int>()
                              % this->vm()["crb-scm-level"].as<int>() ).str().c_str() );

        for ( int k = 0; k < ckconv.size(); ++k )
        {
            osck << k << "  "  << std::setprecision( 16 ) << ckconv[k] << "\n";
        }

        scm_type::parameter_type mu( M_scm->Dmu() );
        std::ofstream ofs( ( boost::format( "eval_K_%1%_Mp_%2%_Ma_%3%_Xi_%4%_L_%5%.dat" )
                             % M_scm->KMax()
                             % M_scm->Mplus()
                             % M_scm->Malpha()
                             % this->vm()["crb-scm-sampling-size"].as<int>()
                             % this->vm()["crb-scm-level"].as<int>() ).str().c_str() );
        ;
#if 0
        mu << 0.4, 10, 1, 1;
        run( ofs, mu, M_scm->KMax() );
        mu << 0.9, 20, 1, 1;
        run( ofs, mu, M_scm->KMax() );
        mu << 50, 20, 1, 1;
        run( ofs, mu, M_scm->KMax() );
        mu << 1, 0.2, 1, 1;
        run( ofs, mu, M_scm->KMax() );
        mu << 0.2, 3, 1, 1;
        run( ofs, mu, M_scm->KMax() );
#endif
        mu << 14.65866591112086, 27.07274082090771, 1, 1;
        run( ofs, mu, M_scm->KMax() );

        std::ofstream ofs2( ( boost::format( "conv_K_%1%_Mp_%2%_Ma_%3%_Xi_%4%_L_%5%.dat" )
                              % M_scm->KMax()
                              % M_scm->Mplus()
                              % M_scm->Malpha()
                              % this->vm()["crb-scm-sampling-size"].as<int>()
                              % this->vm()["crb-scm-level"].as<int>() ).str().c_str() );
        mu << 14.65866591112086, 27.07274082090771, 1, 1;

        //mu << 0.2, 0.2, 1, 1;
        for ( int k = 1; k <= M_scm->KMax(); ++k )
        {
            run( ofs2, mu, k );
        }

#endif

    }

private:

    opusmodel_ptrtype M_opusmodel;
    scm_ptrtype M_scm;
}; // Opus

} // Feel


int
main( int argc, char** argv )
{
    Feel::OpusSCMApp app( argc, argv,
                          Feel::makeAbout(),
                          Feel::makeHeat1DOptions() );

    app.run();
}









