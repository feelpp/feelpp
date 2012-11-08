/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-10

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
   \file opusscm.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-08-10
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/applicationxml.hpp>
#include <feel/options.hpp>

//#include <opusdata.hpp>
//#include <opusmodelrb.hpp>
#include <crbmodel.hpp>
#include <heat1d.hpp>


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
    Feel::AboutData about( "opuseigs" ,
                           "opuseigs" ,
                           "0.1",
                           "2D OPUS/EADS Benchmark (eigenvalue problems)",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 * \class OpusAppEigs
 * \brief Opus application
 *
 * This class implements the Opus application, getting the command
 * line arguments and running the actual code.
 *
 * @author Christophe Prud'homme
 */
class OpusAppEigs   : public ApplicationXML
{
    typedef ApplicationXML super;
public:

    //typedef CRBModel<OpusModelRB<2,1,2> > opusmodel_type;
    typedef CRBModel<Heat1D> opusmodel_type;
    typedef boost::shared_ptr<opusmodel_type> opusmodel_ptrtype;

    OpusAppEigs( int argc, char** argv, AboutData const& ad, po::options_description const& od )
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
        Parameter h( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="2e-4:10:1e-3" );
        this->
        addParameter( Parameter( _name="order",_type=DISC_ATTR,_latex="N_T",_values=boost::lexical_cast<std::string>( 2 ).c_str() ) )
        .addParameter( Parameter( _name="kic",_type=CONT_ATTR,_latex="k_{\\mathrm{IC}}",_values="0.2:10:150" ) )
        .addParameter( Parameter( _name="D",_type=CONT_ATTR,_latex="D",_values="1e-4:10:1e-3" ) )
        .addParameter( h );

        LOG(INFO) << "parameter added\n";
        this->
        addOutput( Output( _name="eigmin",_latex="\\alpha^{\\mathcal{N}}(\\mu)" ) )
        .addOutput( Output( _name="eigmax",_latex="\\gamma^{\\mathcal{N}}(\\mu)" ) );
        LOG(INFO) << "output added\n";

    }

    void run()
    {
        std::cout << "start run()\n";
        std::cout << "hsize=" << this->vm()["hsize"].as<double>() << std::endl;
#if 0
        //std::cout << "D=" << this->vm()["fluid-flow-rate"]. as<double>() << std::endl;
        //std::cout << "kic=" << this->vm()["kic"].as<double>() << std::endl;
        this->addParameterValue( 2 )
        .addParameterValue( this->vm()["hsize"]. as<double>() )
        .addParameterValue( this->vm()["fluid-flow-rate"]. as<double>() )
        .addParameterValue( this->vm()["kic"]. as<double>() );
        std::cout << "parameter defined\n";
        RunStatus ierr = this->preProcessing();

        if ( ierr == RUN_EXIT )
            return;

#endif
        double eigmin, eigmax;
        std::map<double, boost::tuple<double,double,double,int,double,double> > res = M_opusmodel->run();

        boost::tie( boost::tuples::ignore, eigmin, eigmax, boost::tuples::ignore, boost::tuples::ignore, boost::tuples::ignore ) = res.begin()->second;

#if 0
        this->addOutputValue( eigmin ).addOutputValue( eigmax );
        this->postProcessing();
#endif

        std::map<double, boost::tuple<double,double,double,int,double,double> > res2 = M_opusmodel->runq();

        std::map<double, boost::tuple<double,double,double,int,double,double> >::iterator it = res.begin();
        std::map<double, boost::tuple<double,double,double,int,double,double> >::iterator end = res.end();

        for ( ; it != end; ++it )
        {
            std::cout <<"| " << std::setprecision( 4 ) << it->second.get<0>()
                      << " | " << std::setprecision( 16 ) << it->second.get<01>()
                      << " | " << it->second.get<4>()
                      << " | " << it->second.get<2>()
                      << " | " << it->second.get<5>()
                      << " | " << it->second.get<3>()
                      << " |\n";
        }

        std::map<double, boost::tuple<double,double,double,int,double,double> >::iterator it2 = res2.begin();
        std::map<double, boost::tuple<double,double,double,int,double,double> >::iterator end2 = res2.end();

        for ( ; it2 != end2; ++it2 )
        {
            std::cout <<"| " << std::setprecision( 4 ) << it2->second.get<0>()
                      << " | " << std::setprecision( 16 ) << it2->second.get<01>()
                      << " | " << it2->second.get<4>()
                      << " | " << it2->second.get<2>()
                      << " | " << it2->second.get<5>()
                      << " | " << it2->second.get<3>()
                      << " |\n";
        }
    }

private:

    opusmodel_ptrtype M_opusmodel;

}; // Opus

} // Feel


int
main( int argc, char** argv )
{
    Feel::OpusAppEigs app( argc, argv,
                           Feel::makeAboutHeat1D(),
                           Feel::makeOptions() );
    //Feel::OpusData::makeOptions()
    //.add( Feel::opusModelOptions() )
    //.add( Feel::solvereigen_options() ) );

    app.run();
}









