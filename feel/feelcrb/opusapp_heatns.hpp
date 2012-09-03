/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-06-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file opusapp.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-06-18
 */
#ifndef __OpusApp_heatns_H
#define __OpusApp_heatns_H 1

#include <boost/assign/std/vector.hpp>
#include <feel/feelcrb/crb_heatns.hpp>
#include <feel/feelcrb/crbmodel_heatns.hpp>
#include <boost/serialization/version.hpp>

namespace Feel
{
po::options_description opusapp_heatns_options( std::string const& prefix );
std::string _o( std::string const& prefix, std::string const& opt )
{
    std::string o = prefix;

    if ( !o.empty() )
        o += ".";

    return o + opt;
}


enum class SamplingMode_heatns
{
    RANDOM = 0, EQUIDISTRIBUTED = 1
};

#define prec 4
#define Pdim 7
#define fill ' '
#define dmanip std::scientific << std::setprecision( prec )
#define hdrmanip(N) std::setw(N) << std::setfill(fill) << std::right
#define tabmanip(N) std::setw(N) << std::setfill(fill) << std::right << dmanip


/**
 * \class OpusApp
 * \brief Certified Reduced Basis application
 *
 * @author Christophe Prud'homme
 */
template<typename ModelType>
class OpusApp_heatns   : public Application
{
    typedef Application super;
public:

    typedef CRBModel_heatns<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB_heatns<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    OpusApp_heatns( AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_heatns_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode_heatns )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
    {
        this->init();
    }

    OpusApp_heatns( AboutData const& ad, po::options_description const& od, CRBModelMode_heatns mode )
        :
        super( ad, opusapp_heatns_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( mode )
    {
        this->init();
    }

    OpusApp_heatns( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, opusapp_heatns_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode_heatns )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
    {
        this->init();
    }
    OpusApp_heatns( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode_heatns mode )
        :
        super( argc, argv, ad, opusapp_heatns_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( mode )
    {
        this->init();        
    }
    
    void init()
    {
        
        try
        {
            M_current_path = fs::current_path();

            std::srand( static_cast<unsigned>( std::time( 0 ) ) );
            Debug() << "[OpusApp_heatns] constructor " << this->about().appName()  << "\n";
            this->changeRepository( boost::format( "%1%/h_%2%/" )
                                    % this->about().appName()
                                    % this->vm()["hsize"].template as<double>()
                                  );
            Debug() << "[OpusApp_heatns] ch repo" << "\n";
            this->setLogs();
            Debug() << "[OpusApp_heatns] set Logs" << "\n";
            Debug() << "[OpusApp_heatns] mode:" << ( int )M_mode << "\n";
            
            
            model = crbmodel_ptrtype( new crbmodel_type( this->vm(),M_mode ) );
            Debug() << "[OpusApp_heatns] get model done" << "\n";
            
            crb = crb_ptrtype( new crb_type( this->about().appName(),
                                             this->vm() ,
                                             model ) );
            Debug() << "[OpusApp_heatns] get crb done" << "\n";
            //crb->setTruthModel( model );
            Debug() << "[OpusApp_heatns] constructor done" << "\n";
        }

        catch ( boost::bad_any_cast const& e )
        {
            std::cout << "[OpusApp_heatns] a bad any cast occured, probably a nonexistant or invalid  command line/ config options\n";
            std::cout << "[OpusApp_heatns] exception reason: " << e.what() << "\n";
        }
        
    }

    void setMode( std::string const& mode )
    {        
        if ( mode == "pfem" ) M_mode = CRBModelMode_heatns::PFEM;

        if ( mode == "crb" ) M_mode = CRBModelMode_heatns::CRB;

        if ( mode == "scm" ) M_mode = CRBModelMode_heatns::SCM;

        if ( mode == "scm_online" ) M_mode = CRBModelMode_heatns::SCM_ONLINE;

        if ( mode == "crb_online" ) M_mode = CRBModelMode_heatns::CRB_ONLINE;
    }
    void setMode( CRBModelMode_heatns mode )
    {        
        M_mode = mode;
    }

    void loadDB()
    {
        if ( M_mode == CRBModelMode_heatns::PFEM )
            return;

        if ( !crb->scm()->isDBLoaded() || crb->scm()->rebuildDB() )
        {
            if ( M_mode == CRBModelMode_heatns::SCM )
            {
                std::cout << "No SCM DB available, do scm offline computations first...\n";
                if( crb->scm()->doScmForMassMatrix() )
                    crb->scm()->setScmForMassMatrix( true );

                crb->scm()->offline();
            }
        }
        else
        {
            crb->scm()->computeAffineDecomposition();
        }

        if ( !crb->isDBLoaded() || crb->rebuildDB() )
        {
            if ( M_mode == CRBModelMode_heatns::CRB )
                //|| M_mode == CRBModelMode_heatns::SCM )
            {
                std::cout << "No CRB DB available, do crb offline computations...\n";
                crb->offline();
            }

            else if ( M_mode != CRBModelMode_heatns::SCM )
                throw std::logic_error( "CRB/SCM Database could not be loaded" );
        }

        if( crb->isDBLoaded() )
        {
            int current_dimension = crb->dimension();
            int dimension_max = this->vm()["crb.dimension-max"].template as<int>();
            if( current_dimension < dimension_max )
                crb->offline();
        }
    }
    FEELPP_DONT_INLINE
    void run()
    {
        
        std::cout << "\n-------->OpusApp_heatns::run(). \n \n";
        
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        this->loadDB();
        int run_sampling_size = this->vm()[_o( this->about().appName(),"run.sampling.size" )].template as<int>();
        SamplingMode_heatns run_sampling_type = ( SamplingMode_heatns )this->vm()[_o( this->about().appName(),"run.sampling.mode" )].template as<int>();
        int output_index = this->vm()["crb.output-index"].template as<int>();
        //int output_index = this->vm()[_o(this->about().appName(),"output.index")].template as<int>();

        
        run_sampling_size = 1;
        std::cout << "-------->run_sampling_size = " << run_sampling_size << std::endl;
        std::cout << "-------->output_index = " << output_index << std::endl;

        
        
        typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );

        switch ( run_sampling_type )
        {
        default:
        case SamplingMode_heatns::RANDOM:
                std::cout << "-------->SamplingMode_heatns::RANDOM\n";
            Sampling->randomize( run_sampling_size  );
            break;

        case SamplingMode_heatns::EQUIDISTRIBUTED:
                std::cout << "-------->SamplingMode_heatns::EQUIDISTRIBUTED\n";
            Sampling->equidistribute( run_sampling_size  );
            break;
        }

        std::map<CRBModelMode_heatns,std::vector<std::string> > hdrs;
        using namespace boost::assign;
        std::vector<std::string> pfemhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" );
        std::vector<std::string> crbhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" )( "RB Output" )( "Error Bounds" )( "CRB Time" )( "Relative error PFEM/CRB" )( "Conditionning" );
        std::vector<std::string> scmhdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "FEM" )( "FEM Time" )( "Rel.(FEM-Lb)" );
        std::vector<std::string> crbonlinehdrs = boost::assign::list_of( "RB Output" )( "Error Bounds" )( "CRB Time" );
        std::vector<std::string> scmonlinehdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "Rel.(FEM-Lb)" );
        hdrs[CRBModelMode_heatns::PFEM] = pfemhdrs;
        hdrs[CRBModelMode_heatns::CRB] = crbhdrs;
        hdrs[CRBModelMode_heatns::SCM] = scmhdrs;
        hdrs[CRBModelMode_heatns::CRB_ONLINE] = crbonlinehdrs;
        hdrs[CRBModelMode_heatns::SCM_ONLINE] = scmonlinehdrs;
        std::ostringstream ostr;
        
        if( crb->printErrorDuringOfflineStep() )
            crb->printErrorsDuringRbConstruction();
        if ( crb->showMuSelection() )
            crb->printMuSelection();
        
        printParameterHdr( ostr, model->parameterSpace()->dimension(), hdrs[M_mode] );
        BOOST_FOREACH( auto mu, *Sampling )
        {
            int size = mu.size();
            std::cout << "mu = [ ";
            
            for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";
            
            std::cout<< mu[size-1]<<" ] \n";
            
            switch ( M_mode )
            {
            case  CRBModelMode_heatns::PFEM:
            {
                std::cout << "-------->CRBModelMode_heatns::PFEM\n";
                boost::timer ti;
                //model->solve( mu );
                std::vector<double> o = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );
                std::cout << "output=" << o[0] << "\n";
                printEntry( ostr, mu, o );

            }
            break;

            case  CRBModelMode_heatns::CRB:
            {
                std::cout << "-------->CRBModelMode_heatns::CRB\n";
                boost::timer ti;
                
                //model->solve( mu );
                std::vector<double> ofem = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );
                ti.restart();
                auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() );

                double relative_error = std::abs( ofem[0]-o.template get<0>() ) /ofem[0];
                double relative_estimated_error = o.template get<1>() / ofem[0];
                double condition_number = o.template get<3>();
                
                if ( crb->errorType()==2 )
                {

                    std::vector<double> v = boost::assign::list_of( ofem[0] )( ofem[1] )( o.template get<0>() )( relative_estimated_error )( ti.elapsed() )( relative_error )( condition_number );
                    std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions\n";
                    std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) % o.template get<2>() ).str().c_str() ,std::ios::out | std::ios::app );
                    printEntry( file_summary_of_simulations, mu, v );
                    printEntry( ostr, mu, v );
                    file_summary_of_simulations.close();
                }

                else
                {
                    std::vector<double> v = boost::assign::list_of( ofem[0] )( ofem[1] )( o.template get<0>() )( relative_estimated_error )( ti.elapsed() ) ( relative_error )( condition_number ) ;
                    std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions  (relative error estimation on this output : " << relative_estimated_error<<") \n";
                    std::ofstream file_summary_of_simulations( ( boost::format( "summary_of_simulations_%d" ) % o.template get<2>() ).str().c_str() ,std::ios::out | std::ios::app );
                    printEntry( file_summary_of_simulations, mu, v );
                    printEntry( ostr, mu, v );
                    file_summary_of_simulations.close();
                }

            }
            break;

            case  CRBModelMode_heatns::CRB_ONLINE:
            {
                std::cout << "-------->CRBModelMode_heatns::CRB_ONLINE\n";
                boost::timer ti;
                ti.restart();
                auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() );

                if ( crb->errorType()==2 )
                {
                    std::vector<double> v = boost::assign::list_of( o.template get<0>() )( ti.elapsed() );
                    std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions\n";
                    printEntry( ostr, mu, v );
                }

                else
                {
                    std::vector<double> v = boost::assign::list_of( o.template get<0>() )( o.template get<1>() )( ti.elapsed() );
                    std::cout << "output=" << o.template get<0>() << " with " << o.template get<2>() << " basis functions  (error estimation on this output : " << o.template get<1>()<<") \n";
                    printEntry( ostr, mu, v );
                }

            }
            break;

            case  CRBModelMode_heatns::SCM:
            {
                std::cout << "-------->CRBModelMode_heatns::SCM\n";
                auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                printEntry( ostr, mu, o );
            }
            break;

            case  CRBModelMode_heatns::SCM_ONLINE:
            {
                std::cout << "-------->CRBModelMode_heatns::SCM_ONLINE\n";
                auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                printEntry( ostr, mu, o );
            }
            break;

            }

            std::cout << "------------------------------------------------------------\n";
        }
        std::cout << ostr.str() << "\n";

    }
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
    {

        std::cout << "\n-------->OpusApp_heatns::run( const double * X, unsigned long N, double * Y, unsigned long P ). \n \n";
        
        switch ( M_mode )
        {
        case  CRBModelMode_heatns::PFEM:
        {
            model->run( X, N, Y, P );
        }
        break;

        case  CRBModelMode_heatns::SCM:
        case  CRBModelMode_heatns::SCM_ONLINE:
        {
            crb->scm()->run( X, N, Y, P );
        }
        break;

        case  CRBModelMode_heatns::CRB:
        case  CRBModelMode_heatns::CRB_ONLINE:
        {
            crb->run( X, N, Y, P );
        }
        break;
        }

        fs::current_path( M_current_path );
    }

private:
    int printParameterHdr( std::ostream& os, int N, std::vector<std::string> outputhdrs )
    {
        for ( int i = 0; i < N; ++i )
        {
            std::ostringstream s;
            s << "mu" << i;
            os  << hdrmanip( prec+7 ) << s.str();
        }

        BOOST_FOREACH( auto output, outputhdrs )
        {
            os << hdrmanip( 15 ) << output;
        }
        os << "\n";

        return N*( prec+7 )+outputhdrs.size()*15;
    }
    void printEntry( std::ostream& os,
                     typename ModelType::parameter_type const& mu,
                     std::vector<double> const& outputs )
    {
        for ( int i = 0; i < mu.size(); ++i )
            os  << std::right <<std::setw( prec+7 ) << dmanip << mu[i];

        BOOST_FOREACH( auto o, outputs )
        {
            os << tabmanip( 15 ) << o;
        }
        os << "\n";
    }

private:
    CRBModelMode_heatns M_mode;
    crbmodel_ptrtype model;
    crb_ptrtype crb;

    fs::path M_current_path;
}; // OpusApp
} // Feel

#endif /* __OpusApp_heatns_H */

