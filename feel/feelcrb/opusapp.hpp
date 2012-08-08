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
#ifndef __OpusApp_H
#define __OpusApp_H 1

#include <boost/assign/std/vector.hpp>
#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <boost/serialization/version.hpp>

namespace Feel
{
po::options_description opusapp_options( std::string const& prefix );
std::string _o( std::string const& prefix, std::string const& opt )
{
    std::string o = prefix;

    if ( !o.empty() )
        o += ".";

    return o + opt;
}


enum class SamplingMode
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
class OpusApp   : public Application
{
    typedef Application super;
public:

    typedef CRBModel<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    OpusApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
        {
            this->init();
        }

    OpusApp( AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( mode )
        {
            this->init();
        }

    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
        M_mode( ( CRBModelMode )this->vm()[_o( this->about().appName(),"run.mode" )].template as<int>() )
        {
            this->init();
        }
    OpusApp( int argc, char** argv, AboutData const& ad, po::options_description const& od, CRBModelMode mode )
        :
        super( argc, argv, ad, opusapp_options( ad.appName() ).add( od ).add( crbOptions() ).add( feel_options() ) ),
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
                Debug() << "[OpusApp] constructor " << this->about().appName()  << "\n";

                if( vm().count("hsize") && !vm().count("geofile") )
                    {
                        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                                % this->about().appName()
                                                % this->vm()["hsize"].template as<double>()
                                                );
                    }
                if( vm().count("geofile") )
                    {
                        this->changeRepository( boost::format( "%1%/%2%/" )
                                                % this->about().appName()
                                                % this->vm()["geofile"].template as<std::string>()
                                                );
                    }

                Debug() << "[OpusApp] ch repo" << "\n";
                this->setLogs();
                Debug() << "[OpusApp] set Logs" << "\n";
                Debug() << "[OpusApp] mode:" << ( int )M_mode << "\n";
                model = crbmodel_ptrtype( new crbmodel_type( this->vm(),M_mode ) );
                Debug() << "[OpusApp] get model done" << "\n";


                crb = crb_ptrtype( new crb_type( this->about().appName(),
                                                 this->vm() ) );
                Debug() << "[OpusApp] get crb done" << "\n";
                crb->setTruthModel( model );
                Debug() << "[OpusApp] constructor done" << "\n";
            }

            catch ( boost::bad_any_cast const& e )
            {
                std::cout << "[OpusApp] a bad any cast occured, probably a nonexistant or invalid  command line/ config options\n";
                std::cout << "[OpusApp] exception reason: " << e.what() << "\n";
            }

        }

    void setMode( std::string const& mode )
        {
            if ( mode == "pfem" ) M_mode = CRBModelMode::PFEM;

            if ( mode == "crb" ) M_mode = CRBModelMode::CRB;

            if ( mode == "scm" ) M_mode = CRBModelMode::SCM;

            if ( mode == "scm_online" ) M_mode = CRBModelMode::SCM_ONLINE;

            if ( mode == "crb_online" ) M_mode = CRBModelMode::CRB_ONLINE;
        }
    void setMode( CRBModelMode mode )
        {
            M_mode = mode;
        }

    void loadDB()
        {
            if ( M_mode == CRBModelMode::PFEM )
                return;

            if ( !crb->scm()->isDBLoaded() || crb->scm()->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::SCM )
                {
                    std::cout << "No SCM DB available, do scm offline computations first...\n";
                    if( crb->scm()->doScmForMassMatrix() )
                        crb->scm()->setScmForMassMatrix( true );

                    crb->scm()->offline();
                }
            }

            if ( !crb->isDBLoaded() || crb->rebuildDB() )
            {
                if ( M_mode == CRBModelMode::CRB )
                    //|| M_mode == CRBModelMode::SCM )
                {
                    std::cout << "No CRB DB available, do crb offline computations...\n";
                    crb->offline();
                }

                else if ( M_mode != CRBModelMode::SCM )
                    throw std::logic_error( "CRB/SCM Database could not be loaded" );
            }
        }

    FEELPP_DONT_INLINE
    void run()
        {
            if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }

            this->loadDB();
            int run_sampling_size = this->vm()[_o( this->about().appName(),"run.sampling.size" )].template as<int>();
            SamplingMode run_sampling_type = ( SamplingMode )this->vm()[_o( this->about().appName(),"run.sampling.mode" )].template as<int>();
            int output_index = this->vm()["crb.output-index"].template as<int>();
            //int output_index = this->vm()[_o(this->about().appName(),"output.index")].template as<int>();

            typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( model->parameterSpace() ) );

            switch ( run_sampling_type )
            {
            default:
            case SamplingMode::RANDOM:
                Sampling->randomize( run_sampling_size  );
                break;

            case SamplingMode::EQUIDISTRIBUTED:
                Sampling->equidistribute( run_sampling_size  );
                break;
            }

            std::map<CRBModelMode,std::vector<std::string> > hdrs;
            using namespace boost::assign;
            std::vector<std::string> pfemhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" );
            std::vector<std::string> crbhdrs = boost::assign::list_of( "FEM Output" )( "FEM Time" )( "RB Output" )( "Error Bounds" )( "CRB Time" )( "Relative error PFEM/CRB" )( "Conditionning" );
            std::vector<std::string> scmhdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "FEM" )( "FEM Time" )( "Rel.(FEM-Lb)" );
            std::vector<std::string> crbonlinehdrs = boost::assign::list_of( "RB Output" )( "Error Bounds" )( "CRB Time" );
            std::vector<std::string> scmonlinehdrs = boost::assign::list_of( "Lb" )( "Lb Time" )( "Ub" )( "Ub Time" )( "Rel.(FEM-Lb)" );
            hdrs[CRBModelMode::PFEM] = pfemhdrs;
            hdrs[CRBModelMode::CRB] = crbhdrs;
            hdrs[CRBModelMode::SCM] = scmhdrs;
            hdrs[CRBModelMode::CRB_ONLINE] = crbonlinehdrs;
            hdrs[CRBModelMode::SCM_ONLINE] = scmonlinehdrs;
            std::ostringstream ostr;

            if( crb->printErrorDuringOfflineStep() )
                crb->printErrorsDuringRbConstruction();
            if ( crb->showMuSelection() )
                crb->printMuSelection();

            auto exporter = Exporter<typename crbmodel_type::mesh_type>::New( "ensight" );
            exporter->step( 0 )->setMesh( model->functionSpace()->mesh() );

            printParameterHdr( ostr, model->parameterSpace()->dimension(), hdrs[M_mode] );
            BOOST_FOREACH( auto mu, *Sampling )
            {

                int size = mu.size();
                std::cout << "mu = [ ";

                for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";

                std::cout<< mu[size-1]<<" ] ";

                std::ostringstream mu_str;
                for ( int i=0; i<size-1; i++ ) mu_str << std::scientific << std::setprecision( 5 ) << mu[i] <<",";
                mu_str << std::scientific << std::setprecision( 5 ) << mu[size-1];

                switch ( M_mode )
                {
                case  CRBModelMode::PFEM:
                {
                    boost::timer ti;
                    model->solve( mu );
                    std::vector<double> o = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );
                    std::cout << "output=" << o[0] << "\n";
                    printEntry( ostr, mu, o );

                }
                break;

                case  CRBModelMode::CRB:
                {
                    std::cout << "CRB mode\n";
                    LOG(INFO) << "solve u_fem\n";
                    google::FlushLogFiles(google::GLOG_INFO);
                    boost::timer ti;
                    auto u_fem = model->solve( mu );
                    std::ostringstream u_fem_str;
                    u_fem_str << "u_fem(" << mu_str.str() << ")";
                    u_fem.setName( u_fem_str.str()  );

                    LOG(INFO) << "compute output\n";
                    google::FlushLogFiles(google::GLOG_INFO);

                    exporter->step(0)->add( u_fem.name(), u_fem );

                    std::vector<double> ofem = boost::assign::list_of( model->output( output_index,mu ) )( ti.elapsed() );
                    ti.restart();
                    LOG(INFO) << "solve crb\n";
                    google::FlushLogFiles(google::GLOG_INFO);
                    auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() );

                    if( this->vm()["crb.rebuild-database"].template as<bool>() )
                        {
                            auto u_crb = crb->expansion( mu );
                            std::ostringstream u_crb_str;
                            u_crb_str << "u_crb(" << mu_str.str() << ")";
                            u_crb.setName( u_crb_str.str()  );
                            exporter->step(0)->add( u_crb.name(), u_crb );
                        }

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

                case  CRBModelMode::CRB_ONLINE:
                {
                    std::cout << "CRB Online mode\n";
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

                case  CRBModelMode::SCM:
                {
                    std::cout << "SCM mode\n";
                    auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                    printEntry( ostr, mu, o );
                }
                break;

                case  CRBModelMode::SCM_ONLINE:
                {
                    std::cout << "SCM online mode\n";
                    auto o = crb->scm()->run( mu, crb->scm()->KMax() );
                    printEntry( ostr, mu, o );
                }
                break;

                }

                std::cout << "------------------------------------------------------------\n";
            }
            exporter->save();
            std::cout << ostr.str() << "\n";

        }
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
        {

            switch ( M_mode )
            {
            case  CRBModelMode::PFEM:
            {
                model->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::SCM:
            case  CRBModelMode::SCM_ONLINE:
            {
                crb->scm()->run( X, N, Y, P );
            }
            break;

            case  CRBModelMode::CRB:
            case  CRBModelMode::CRB_ONLINE:
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
    CRBModelMode M_mode;
    crbmodel_ptrtype model;
    crb_ptrtype crb;

    fs::path M_current_path;
}; // OpusApp
} // Feel

#endif /* __OpusApp_H */

