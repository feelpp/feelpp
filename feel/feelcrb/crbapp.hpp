/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-03-22

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file opuseadscrbapp.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-22
 */
#ifndef __OpusEadsCrbApp_H
#define __OpusEadsCrbApp_H 1

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>


namespace Feel
{
/**
 * \class CRBApp
 * \brief Certified Reduced Basis application
 *
 * @author Christophe Prud'homme
 */
template<typename ModelType>
class CRBApp   : public Application
{
    typedef Application super;
public:

    typedef CRBModel<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;
    typedef CRB<crbmodel_type> crb_type;
    typedef boost::shared_ptr<crb_type> crb_ptrtype;

    CRBApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, crbOptions().add( od ) )
    {
        this->init();
    }

    CRBApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, crbOptions().add( od ) )
    {
        this->init();
    }
    void init()
    {
        std::srand( static_cast<unsigned>( std::time( 0 ) ) );
        std::cerr << "[CRBApp] constructor " << this->about().appName()  << std::endl;

        if ( this->vm().count( "crb.output-index" ) )
            this->changeRepository( boost::format( "%1%/h_%2%/s_%3%" )
                                    % this->about().appName()
                                    % this->vm()["hsize"].template as<double>()
                                    % this->vm()["crb.output-index"].template as<int>()
                                  );

        else
            this->changeRepository( boost::format( "%1%/h_%2%/" )
                                    % this->about().appName()
                                    % this->vm()["hsize"].template as<double>()
                                  );

        std::cerr << "[CRBApp] ch repo" << std::endl;
        this->setLogs();
        std::cerr << "[CRBApp] set Logs" << std::endl;
        opus = crbmodel_ptrtype( new crbmodel_type( this->vm() ) );
        std::cerr << "[CRBApp] get model done" << std::endl;
        crb = crb_ptrtype( new crb_type( this->about().appName(),
                                         this->vm() ) );

        std::cerr << "[CRBApp] get crb done" << std::endl;
        crb->setTruthModel( opus );
        std::cerr << "[CRBApp] constructor done" << std::endl;

    }
    void setOutput( int i = 0, CRBErrorType error_type = ( int )CRB_RESIDUAL , int maxiter = 10 )
    {
        auto  ckconv = crb->offline();

        if ( ckconv.size() )
        {
#if 0
            double max_ei=0;
            double min_ei=0;
            crb->computeErrorEstimationEfficiencyIndicator ( opus->parameterSpace(), max_ei, min_ei,10 );
#endif

            for ( auto it = ckconv.left.begin(); it != ckconv.left.end(); ++it )
            {
                LOG(INFO) << "ckconv[" << it->first <<"]=" << it->second << "\n";
            }
        }
    }
    void run()
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        if ( !crb->isDBLoaded() )
        {
            std::cout << "No DB available, do offline computations first...\n";
            crb->offline();
        }

        typename crb_type::sampling_ptrtype Sampling( new typename crb_type::sampling_type( opus->parameterSpace() ) );
        Sampling->randomize( 10 );
        int crb_error_type = crb->errorType();
        int output_index = crb->outputIndex();
        BOOST_FOREACH( auto mu, *Sampling )
        {
            double sfem = opus->output( output_index, mu );
            int size = mu.size();
            std::cout << "------------------------------------------------------------\n";
            std::cout << "tolerance : " << this->vm()["crb.online-tolerance"].template as<double>() << "\n";
            std::cout << "mu = [ ";

            for ( int i=0; i<size-1; i++ ) std::cout<< mu[i] <<" , ";

            std::cout<< mu[size-1]<<" ] \n";
            auto o = crb->run( mu,  this->vm()["crb.online-tolerance"].template as<double>() );

            if ( crb_error_type==2 )
            {
                std::cout << "output=" << o.get<0>() << " with " << o.get<2>() << " basis functions\n";
                std::cout << "output obtained using FEM : "<<sfem<<std::endl;
            }

            else
            {
                std::cout << "output=" << o.get<0>() << " with " << o.get<2>() << " basis functions  (error estimation on this output : " << o.get<1>()<<") \n";
                std::cout << "output obtained using FEM : "<<sfem<<std::endl;
            }

            std::cout << "------------------------------------------------------------\n";
        }
    }
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
    {
        crb->run( X, N, Y, P );
    }

    crbmodel_ptrtype opus;
    crb_ptrtype crb;
}; // CRBApp
} // Feel

#endif /* __CrbApp_H */
