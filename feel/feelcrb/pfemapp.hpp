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
   \file pfemapp.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-22
 */
#ifndef PFEMAPP_HPP
#define PFEMAPP_HPP

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>


#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feelcrb/pfemapp.hpp>


namespace Feel
{
po::options_description pfemapp_options( std::string const& prefix );
std::string _o( std::string const& prefix, std::string const& opt )
{
    std::string o = prefix;

    if ( !o.empty() )
        o += "-";

    return o + opt;
}
/**
 * \class PFemApp
 * \brief Parametrized Finite Element Method Application
 *
 * @author Christophe Prud'homme
 */
template<typename ModelType>
class PFemApp   : public Application
{
    typedef Application super;
public:

    typedef CRBModel<ModelType> crbmodel_type;
    typedef boost::shared_ptr<crbmodel_type> crbmodel_ptrtype;


    PFemApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, pfemapp_options( ad.appName() ).add( od ) )
    {
        this->init();
    }

    PFemApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, pfemapp_options( ad.appName() ).add( od ) )
    {
        this->init();
    }
    void init()
    {
        std::srand( static_cast<unsigned>( std::time( 0 ) ) );
        this->setLogs();
        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].template as<double>()
                              );
        std::cout << "[PFemApp] build model " << this->about().appName() << "\n";
        model = crbmodel_ptrtype( new crbmodel_type( this->vm() ) );

        std::cout << "build model " << this->about().appName() << " done\n";
    }
    void run()
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        typename crbmodel_type::parameter_type mu( model->parameterSpace() );
        int mutype = this->vm()[_o( this->about().appName(),"mu-type" )].template as<int>();

        if ( mutype <= 0 )
        {
            bool broadcast = true;
            mu = crbmodel_type::parameterspace_type::logRandom( model->parameterSpace(), broadcast );
        }

        else if ( mutype == 1 )
        {
            mu = model->parameterSpace()->min();
        }

        else if ( mutype == 2 )
        {
            //mu = model->parameterSpace()->middle();
            mu = model->parameterSpace()->max();
        }

        else
        {
            mu = model->parameterSpace()->max();
        }

        std::cout << "[PFemApp] running " << this->about().appName() << " with mu = [";
        int size = mu.size();

        for ( int i=0; i<size-1; i++ ) std::cout<<mu( i )<<" , ";

        std::cout<< mu( size-1 )<<"] "<<std::endl;

        model->solve( mu );

        auto Xh = model->functionSpace();
        auto u = Xh->element();
        bool need_to_solve=true;

        for ( int l =0; l < model->Nl(); ++l )
        {
            double o = model->output( l,mu , u , need_to_solve );
            std::cout << "[PFemApp] output " << l << " of  " << this->about().appName() << " = " << o << "\n";
        }
    }

    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
    {

        model->run( X, N, Y, P );

    }

    crbmodel_ptrtype model;


}; // Opus



} // Feel

#endif
