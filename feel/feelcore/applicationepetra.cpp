/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-04-25

  Copyright (C) 2006, 2009 EPFL

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
   \file application.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-04-25
 */
#include <boost/concept_check.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

namespace Feel
{
#if defined( FEELPP_HAS_TRILINOS_EPETRA )

po::options_description
makeMLEpetraOptions()
{
    po::options_description mlEpetraOptions( "ML Epetra options" );
    mlEpetraOptions.add_options()
    ( "max-levels", po::value<int>()->default_value( 6 ), "max levels number" )
    ( "aggregation-type", po::value<std::string>()->default_value( "MIS" ), "aggregation type" )
    ( "coarse-type", po::value<std::string>()->default_value( "Amesos_KLU" ), "coarse type" )
    ( "increasing-or-decreasing", po::value<std::string>()->default_value( "decreasing" ), "increasing or decreasing level indices" );
    return mlEpetraOptions;
}
FEELPP_NO_EXPORT
po::options_description
epetraOptions()
{
    po::options_description epetra( "EPETRA options" );
    epetra.add_options()
    ( "disable-epetra", "disable epetra" )
    ;
    return epetra;
}

void
Application::init( MPI_Comm& _comm )
{
    if ( _S_is_Initialized == false )
    {
        _S_comm = boost::shared_ptr<comm_type>( new Epetra_MpiComm( _comm ) );

        _S_is_Initialized = true;
    }
}


#if defined(FEELPP_HAS_MPI)
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          MPI_Comm comm )
    :
    super( argc, argv, ad, epetraOptions().add( makeMLEpetraOptions() ), comm )
{
    init( comm );
}


#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad )
    :
    super( argc, argv, ad, epetraOptions().add( makeMLEpetraOptions() ) )
{
    init( comm );
}
#endif /* FEELPP_HAS_MPI */


#if defined(FEELPP_HAS_MPI)
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od,
                          MPI_Comm comm )
    :
    super( argc, argv, ad, epetraOptions().add( makeMLEpetraOptions() ).add( od ), comm )
{
    //cout << "hallo FEELPP_HAS_MPI2\n" << endl;
    init( comm );
}

#else
Application::Application( int argc,
                          char** argv,
                          AboutData const& ad,
                          po::options_description const& od )
    :
    super( argc, argv, ad, epetraOptions().add( makeMLEpetraOptions() ).add( od ) )
{
    init( comm );
}



#endif /* FEELPP_HAS_MPI */


Application::~Application()
{
}

boost::shared_ptr<Epetra_MpiComm> Application::_S_comm;

bool Application::_S_is_Initialized = false;

#endif // FEELPP_HAS_TRILINOS_EPETRA



}

