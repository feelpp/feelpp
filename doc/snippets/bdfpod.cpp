/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-03-15

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
   \file dofpoints.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-03-15
 */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

}
int main( int argc, char** argv )
{
    double hsize = 2;
    double T=2;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "hsize", po::value<double>( &hsize )->default_value( 2 ), "h size" )
    ( "Tfinal", po::value<double>( &T )->default_value( 2 ), "final time" )
    ;
    desc.add( Feel::feel_options() );
    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );

    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<2,4> > mesh_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name="ellipsoid-2",
                                        _usenames=true,
                                        _shape="ellipsoid",
                                        _dim=2,
                                        _order=4,
                                        _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    typedef FunctionSpace<mesh_type,bases<Lagrange<1,Scalar> > > space_type;
    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto mybdf = bdf( _space=Xh, _vm=vm, _name="bdfpod", _initial_time=0., _final_time=T, _time_step=1. );

    for ( mybdf->start(); mybdf->isFinished() == false; mybdf->next() )
    {
        std::cout << "   o t=" << mybdf->time() << "\n";
        u = vf::project( _space=Xh, _range=elements( mesh ), _expr=cst( mybdf->time() ) );
        mybdf->shiftRight( u );
        std::cout << "element iter="  << mybdf->iteration() << " :"  << mybdf->unknown( 0 ) << "\n";
    }

    std::cout << " -- Ndof = "  << Xh->nLocalDof() << "\n";
    std::cout << " -- K = "  << mybdf->timeValues().size() << "\n";
    mybdf->setRestart( true );
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pod( mybdf->timeValues().size(), mybdf->timeValues().size() );
    auto bdfi = mybdf->deepCopy();
    auto bdfj = mybdf->deepCopy();

    for ( bdfi->start(); !bdfi->isFinished(); bdfi->next() )
    {
        int i = bdfi->iteration()-1;
        bdfi->loadCurrent();
        std::cout << "====================\n";
        std::cout << "element i="  << i << ":" << bdfi->unknown( 0 ) << "\n";
        std::cout << "====================\n";
        std::cout << "--------------------\n";

        for ( bdfj->start(); !bdfj->isFinished() && ( bdfj->iteration() < bdfi->iteration() ); bdfj->next() )
        {
            int j = bdfj->iteration()-1;
            bdfj->loadCurrent();
            std::cout << "element i="  << i << ":" << bdfi->unknown( 0 ) << "\n";
            std::cout << "element j="  << j << ":\n"  << bdfj->unknown( 0 ) << "\n";
            pod( i,j ) = inner_product( bdfj->unknown( 0 ), bdfi->unknown( 0 ) );
            pod( j,i ) = pod( i,j );

        }

        pod( i,i ) = inner_product( bdfi->unknown( 0 ), bdfi->unknown( 0 ) );
    }

    std::cout << "pod=" << pod << "\n";

}
