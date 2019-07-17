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

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

}
int main( int argc, char** argv )
{
    double T=5;

    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     //_desc=laplacianoptions,
                     _about=about(_name="bdfpod",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));



    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    auto Xh = Pch<1>( mesh );
    auto u = Xh->element();
    auto mybdf = bdf( _space=Xh, _name="bdfpod", _initial_time=0., _final_time=T, _time_step=1.,_rank_proc_in_files_name=1 );

    for ( mybdf->start(); mybdf->isFinished() == false; mybdf->next() )
    {
        if ( mesh->worldComm().isMasterRank() )
            std::cout << "mybdf iteration " << mybdf->iteration() << " with time " << mybdf->time() << "\n";
        u.on( _range=elements( mesh ), _expr=cst( mybdf->time() ) );
        mybdf->shiftRight( u );
    }

    if ( mesh->worldComm().isMasterRank() )
    {
        std::cout << " -- Ndof = "  << Xh->nDof() << "\n";
        std::cout << " -- K = "  << mybdf->timeValues().size() << "\n";
    }
    mybdf->setRestart( true );


    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> pod( mybdf->timeValues().size()-1, mybdf->timeValues().size()-1 );
    auto bdfi = mybdf->deepCopy();
    auto bdfj = mybdf->deepCopy();

    for ( bdfi->start(); !bdfi->isFinished(); bdfi->next() )
    {
        int i = bdfi->iteration() - 1;
        bdfi->loadCurrent();
        if ( mesh->worldComm().isMasterRank() )
            std::cout << "====================\n"
                      << "element i="  << i+1 << " with time " << bdfi->time() << "\n";

        for ( bdfj->start(); !bdfj->isFinished() && ( bdfj->iteration() < bdfi->iteration() ); bdfj->next() )
        {
            int j = bdfj->iteration() -1;
            bdfj->loadCurrent();

            if ( mesh->worldComm().isMasterRank() )
                std::cout << "element j="  << j+1 <<  "\n";
            pod( i,j ) = inner_product( bdfj->unknown( 0 ), bdfi->unknown( 0 ) );
            pod( j,i ) = pod( i,j );
        }
        pod( i,i ) = inner_product( bdfi->unknown( 0 ), bdfi->unknown( 0 ) );
    }


    if ( mesh->worldComm().isMasterRank() )
        std::cout << "pod=\n" << pod << "\n";

    Eigen::EigenSolver< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > eigenSolver;
    eigenSolver.compute( pod );
    auto eigenValues = eigenSolver.eigenvalues();
    //double myeigen = real(eigenValues[k]);
    auto eigenVectors = eigenSolver.eigenvectors();
    if ( mesh->worldComm().isMasterRank() )
        std::cout << "eigenValues=\n" << eigenValues << "\n";

    return 0;
}
