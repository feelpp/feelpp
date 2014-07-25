/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@imag.fr>
       Date: 2014-07-22

  Copyright (C) 2014 Université de Strasbourg

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
   \file explicit.cpp
   \author Abdoulaye Samake <abdoulaye.samake@imag.fr>
   \date 2014-07-22
 */
#include <feel/feel.hpp>

namespace Feel
{
template<int Dim, int Order>
class Explicit : public Simget
{
    typedef Simget super;
public:
    /**
     * Constructor
     */
    Explicit()
        :
        super()
    {
    }

    void run();
}; // Explicit

template<int Dim,int Order>
void
Explicit<Dim,Order>::run()
{
    Environment::changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                                   % this->about().appName()
                                   % Dim
                                   % Order
                                   % doption("gmsh.hsize") );

    auto commWorld = Environment::worldComm();
    auto commSelf = Environment::worldComm().subWorldCommSeq();

    auto mesh = createGMSHMesh( _mesh=new Mesh<Simplex<Dim>>(commSelf),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER|MESH_PROPAGATE_MARKERS,
                                _desc=domain( _name=(boost::format( "domain-%1%" ) % commWorld.globalRank() ).str() ,
                                              _addmidpoint=false,
                                              _usenames=false,
                                              _shape="hypercube",
                                              _dim=Dim,
                                              _h=doption("gmsh.hsize"),
                                              _convex="Simplex",
                                              _substructuring=true
                                              ),
                                _structured=2,
                                _partitions=commSelf.localSize(),
                                _worldcomm=commSelf );



    auto Xh  = FunctionSpace<Mesh<Simplex<Dim>>, bases<Lagrange<Order,Scalar>>>::New( _mesh=mesh,
                                                                                      _worldscomm=Environment::worldsCommSeq(1) );

    /* CrossPoints */
    std::string str = "CrossPoints";
    auto dftcrp = Xh->dof()->markerToDof(str);
    std::cout<<"Number of degrees of freedom for CrossPoints= "<< std::distance(dftcrp.first,dftcrp.second) <<"\n";

    /* WireBasket */
    str = "WireBasket";
    auto dftwrb = Xh->dof()->markerToDof(str);
    std::cout<<"Number of degrees of freedom for WireBasket= "<< std::distance(dftwrb.first,dftwrb.second) <<"\n";

} // Explicit::run

} // Feel

int main(int argc, char** argv) {
    using namespace Feel;
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="doc_explicit",
                                  _author="Abdoulaye Samake",
                                  _email="abdoulaye.samake@imag.fr") );
    Application app;
    app.add(new Explicit<3,1>());
    app.run();
}
