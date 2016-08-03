/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 26 Feb 2016

 Copyright (C) 2016 Feel++ Consortium

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

#include <lc_model2.hpp>

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "laminacribrosa" ,
                     "laminacribrosa" ,
                     "0.1",
                     "Lamina cribrosa Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "" );
    about.addAuthor( "Daniele Prada", "developer", "", "" );
    return about;

}

int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=FeelModels::makeMixedPoissonOptions("mixedpoisson"),             
                           _desc_lib=FeelModels::makeMixedPoissonLibOptions("mixedpoisson").add(feel_options())
                           /*_desc=FeelModels::makeLCHDGOptions(),
                           _desc_lib=FeelModels::makeLCHDGLibOptions().add(feel_options())*/ );

    typedef FeelModels::LaminaCribrosa<FEELPP_DIM,FEELPP_ORDER> lc_type;

    lc_type LC ;
    auto mesh = loadMesh( _mesh=new lc_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<1>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<1>(mesh) ) ) Idhv;
    if ( soption( "mixedpoisson.gmsh.submesh" ).empty() )
        LC.init(mesh, 1, 1);
    else
    {
        auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("mixedpoisson.gmsh.submesh")), Environment::worldComm() );
        Idh = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmesh), _imageSpace=Pdh<1>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<1>(mesh) );
        LC.init( cmesh, 1, 1, mesh );
    }
    
    if ( LC.isStationary() )
    {
        LC.solve();
	LC.second_step();    
        LC.exportResults( mesh );
    }
    else
    {
        for ( ; !( LC.timeStepBase()->isFinished() && LC.timeStepBase_statevar()->isFinished() ) ; LC.updateTimeStepBDF() )
        {
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << LC.time() << "s \n";
            Feel::cout << "============================================================\n";
		    // 3D model
            LC.solve();
		    // OD model
		    LC.second_step();
            LC.exportResults( mesh, Idh, Idhv );
        }
    }

 
    return 0;
}



