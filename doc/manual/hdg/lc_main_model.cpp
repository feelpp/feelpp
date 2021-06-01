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

#include "lc_model.hpp"

using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "laminacribrosa" ,
                     "laminacribrosa" ,
                     "0.1",
                     "Lamina-Cribrosa-Model",
                     AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Christophe Trophime", "developer", "christophe.trophime@lncmi.cnrs.fr", "" );
    about.addAuthor( "Lorenzo Sala", "developer", "", "" );
    return about;

}

int main(int argc, char *argv[])
{
    using namespace Feel;
    Feel::Environment env( _argc=argc,
                           _argv=argv,
                           _about=makeAbout(),
                           _desc=FeelModels::makeMixedPoissonOptions("", "mixedpoisson"),             
                           _desc_lib=FeelModels::makeMixedPoissonLibOptions("", "mixedpoisson").add(feel_options())
    					   );

    typedef FeelModels::LaminaCribrosa<FEELPP_DIM,FEELPP_ORDER, FEELPP_G_ORDER> lc_type;
	// typedef FeelModels::MixedPoisson<FEELPP_DIM,FEELPP_ORDER, FEELPP_G_ORDER> mp_type;
	
	auto LC = lc_type::New("mixedpoisson");
	
    auto mesh = loadMesh( _mesh=new lc_type::mesh_type );
    decltype( IPtr( _domainSpace=Pdh<FEELPP_ORDER>(mesh), _imageSpace=Pdh<FEELPP_ORDER>(mesh) ) ) Idh ;
    decltype( IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(mesh), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) ) ) Idhv;
    
	if ( soption( "gmsh.submesh" ).empty() )
    {
        LC -> init( mesh );
    }
    else
    {
		Feel::cout << "Using submesh: " << soption("gmsh.submesh") << std::endl;
        auto cmesh = createSubmesh( mesh, markedelements(mesh,soption("gmsh.submesh")) );
        Idh = IPtr( _domainSpace=Pdh<FEELPP_ORDER>(cmesh), _imageSpace=Pdh<FEELPP_ORDER>(mesh) );
        Idhv = IPtr( _domainSpace=Pdhv<FEELPP_ORDER>(cmesh), _imageSpace=Pdhv<FEELPP_ORDER>(mesh) );
        LC -> init( cmesh, mesh );
    }
    

    if ( LC -> isStationary() )
    {
        Feel::cout << " This model could not be stationary! " << std::endl;
        return -1;
    }
    else
    {
#ifndef USE_SAME_MAT
        LC -> assembleCstPart();
#endif
        for ( ; !(LC->timeStepBase()->isFinished()) ; LC->updateTimeStepBDF() )
        {
#ifdef USE_SAME_MAT
			LC->assembleCstPart();
#endif		
            Feel::cout << "============================================================\n";
            Feel::cout << "time simulation: " << LC->time() << "s \n";
            Feel::cout << "============================================================\n";
	        // 3D model
            LC->assembleNonCstPart();
            LC->solve();
	        // OD model
	        LC->second_step();

            LC->exportResults( mesh , Idh, Idhv );
        }
    }

    return 0;
}



