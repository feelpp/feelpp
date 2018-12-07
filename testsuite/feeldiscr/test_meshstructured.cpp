/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

   Author(s): Thomas Lantz
   Date: 2015-04-27

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
   \file test_multiscaleimage.cpp
   \author Thomas Lantz
   \date 2015-04-27
 */

#define USE_BOOST_TEST 1
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_meshstructured
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelvf/msi.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/norml2.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/operations.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelpoly/multiscalequadrature.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/hbf.hpp>
#if defined(FEELPP_HAS_FFTW)
#include <fftw3.h>
#endif
#include <feel/feeldiscr/meshstructured.hpp>

namespace Feel
{

inline
AboutData
makeAbout()
{
    AboutData about( "test_meshstructured" ,
                     "test_meshstructured" ,
                     "8.9e-3",
                     "Test structured mesh with image",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2014 Feel++ Consortium" );

    about.addAuthor( "Thomas Lantz", "student", "", "" );
    return about;
}

template < typename T  >
using holo3_image = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> ;

class Test_MeshStructured
{
    public :


    Test_MeshStructured(int nx,int ny)
    {
        ima = holo3_image<float>(ny,nx);
        ima2 = holo3_image<float>(ny,nx);
        for (int i=0;i<nx;i++)
        {
            for (int j=0;j<ny;j++)
            {
                ima(j,i)=1;
                ima2(j,i)=0;
            }
        }
    }
    
    
    void run (double pixelsize) 
        {
            tic();
            auto mesh = std::make_shared<MeshStructured>( ima.rows(), ima.cols(), pixelsize, ima,ima,
                                                          Environment::worldCommPtr(),"", false, false);
            
            mesh->components().reset();
            mesh->components().set( size_type(MESH_NO_UPDATE_MEASURES) );
            mesh->updateForUse();
            toc("mesh");
    
            auto Vh = Pch<1> ( mesh ) ;
            auto u = Vh->element () ;
            auto v = Vh->element () ;
     
            auto imaP = vf::project(Vh,elements(mesh),Px()); 
            auto l = form1( _test=Vh );
            l = integrate(_range=elements(mesh),
                          _expr=id(v)*idv(imaP));
                     
            auto a = form2( _trial=Vh, _test=Vh);
            a = integrate( _range=elements(mesh),
                           _expr=idt(u)*id(v));
                      
            a.solve(_rhs=l,_solution=u) ;

            auto e = exporter( _mesh=mesh );
            e->add( "u", u );
            e->addRegions(); 
            e->save();

            auto n = normL2 ( _range=elements( mesh ), _expr=idv(u)-idv(imaP));
            std::cout << " Norm2 " << n << std::endl;

        }

    void runImage  (double pixelsize) 
        {

            auto mesh = std::make_shared<MeshStructured>( ima.rows(), ima.cols(), pixelsize, ima,ima,
                                                          Environment::worldCommPtr(),"", false, false);
            mesh->components().reset();
            mesh->components().set( size_type(MESH_NO_UPDATE_MEASURES) );
            mesh->updateForUse();
    
            auto Vh = Pch<1> ( mesh ) ;
            auto u = Vh->element () ;
            auto v = Vh->element () ;
     
            auto imaP = vf::project(Vh,elements(mesh),cst(1.)); 
            auto imaP2 = vf::project(Vh,elements(mesh),cst(0.)); 
            auto l = form1( _test=Vh );
            l = integrate(_range=elements(mesh),
                          //_expr=grad(v)*vec(idv(ima),idv(ima)));
                          _expr=grad(v)*vec(idv(imaP),idv(imaP2)));
    
            auto a = form2( _trial=Vh, _test=Vh);
            a = integrate( _range=elements(mesh),
                           //_expr=gradt(u)*trans(grad(v)));
                           _expr=gradt(u)*trans(grad(v)));
    
            a.solve(_rhs=l,_solution=u) ;
    
            auto value  = vf::project(Vh,elements(mesh),Py()); 
            auto n = normL2 ( _range=elements( mesh ), _expr=idv(u)-idv(value)); 
            std::cout << " Norm2 " << n << std::endl;

        }
  
    void runImageHbf  (std::string hbf1, std::string hbf2, std::string solution, double pixelsize) 
        {
#if 0
            auto x = readHBF( hbf1 );
            auto y = readHBF( hbf2 );
            auto z = readHBF( solution );

            auto mesh = std::make_shared<MeshStructured>( x.rows(),
                                                            x.cols(),
                                                            pixelsize,
                                                            cx,cy,
                                                            /*NULL,*/
                                                            Environment::worldCommPtr(),
                                                            "",
                                                            false,
                                                            false);
            mesh->components().reset();
            mesh->components().set( size_type(MESH_NO_UPDATE_MEASURES) );
            mesh->updateForUse();
    
            auto Vh = Pch<1> ( mesh ) ;
            auto u = Vh->element () ;
            auto v = Vh->element () ;
            auto px = Vh->element () ;
            auto py = Vh->element () ;
            Hbf2FeelppStruc h2f( nx, ny, Vh );
            Hbf2FeelppStruc h2fS( z.cols(), z.rows(), Vh );
            px = h2f( x);
            py = h2f( y );
            auto s = h2fS ( z);
      
            auto l = form1( _test=Vh );
            l = integrate(_range=elements(mesh),
                          _expr=grad(v)*vec(idv(px),idv(py)));
    
            auto a = form2( _trial=Vh, _test=Vh);
            a = integrate( _range=elements(mesh),
                           //_expr=gradt(u)*trans(grad(v)));
                           _expr=gradt(u)*trans(grad(v)));
    
            a.solve(_rhs=l,_solution=u) ;
    
            std::cout << "L2 error norm : " << normL2( _range=elements(mesh), _expr=idv(u)-idv(s) ) << "\n"; 
#endif   

        }


    void runParallel ( double pixelsize )
        {
            tic();
            Feel::cout << "[nbPt X] x [nbPt Y] = " << ima.rows() << "\t" << ima.cols() << std::endl;
            auto mesh = std::make_shared<MeshStructured>( ima.rows(), ima.cols(), pixelsize, ima,ima,
                                                          Environment::worldCommPtr(),"", false, false);
            Feel::cout <<"  mesh->numGlobalElements() " << mesh->numGlobalElements() << std::endl;
            Feel::cout <<"  mesh->numGlobalPoints()   " << mesh->numGlobalPoints()   << std::endl;
            mesh->components().reset();
            Feel::cout <<"  mesh->numGlobalElements() " << mesh->numGlobalElements() << std::endl;
            Feel::cout <<"  mesh->numGlobalPoints()   " << mesh->numGlobalPoints()   << std::endl;
            //mesh->components().set( size_type(MESH_NO_UPDATE_MEASURES) );
            mesh->components().set( size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES) );
            Feel::cout <<"  mesh->numGlobalElements() " << mesh->numGlobalElements() << std::endl;
            Feel::cout <<"  mesh->numGlobalPoints()   " << mesh->numGlobalPoints()   << std::endl;
            mesh->updateForUse();
            Feel::cout <<"  mesh->numGlobalElements() " << mesh->numGlobalElements() << std::endl;
            Feel::cout <<"  mesh->numGlobalPoints()   " << mesh->numGlobalPoints()   << std::endl;
 

            auto Vh = Pch<1>(mesh,true);
            std::cout << "Rank - nDof \t" << Environment::worldComm().localRank() << "\t" <<  Vh->nLocalDof() << std::endl;
            //Vh->map().showMe(true,std::cout);
            std::cout << "Rank - Ghost Local Dof :" << Environment::worldComm().localRank() << "\t"<< Vh->dof()->nLocalGhosts() << std::endl;
            //auto a = form2(Vh,Vh);
            //a.matrixPtr()->printMatlab("theMat");

            auto e = exporter( _mesh=mesh );
            e->addRegions();
            e->save(); 

        }



    private :

    holo3_image<float> ima;
    holo3_image<float> ima2;
    
};

#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( multiscale_suite )

double ps=1e-2;

BOOST_AUTO_TEST_CASE( test_run0 )
{
    Test_MeshStructured tms0(13,4);
    tms0.runImage(ps);
    //tms0.runParallel(ps);
}


BOOST_AUTO_TEST_CASE( test_run2 )
{
    Test_MeshStructured tms2(12,24);
    tms2.runParallel(ps);
}

BOOST_AUTO_TEST_CASE( test_run3 )
{
    Test_MeshStructured tms3(42,18);
    tms3.runParallel(ps);
}

/*
BOOST_AUTO_TEST_CASE( test_run4 )
{
    Test_MeshStructured tms4(43,9);
    tms4.runParallel(ps);
}
*/


/*
BOOST_AUTO_TEST_CASE( test_run1 )
{
    Test_MeshStructured tms1(18,22);
    tms1.runImage(8.9e-3);
}
*/





BOOST_AUTO_TEST_SUITE_END()
#else
std::cout << "USE_BOOST_TEST non define" << std::endl;
#endif


}
