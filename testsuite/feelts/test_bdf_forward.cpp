/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-28

  Copyright (C) 2011-2014 Feel++ Consortium

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_bdf_forward
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using Feel::project;

inline
AboutData
makeAbout()
{
    AboutData about( "test_bdf_forward" ,
                     "test_bdf_forward" ,
                     "0.2",
                     "nD(n=2,3) test forward",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dolle", "developer", "gdolle@unistra.fr", "" );
    return about;

}

template<int Dim>
class Test
{
public :
    using mesh_type = Mesh<Simplex<Dim,1>>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    Test():
        M_mesh( createGMSHMesh( _mesh=new mesh_type,
                    _desc=domain( _name=Dim > 2 ? "sphere" : "circle",
                        _shape="ellipsoid",
                        _dim=Dim,
                        _xmin=-1,
                        _ymin=-1,
                        _zmin=-1,
                        _h=doption("gmsh.hsize") ) ) )
                {}

    // Run BDF forward one or several times.
    void forward( int N, bool restart = true, std::string fileformat="hdf5")
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );

        for(int i=0;i<N;i++)
        {
            BOOST_TEST_MESSAGE( "Test bdf initial settings" );
            ts->setfileFormat(fileformat);
            ts->setRestart(restart);
            ts->setTimeInitial(20);
            ts->setTimeFinal(100);
            ts->setTimeStep(10);

            ts->start();

            BOOST_CHECK( ts->iteration() < ts->iterationNumber() );
            BOOST_CHECK( ts->timeInitial() == 20 );
            BOOST_CHECK( ts->timeFinal() == 100 );
            BOOST_CHECK( ts->timeStep() == 10 );

            for ( ; ts->isFinished() == false; ts->next() )
                ;// Go to end
        }
    }

    void forwardLoad( int N, bool restart = true, bool shift=false, std::string fileformat="hdf5")
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );
        auto u = Xh->element();
        u.zero();

        for(int i=0;i<N;i++)
        {
            BOOST_TEST_MESSAGE( "Test bdf initial settings" );
            ts->setfileFormat(fileformat);
            ts->setRestart(restart);
            ts->setTimeInitial(20);
            ts->setTimeFinal(100);
            ts->setTimeStep(10);

            ts->start(u);

            BOOST_CHECK( ts->iteration() < ts->iterationNumber() );
            BOOST_CHECK( ts->timeInitial() == 20 );
            BOOST_CHECK( ts->timeFinal() == 100 );
            BOOST_CHECK( ts->timeStep() == 10 );

            for ( ; ts->isFinished() == false; )
            {
                LOG(INFO) << "BDF forward iteration: " << ts->iteration() << ", time:" << ts->time();
                const int iter = ts->iteration();
                u.on(_range=elements(M_mesh), _expr=cst(iter) );
                ( shift ) ? ts->next(u) : ts->next();
            }
        }
    }

    private:
        mesh_ptrtype M_mesh;
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( bdf_forward )

BOOST_AUTO_TEST_CASE( test_forward_single )
{
    Test<2> t;
    t.forward(1,false);
}

BOOST_AUTO_TEST_CASE( test_forward_single_restart)
{
    Test<2> t;
    t.forward(1,true);
}

BOOST_AUTO_TEST_CASE( test_forward_multi_restart )
{
    Test<2> t;
    t.forward(2,true);
}

BOOST_AUTO_TEST_CASE( test_forward_load_single )
{
    Test<2> t;
    t.forwardLoad(1,false);
}

BOOST_AUTO_TEST_CASE( test_forward_load_single_restart )
{
    Test<2> t;
    t.forwardLoad(1,true);
}

BOOST_AUTO_TEST_CASE( test_forward_load_multi_restart )
{
    Test<2> t;
    t.forwardLoad(2,true);
}

BOOST_AUTO_TEST_SUITE_END()
