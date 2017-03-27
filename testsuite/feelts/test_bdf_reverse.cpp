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
#define BOOST_TEST_MODULE test_bdf_reverse
#include <testsuite/testsuite.hpp>

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
    AboutData about( "test_bdf_reverse" ,
                     "test_bdf_reverse" ,
                     "0.2",
                     "nD(n=2,3) test reverse",
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
    using mesh_ptrtype = boost::shared_ptr<mesh_type>;

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

    void reverse()
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );

        BOOST_TEST_MESSAGE( "Test bdf reverse settings" );
        ts->setTimeInitial(20);
        ts->setTimeFinal(100);
        ts->setTimeStep(10);
        ts->setReverse(true);
        ts->start();
        CHECK( ts->isReverse() == true );
        BOOST_CHECK_MESSAGE( ts->timeInitial() == 100, "current Ti=" << ts->timeInitial() << "\n");
        BOOST_CHECK_MESSAGE( ts->timeFinal() == 20, "current Tf=" << ts->timeFinal() << "\n");
        BOOST_CHECK_MESSAGE( ts->timeStep() == -10, "current dt=" << ts->timeStep() << "\n");
        for ( ; ts->isFinished() == false; ts->next() )
            ;// Go to end

        // BDF State is STOPPED here.
    }

    void reverseInit()
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );
        auto u = Xh->element();
        u.zero();

        BOOST_TEST_MESSAGE( "Test bdf reverse settings" );
        ts->setTimeInitial(20);
        ts->setTimeFinal(100);
        ts->setTimeStep(10);
        ts->setReverse(true);
        ts->start(u);
        CHECK( ts->isReverse() == true );
        BOOST_CHECK_MESSAGE( ts->timeInitial() == 100, "current Ti=" << ts->timeInitial() << "\n");
        BOOST_CHECK_MESSAGE( ts->timeFinal() == 20, "current Tf=" << ts->timeFinal() << "\n");
        BOOST_CHECK_MESSAGE( ts->timeStep() == -10, "current dt=" << ts->timeStep() << "\n");
        for ( ; ts->isFinished() == false; ts->next(u) )
        {
            LOG(INFO) << "BDF reverse iteration: " << ts->iteration() << ", time:" << ts->time();
            const int iter = ts->iteration();
            u.on(_range=elements(M_mesh), _expr=cst(iter) );
        }

        // BDF State is STOPPED here.
    }

    // Test several successive reverse
    // Bdf goes from 0->T then T->0, N times.
    void forwardReverse( int n=1 )
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );

        for(int i=0;i<n;i++)
        {
            BOOST_TEST_MESSAGE( "Test bdf initial settings" );
            ts->setReverse(false);
            ts->setTimeInitial(20);
            ts->setTimeFinal(100);
            ts->setTimeStep(10);

            ts->start();

            BOOST_CHECK( ts->isReverse() == false );
            BOOST_CHECK( ts->timeInitial() == 20 );
            BOOST_CHECK( ts->timeFinal() == 100 );
            BOOST_CHECK( ts->timeStep() == 10 );

            for ( ; ts->isFinished() == false; ts->next() )
                ;// Go to end

            // BDF State is STOPPED here.

            BOOST_TEST_MESSAGE( "Test bdf reverse settings" );
            ts->setReverse(true);
            ts->start();

            BOOST_CHECK( ts->isReverse() == true );
            BOOST_CHECK( ts->timeInitial() == 100 );
            BOOST_CHECK( ts->timeFinal() == 20 );
            BOOST_CHECK( ts->timeStep() == -10 );

            for ( ; ts->isFinished() == false; ts->next() )
                ;

            // BDF State is STOPPED here.
        }
    }

    // Test BDF load solution from HDF5 file, when reverse is called.
    // BDF run from time T to 0 but iter 0 to N. We have to be careful with
    // HDF5 files name indexed by iteration number.
    void forwardReverseLoad( std::string fileformat="hdf5" )
    {
        auto Xh = Pch<1>( M_mesh );
        auto ts = bdf( _space=Xh, _name="ts" );
        auto u = Xh->element();
        u.zero();

        BOOST_TEST_MESSAGE( "Test bdf initial settings" );
        ts->setfileFormat(fileformat);
        ts->setTimeInitial(20);
        ts->setTimeFinal(100);
        ts->setTimeStep(10);
        // Restart from iteration 0
        ts->start(u);

        BOOST_CHECK( ts->isReverse() == false );
        BOOST_CHECK( ts->timeInitial() == 20 );
        BOOST_CHECK( ts->timeFinal() == 100 );
        BOOST_CHECK( ts->timeStep() == 10 );

        std::vector<double> v;
        v.reserve((100-20)/10);
        for ( ; ts->isFinished() == false; ts->next(u) )
        {
            LOG(INFO) << "BDF forward iteration: " << ts->iteration() << ", time:" << ts->time();
            const int iter = ts->iteration();
            u.on(_range=elements(M_mesh), _expr=cst(iter) );
            double val = integrate( _range=elements(M_mesh),
                                    _expr=idv(u) ).evaluate()(0,0);
            v.push_back(val);
        }

        // BDF State is STOPPED (reinit) here. u(x,t) saved in HDF5 files.

        u.zero();
        BOOST_TEST_MESSAGE( "Test bdf reverse load solution" );
        ts->setReverse(true);
        // This change the order to read files from filename-N to filename-0.
        ts->setReverseLoad(true);
        ts->setRestart(true);
        ts->setRestartAtLastSave(false);
        ts->start();

        BOOST_CHECK( ts->isReverse() == true );
        BOOST_CHECK( ts->timeInitial() == 100 );
        BOOST_CHECK( ts->timeFinal() == 20 );
        BOOST_CHECK( ts->timeStep() == -10 );
        int k=v.size()-1;
        for ( ; ts->isFinished() == false; ts->next(), --k )
        {
             LOG(INFO) << "BDF reversed iteration:" << ts->iteration() << ", time:" << ts->time();
             ts->loadCurrent();
             auto w = ts->unknown(0);
             double val = integrate( _range=elements(M_mesh),
                                     _expr=idv(w) ).evaluate()(0,0);
             double err= (val-v[k])*(val-v[k]);
             LOG(INFO) << "err = " << err << " val=" << val << "v[" << k << "]=" << v[k];
             BOOST_TEST_MESSAGE( "err = " << err << " val=" << val << "v[" << k << "]=" << v[k] );
             BOOST_CHECK(  ts->time() <= ts->timeInitial() );
             BOOST_CHECK(  ts->time() >= ts->timeFinal() );
             BOOST_CHECK(  err < 1e-10 );
        }

        // BDF State is STOPPED here.
    }

    private:
        mesh_ptrtype M_mesh;
};


FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( bdf_reverse )

BOOST_AUTO_TEST_CASE( test_reverse )
{
    Test<2> t;
    t.reverse();
}

BOOST_AUTO_TEST_CASE( test_reverse_init )
{
    Test<2> t;
    t.reverseInit();
}

BOOST_AUTO_TEST_CASE( test_forward_reverse )
{
    Test<2> t;
    t.forwardReverse();
}

BOOST_AUTO_TEST_CASE( test_forward_reverse_n )
{
    Test<2> t;
    t.forwardReverse(2);
}

BOOST_AUTO_TEST_CASE( test_forward_reverse_h5load )
{
    Test<2> t;
    t.forwardReverseLoad();
}

BOOST_AUTO_TEST_SUITE_END()
