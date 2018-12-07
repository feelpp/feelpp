/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

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
#if defined(USE_BOOST_TEST)
#define BOOST_TEST_MODULE test_bdf3
#include <feel/feelcore/testsuite.hpp>
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using Feel::project;

inline
AboutData
makeAbout()
{
    AboutData about( "test_bdf3" ,
                     "test_bdf3" ,
                     "0.2",
                     "nD(n=2,3) test bdf3",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2015 Feel++ Consortium" );

    about.addAuthor( "Guillaume Dollé", "developer", "gdolle@unistra.fr", "" );
    return about;

}

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description opts( "Test Environment options" );
    opts.add_options()
        ( "N", Feel::po::value<int>()->default_value( 2 ), "Number of problems to solve successively" )
    ;
    return opts.add( Feel::feel_options() );
}

template<int Dim>
class Test
{
public :
    using mesh_type = Mesh<Simplex<Dim,1>>;
    using element_type = typename Pch_type<mesh_type,1>::element_type;
    using timeset_type = typename Exporter<mesh_type>::timeset_type;
    using timeset_ptrtype = std::shared_ptr<timeset_type>;
    using space_type = Pch_type<mesh_type,1>;
    using bdf_type = Bdf<space_type>;
    using bdf_ptrtype = std::shared_ptr<bdf_type>;

    // Test creation of new timeSet and exporter with bdf.
    // It solves successively three laplacian problem
    // $$ -\Delta u = \{f1,f2,f3\} $$
    // (f1 = 100, f2 = 200, f3 = 300, ... , fN = N*100)
    //
    // Remark: Loop first on sources, then on time!
    //
    // \param resetTimeSet clear timeset and reuse the default one (default: false).
    // \param withPrefix change the first timeset (default) prefix (default:
    // false).
    void runtest( bool resetTimeSet=false,
                  bool withPrefix = false)
    {
        auto mesh = loadMesh( _mesh=new mesh_type );
        const int N = ioption("N");
        const double c = -0.2;
        auto Xh = space_type::New( mesh );
        auto u = Xh->element();
        auto v = Xh->element();

        std::vector<element_type> U;
        for( int i =0; i<N; i++ )
            U.push_back( Xh->element() );

        auto ts = bdf( _space=Xh, _name="mybdf" );

        //stiffness matrix
        auto a = form2( _test=Xh, _trial=Xh );
        auto at = form2( _test=Xh, _trial=Xh ); // Time dependent.
        auto e = exporter( _mesh=mesh, _name="test_bdf" );

        a = integrate( _range = elements( mesh ),
                       _expr = gradt( u )*trans( grad( v ) )
                               + c * ts->polyDerivCoefficient(0)*idt(u)*id(v) );

        auto l = form1(_test=Xh);
        auto lt = form1(_test=Xh); // Time dependent.

        a += on( _range=boundaryfaces(mesh),
                 _rhs=l,
                 _element=u,
                 _expr=cst(0.) );

        // Initialize bdf unknowns.

        for( int i=0; i<N; i++ )
        {
            // ----------------------------------------------------------------
            LOG(INFO) << "\n" << std::string(60,'*') << "\n"
                      << "RHS SOLVE SOURCE " << i << "\n"
                      << std::string(60,'*');
            // ----------------------------------------------------------------
            std::string prefix = ( boost::format( "rhs_%1%" ) % i ).str();
            std::string ui_name = ( boost::format( "u%1%" ) % i ).str();
            std::string bdf_name = ( boost::format( "mybdf%1%" ) % i ).str();

            auto tsi = bdf( _space=Xh, _name=bdf_name );
            tsi->initialize(U[i]);

            l.zero();
            l = integrate( _range=elements( mesh ),
                           _expr=(i+1)*cst(100.)  );

            if( e->doExport() )
            {
                // Create a new timeset for each rhs.
                if( i > 0 )
                {
                    if( resetTimeSet==true )
                    {
                        e->timeSet(0)->clear();
                    }
                    else
                    {
                        auto timeset = std::make_shared<timeset_type>(
                            timeset_type( prefix ) );
                        e->addTimeSet( timeset );
                    }
                }
                // Set all these for the last timeset.
                if( withPrefix )
                    e->setPrefix( prefix );
                e->setMesh( mesh, EXPORTER_GEOMETRY_STATIC );
                e->step(0)->add( ui_name, U[i] );
                e->save();
            }

            for( tsi->start();  tsi->isFinished() == false; tsi->next(U[i]) )
            {
                at=a;
                lt=l;

                auto bdf_poly = tsi->polyDeriv();

                lt += integrate( _range=elements(mesh),
                                 _expr=c*idv(bdf_poly)*id(v) );

                at += on( _range=boundaryfaces(mesh),
                          _rhs=lt,
                          _element=U[i],
                          _expr=cst(0.) );

                at.solve( _rhs=lt,
                          _solution=U[i] );

                // -------------------------------------------------------------
                LOG(INFO) << "Export source: " << i
                          << " time iteration" << tsi->iteration();
                // -------------------------------------------------------------
                auto time = tsi->time();
                if( e->doExport() )
                {
                    e->step(time)->add( ui_name, U[i] );
                }

                if( e->doExport() )
                {
                    e->save();
                }
            } // Bdf loop.
        } // Source loop.
    } // test2.
};


#if defined(USE_BOOST_TEST)
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() );
BOOST_AUTO_TEST_SUITE( bdf3 )

// Test creation of different timeset (N timesets), keep default prefix.
BOOST_AUTO_TEST_CASE( test_1 )
{
    Test<2> test;
    test.runtest();
}

// Test reset default timeset (only 1 timeset), keep default prefix.
BOOST_AUTO_TEST_CASE( test_2 )
{
    Test<2> test;
    test.runtest(true);
}

// Test reset default timeset, rename default prefix.
BOOST_AUTO_TEST_CASE( test_3 )
{
    Test<2> test;
    test.runtest(true, true);
}

BOOST_AUTO_TEST_SUITE_END()
#else
int main(int argc, char** argv )
{
    Feel::Environment env( _argc=argc, _argv=argv,
                           _desc=makeOptions(),
                           _about=makeAbout() );
    Test<2> test;

    // Test creation of different timeset (N timesets), keep default prefix.
    test.runtest();

    // Test reset default timeset (only 1 timeset), keep default prefix.
    //test.runtest(true);

    // Test reset default timeset, rename default prefix.
    //test.runtest(true, true);
}
#endif

