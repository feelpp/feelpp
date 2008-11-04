/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-02-01

  Copyright (C) 2006 EPFL

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-02-01
 */
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <iostream>

#include <boost/foreach.hpp>
#include <boost/plugin.hpp>

#include <life/lifecore/application.hpp>




#include <polyvis.hpp>


LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "life_polyvis" ,
                            "life_polyvis" ,
                            "0.1",
                            "Generate ensight files to visualize polynomials",
                            LifeV::AboutData::License_LGPL,
                            "Copyright (c) 2006 EPFL");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@epfl.ch", "");
    return about;

}

class PolyvisApp : public LifeV::Application
{
public:
    typedef LifeV::Application super;

    PolyvisApp( int argc,  char** argv, LifeV::AboutData const& ad, LifeV::po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {}

    void run() ;
};
void
PolyvisApp::run()
{
    try {
        std::string poly = vm()["poly"].as<std::string>();
        std::string lib = vm()["lib"].as<std::string>();
        /* get the handle of the library */
        boost::plugin::dll d ( lib );

        boost::plugin::plugin_factory <Polyvis> pf (d);


        std::cout << "*** Creating an instance of plugin class " << poly << " out of lib " << lib << "\n";
        std::auto_ptr <Polyvis> p (pf.create ( poly,
                                               "Lagrange",
                                               vm()["dim"].as<int>(),
                                               vm()["order"].as<int>()
                                               ));

        std::cout << "*** Calling method of the created instance\n";
    }
    catch ( std::logic_error const & e )
        {
            /* report error, and skip the library */
            std::cerr << "Could not load polynomial family: " << e.what () << std::endl;
        }

}

int
main( int argc, char** argv )
{
    LifeV::po::options_description desc("Specific options");
    desc.add_options()
        ("lib", LifeV::po::value<std::string>()->default_value("stdpoly.so"), "library of polynomials")
        //("list", "list components in pluging")
        ("poly,p", LifeV::po::value<std::string>()->default_value("Lagrange"), "polynomials to display")
        ("dim,d", LifeV::po::value<int>()->default_value(2), "dimension")
        ("order,o", LifeV::po::value<int>()->default_value(1), "polynomials order")
        ("convex,c", LifeV::po::value<int>()->default_value(1), "convex")
        ;

    PolyvisApp app( argc, argv, makeAbout(), desc );

    app.run();



}

#if 0
const int OrderPts = 10;


template<typename P, template<LifeV::uint16_type> class PS >
void
visu( std::string const& s, LifeV::Polynomial<P,PS> const& poly )
{
    using namespace LifeV;
    const uint16_type nDim = P::nDim;
    PointSetToMesh<LifeV::Simplex<nDim,OrderPts>, double> p2m;
    typedef typename PointSetToMesh<LifeV::Simplex<nDim,1>, double>::mesh_type mesh_type;

    //GaussLobatto<<Simplex<nDim,OrderPts>, OrderPts, double> GS;
    PointSetEquiSpaced<Simplex<nDim,OrderPts>, double> GS;
    p2m.addBoundaryPoints( false );
    p2m.visit( &GS );

    typedef Space<mesh_type, FEM_PK<nDim, 1, PS, Continuous> > P1_type;
    P1_type P1( p2m.mesh() );

    std::vector<typename P1_type::element_type > u( 1 );
    ublas::matrix<double> evalpoly( poly.evaluate( GS.points() ) );
    // evaluate the dubiner polynomials at the lobatto nodes
    for ( int i = 0;i < 1; ++i )
    {
        u[i] = P1.newElement( "u" );
        if ( poly.is_scalar )
            {
                //u[i].resize( evalpoly.size2() );
                u[i] = ublas::row( evalpoly, 0 );
            }
        else
            {
                int nC = poly.nComponents;
                ublas::matrix<double> m( ublas::project( evalpoly,
                                                         ublas::range( i*nC, (i+1)*nC ),
                                                         ublas::range( 0, evalpoly.size2() ) ) );
                std::cout << "m = " << m << "\n"
                          << "v2m(m) = " << vectorialToMatrix( m, nC ) << "\n";
                u[i] = ublas::row( vectorialToMatrix( m, nC ), 0 );
            }

        std::cout << s << "_u_" << i << u[i] << "\n";
        std::cout << s << "_p_" << i << poly.evaluate( GS.points() ) << "\n";
    }


    typedef ExporterEnsight<mesh_type> export_type;
    typedef typename ExporterEnsight<mesh_type>::timeset_type timeset_type;
    export_type __ensight( s );
    typename export_type::timeset_ptrtype __ts( new timeset_type( s ) );
    __ts->setTimeIncrement( 1.0);
    __ensight.addTimeSet( __ts );

    typename timeset_type::step_ptrtype __step = __ts->step( 1.0 );

    __step->setMesh( p2m.mesh() );
    for ( int i = 0;i < 1; ++i )
    {
        std::ostringstream str;
        str << s << "_p_" << i;
        if ( poly.is_scalar )
            __step->addNodalScalar( str.str(), u[i].size(), u[i].begin(), u[i].end()  );
        else
            __step->addNodalVector( str.str(), u[i].size(), u[i].begin(), u[i].end()  );
    }
    __ensight.save();
}
#endif


#if 0
int main()
{
    using namespace LifeV;

    //fem::Lagrange<2,1,Scalar> lag2;
    //visu( "lag_2_1", lag2.functionShape() );

    fem::Lagrange<2,5,Scalar> lag25;
    visu( "lag_2_5", lag25.functionShape() );

    //fem::Lagrange<3,1,Scalar> lag3;
    //visu( "lag_3_1", lag3.functionShape() );

    //fem::RaviartThomas<3,1> fert;
    //visu( "fert_3_1", fert.functionShape() );
#if 0
    DivFreePolynomialSet<2,1> dfpset;
    visu( "divfree_2_1", dfpset );
    std::cout << "div(divfree)=" << ublas::norm_inf( divergence( dfpset ).coeff() ) << "\n";
    std::cout << "div(divfree 0)=" << ublas::norm_inf( divergence( dfpset.polynomial( 0 ) ).coeff() ) << "\n";
#endif

#if 0
    // rt
    RaviartThomasPolynomialSet<2,0> rt;
    std::cout << "rt.coeff() " << rt.coeff() << "\n";
    visu( "rt_2_0", rt );

    fem::RaviartThomas<2,0> fert;
    visu( "fert_2_0", fert.functionShape() );

    fem::RaviartThomas<2,1> fert1;
    visu( "fert_2_1", fert1.functionShape() );

    fem::RaviartThomas<2,2> fert2;
    visu( "fert_2_2", fert2.functionShape() );
#endif


    // bubble
#if 0
    fem::Lagrange<2,1,Scalar> lag;
    std::cout << lag.evaluate( lag.points() ) << "\n";
    //Bubble<2,Vectorial> bubble_v;
    //LifeV::P1BubblePolynomialSet<2, Scalar> p1bubbleset;

    LifeV::node<double>::type x(2);
    x(0) = 0;
    x(1) = 0;
    fem::Lagrange<2,1,Vectorial> lag_v;
    std::cout << "lag_v([0,0]) = " << lag_v.evaluate( x ) << "\n";
    std::cout << "lag_v(0)[0]([0,0]) = " << lag_v.functionShape().polynomial(0)[0].evaluate( x ) << "\n";

    std::cout << "lag_v(pts) = " << lag_v.evaluate( lag.points() ) << "\n";

    //visu( "lag_v", lag_v.functionShape() );

    std::cout << "lag_v.p(0) = " << lag_v.functionShape().polynomial( 0 ).evaluate( x ) << "\n";

    Bubble<2,Vectorial> bubble_v;
    std::cout << "bubble_v.coeff() = " << bubble_v.coeff() << "\n";
    std::cout << "bubble_v(pts) = " << bubble_v.evaluate( lag.points() ) << "\n";
    visu( "bubble_v", bubble_v );
#endif


#if 0
    typedef OrthonormalPolynomialSet<2,2,Vectorial> space_type;
    space_type pset;

    ublas::vector<double> n(2);
    n(0)=1;
    n(1)=1;

    ublas::vector<double> p(2);
    p(0)=0.5;
    p(1)=-0.5;


    functional::DirectionalComponentPointEvaluation<space_type> f( pset, n, p );

    std::cout << f.coeff() << "\n";



#if 0
    LifeV::OrthogonalPolynomialSet<2, 2, Scalar> ortho;
    visu( "ortho", ortho );

    LifeV::P1BubblePolynomialSet<2, Scalar> p1bubbleset;
    visu( "p1b", p1bubbleset );

    FEM_PK<2,1,Scalar> p1;
    visu( "p1", p1.functionShape() );
    FEM_PK<2,2,Scalar> p2;
    visu( "p2", p2.functionShape() );
#endif
    P1Bubble<2,Scalar> p1b;
    //visu( "p1bubble", p1b.functionShape() );
#endif
}
#endif // 0


