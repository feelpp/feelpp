/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
      Date: 2011-08-08

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
/**
   \file test_lift.cpp
   \author Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
   \date 2011-08-08
*/

#include <feel/options.hpp>

#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feeldiscr/operatorlift.hpp>

#include <feel/feeldiscr/region.hpp>

#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelfilters/exporter.hpp>

#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description testLiftoptions("TestLift options");
    testLiftoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ("nu", po::value<double>()->default_value( 1 ), "grad.grad coefficient")
        ("weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
        ("penaldir", Feel::po::value<double>()->default_value( 10 ),
         "penalisation parameter for the weak boundary Dirichlet formulation")
        ;
    return testLiftoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "testlift" ,
                     "testlift" ,
                     "0.2",
                     "nD(n=1,2,3) TestLift on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

template<int Dim>
class TestLift
    :
    public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 2;

    typedef double value_type;

    typedef Backend<value_type> backend_type;

    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;

    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;

    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef Simplex<Dim> convex_type;

    typedef Mesh<convex_type> mesh_type;

    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef bases<Lagrange<Order,Vectorial> > basis_type;

    typedef FunctionSpace<mesh_type, basis_type> space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::element_type element_type;

    typedef Exporter<mesh_type> export_type;

    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestLift( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:


    backend_ptrtype M_backend;

    double meshSize;

    std::string shape;
}; // TestLift

template<int Dim> const uint16_type TestLift<Dim>::Order;

template<int Dim>
void
TestLift<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestLift<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
TestLift<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );


    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                        _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % 1 ).str() ,
                                                      _usenames=false,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=X[0],
                                                       _xmin=-1.0,
                                                      _xmax=1.0,
                                                      _ymin=-1.0,
                                                      _ymax=1.0 ) );



    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );

    using namespace Feel::vf;

    value_type pi = M_PI;

    // auto g = cos(pi*Py());
    auto g = cos(pi*Py())*unitX() + cos(pi*Py())*unitY() ;

    WeakDirichlet dir_type = (WeakDirichlet)this->vm()["weakdir"].template as<int>();

    auto op_lift = operatorLift(Xh, M_backend, 20.0, dir_type);

    auto glift1 = op_lift->lift( _range=markedfaces(mesh,1),_expr=trans(g));

    auto glift2 = (*op_lift)(markedfaces(mesh,1),trans(g));

    auto ux = idv(glift1)(0,0);

    auto uy = idv(glift1)(1,0);

    auto tempx = trace(hessv(glift1)(0,0));

    auto tempy = trace(hessv(glift1)(1,0));

    auto temp = tempx*oneX() + tempy*oneY();

    // auto verif1 = vf::project(Xh, elements(mesh), idv(temp));
    //auto verif1 = vf::project(Xh, elements(mesh), trace(hessv(glift1)));
    //auto verif2 = vf::project(Xh, elements(mesh), trace(hessv(glift2)));

    auto gproj = vf::project(Xh, elements(mesh), g);


    export_ptrtype exporter( export_type::New( this->vm(),
                                               (boost::format( "%1%-%2%-%3%" )
                                                % this->about().appName()
                                                % shape
                                                % Dim).str() ) );
    if ( exporter->doExport() )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( mesh );

        exporter->step(0)->add( "glift1", glift1 );
        exporter->step(0)->add( "glift2", glift2 );

        // exporter->step(0)->add( "verif1", verif1 );
        // exporter->step(0)->add( "verif2", verif2 );
        exporter->step(0)->add( "g", gproj );

        exporter->save();
        Log() << "exportResults done\n";
    }
} // TestLift::run

int
main( int argc, char** argv )
{
    Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    app.add( new TestLift<2>( app.vm(), app.about() ) );

    app.run();

}





