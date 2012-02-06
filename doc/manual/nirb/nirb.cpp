/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-02-01

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file nirb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-02-01
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
	// meshs parameters
        ("hfinsize", po::value<double>()->default_value( 0.1 ), "fine mesh size")
        ("hcoarsesize", po::value<double>()->default_value( 0.5 ), "coarse mesh size")
	
	// Reduced basis parameters        
	("NbSnapshot",po::value<int>()->default_value(100),"numbers of snapshot computed")
        ("sizeBR",po::value<int>()->default_value(15),"size of reduced basis")
	
	("mu", po::value<double>()->default_value( 0 ), "angle in [0,pi/2]")
        ;
    return laplacianoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "nirb-test" ,
                     "nirb-test" ,
                     "0.2",
                     "Non intrusive reduced basis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


class NIRBTEST
    :
    public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 3;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;


    /**
     * Constructor
     */
    NIRBTEST( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        FineMeshSize( this->vm()["hfinsize"].as<double>() ),
        CoarseMeshSize( this->vm()["hcoarsesize"].as<double>() ),
	NbSnapshot( this->vm()["NbSnapshot"].as<int>() ),
	sizeBR(this->vm()["sizeBR"].as<int>()),
        mu( this->vm()["mu"].as<double>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );
    element_type blackbox( space_ptrtype Xh, double param );
    
private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double FineMeshSize;
    double CoarseMeshSize;

    // Reduced basis parameter
    int NbSnapshot,sizeBR;
    double mu;
}; // NIRBTEST
const uint16_type NIRBTEST::Order;
void
NIRBTEST::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute NIRBTEST<" << 2 << ">\n";
    std::vector<double> X(5);
    X[0] = FineMeshSize;
    X[1] = CoarseMeshSize;
    X[2] = NbSnapshot;
    X[3] = sizeBR;
    X[4] = mu;
    std::vector<double> Y( 1 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
void
NIRBTEST::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "%1%/P%2%/" )
                                       % this->about().appName()
                                       % Order );


    mesh_ptrtype meshFine = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name="Omega1",
                                                       _usenames=true,
                                                       _shape="hypercube",
                                                       _dim=2,
                                                       _h=X[0] ) );
    mesh_ptrtype meshCoarse = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name="Omega2",
                                                       _usenames=true,
                                                       _shape="hypercube",
                                                       _dim=2,
                                                       _h=X[1] ));

    space_ptrtype XhFine = space_type::New( meshFine );
    space_ptrtype XhCoarse = space_type::New( meshCoarse );

    export_ptrtype exporterFine( export_type::New( this->vm(), "nirbInFine" ) );
    export_ptrtype exporterCoarse( export_type::New( this->vm(), "nirbInCoarse" ) );
    export_ptrtype exporter2Grid(export_type::New( this->vm(), "nirbOut"));

    exporterFine->step(0)->setMesh( meshFine );
    exporterCoarse->step(0)->setMesh( meshCoarse );
    exporter2Grid->step(0)->setMesh( meshFine );


    //STEP ONE : Construction of the "non intruisive reduced basis (nirb) functions"
    //Computation of the X[3] snapshots solution on Xhfine
    //Selection of X[2] fonctions to build the "reduced basis" using a POD technique
    //Orthogonalisation in L2 and H1 norm of "reduced basis function"
    //Save this final functions
     

   //STEP TWO : Approximation of the solution using the "nirb" functions for a choosen mu 
    	



    for( int i =0; i < sizeBR; ++i )
    {
        double p = i*M_PI/(2*9);
        auto uFine = blackbox( XhFine, p );
        auto uCoarse = blackbox( XhCoarse, p );

        exporterFine->step(0)->add( (boost::format("uFine-%1%")%i).str(), uFine );
        exporterCoarse->step(0)->add( (boost::format("uCoarse-%1%")%i).str(), uCoarse );
    }

    exporterFine->save();
    exporterCoarse->save(); 
    exporter2Grid->save();


}
NIRBTEST::element_type
NIRBTEST::blackbox( space_ptrtype Xh, double param )
{

    auto u = Xh->element();
    auto v = Xh->element();
    std::cout << "Xh ndof=" << Xh->nLocalDof() << "\n";
    value_type pi = M_PI;

    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto velocity = vec(cst(cos(param)),cst(sin(param)));
    auto g = ( chi(abs(Px()-1) < 1e-10 )*Px()*Px()+
               chi(abs(Py()-1) < 1e-10 )*Py()*Py() );
    auto F = M_backend->newVector( Xh );

    form1( _test=Xh, _vector=F ) =
        integrate( _range=boundaryfaces(Xh->mesh()),
                   _expr=g*(-0.01*dn(v)+30*id(v)/hFace())  );

    auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );
    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=0.01*gradt(u)*trans(grad(v))+ (gradt(u)*velocity)*id(v) );

    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate( boundaryfaces(Xh->mesh()),
                   -0.01*dnt(u)*id(v)-0.01*dn(v)*idt(u)+30*id(v)*idt(u)/hFace());

    backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );

    return u;
} // NIRBTEST::run

/**
 * main function: entry point of the program
 */
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
    /**
     * register the simgets
     */
    app.add( new NIRBTEST( app.vm(), app.about() ) );

    /**
     * run the application
     */
    app.run();
}






