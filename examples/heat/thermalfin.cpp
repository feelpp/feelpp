/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-06-11

  Copyright (C) 2007-2008 Université Joseph Fourier (Grenoble I)

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
   \file thermalfin.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-06-11
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backendgmm.hpp>
#include <life/lifealg/backendpetsc.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>



std::pair<std::string, std::string> makefin( double hsize );

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description thermalfinoptions("Thermalfin options");
    thermalfinoptions.add_options()
        // mesh parameters
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")

        // physical coeff
        ("k0", Life::po::value<double>()->default_value( 1 ), "k0 diffusion parameter")
        ("k1", Life::po::value<double>()->default_value( 1 ), "k1 diffusion parameter")
        ("k2", Life::po::value<double>()->default_value( 1 ), "k2 diffusion parameter")
        ("k3", Life::po::value<double>()->default_value( 1 ), "k3 diffusion parameter")
        ("Bimin", Life::po::value<double>()->default_value( 0.01 ), "minimum value of Biot number")
        ("Bimax", Life::po::value<double>()->default_value( 1 ), "maximum value of Biot number")

        ("N", Life::po::value<int>()->default_value( 1 ), "number of samples withing parameter space")

        // export
        ("export", "export results(ensight, data file(1D)")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return thermalfinoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "thermalfin" ,
                            "thermalfin" ,
                            "0.2",
                            "nD(n=1,2,3) Thermalfin on simplices or simplex products",
                            Life::AboutData::License_GPL,
                            "Copyright (c) 2006, 2007 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
/**
 * Thermal fin application
 *
 */
class ThermalFin
    :
    public Application
{
    typedef Application super;
public:

#define Entity Simplex

    // -- TYPEDEFS --
    static const uint16_type Dim = 2;
    static const uint16_type Order = 1;
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<2> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous > p0_space_type;
    typedef p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;


    ThermalFin( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( export_type::New( this->vm(), this->about().appName() ) )
    {
        this->changeRepository( boost::format( "%1%/%2%/%3%/" )
                                % this->about().appName()
                                % entity_type::name()
                                % this->vm()["hsize"].as<double>()
                                );
        using namespace Life::vf;


        /*
         * First we create the mesh
         */
        mesh = createMesh( meshSize );

        /*
         * The function space and some associate elements are then defined
         */
        Xh = space_type::New( mesh );
        element_type u( Xh, "u" );
        element_type v( Xh, "v" );

        Fcst = M_backend->newVector( Xh );


        form1( Xh, Fcst, _init=true )  = integrate( markedfaces(mesh,1), _Q<imOrder>(), id(v) );

        if ( this->vm().count( "export-matlab" ) )
            Fcst->printMatlab( "F.m" );


        /*
         * Construction of the left hand side
         */
        Dcst = M_backend->newMatrix( Xh, Xh );

        form2( Xh, Xh, Dcst, _init=true ) = integrate( markedelements(mesh,1), _Q<imOrder>(), ( gradt(u)*trans(grad(v))) );

        Dcst->close();
        if ( this->vm().count( "export-matlab" ) )
            Dcst->printMatlab( "Dcst" );

    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;

    mesh_ptrtype mesh;
    space_ptrtype Xh;

    sparse_matrix_ptrtype Dcst;
    vector_ptrtype Fcst;

    boost::shared_ptr<export_type> exporter;

}; // ThermalFin

ThermalFin::mesh_ptrtype
ThermalFin::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    Gmsh gmsh;
    std::string mesh_name, mesh_desc;
    boost::tie( mesh_name, mesh_desc ) = ::makefin(meshSize);
    std::string fname = gmsh.generate( mesh_name, mesh_desc );

    ImporterGmsh<mesh_type> import( "Mesh.0.1.msh" );
    mesh->accept( import );
    return mesh;
} // ThermalFin::createMesh


void
ThermalFin::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    using namespace Life::vf;

    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    double Bimin = this->vm()["Bimin"].as<double>();
    double Bimax = this->vm()["Bimax"].as<double>();
    int N = this->vm()["N"].as<int>();

    for( int i = 0; i < N; ++i )
        {
            /*
             * Construction of the left hand side
             */
            sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );


            //D->addMatrix( 1.0, *Dcst );

            int Nb = 2;
            ublas::matrix<double> k( Nb, Nb );
            k( 0 , 0 ) = this->vm()["k0"].as<double>();
            k( 1 , 0 ) = this->vm()["k1"].as<double>();
            k( 0 , 1 ) = this->vm()["k2"].as<double>();
            k( 1 , 1 ) = this->vm()["k3"].as<double>();
            for( int r = 0; r < Nb; ++r )
                for( int c = 0; c < Nb; ++c )
                    {

                        form2( Xh, Xh, D ) += integrate( markedelements(mesh,Nb*c+r+1), _Q<imOrder>(),
                                                         k(r,c)*gradt(u)*trans(grad(v)) );
                    }

            double Bi;
            if ( N == 1 )
                Bi = 0.1;
            else
                Bi = math::exp( math::log(Bimin)+value_type(i)*(math::log(Bimax)-math::log(Bimin))/value_type(N-1) );
            Log() << "Bi = " << Bi << "\n";

            form2( Xh, Xh, D ) += integrate( markedfaces(mesh,2), _Q<imOrder>(), Bi*idt(u)*id(v) );

            D->close();

            if ( this->vm().count( "export-matlab" ) )
                D->printMatlab( "D" );

            this->solve( D, u, Fcst );

            double moy_u = ( integrate( markedfaces(mesh,1), _Q<imOrder>(), idv(u) ).evaluate()(0,0) /
                             integrate( markedfaces(mesh,1), _Q<imOrder>(), constant(1.0) ).evaluate()(0,0) );
            std::cout.precision( 5 );
            std::cout << std::setw( 5 ) << k(0,0) << " "
                      << std::setw( 5 ) << k(1,0) << " "
                      << std::setw( 5 ) << k(0,1) << " "
                      << std::setw( 5 ) << k(1,1) << " "
                      << std::setw( 5 ) << Bi << " "
                      << std::setw( 10 ) << moy_u << "\n";
            this->exportResults( u );



        }
} // ThermalFin::run

void
ThermalFin::solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F  )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace()->map() ) );
    M_backend->solve( D, D, U, F );
    u = *U;
} // ThermalFin::solve


void
ThermalFin::exportResults( element_type& U )
{
    if ( this->vm().count( "export" ) )
        {
            exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
            exporter->step(1.)->add( "ProcessId",
                           regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );
            exporter->step(1.)->add( "Temperature", U );
            exporter->save();

        }
} // ThermalFin::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* define and run application */
    ThermalFin thermalfin( argc, argv, makeAbout(), makeOptions() );

    thermalfin.run();
}





