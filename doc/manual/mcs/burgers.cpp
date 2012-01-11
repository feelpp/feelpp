/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-01-09

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \file Burgers.cpp
   \author Mehdi DAOU, Mamadou BATHILY
   \date 2011-02-11
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

// here we define all parameter that we can to re-define when we execute this application
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description burgersoptions("Burgers problem options");
    burgersoptions.add_options()
        ("nu", Feel::po::value<double>()->default_value( 0.01 ), "value of viscosity")
        ("dt", Feel::po::value<double>()->default_value( 1.0 ), "step time")
        ("penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Feel::po::value<double>()->default_value( 0.05 ), "first h value to start convergence")
        ;
    return burgersoptions.add( Feel::feel_options() );
}
// information about this class
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "Burgers" ,
                           "Burgers" ,
                           "0.1",
                           "nD(n=1) Burgers problem",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2009 Université Joseph Fourier");

    about.addAuthor("Mehdi DAOU, Mamadou Bathily", "developer", "", "");
    return about;

}

namespace Feel
{
    using namespace Feel::vf;
    /**
     * Burgers Problem
     *
     * solve \f$\frac{\partial t}{\partial u}-\nu \Delta u + u*\nabla u = 0, u_I = 0\f$ on \f$ I \f$ , \f$ u_\Omega = \sin (x) \f$ on \f$\Omega\f$
     */
    template<int Dim,
             int Order = 1,
             template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
    class Burgers
        :
            public Application
    {
        typedef Application super;

        public:

            // -- TYPEDEFS --
            typedef Burgers<Dim,Order, Entity> self_type;

            typedef double value_type;

            typedef Backend<value_type> backend_type;
            typedef boost::shared_ptr<backend_type> backend_ptrtype;

            /*matrix*/
            typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
            typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
            typedef typename backend_type::vector_type vector_type;
            typedef typename backend_type::vector_ptrtype vector_ptrtype;

            /*mesh*/
            typedef Entity<Dim, 1,Dim> entity_type;
            typedef Mesh<entity_type> mesh_type;
            typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

            /*basis*/
            typedef bases<Lagrange<Order, Scalar> > basis_type;

            /* number of dofs per element */
            static const uint16_type nLocalDof =
                                             boost::remove_reference<typename fusion::result_of::at<basis_type,mpl::int_<0> >::type>::type::nLocalDof;

            /*space*/
            typedef FunctionSpace<mesh_type, basis_type> functionspace_type;
            typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
            typedef typename functionspace_type::element_type element_type;
            typedef typename element_type::template sub_element<0>::type element_0_type;
            typedef typename element_type::template sub_element<1>::type element_1_type;

            typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
            typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
            typedef FsFunctionalLinear<functionspace_type> funlin_type;
            typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

            /* export */
            typedef Exporter<mesh_type> export_type;
            typedef boost::shared_ptr<export_type> export_ptrtype;

            /**
             * Constructor
             */
            Burgers( int argc, char** argv, AboutData const& ad, po::options_description const& od );

            /**
             * run the convergence test
             */
            void run();

    private:

        /**
         * export results to ensight format (enabled by  --export cmd line options)
         */
        void exportResults( element_type& u, double t );

    private:

        backend_ptrtype M_backend;
        // mesh size
        double meshSize;
        // viscosity
        double M_nu;
        // step time
        double dt;

        functionspace_ptrtype M_Xh;

        export_ptrtype exporter;
    }; // Burgers

    template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
    Burgers<Dim,Order,Entity>::Burgers( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        M_nu( this->vm()["nu"].template as<double>() ),
        dt(this->vm()["dt"].template as<double>() ),
        M_Xh(),
        exporter()
    {
        if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

        this->changeRepository( boost::format( "doc/tutorial/mcs/%1%/%2%/P%3%/h_%4%/nu_%5%" )
                                % this->about().appName()
                                % entity_type::name()
                                % Order
                                % this->vm()["hsize"].template as<double>()
                                % this->vm()["nu"].template as<double>()
                                );
        //creation of the mesh
        value_type pi = M_PI;
        mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                            _desc=domain( _name= (boost::format( "%1%-%2%-%3%" ) % "simplex" % Dim % 1).str() ,
                                                          _shape="simplex",
                                                          _dim=Dim,
                                                          _order=1,
                                                          _xmin=0.0,
                                                          _xmax=2*pi,
                                                          _h=meshSize ) );

        M_Xh = functionspace_ptrtype( functionspace_type::New( mesh ) );
        exporter = export_ptrtype( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) );
    }


    template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
    void Burgers<Dim, Order, Entity>::run()
    {
        using namespace Feel::vf;
        //mesh
        mesh_ptrtype mesh = M_Xh->mesh();
        //creation of vector and matrix
        element_type u( M_Xh, "u" );
        element_type v( M_Xh, "v" );
        element_type vv( M_Xh, "vv" );

        // u0's value
        auto M_c = sin(Px());
        auto uo=vf::project(M_Xh,elements(mesh),M_c);

        value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
        auto t = 0.0;
        auto f = 0;
        auto g = 0;
        //export result
        exportResults( uo, t);
        // definition of our vector and matrix that we use in our
        auto F = M_backend->newVector(M_Xh);
        auto D = M_backend->newMatrix(M_Xh,M_Xh);

        // we begin resolution with t=dt
        t=dt;


       //velocity at t=0
        vv = uo;
        //time loop
        for(t = dt;t<=100.0;t+=dt)
        {
            //creation of the matrix D
            form2( _test=M_Xh, _trial=M_Xh, _matrix=D , _init=true) = integrate( elements( mesh ), M_nu*gradt(u)*trans(grad(v)) );
            form2( M_Xh, M_Xh, D ) +=  integrate( elements( mesh ),  ( ( idt(u)/dt) + idv(vv)*gradt(u) )*id(v) );
            form2( M_Xh, M_Xh, D ) +=  integrate( boundaryfaces(mesh),
                                   ( - trans(id(v))*(gradt(u)*N())));
            D->close();

            //bondary conditions
            form2( M_Xh, M_Xh, D ) +=
                on( markedfaces(mesh,1), u, F, cst(0.) );
            form2( M_Xh, M_Xh, D ) +=
                on( markedfaces(mesh,3), u, F, cst(0.) );

            //creation of th risidual vector
            form1( M_Xh, F , _init=true) = integrate( elements( mesh ),  f*id(v)+idv(vv)*id(v)/dt);
            form1( M_Xh, F ) +=
            integrate( boundaryfaces(mesh), -(grad(v)*N())*g );
            F->close();

            //resolution
            backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );

            //we save the u^n-1
            vv = u;
            exportResults( u, t);
        }
    }

    template<int Dim, int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
    void Burgers<Dim, Order, Entity>::exportResults( element_type& U, double t)
    {
        if ( exporter->doExport() )
        {
            Log() << "exportResults starts\n";
            exporter->step(t)->setMesh( U.functionSpace()->mesh() );
            exporter->step(t)->add( "u", U );
            exporter->save();
        }
    }
} // Feel

int main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nDim = 1;
    const int nOrder = 1;
    typedef Feel::Burgers<nDim, nOrder> burgers_app_type;

    /* instantiate application */
    burgers_app_type burgers( argc, argv, makeAbout(), makeOptions() );

    /* run application */
    burgers.run();
}

