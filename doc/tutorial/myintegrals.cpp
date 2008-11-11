/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-07

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file myintegrals.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-07
 */
#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>


using namespace Life;

inline
po::options_description
makeOptions()
{
    po::options_description myintegralsoptions("MyintegralsMyintegrals options");
    myintegralsoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ;
    return myintegralsoptions.add( Life::life_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "myintegrals" ,
                     "myintegrals" ,
                     "0.2",
                     "nD(n=1,2,3) MyIntegrals on simplices or simplex products",
                     Life::AboutData::License_GPL,
                     "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


/**
 * MyIntegrals: compute integrals over a domain of \f$\mathbb{R}^d,\ d=-1,2,3\f$
 *
 * @author Christophe Prud'homme
 */
template<int Dim>
class MyIntegrals
    :
    public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2;

    typedef double value_type;
    typedef Application application_type;

    /*mesh*/
    typedef Simplex<Dim, 1,Dim> entity_type;
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*quadrature*/
    typedef IM<Dim, imOrder, value_type, Simplex> im_type;

    MyIntegrals( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        meshSize( this->vm()["hsize"].template as<double>() )
    {
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

    double meshSize;

}; // MyIntegrals

template<int Dim> const uint16_type MyIntegrals<Dim>::imOrder;

template<int Dim>
typename MyIntegrals<Dim>::mesh_ptrtype
MyIntegrals<Dim>::createMesh( double meshSize )
{
    mesh_ptrtype mesh( new mesh_type );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,Simplex> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    mesh->setComponents( MESH_PARTITION| MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
    mesh->updateForUse();
    return mesh;
} // MyIntegrals::createMesh


template<int Dim>
void
MyIntegrals<Dim>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    //    int maxIter = 10.0/meshSize;
    using namespace Life::vf;


    this->changeRepository( boost::format( "%1%/%2%/h_%3%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % this->vm()["hsize"].template as<double>()
                            );
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );

    im_type im;

    /*
     * Compute domain Area
     */
    double local_domain_area = integrate( elements(mesh), im,
                                          constant(1.0)).evaluate()(0,0);
    double global_domain_area=local_domain_area;
    if ( Application::nProcess() > 1 )
        mpi::all_reduce( Application::comm(),
                         local_domain_area,
                         global_domain_area,
                         std::plus<double>() );
    Log() << "int_Omega = " << global_domain_area
          << "[ " << local_domain_area << " ]\n";

    /*
     * Compute domain perimeter
     */
    double local_boundary_length = integrate( boundaryfaces(mesh), im,
                                            constant(1.0)).evaluate()(0,0);
    double global_boundary_length = local_boundary_length;
    if ( Application::nProcess() > 1 )
        mpi::all_reduce( Application::comm(),
                         local_boundary_length,
                         global_boundary_length,
                         std::plus<double>() );
    Log() << "int_Omega = " << global_boundary_length
          << "[ " << local_boundary_length << " ]\n";

    /*
     * Compute \int f where f= x^2 + y^2 + z^2
     */
    double local_intf = integrate( elements(mesh), im,
                                   Px()*Px() + Py()*Py() + Pz()*Pz()
                                   ).evaluate()(0,0);
    double global_intf = local_intf;
    if ( Application::nProcess() > 1 )
        mpi::all_reduce( Application::comm(),
                         local_intf,
                         global_intf,
                         std::plus<double>() );
    Log() << "int_Omega = " << global_intf
          << "[ " << local_intf << " ]\n";


} // MyIntegrals::run

int
main( int argc, char** argv )
{
    MyIntegrals<2> myintegrals( argc, argv, makeAbout(), makeOptions() );

    /* assertions handling */
    Life::Assert::setLog( "myintegrals.assert");

    myintegrals.run();
}





