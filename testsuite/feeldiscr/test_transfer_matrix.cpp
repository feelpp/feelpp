/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake<abdoulaye.samake@imag.fr>
       Date: 2013-04-17

  Copyright (C) 2013 Universite de Strasbourg

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
   \file test_transfer_matrix.cpp
   \author Abdoulaye Samake<abdoulaye.samake@imag.fr>
   \date 2013-04-17
*/

#define BOOST_TEST_MODULE test_transfer_matrix
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testtransferoptions( "TransferMatrix options" );
    testtransferoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
        ;
    return testtransferoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_transfer" ,
                     "test_transfer" ,
                     "0.2",
                     "nD(n=1,2,3) TransferMatrix on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@imag.fr", "" );
    return about;

}


template<int Dim>
class TransferMatrix
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef Mesh< Simplex<Dim,1,Dim>,value_type,0> sub1_mesh_type;
    typedef Mesh< Simplex<Dim,1,Dim>,value_type,1> sub2_mesh_type;
    typedef meshes<sub1_mesh_type,sub2_mesh_type> mesh_type;
    typedef bases<Lagrange<Order,Scalar>,Lagrange<Order,Scalar>> basis_type;
    typedef FunctionSpace< mesh_type, basis_type > space_type;

    /**
     * Constructor
     */
    TransferMatrix()
        :
        super(),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {}

    void run();

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
}; // TransferMatrix

template<int Dim> const uint16_type TransferMatrix<Dim>::Order;

template<int Dim>
void
TransferMatrix<Dim>::run()
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );


    auto mesh1 = createGMSHMesh( _mesh=new sub1_mesh_type,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % int(1)).str() ,
                                               _addmidpoint=false, _usenames=false, _shape=shape, _dim=Dim,
                                               _h=meshSize, _xmin=0., _xmax=0.5, _ymin=0., _ymax=1.,
                                               _zmin=0., _zmax=1., _substructuring=true) );

    auto mesh2 = createGMSHMesh( _mesh=new sub2_mesh_type,
                                 _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                 _desc=domain( _name=(boost::format( "%1%-%2%-%3%" ) % shape % Dim % int(2)).str() ,
                                               _addmidpoint=false, _usenames=false, _shape=shape, _dim=Dim,
                                               _h=meshSize, _xmin=0.5, _xmax=1., _ymin=0., _ymax=1.,
                                               _zmin=0., _zmax=1., _substructuring=true) );

    auto mesh = fusion::make_vector(mesh1, mesh2);
    auto Xh = space_type::New( _mesh=mesh );
    auto u = Xh->element();
    auto u1 = u.template element<0>();
    auto u2 = u.template element<1>();

    auto Xh1 = Xh->template functionSpace<0>();
    Xh1->dof()->printDofMarker("dof1.dat");
    auto Xh2 = Xh->template functionSpace<1>();
    Xh2->dof()->printDofMarker("dof2.dat");

    // Slave and Master function spaces
    auto Xh_master = Xh1->trace(markedfaces(mesh1,"EAST"));
    Xh_master->dof()->printDofMarker("master.dat");
    auto Xh_slave = Xh2->trace(markedfaces(mesh2,"WEST"));
    Xh_slave->dof()->printDofMarker("slave.dat");

    saveGMSHMesh(_mesh=Xh_slave->mesh(),_filename="trace_slave.msh");
    saveGMSHMesh(_mesh=Xh_master->mesh(),_filename="trace_master.msh");

    // assemble C_s and C_m
    auto us = Xh_slave->element();
    auto vs = Xh_slave->element();
    auto Cs = M_backend->newMatrix( _trial=Xh_slave, _test=Xh_slave );
    form2( _trial=Xh_slave, _test=Xh_slave, _matrix=Cs ) =
        integrate(_range=elements(Xh_slave->mesh()), _expr=-idt(us)*id(vs) );

    auto um = Xh_master->element();
    auto Cm = M_backend->newMatrix( _trial=Xh_master, _test=Xh_slave );
    form2( _trial=Xh_master, _test=Xh_slave, _matrix=Cm ) =
        integrate( _range=elements(Xh_slave->mesh()), _expr=idt(um)*id(vs) );


    auto g = sin(pi*Px())*cos(pi*Py())*cos(pi*Pz());
    auto g1 = vf::project( Xh1, elements(Xh1->mesh()), g );
    auto g2 = vf::project( Xh2, elements(Xh2->mesh()), g );

    int rsize = Xh->nDof()-u2.extractValuesWithMarker( "WEST", M_backend )->size();

    auto _vec = M_backend->newVector( Xh1 );
    for (size_type i = 0 ; i < g1.size(); ++i)
        _vec->set(i, g1(i) );

    for (size_type i = g1.size() ; i < _vec->size(); ++i)
        _vec->set(i, g2(i) );

    auto u_m = M_backend->newVector( Xh_master );
    auto u_s = M_backend->newVector( Xh_slave );
    auto res = M_backend->newVector( Xh_slave );

    auto Ihm = opInterpolation( _domainSpace=Xh1,
                                _imageSpace=Xh_master,
                                _range=elements(Xh_master->mesh()),
                                _backend=M_backend );

    Ihm->matPtr()->multVector(_vec,u_m);
    Cm->multVector( u_m, res );
    M_backend->solve(_matrix=Cs, _solution=u_s, _rhs=res,_pcfactormatsolverpackage="mumps");
    res->zero();
    res->add(1.,u_m);
    res->add(-1.,u_s);
    std::cout<<"Norm= "<< res->l2Norm() <<"\n";
    BOOST_CHECK_SMALL( res->l2Norm(),1e-11 );

} // TransferMatrix::run


/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( transfer_matrix )
BOOST_AUTO_TEST_CASE( TransferCase )
{
    Application app;

    app.add( new TransferMatrix<2>() );
    app.run();
}

BOOST_AUTO_TEST_SUITE_END()
