/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@imag.fr>
       Date: 2013-04-24

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
   \file test_sub_matrix.cpp
   \author Abdoulaye Samake<abdoulaye.samake@imag.fr>
   \date 2013-04-24
*/

#define BOOST_TEST_MODULE test_sub_matrix
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>


/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testsubmatrixoptions( "TestSubmatrix options" );
    testsubmatrixoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain" )
        ;
    return testsubmatrixoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_sub_matrix" ,
                     "test_sub_matrix" ,
                     "0.2",
                     "nD(n=1,2,3) TestSubmatrix on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Universite de Strasbourg" );

    about.addAuthor( "Abdoulaye Samake", "developer", "abdoulaye.samake@imag.fr", "" );
    return about;

}

template<int Dim>
class TestSubMatrix
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 2;
    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef Simplex<Dim,1,Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;
    typedef Exporter<mesh_type> export_type;

    /**
     * Constructor
     */
    TestSubMatrix()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
    {
    }

    void run();

private:

    backend_ptrtype M_backend;
    double meshSize;
    std::string shape;
}; // TestSubMatrix

template<int Dim> const uint16_type TestSubMatrix<Dim>::Order;


template<int Dim>
void
TestSubMatrix<Dim>::run()
{

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%/h_%5%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % Order
                                       % meshSize );


    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                              _addmidpoint=false,
                                              _usenames=false,
                                              _shape=shape,
                                              _dim=Dim,
                                              _h=meshSize,
                                              _xmin=0.,
                                              _xmax=1.,
                                              _ymin=0.,
                                              _ymax=1.,
                                              _zmin=0.,
                                              _zmax=1.,
                                              _substructuring=true) );

    auto Xh = space_type::New( _mesh=mesh );
    auto u = Xh->element();
    auto v = Xh->element();

    auto A = M_backend->newMatrix( _trial=Xh, _test=Xh ) ;

    form2( _trial=Xh, _test=Xh, _matrix=A ) =
        integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );

    A->close();

    int nr = Xh->nDof()/3;
    int nc = Xh->nDof()/2;

    std::vector<size_type> rows(nr);
    std::vector<size_type> cols(nc);

    for (size_type i = 0; i < nr; ++i)
        rows[i] = i;

    for (size_type i = 0; i < nc; ++i)
        cols[i] = i;

    auto SubA = M_backend->newMatrix(nr,nc,nr,nc);
    A->createSubmatrix(*SubA,rows,cols);

    for (size_type i= 0; i < nr; ++i)
    {
        for (size_type j = 0; j < nc; ++j)
        {
            BOOST_CHECK_EQUAL( SubA->operator()(i,j), A->operator()(i,j) );
        }
    }

} // TestSubMatrix::run


/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );
BOOST_AUTO_TEST_SUITE( submatrix )
BOOST_AUTO_TEST_CASE( SubmatrixCase )
{
    Application app;

    //app.add( new TestSubMatrix<1>() );
    app.add( new TestSubMatrix<2>() );
    //app.add( new TestSubMatrix<3>() );
    app.run();
}

BOOST_AUTO_TEST_SUITE_END()
