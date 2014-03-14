/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-06-04

  Copyright (C) 2008-2012 Université Joseph Fourier (Grenoble I)

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
   \file laplacian_multiplicateur_lagrange.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-06-04
 */
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description laplacian_mloptions( "LaplacianLM options" );
    laplacian_mloptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size in domain" )

    ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )

    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return laplacian_mloptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "laplacian_ml" ,
                           "laplacian_ml" ,
                           "0.1",
                           "nD(n=1,2,3) LaplacianLM on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Mourad Ismail", "developer", "mourad.ismail@ujf-grenoble.fr", "" );
    return about;

}

std::pair<std::string,std::string> createRing( int Dim, double h, double rmin, double rmax );

namespace Feel
{
using namespace vf;
/**
 * Fat boundary method for the laplacian
 *
 */
template<int Dim, int Order>
class LaplacianLM
    :
public Application
{
    typedef Application super;
public:

#define Entity Simplex

    // -- TYPEDEFS --
    static const uint16_type imOrder = 2*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<Order, Scalar>,
            Lagrange<0, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;

    typedef typename functionspace_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;

    /*quadrature*/
    template<int IMORDER> struct MyIM : public IM<Dim, IMORDER, value_type, Entity> {};

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /** constructor */
    LaplacianLM();

    /** mesh generation */
    mesh_ptrtype createMesh();

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u, element_type&, element_type& );

private:

    backend_ptrtype M_backend;

    double h;
    double penalisation_bc;

    mesh_ptrtype mesh;

    functionspace_ptrtype Xh;

    export_ptrtype exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // LaplacianLM

template<int Dim, int Order>
LaplacianLM<Dim,Order>::LaplacianLM()
    :
    super(),
    M_backend( backend_type::build( this->vm() ) ),

    // Data
    h( this->vm()["hsize"].template as<double>() ),
    penalisation_bc( this->vm()["penalbc"].template as<value_type>() ),

    // spaces
    Xh(),

    // export
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),

    //
    timers()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }



    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % h
                          );

    LOG(INFO) << "create mesh\n";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "hypercube-%1%" )  % Dim ).str() ,
                                         _usenames=true,
                                         _shape="hypercube",
                                         _h=h ),
                           _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                           _partitions=this->comm().size()  );

    LOG(INFO) << "create space\n";
    Xh = functionspace_type::New( mesh );

    LOG(INFO) << "print space info\n";
    Xh->printInfo();
    LOG(INFO) << "print space 0 info\n";
    Xh->template functionSpace<0>()->printInfo();
    LOG(INFO) << "print space 1 info\n";
    Xh->template functionSpace<1>()->printInfo();

    LOG(INFO) << "Constructor done\n";
}
template<int Dim, int Order>
void
LaplacianLM<Dim, Order>::run()
{
    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    timers["init"].first.restart();

    element_type U( Xh, "u" );
    element_type V( Xh, "v" );
    // trial functions
    element_0_type u = U.template element<0>() ;
    element_1_type lambda = U.template element<1>() ;
    // test functions
    element_0_type v = V.template element<0>() ;
    element_1_type nu = V.template element<1>() ;

    auto g= sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );;
    auto grad_g = vec(
                      +pi*cos( pi*Px() )*cos( pi*Py() ),
                      -pi*sin( pi*Px() )*sin( pi*Py() )
                  );
    auto f = 2*pi*pi*g;

    auto M = M_backend->newMatrix( Xh, Xh );

    form2( _test=Xh, _trial=Xh, _matrix=M ) = integrate( _range=elements( mesh ),
            _expr=gradt( u )*trans( grad( v ) ) + id( u )*idt( lambda ) + idt( u )*id( nu ) + 0*idt( lambda )*id( nu ) );

    double area = integrate( _range=elements( mesh ), _expr=constant( 1.0 ) ).evaluate()( 0,0 );
    double mean = integrate( _range=elements( mesh ), _expr=g ).evaluate()( 0, 0 )/area;
    LOG(INFO) << "int g  = " << mean << "\n";
    vector_ptrtype F( M_backend->newVector( Xh ) );
    form1( Xh, F ) = ( integrate( _range=elements( mesh ), _expr=f*id( v ) )+
                       integrate( _range=boundaryfaces( mesh ), _expr=( trans( grad_g )*N() )*id( v ) ) +
                       integrate( _range=elements( mesh ), _expr=mean*id( nu ) )
                     );

    if ( this->vm().count( "export-matlab" ) && Environment::numberOfProcessors() == 1 )
    {
        M->printMatlab( "M.m" );
        F->printMatlab( "F.m" );
    }

    M_backend->solve( _matrix=M, _solution=U, _rhs=F );

    LOG(INFO) << "lambda = " << lambda( 0 ) << "\n";
    LOG(INFO) << "area   = " << area << "\n";
    LOG(INFO) << "int g  = " << integrate( _range=elements( mesh ), _expr=g ).evaluate()( 0, 0 )/area << "\n";
    LOG(INFO) << "int u  = " << integrate( _range=elements( mesh ), _expr=idv( u ) ).evaluate()( 0, 0 )/area << "\n";
    LOG(INFO) << "error  = " << math::sqrt( integrate( _range=elements( mesh ), _expr=( idv( u )-g )*( idv( u )-g ) ).evaluate()( 0, 0 ) ) << "\n";

    v = vf::project( _space=Xh->template functionSpace<0>(), _range=elements( mesh ), _expr=g );

    element_type E( Xh, "e" );
    element_0_type e = E.template element<0>();
    e = vf::project( _space=Xh->template functionSpace<0>(), _range=elements( mesh ), _expr=idv( u )-g );

    exportResults( U, V, E );

} // LaplacianLM::run


template<int Dim, int Order>
void
LaplacianLM<Dim, Order>::exportResults( element_type& U, element_type& V, element_type& E )

{
    timers["export"].first.restart();

    LOG(INFO) << "exportResults starts\n";
    exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 0 )->add( "u", U.template element<0>() );
    exporter->step( 0 )->add( "exact", V.template element<0>() );
    exporter->step( 0 )->add( "error", E.template element<0>() );
    //exporter->step(0)->add( "l", U.template element<1>() );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // LaplacianLM::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 4;

    typedef Feel::LaplacianLM<nDim, nOrder> laplacian_ml_type;

    /* define and run application */
    laplacian_ml_type laplacian_ml;

    laplacian_ml.run();
}








