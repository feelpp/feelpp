/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@math.unistra.fr>
       Date: 2012-09-26

  Copyright (C) 2012 Université de Strasbourg

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
   \file stokes_periodic.cpp
   \author Christophe Prud'homme <christophe.prudhomme@math.unistra.fr>
   \date 2012-09-26
 */
#include <boost/foreach.hpp>

#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description periodicoptions( "Periodic Stokes options" );
    periodicoptions.add_options()
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size in domain" )
        ( "wp", Feel::po::value<int>()->default_value( 1 ), "with particules or not" )
        ( "px", Feel::po::value<double>()->default_value( 0. ), "x coordinates of particules" )
        ( "py", Feel::po::value<double>()->default_value( 0. ), "y coordinates of particules" )
        ( "penalbc", Feel::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions" )

        ( "export-matlab", "export matrix and vectors in matlab" )
        ;
    return periodicoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokes_periodic" ,
                           "stokes_periodic" ,
                           "0.1",
                           "2D Periodic Stokes",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2012 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

namespace Feel
{
/**
 * Fat boundary method for the laplacian
 *
 */
template<int Dim, int Order>
class PeriodicStokes
    :
public Simget
{
    typedef Simget super;
public:

#define Entity Simplex

    typedef double value_type;

    /*mesh*/
    typedef Entity<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef bases<Lagrange<Order+1, Vectorial>, Lagrange<Order, Scalar> > basis_composite_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_composite_type, Periodicity<Periodic<>, NoPeriodicity > > functionspace_composite_type;
    typedef boost::shared_ptr<functionspace_composite_type> functionspace_composite_ptrtype;
    //BOOST_MPL_ASSERT( ( boost::is_same<mpl::bool_<functionspace_composite_type::template sub_functionspace<0>::type::element_type::is_periodic>,mpl::bool_<true> > ) );
    //BOOST_MPL_ASSERT( ( boost::is_same<mpl::bool_<functionspace_composite_type::template sub_functionspace<1>::type::element_type::is_periodic>,mpl::bool_<false> > ) );


    typedef typename functionspace_composite_type::element_type element_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /** constructor */
    PeriodicStokes();

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& U );

private:


    double h;
    double penalisation_bc;

    mesh_ptrtype mesh;

    functionspace_composite_ptrtype Xh;

    node_type translat;

    export_ptrtype exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Periodic

template<int Dim, int Order>
PeriodicStokes<Dim,Order>::PeriodicStokes()
    :
    super(),

    // Data
    h( doption("hsize") ),
    penalisation_bc( this->vm()["penalbc"].template as<value_type>() ),

    // spaces
    Xh(),

    // translation for pertiodicity
    translat(2),

    // export
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),

    //
    timers()
{
    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % h
                          );

    LOG(INFO) << "create mesh\n";
    const std::string shape = "hypercube";
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                         _usenames=false,
                                         _shape=shape,
                                         _dim=Dim,
                                         _h=h,
                                         _xmin=-1, _ymin=-1, _xmax=1, _ymax=1 ) );


    LOG(INFO) << "create space\n";
    // node_type trans(2);
    translat[0]=0;
    translat[1]=2;
    Xh = functionspace_composite_type::New( _mesh=mesh, _periodicity=periodicity( Periodic<>(2, 4, translat), NoPeriodicity() ) );

    LOG(INFO) << "Xh print space info\n";
    Xh->printInfo();

    LOG(INFO) << "Constructor done\n";
}

template<int Dim, int Order>
void
PeriodicStokes<Dim, Order>::run()
{
    timers["init"].first.restart();

    auto U = Xh->element( "u" );
    auto u = U.template element<0>();
    auto p = U.template element<1>();
    auto V = Xh->element( "v" );
    auto v = V.template element<0>();
    auto q = V.template element<1>();


    auto M = backend()->newMatrix( Xh, Xh );
    auto F = backend()->newVector( Xh );

    auto a = form2( Xh, Xh, _matrix=M );
    a = integrate( _range=elements( mesh ), _expr=trace(gradt( u )*trans( grad( v ) ) ));
    if ( ioption("wp") )
    {
        double px = doption("px");
        double py = doption("py");

        //auto dist2center = norm2(P()-cst(px)*oneX()-cst(py)*oneY());
        auto dist= norm2(P()-cst(px)*oneX()-cst(py)*oneY());
        auto distMinusTrans = norm2(P()-(cst(px)-translat[0])*oneX()-(cst(py)-translat[1])*oneY());
        auto distPlusTrans  = norm2(P()-(cst(px)+translat[0])*oneX()-(cst(py)+translat[1])*oneY());

        auto dist2center = min( dist, min(distMinusTrans, distPlusTrans) );

        a += integrate( _range=elements( mesh ), _expr=chi(dist2center < 0.25)*1e5*trace(gradt( u )*trans( grad( v ) ) ));
    }
    a+= integrate( _range=elements( mesh ), _expr=-idt(p)*div(v)+id(q)*divt(u) );
    a+= integrate( _range=elements( mesh ), _expr=1e-6*idt(p)*id(q) );

    auto b = form1( Xh, _vector=F );

    std::list<int> bdys = {1,3};
    BOOST_FOREACH( auto bdy, bdys )
    {
        a += integrate( _range=markedfaces( mesh, bdy ),
                        _expr=( trans(+idt(p)*N()-gradt( u )*N())*id( v )
                                +trans(+id(q)*N()-grad( v )*N())*idt( u )
                                + penalisation_bc*trans(id( u ))*idt( v )/hFace() ) );

        b += integrate( _range=markedfaces( mesh, bdy ),
                        _expr=((bdy==1)?1:-1)*trans(oneY())*(id(q)*N()-grad( v )*N()+ penalisation_bc*id( u )/hFace() ) );
    }
    backend()->solve( _matrix=M, _solution=U, _rhs=F );

    exportResults( U );

} // PeriodicStokes::run

template<int Dim, int Order>
void
PeriodicStokes<Dim, Order>::exportResults( element_type& U )
{
    timers["export"].first.restart();

    LOG(INFO) << "exportResults starts\n";

    exporter->step( 1. )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 1. )->add( {"u","p"}, U );

    // position of the particle
    auto u = U.template element<1>();
    double px = doption("px");
    double py = doption("py");
    auto dist= norm2(P()-cst(px)*oneX()-cst(py)*oneY());
    auto distMinusTrans = norm2(P()-(cst(px)-translat[0])*oneX()-(cst(py)-translat[1])*oneY());
    auto distPlusTrans  = norm2(P()-(cst(px)+translat[0])*oneX()-(cst(py)+translat[1])*oneY());
    auto dist2center = min( dist, min(distMinusTrans, distPlusTrans) );
    u = project( _space=u.functionSpace(), _range=elements(u.functionSpace()->mesh()), _expr= chi(dist2center < 0.25));
    exporter->step( 1. )->add( "particule", u );

    exporter->save();
    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // PeriodicStokes::export
} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Application app;
    app.add( new PeriodicStokes<2,1>() );
    app.run();
}










