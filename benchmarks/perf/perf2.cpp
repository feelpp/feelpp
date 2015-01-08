/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-07

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)

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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-07
 */
//#include <google/profiler.h>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelpoly/polynomialset.hpp>


#include <feel/feelvf/vf.hpp>


using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description myintegralsoptions( "MyIntegrals options" );
    myintegralsoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "nthreads", po::value<int>()->default_value( 2 ), "nthreads" )
    ;
    return myintegralsoptions.add( Feel::feel_options() );
}
inline
AboutData
makeAbout()
{
    AboutData about( "perf2" ,
                     "perf2" ,
                     "0.3",
                     "nD(n=1,2,3) MyIntegrals on simplices or simplex products",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


/**
 * MyIntegrals: compute integrals over a domain
 * \see the \ref ComputingIntegrals section in the tutorial
 * @author Christophe Prud'homme
 */
template<int Dim>
class MyIntegrals
    :
public Simget
{
    typedef Simget super;
public:
    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef FunctionSpace<mesh_type, bases<Lagrange<2, Scalar> >, double> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    MyIntegrals( po::variables_map const& vm, AboutData const& about )
        :
        super(),
        meshSize( doption("hsize") ),
        shape( soption("shape")  ),
        nthreads( ioption("nthreads")  ),
        backend( Backend<double>::build( soption("backend") ) )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    double meshSize;
    std::string shape;
    int nthreads;
    boost::shared_ptr<Backend<double> > backend;
}; // MyIntegrals


template<int Dim>
void
MyIntegrals<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute MyIntegrals<" << Dim << ">\n";
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
MyIntegrals<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    using namespace Feel::vf;

    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "benchmarks/perf2/%1%/%2%/h_%3%/" )
                                       % this->about().appName()
                                       % shape
                                       % meshSize );

#if defined(FEELPP_HAS_TBB)
    /*
     * First we create the mesh
     */
    tbb::tick_count t0 = tbb::tick_count::now();
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name= ( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                _shape=shape,
                                                _dim=Dim,
                                                _h=X[0] ) );

    tbb::tick_count t1 = tbb::tick_count::now();
    double t = ( t1-t0 ).seconds();
    std::cout << "mesh: " << t << "s\n";
    t0 = tbb::tick_count::now();
    mesh->setComponents( MESH_PARTITION| MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    //mesh->setComponents( 0 );
    //ProfilerStart( "/tmp/updateforuse.prof" );
    mesh->updateForUse();
    //ProfilerStop();
    t1 = tbb::tick_count::now();
    t = ( t1-t0 ).seconds();
    std::cout << "update mesh: " << t << "s\n";

    auto Xh = space_type::New( mesh );
    auto u = Xh->element();
    auto M = backend->newMatrix( Xh, Xh );
    u = vf::project( Xh, elements( mesh ), Px()*Px()+Py()*Py()+Pz()*Pz() );

    double overhead = 0;//4.5e-2;
    auto myformexpr = idt( u )*id( u )*Px()*cos( Py() );
    auto myexpr = vf::P()*trans( vf::P() )*sin( Px() )*cos( Py() )*cos( Pz() );
    //auto myexpr = vf::P()*trans(vf::P())*idv(u);
    std::cout << "nb elts in mesh =" << mesh->numElements() << std::endl;
    /*
     * Compute domain Area
     */
    //# marker1 #
    double local_domain_area;

    form2( Xh,Xh,_matrix=M,_init=true );
#if 1
    int n = tbb::task_scheduler_init::default_num_threads();
    double initt;
    std::vector<double> speedup( n );
    {
        std::cout << 1 << " thread" << std::endl;
        tbb::task_scheduler_init init( 1 );
        std::cout << "is_active: " << init.is_active() << "\n";
        tbb::task_scheduler_init init2( 2 );
        std::cout << "is_activ2e: " << init2.is_active() << "\n";
        t0 = tbb::tick_count::now();
        //local_domain_area = integrate( elements(mesh), trace(vf::P()*trans(vf::P()))*idv(u)).evaluate()(0,0);
        auto res  = integrate( elements( mesh ), myexpr, _Q<10>() ).evaluate();
        local_domain_area = res( 0,0 );
        //form2(Xh,Xh,M)=integrate( elements(mesh), myformexpr);
        t1 = tbb::tick_count::now();
        initt = ( t1-t0 ).seconds();
        std::cout << "time: " << initt << " for " << "1 thread" << std::endl;
        speedup[0] = 1;
    }
    std::cout << "------------------------------------------------------------\n";

    for ( int p=ioption("nthreads"); p<=n; ++p )
    {
        std::cout << p << " threads" << std::endl;
        tbb::task_scheduler_init init( p );
        std::cout << "is_active: " << init.is_active() << "\n";
        t0 = tbb::tick_count::now();
        //form2(Xh,Xh,M)=integrate( elements(mesh),myformexpr );
        auto res = integrate( elements( mesh ), myexpr,_Q<10>() ).evaluate();
        local_domain_area = res( 0,0 );
        t1 = tbb::tick_count::now();
        double t = ( t1-t0 ).seconds();
        speedup[p-1] = initt/t;

        std::cout << "time: " << t << " for " << p << " threads speedup=" << speedup[p-1] << std::endl;
        std::cout << "area = " << local_domain_area << std::endl;
        std::cout << "------------------------------------------------------------\n";
    }

#endif
#endif // FEELPP_HAS_TBB
} // MyIntegrals::run

int
main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    //app.add( new MyIntegrals<1>( app.vm(), app.about() ) );
    //app.add( new MyIntegrals<2>( app.vm(), app.about() ) );
    app.add( new MyIntegrals<3>( app.vm(), app.about() ) );
    app.run();
}





