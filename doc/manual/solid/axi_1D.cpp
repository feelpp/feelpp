/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-04

  Copyright (C) 2008 Christophe Prud'homme
  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
   \file elaxi.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-01-04
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelvf/vf.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description elaxioptions( "Elaxi options" );
    elaxioptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.25 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 1.0e+5 ), "coeff for weak Dirichlet conditions" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return elaxioptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "elaxi" ,
                           "elaxi" ,
                           "0.1",
                           "Elasticity axisym  on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007 University Joseph Fourier Grenoble 1" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Vuk Milisic", "developer", "vuk.milisic@imag.fr", "" );
    return about;

}


namespace Feel
{
template<typename A, uint16_type i>
class mytag : public A
{
public:
    static const uint16_type TAG = i;

};
/**
 * Diffussion Advection Reaction Solver
 *
 * solve \f$-\epsilon \Delta u -\beta\cdot\nabla u + \mu u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma_{in}\f$
 */
template<int Order,
         template<uint16_type,uint16_type,uint16_type> class Entity = Simplex>
class Elaxi
    :
public Application
{
    typedef Application super;
public:
    static const uint16_type Dim=2;
    // -- TYPEDEFS --
    static const uint16_type imOrder = 4*Order;

    typedef double value_type;

    typedef Application application_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Entity<Dim,1,Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //typedef mytag<fem::Lagrange<Dim, Order, Vectorial, Continuous, double, Entity>,0> basis_u_type;
    typedef Lagrange<Order, Scalar> basis_scalar_type;
    typedef bases<basis_scalar_type,basis_scalar_type> basis_type;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename element_type::template sub_element<0>::type element_0_type;
    typedef typename element_type::template sub_element<1>::type element_1_type;


    /* export */
    typedef Exporter<mesh_type> export_type;


    Elaxi()
        :
        super(),
        M_backend( backend_type::build( soption("backend") ) ),
        meshSize( doption("hsize") ),
        bcCoeff(  doption("bccoeff") ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
        LOG(INFO) << "[Elaxi] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Elaxi] bccoeff = " << bcCoeff << "\n";
        LOG(INFO) << "[Elaxi] export = " << this->vm().count( "export" ) << "\n";

    }

    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double ,element_type& u );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double bcCoeff;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Elaxi


template<int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Elaxi<Order, Entity>::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    this->changeRepository( boost::format( "doc/manual/solid/%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % entity_type::name()
                            % Order
                            % doption("hsize")
                          );
    /*
     * logs will be in <feel repo>/<app name>/<entity>/P<p>/h_<h>
     */
    this->setLogs();

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="axi_1D",
                                                _shape="hypercube",
                                                _usenames=true,
                                                _ymin=1, _ymax=2,
                                                _h=meshSize ),
                                        _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                        _partitions=this->comm().size()  );

    stats["nelt"] = mesh->elements().size();


    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    space_ptrtype Xh = space_type::New( mesh );

    //Xh->dof()->showMe();
    element_type U( Xh, "u" );
    element_type Uk( Xh, "uk" );
    element_type V( Xh, "v" );
    element_type Phi( Xh, "Phi" );
    Uk.zero();
    element_type poids( Xh, "p" );

    element_0_type u0 = U.template element<0>();
    element_0_type v0 = V.template element<0>();
    element_0_type phi0 = Phi.template element<0>();

    element_1_type u1 = U.template element<1>();
    element_1_type v1 = V.template element<1>();
    element_1_type phi1 = Phi.template element<1>();


    phi0=vf::project( Xh->template functionSpace<0>(),elements( mesh ), Py() );
    phi1.zero();




    uint16_type counter_it_newt=0;
    uint16_type max_it_newt=10000;
    Vector<double>::real_type error=1.0e+10;


    /*
     * Data associated with the simulation
     */
    const double tol = 1e-5;
    const double E = 21*1e5;
    const double sigma = 0.28;
    const double mu = E/( 2*( 1+sigma ) );
    const double lambda = E*sigma/( ( 1+sigma )*( 1-2*sigma ) );
    const double density = 50;
    //    const double gravity = -density*9.81;
    const double gravity = -1.0;
    LOG(INFO) << "lambda = " << lambda << "\n"
          << "mu     = " << mu << "\n"
          << "gravity= " << gravity << "\n";
    std::cout << "lambda = " << lambda << "\n"
              << "mu     = " << mu << "\n"
              << "gravity= " << gravity << "\n";
    /*
     * Construction of the constant right hand side
     *
     * \f$ f = \int_\Omega g * v \f$ where \f$ g \f$ is a vector
     * directed in the \f$ z \f$ direction.
     */

    vector_ptrtype rhs( M_backend->newVector( Xh ) );
    vector_ptrtype vec_gravity( M_backend->newVector( Xh ) );
    vector_ptrtype newt_nl_source_term( M_backend->newVector( Xh ) );
    newt_nl_source_term->zero();


    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();


    LOG(INFO) << "Data Summary:\n";
    size_type pattern = Pattern::COUPLED;

    timers["assembly"].first.restart();
    auto D = M_backend->newMatrix( Xh, Xh );
    std::cout << "====================Newton========================\n---->Start\n";
    LOG(INFO) << "====================Newton========================\n---->Start\n";

    while ( ( error>tol ) && ( counter_it_newt <max_it_newt ) )
    {

        std::cout << "iteration #" << counter_it_newt << "\n";
        std::cout << "error=     " << error << "\n";
        std::cout << "==================================================\n\n";

        LOG(INFO)<<  "iteration #" << counter_it_newt << "\n";
        LOG(INFO)<<  "error=     " << error << "\n";
        LOG(INFO)<<  "==================================================\n\n";




        timers["assembly"].first.restart();

        form2( _test=Xh, _trial=Xh, _matrix=D ) =
            integrate( elements( mesh ), 2.0*(
                           //idt(u1)*id(v1)/Py()
                           idt( u0 )*id( v0 )/Py()
                           +dyt( u0 )*dy( v0 )*Py()
                           +dxt( u0 )*dx( v0 )*Py()
                           +dyt( u1 )*dy( v1 )*Py()
                           +dxt( u1 )*dx( v1 )*Py()
                       ) );


        LOG(INFO) << "[elaxi] matrix local assembly done\n";
        D->close();
        LOG(INFO) << "[elaxi] vector/matrix global assembly done\n";


        /*
         * Construction of the variable right hand side
         *  (depending on the newton iterations)
         */
        //(*rhs)=(*vec_gravity);
        ( *rhs ).zero();
#if 1
        std::cout << "phi0"<< phi0.l2Norm()<<"\n";
        //      phi0.print();
        std::cout << "phi1"<< phi1.l2Norm()<<"\n";
        //      Phi.print();
#endif

        form1( _test=Xh, _vector=newt_nl_source_term ) =
            integrate( elements( mesh ), 2.0*(
                           //idv(phi1)*id(v1)/Py()
                           dyv( phi1 )*dy( v1 )*val( Py() )
                       ) );

#if 0
        newt_nl_source_term->print();
#endif
        rhs->add( -1.0,newt_nl_source_term );
        rhs->close();

        std::cout << "rhs->l2Norm= " << rhs->l2Norm() << "\n";

        std::cout << "----> Block marked dofs\n";

        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( _range=boundaryfaces( mesh ), _element=u0, _rhs=rhs, _expr=constant( 0. ) );
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( _range=boundaryfaces( mesh ), _element=u1, _rhs=rhs, _expr=constant( 0. ) );

        std::cout << "rhs->l2Norm= " << rhs->l2Norm() << "\n";

        error=rhs->l2Norm();


        LOG(INFO) << "[elaxi] starting solve for D\n";
        M_backend->solve( _matrix=D, _solution=U, _rhs=rhs );
        std::cout << "rhs->l2Norm= " << rhs->l2Norm() << "\n";

        u1.zero();

        LOG(INFO) << "[elaxi] solve for D done\n";


        error=rhs->l2Norm();


        /*
          Preparing data for the next step
        */
        std::cout << "U->l2Norm= " << U.l2Norm() << "\n";
        std::cout << "Uk->l2Norm= " << Uk.l2Norm() << "\n";
        Phi+=U;
        Uk+=U;
        counter_it_newt++;
        this->exportResults( counter_it_newt, Uk );
    } //while



} //run


template<int Order, template<uint16_type,uint16_type,uint16_type> class Entity>
void
Elaxi<Order, Entity>::exportResults( double time, element_type& U )
{
    timers["export"].first.restart();


    exporter->step( time )->setMesh( U.functionSpace()->mesh() );
    exporter->step( time )->add( "u0", U.template element<0>() );
    exporter->step( time )->add( "u1", U.template element<1>() );
    exporter->save();

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Elaxi::export


} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    /* change parameters below */
    const int nOrder = 2;

    typedef Feel::Elaxi<nOrder, Simplex> elaxi_type;

    /* assertions handling */
    Feel::Assert::setLog( "elaxi.assert" );

    /* define and run application */
    elaxi_type elaxi;
    elaxi.run();
}


