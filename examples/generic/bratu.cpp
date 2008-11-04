/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-04-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file bratu.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-04-14
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/backend.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifediscr/operatorlinear.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>
#include <life/lifepoly/polynomialset.hpp>


#include <life/lifevf/vf.hpp>




inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description bratuoptions("Bratu problem options");
    bratuoptions.add_options()
        ("dt", Life::po::value<double>()->default_value( 1 ), "time step value")
        ("ft", Life::po::value<double>()->default_value( 1 ), "final time value")
        ("lambda", Life::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Bratu problem")

        ("order", Life::po::value<int>()->default_value( 2 ), "order of time discretisation")
        ("diff", Life::po::value<double>()->default_value( 1 ), "diffusion parameter")
        ("penal", Life::po::value<double>()->default_value( 10 ), "penalisation parameter")
        ("penalbc", Life::po::value<double>()->default_value( 10 ), "penalisation parameter for the weak boundary conditions")
        ("hsize", Life::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("bctype", Life::po::value<int>()->default_value( 1 ), "0 = strong Dirichlet, 1 = weak Dirichlet")
        ("export", "export results(ensight, data file(1D)")
        ("export-mesh-only", "export mesh only in ensight format")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return bratuoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "bratu" ,
                           "bratu" ,
                           "0.1",
                           "nD(n=1,2,3) Bratu problem on simplices or simplex products",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008 Université Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Life
{
using namespace Life::vf;
/**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
template<int Dim,
         int Order,
         typename Cont,
         template<uint16_type,uint16_type,uint16_type> class Entity,
         template<uint16_type> class FType>
class Bratu
    :
    public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const uint16_type imOrder = 4*Order;

    typedef Bratu<Dim,Order, Cont, Entity, FType> self_type;

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
    typedef Mesh<GeoEntity<entity_type> > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<fem::Lagrange<Dim, 0, Scalar, Discontinuous> > > p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    template<typename Conti = Cont>
    struct space
    {
        /*basis*/
#if 0
        typedef typename mpl::if_<mpl::bool_<Conti::is_continuous>,
                                  mpl::identity<fusion::vector<fem::Lagrange<Dim, Order, FType, Conti, double, Entity> > >,
                                  mpl::identity<fusion::vector<OrthonormalPolynomialSet<Dim, Order, FType, double, Entity> > > >::type::type basis_type;
#else
        typedef fusion::vector<fem::Lagrange<Dim, Order, FType, Conti, double, Entity> > basis_type;

#endif
        /* number of dofs per element */
        static const uint16_type nLocalDof = boost::remove_reference<typename fusion::result_of::at<basis_type,mpl::int_<0> >::type>::type::nLocalDof;

        /*space*/
        typedef FunctionSpace<mesh_type, basis_type, value_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type element_type;
        typedef typename element_type::template sub_element<0>::type element_0_type;
        typedef typename element_type::template sub_element<1>::type element_1_type;
    };
    typedef typename space<Cont>::type functionspace_type;
    typedef typename space<Cont>::ptrtype functionspace_ptrtype;
    typedef typename space<Cont>::element_type element_type;

    typedef OperatorLinear<functionspace_type,functionspace_type> oplin_type;
    typedef boost::shared_ptr<oplin_type> oplin_ptrtype;
    typedef FsFunctionalLinear<functionspace_type> funlin_type;
    typedef boost::shared_ptr<funlin_type> funlin_ptrtype;

    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Entity> im_type;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;

    Bratu( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        M_lambda( this->vm()["lambda"].template as<double>() ),
        M_Xh(),
        exporter( Exporter<mesh_type>::New( this->vm()["exporter"].template as<std::string>() )->setOptions( this->vm() ) ),
        timeSet( new timeset_type( "bratu" ) )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
        exporter->setPrefix( "bratu" );

        if ( this->vm().count( "help" ) )
            {
                std::cout << this->optionsDescription() << "\n";
                return;
            }




        this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                                % this->about().appName()
                                % entity_type::name()
                                % Order
                                % this->vm()["hsize"].template as<double>()
                            );

        mesh_ptrtype mesh = createMesh( meshSize );

        M_Xh = functionspace_ptrtype( space<Cont>::type::New( mesh ) );
    }

    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize, double ymin = 0, double ymax = 1 );

    /**
     * run the convergence test
     */
    void run();


    void updateResidual( const vector_ptrtype& X, vector_ptrtype& R );
    void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    void updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J);

private:



    /**
     * solve the system
     */
    void solve( sparse_matrix_ptrtype& D, element_type& u, vector_ptrtype& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename f1_type, typename f2_type, typename f3_type>
    void exportResults( double time,
                        f1_type& u,
                        f2_type& v,
                        f3_type& e );

private:

    backend_ptrtype M_backend;

    double meshSize;
    double M_lambda;

    functionspace_ptrtype M_Xh;
    oplin_ptrtype M_oplin;
    oplin_ptrtype M_jac;
    funlin_ptrtype M_residual;

    export_ptrtype exporter;
    typename export_type::timeset_ptrtype timeSet;


}; // Bratu

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
typename Bratu<Dim,Order,Cont,Entity,FType>::mesh_ptrtype
Bratu<Dim,Order,Cont,Entity,FType>::createMesh( double meshSize, double ymin, double ymax )
{
    mesh_ptrtype mesh( new mesh_type );
    //mesh->setRenumber( false );

    GmshTensorizedDomain<entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,Entity> td;
    td.setCharacteristicLength( meshSize );
    if ( Dim >=2 )
        td.setY( std::make_pair( ymin, ymax ) );
    std::string fname = td.generate( entity_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );

    return mesh;
} // Bratu::createMesh


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateResidual( const vector_ptrtype& X, vector_ptrtype& R )
{
    boost::timer ti;
    Debug() << "[updateResidual] start\n";
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    im_type im;
    u = *X;
    AUTO( g, constant(0.0) );
    //std::cout << "u = " << u << "\n";
    *M_residual =
        integrate( elements( mesh ), im, + gradv(u)*trans(grad(v)) + M_lambda*exp(idv(u))*id(v) ) +
        integrate( boundaryfaces(mesh), im,
                   //integrate( markedfaces(mesh,1), im,
                   ( - trans(id(v))*(gradv(u)*N())
                     - trans(idv(u))*(grad(v)*N())
                     + penalisation_bc*trans(idv(u))*id(v)/hFace())-
                   g*( - grad(v)*N() + penalisation_bc*id(v)/hFace() )
                   );

    M_residual->close();
    *R = M_residual->container();
    Debug() << "[updateResidual] done in " << ti.elapsed() << "s\n";
                   }
template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    boost::timer ti;
    Debug() << "[updateJacobian] start\n";
    static bool is_init = false;
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    mesh_ptrtype mesh = M_Xh->mesh();
    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );
    u = *X;
    im_type im;
    if ( is_init == false )
        {
            *M_jac = integrate( elements( mesh ), im, M_lambda*(exp(idv(u)))*idt(u)*id(v) );
            is_init = true;
        }
    else
        {
            M_jac->matPtr()->zero();
            *M_jac += integrate( elements( mesh ), im, M_lambda*(exp(idv(u)))*idt(u)*id(v) );
        }
    M_jac->close();
    M_jac->matPtr()->addMatrix( 1.0, M_oplin->mat() );
    J = M_jac->matPtr();
    Debug() << "[updateJacobian] done in " << ti.elapsed() << "s\n";
}
template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::updateResidualJacobian( const vector_ptrtype& X, vector_ptrtype& R, sparse_matrix_ptrtype& J)
{
}

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::run()
{
    using namespace Life::vf;
    mesh_ptrtype mesh = M_Xh->mesh();

    element_type u( M_Xh, "u" );
    element_type v( M_Xh, "v" );

    im_type im;

    value_type penalisation = this->vm()["penal"].template as<value_type>();
    value_type penalisation_bc = this->vm()["penalbc"].template as<value_type>();
    int bctype = this->vm()["bctype"].template as<int>();
    value_type dt = this->vm()["dt"].template as<value_type>();
    value_type ft = this->vm()["ft"].template as<value_type>();
    value_type order = this->vm()["order"].template as<int>();

    M_oplin = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    *M_oplin =
        integrate( elements( mesh ), im, gradt(u)*trans(grad(v)) ) +
        integrate( boundaryfaces(mesh), im,
                   ( - trans(id(v))*(gradt(u)*N())
                     - trans(idt(u))*(grad(v)*N())
                     + penalisation_bc*trans(idt(u))*id(v)/hFace()) );
    M_oplin->close();

    M_jac = oplin_ptrtype( new oplin_type( M_Xh, M_Xh, M_backend ) );
    M_residual = funlin_ptrtype( new funlin_type( M_Xh, M_backend ) );



    M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual, boost::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian, boost::ref( *this ), _1, _2 );

    u = project( M_Xh, elements(mesh), constant(0.) );

    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    vector_ptrtype R( M_backend->newVector( u.functionSpace() ) );
    this->updateResidual( U, R );
    sparse_matrix_ptrtype J;
    this->updateJacobian( U, J );
    solve( J, u, R );

    *U = u;
    this->updateResidual( U, R );
    std::cout << "R( u ) = " << M_backend->dot( U, R ) << "\n";

    exportResults( 0.0, u, u, u );

} // Bratu::run

template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
void
Bratu<Dim, Order, Cont, Entity, FType>::solve( sparse_matrix_ptrtype& D,
                                                   element_type& u,
                                                   vector_ptrtype& F )
{
    vector_ptrtype U( M_backend->newVector( u.functionSpace() ) );
    *U = u;
    M_backend->nlSolve( D, U, F, 1e-10, 10 );
    u = *U;


} // Bratu::solve


template<int Dim, int Order, typename Cont, template<uint16_type,uint16_type,uint16_type> class Entity, template<uint16_type> class FType>
template<typename f1_type, typename f2_type, typename f3_type>
void
Bratu<Dim, Order, Cont, Entity,FType>::exportResults( double time,
                                                     f1_type& U,
                                                     f2_type& V,
                                                     f3_type& E )
{

    Log() << "exportResults starts\n";
    typename timeset_type::step_ptrtype timeStep = timeSet->step( time );
    timeStep->setMesh( U.functionSpace()->mesh() );
    //timeStep->setMesh( this->createMesh( meshSize/2, 0.5, 1 ) );
    //timeStep->setMesh( this->createMesh( meshSize/Order, 0, 1 ) );
    //timeStep->setMesh( this->createMesh( meshSize, 0, 1 ) );
    if ( !this->vm().count( "export-mesh-only" ) )
        {
            timeStep->add( "pid",
                           regionProcess( boost::shared_ptr<p0_space_type>( new p0_space_type( U.functionSpace()->mesh() ) ) ) );


            timeStep->add( "u", U );
            timeStep->add( "v", V );
            timeStep->add( "e", E );
        }
    exporter->save();


    if ( Dim == 1 )
        {
            std::ostringstream fname_u;
            fname_u << "u-" << Application::processId() << "." << boost::format( "%.2f" ) % time << ".dat";
            std::ofstream ofs3( fname_u.str().c_str() );
            typename mesh_type::element_iterator it = U.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            typename mesh_type::element_iterator en = U.functionSpace()->mesh()->endElementWithProcessId( Application::processId() );
            if ( !U.areGlobalValuesUpdated() )
                U.updateGlobalValues();
            for( ; it!=en; ++it )
                {
                    for( size_type i = 0; i < space<Cont>::nLocalDof; ++i )
                        {
                            size_type dof0 = boost::get<0>( U.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                            ofs3 << std::setw( 5 ) << it->id() << " "
                                 << std::setw( 5 ) << i << " "
                                 << std::setw( 5 ) << dof0 << " "
                                 << std::setw( 15 ) << U.globalValue( dof0 ) << " ";
                            value_type a = it->point(0).node()[0];
                            value_type b = it->point(1).node()[0];
                            if ( i == 0 )
                                ofs3 << a;
                            else if ( i == 1 )
                                ofs3 <<  b;
                            else
                                ofs3 <<  a + (i-1)*(b-a)/(space<Continuous>::nLocalDof-1);

                            ofs3 << "\n";

                        }
                }
            ofs3.close();

            std::ostringstream fname_v;
            fname_v << "values-" << Application::processId() << "." << boost::format( "%.2f" ) % time << ".dat";
            std::ofstream ofs2( fname_v.str().c_str() );
            it = V.functionSpace()->mesh()->beginElementWithProcessId( Application::processId() );
            en = V.functionSpace()->mesh()->endElementWithProcessId(  Application::processId() );
            if ( !V.areGlobalValuesUpdated() ) V.updateGlobalValues();
            if ( !E.areGlobalValuesUpdated() ) E.updateGlobalValues();
            for( ; it!=en; ++it )
                {
                    for( size_type i = 0; i < space<Continuous>::nLocalDof; ++i )
                        {
                            size_type dof0 = boost::get<0>( V.functionSpace()->dof()->localToGlobal( it->id(), i ) );
                            ofs2 << std::setw( 5 ) << it->id() << " "
                                 << std::setw( 5 ) << i << " "
                                 << std::setw( 5 ) << dof0 << " "
                                 << std::setw( 15 ) << V.globalValue( dof0 ) << " "
                                 << std::setw( 15 ) << E.globalValue( dof0 ) << " ";
                            value_type a = it->point(0).node()[0];
                            value_type b = it->point(1).node()[0];
                            if ( i == 0 )
                                ofs2 << a;
                            else if ( i == 1 )
                                ofs2 <<  b;
                            else
                                ofs2 <<  a + (i-1)*(b-a)/(space<Continuous>::nLocalDof-1);
                            ofs2 << "\n";

                        }
                }

        }


} // Bratu::export
} // Life




int
main( int argc, char** argv )
{
    using namespace Life;

    /* change parameters below */
    const int nDim = 2;
    const int nOrder = 2;
    typedef Continuous MyContinuity;
    //typedef Discontinuous MyContinuity;
    //typedef Life::Bratu<nDim, nOrder, MyContinuity, SimplexProduct, Scalar> bratu_type;

    typedef Life::Bratu<nDim, nOrder, MyContinuity, Simplex, Scalar> bratu_type;

    /* define and run application */
    bratu_type bratu( argc, argv, makeAbout(), makeOptions() );

    bratu.run();
}





