/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2006-07-04

  Copyright (C) 2006 EPFL
  Copyright (C) 2007,2008 Université Joseph Fourier (Grenoble 1)


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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file bench1.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2006-07-04
 */
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

//#include <boost/test/unit_test.hpp>
//using boost::unit_test::test_suite;

#include <boost/program_options.hpp>
#include <boost/lambda/bind.hpp>

#include <lifeconfig.h>


#include <life/options.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifediscr/functionspace.hpp>
#include <life/lifepoly/im.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>

#include <life/lifealg/vectorublas.hpp>
#include <life/lifealg/matrixublas.hpp>
#include <life/lifealg/backend.hpp>
#include <life/lifevf/vf.hpp>

#include <benchmarks/logs.hpp>

using namespace Life;
using namespace Life::vf;

Life::AboutData
makeAbout()
{
    Life::AboutData about( "bench1" ,
                           "bench1" ,
                           "0.1",
                           "assembly performance",
                           Life::AboutData::License_LGPL,
                           "Copyright (c) 2005,2006 EPFL"
                           "Copyright (c) 2007,2008 Université Joseph Fourier (Grenoble 1)"
                           );

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

Life::po::options_description
makeOptions()
{
    Life::po::options_description desc("Specific options");
    desc.add_options()
        ("dim", Life::po::value<int>()->default_value( 1 ), "dimension (1,2,3)")
        ("hsize", Life::po::value<double>()->default_value( 0.1 ), "element size")
        ;
    return desc.add( Life::life_options() );
}



/*!
 * \class Bench1
 * \brief Benchmark for assembly performance in 1D, 2D and 3D
 *
 * The benchmark is called as follows:
 * \code
 * bench1 --dim=1 --hsize=0.1
 * bench1 --dim=2 --hsize=0.1
 * bench1 --dim=3 --hsize=0.1
 * \endcode
 *
 * For a fixed \p hsize the bench is run in 1D, 2D or 3D (by default
 * 1D)
 */
class Bench1
    :
    public Life::Application
{
public:


    /** @name Typedefs
     */
    //@{

    typedef Life::Application super;

    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    Bench1( int argc,
            char** argv,
            Life::AboutData const& ad,
            Life::po::options_description const& od );

    ~Bench1()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void run();

    //@}



protected:

private:

    template<typename FSType, typename IMType> void R( boost::shared_ptr<FSType> const& Xh, IMType const& );
    template<typename FSType, typename IMType> void D( boost::shared_ptr<FSType> const& Xh, IMType const& );
    template<typename FSType, typename IMType> void DR( boost::shared_ptr<FSType> const& Xh, IMType const& );
    template<typename FSType, typename IMType> void ADR( boost::shared_ptr<FSType> const& Xh, IMType const&, mpl::int_<2> );
    template<typename FSType, typename IMType> void ADR( boost::shared_ptr<FSType> const& Xh, IMType const&, mpl::int_<3> );
    /**
     * 1D performance test
     */
    void run1d();

    /**
     * 2D performance test
     */
    void run2d();

    /**
     * 3D performance test
     */
    void run3d();

    /**
     * dimension independant code
     */
    template<typename MeshType, int Order>
    void bench1( boost::shared_ptr<MeshType> & mesh );

private:

    backend_ptrtype M_backend;
    double meshSize;
};

void
clear( boost::shared_ptr<Mesh<LinearTriangle> >& mesh )
{
    //mesh->cleanElementFaces();
}

void
clear( boost::shared_ptr<Mesh<LinearTetra> >& mesh )
{
    //mesh->cleanElementFaces();
    //mesh->cleanElementEdges();
}

Bench1::Bench1( int argc,
            char** argv,
            Life::AboutData const& ad,
            Life::po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( vm()["hsize"].as<double>() )
{
}

void
Bench1::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    this->changeRepository( boost::format( "/benchmarks/%1%/%2$dD/%3$.3f" )
                             % this->about().appName()
                             % this->vm()["dim"].as<int>()
                             % this->vm()["hsize"].as<double>() );
    this->setLogs();
    switch(  vm()["dim"].as<int>() )
        {
        case 1:
            run1d();
            break;
        case 2:
            run2d();
            break;
        case 3:
            run3d();
            break;
        default:
            std::cout << this->optionsDescription() << "\n";
            return;
        }
}

void
Bench1::run1d()
{
    using namespace Life;

    typedef Mesh<LinearTriangle> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    GmshTensorizedDomain<1,1,1,Simplex> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( "bench11d" );

    ImporterGmsh<mesh_type> import( fname );
    aMesh->accept( import );


    BOOST_LOG( app ) << "run2d starts" << std::endl;
    bench1<mesh_type, 1>( aMesh );
    bench1<mesh_type, 2>( aMesh );
    //bench1<mesh_type, 5>( aMesh );
    //bench1<mesh_type, 8>( aMesh );
    BOOST_LOG( app ) << "run2d ends" << std::endl;

}
void
Bench1::run2d()
{
    using namespace Life;

    typedef Mesh<LinearTriangle> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    GmshTensorizedDomain<2,1,2,Simplex> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( "bench12d" );

    ImporterGmsh<mesh_type> import( fname );
    aMesh->accept( import );


    BOOST_LOG( app ) << "run2d starts" << std::endl;
    bench1<mesh_type, 1>( aMesh );
    bench1<mesh_type, 2>( aMesh );
    //bench1<mesh_type, 5>( aMesh );
    //bench1<mesh_type, 8>( aMesh );
    BOOST_LOG( app ) << "run2d ends" << std::endl;

}
void
Bench1::run3d()
{
    using namespace Life;

    typedef Mesh<LinearTetra> mesh_type;
    boost::shared_ptr<mesh_type> aMesh( new mesh_type );

    GmshTensorizedDomain<3,1,3,Simplex> td;
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( "bench13d" );

    ImporterGmsh<mesh_type> import( fname );
    aMesh->accept( import );

    BOOST_LOG( app ) << "run3d starts" << std::endl;
    bench1<mesh_type, 1>( aMesh );
    bench1<mesh_type, 2>( aMesh );
    BOOST_LOG( app ) << "run3d ends" << std::endl;
}

template<typename FSType, typename IMType>
void
Bench1::R( boost::shared_ptr<FSType> const& Xh, IMType const& im  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    typedef fusion::vector<fem::Lagrange<FSType::nDim, 0, Scalar, Discontinuous> > dp0_basis_type;
    typedef FunctionSpace<typename FSType::mesh_type, dp0_basis_type> dp0_space_type;
    typename dp0_space_type::pointer_type P0h = dp0_space_type::New( Xh->mesh() );
    typename dp0_space_type::element_type w( P0h );
    w = vf::project( P0h, elements(Xh->mesh()), Px() );
    BOOST_LOG( app )<< "quad npts: " << im.nPoints() << std::endl;
    form2(Xh,Xh,M,_init=true);
    boost::timer timer;


    //
    // R
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, idt(u)*id(v) );
    BOOST_LOG( app ) << " o- R<const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, idv(w)*idt(u)*id(v) );
    BOOST_LOG( app ) << " o- R<p0 const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, ((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    BOOST_LOG( app ) << " o-   R<xyz> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, idv(w)*((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    BOOST_LOG( app ) << " o-   R<dp0 xyz> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, val((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    BOOST_LOG( app ) << " o-   R<val xyz> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, idv(w)*val((Px()^(3))+(Py()^(2))*Pz())*idt(u)*id(v) );
    BOOST_LOG( app ) << " o-   R<dp0 val xyz> time : " << timer.elapsed() << std::endl;
}

template<typename FSType, typename IMType>
void
Bench1::D( boost::shared_ptr<FSType> const& Xh, IMType const& im  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    form2(Xh,Xh,M,_init=true);
    BOOST_LOG( app )<< "quad npts: " << im.nPoints() << std::endl;
    boost::timer timer;


    //
    // D
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, gradt(u)*trans(grad(v)) );
    BOOST_LOG( app ) << " o- D<const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, dxt(u)*dx(v)+dyt(u)*dy(v) );
    BOOST_LOG( app ) << " o- D<const2> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, dxt(u)*dx(v)+dyt(u)*dy(v)+dzt(u)*dz(v) );
    BOOST_LOG( app ) << " o- D<const3> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, val((Px()^(3))+(Py()^(2))*Pz())*gradt(u)*trans(grad(v)) );
    BOOST_LOG( app ) << " o-   D<xyz> time : " << timer.elapsed() << std::endl;
}
template<typename FSType, typename IMType>
void
Bench1::DR( boost::shared_ptr<FSType> const& Xh, IMType const& im  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    BOOST_LOG( app )<< "quad npts: " << im.nPoints() << std::endl;
    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

    //
    // DR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, gradt(u)*trans(grad(v))+idt( u )*id( v ));
    BOOST_LOG( app ) << " o- DR<const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im, val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v ) ));
    BOOST_LOG( app ) << " o-   DR<xyz> time : " << timer.elapsed() << std::endl;
}
template<typename FSType, typename IMType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, IMType const& im, mpl::int_<2>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype M( M_backend->newMatrix( Xh, Xh ) );

    BOOST_LOG( app )<< "quad npts: " << im.nPoints() << std::endl;
    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

    //
    // ADR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im,
                                           gradt(u)*trans(grad(v))+idt( u )*id( v ) +
                                           (gradt(u)*vec(constant(1.0),constant(1.0)))*id(v));
    BOOST_LOG( app ) << " o- ADR<const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im,
                                           val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v )) +
                                           (gradt(u)*vec(val(Px()^(3)+Py()^(2)*Pz()),val(Px()^(3)+Py()^(2))))*id(v));
    BOOST_LOG( app ) << " o-   ADR<xyz> time : " << timer.elapsed() << std::endl;
}
template<typename FSType, typename IMType>
void
Bench1::ADR( boost::shared_ptr<FSType> const& Xh, IMType const& im, mpl::int_<3>  )
{
    typename FSType::element_type u( Xh );
    typename FSType::element_type v( Xh );
    sparse_matrix_ptrtype  M( M_backend->newMatrix( Xh, Xh ) );

    BOOST_LOG( app )<< "quad npts: " << im.nPoints() << std::endl;
    form2(Xh,Xh,M,_init=true);
    boost::timer timer;

    //
    // ADR
    //
    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im,
                                           gradt(u)*trans(grad(v))+idt( u )*id( v ) +
                                           (gradt(u)*vec(constant(1.0),constant(1.0),constant(1.0)))*id(v));
    BOOST_LOG( app ) << " o- ADR<const> time : " << timer.elapsed() << std::endl;

    timer.restart();
    form2(Xh,Xh,M) += integrate( elements(Xh->mesh()),  im,
                                           val((Px()^(3))+(Py()^(2))*Pz())*(gradt(u)*trans(grad(v))+idt( u )*id( v )) +
                                           (gradt(u)*vec(val(Px()^(3)+Py()^(2)*Pz()),val(Px()^(3)+Py()^(2)),val(Px()^(3))))*id(v));
    BOOST_LOG( app ) << " o-   ADR<xyz> time : " << timer.elapsed() << std::endl;
}
template<typename MeshType, int Order>
void
Bench1::bench1( boost::shared_ptr<MeshType> & mesh )
{


    const int nDim = MeshType::nDim;
    std::ostringstream femstr;
    femstr << "FEM_PK(" << nDim << "," << Order << ")";
    std::ostringstream geostr;
    geostr << "GT_PK(" << nDim << "," << 1 << ")";

    BOOST_LOG( app ) << "------------------------------------------------------------" << std::endl;
    BOOST_LOG( app ) << "dimension : " << nDim << std::endl;
    BOOST_LOG( app ) << "      fem : " << femstr.str() << std::endl;
    BOOST_LOG( app ) << "      geo : " << geostr.str() << std::endl;
    BOOST_LOG( app ) << "++++++++++++++++++++++++++++++" << std::endl;

    typedef fusion::vector<fem::Lagrange<nDim, Order, Scalar, Continuous> > basis_type;
    typedef FunctionSpace<MeshType, basis_type> space_type;
    typename space_type::pointer_type Xh = space_type::New( mesh );



    //IM<nDim, 2*Order, double, Simplex> im;
    IMSimplex<nDim, 2*Order,double> im;
    //IM_PK<nDim, 4> im;
    //IMSimplex<2,2*Order> im;

    //v.space()->gm()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);
    //v.space()->fe()->setCacheInformation( QDR, mesh->numElements(), QDR_NPTS);

    BOOST_LOG( app ) << "dof : " << Xh->nDof() << std::endl
                     << "elt : " << Xh->mesh()->numElements() << std::endl;


    boost::timer timer;

    R( Xh,  im );
    D( Xh,  IMSimplex<nDim, 2*(Order-1),double>() );
    DR( Xh, im );
    ADR( Xh, im, mpl::int_<nDim>() );
    BOOST_LOG( app ) << "------------------------------------------------------------" << std::endl;

}


int
main( int argc, char** argv )
{
    void init_logs();
    init_logs();

    Bench1 bench1(argc, argv, makeAbout(), makeOptions() );

    bench1.run();
}
