/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-11-23

  Copyright (C) 2006,2007 Universit� Joseph Fourier (Grenoble I)

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file laplacian.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-11-23
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/gausslobatto.hpp>


#include <feel/feelvf/vf.hpp>


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
        ("diff", Feel::po::value<double>()->default_value( 1 ), "diffusion parameter")
        ("gamma", Feel::po::value<double>()->default_value( 10 ), "jump penalisation parameter")
        ("delta", Feel::po::value<double>()->default_value( 0 ), "lifting operator penalisation parameter")
        ("theta", Feel::po::value<double>()->default_value( 1 ), "theta=1: symmetric, theta=-1: non-symmetric")
        ("alpha", Feel::po::value<double>()->default_value( 3 ), "Regularity coefficient for example 2")
        ("example", Feel::po::value<int>()->default_value( 1 ), "Example number")
        ("anisomesh", Feel::po::value<int>()->default_value( 0 ), "0: using normal mesh generation, 1: using anisotropic mesh containing 2 elements")
//        ("f", Feel::po::value<double>()->default_value( 1 ), "forcing term")
//        ("g", Feel::po::value<double>()->default_value( 0 ), "boundary term")
        ("hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence")
        ("hvisu", Feel::po::value<double>()->default_value( 0.05 ), "first h value to start convergence")
        ("export", "export results(ensight, data file(1D)")
        ("export-matlab", "export matrix and vectors in matlab" )
        ;
    return laplacianoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "laplacian_hpOPT" ,
                            "laplacian_hpOPT" ,
                            "0.2",
                            "nD(n=1,2,3) Laplacian on simplices or simplex products",
                            Feel::AboutData::License_GPL,
                            "Copyright (c) 2006, 2007 Universit� Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


namespace Feel
{
using namespace vf;

/**
 * Laplacian Solver using discontinous approximation spaces
 *
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 */
template<int Order>
class Laplacian
    :
        public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
	// Define also more accurate for norm computations since a singularity is in the domain
    static const uint16_type Dim = 2;
    static const uint16_type imOrder = 2*Order;
	static const uint16_type imOrder_norm = 5*Order;

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    typedef Hypercube<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*discontinuous basis*/
    typedef fusion::vector<Lagrange<Order, Scalar, Discontinuous, PointSetGaussLobatto> > basis_type;
    //    typedef fusion::vector<fem::Dubiner<Dim, Order, Scalar, Discontinuous, value_type, Hypercube>,
    //						   fem::Dubiner<Dim, Order, Vectorial, Discontinuous, value_type, Hypercube>
    //						   > basis_type;
	/*discontinuous space*/
	typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
	typedef boost::shared_ptr<space_type> space_ptrtype;
	typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    typedef fusion::vector<Lagrange<Order, Vectorial, Discontinuous, PointSetGaussLobatto> > vectorial_basis_type;
    typedef FunctionSpace<mesh_type, vectorial_basis_type> vectorial_space_type;
	typedef boost::shared_ptr<vectorial_space_type> vectorial_space_ptrtype;
	typedef typename vectorial_space_type::element_type vectorial_element_type;
    typedef boost::shared_ptr<vectorial_element_type> vectorial_element_ptrtype;
//	typedef typename element_type::template sub_element<0>::type element_0_type;
//	typedef typename element_type::template sub_element<1>::type element_1_type;

#if 0
        typedef typename mpl::if_<mpl::bool_<Conti::is_continuous>,
                                  mpl::identity<fusion::vector<Lagrange<Order, FType> > >,
                                  mpl::identity<fusion::vector<OrthonormalPolynomialSet<Order, FType> > > >::type::type basis_type;
#endif
    /*continuous basis*/
	typedef fusion::vector<Lagrange<Order, Scalar> > basis_type_cont;

	/*continuous space*/
    typedef FunctionSpace<mesh_type, basis_type_cont, value_type> type_cont;
    typedef boost::shared_ptr<type_cont> ptrtype_cont;
    typedef typename type_cont::element_type element_type_cont;


    /*quadrature*/
    //typedef IM_PK<Dim, imOrder, value_type> im_type;
    typedef IM<Dim, imOrder, value_type, Hypercube> im_type;
	typedef IM<Dim, imOrder_norm, value_type, Hypercube> im_type_norm;

    /* export */
    typedef Exporter<mesh_type> export_type;

    /** constructor */
    Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    /** create the mesh using mesh size \c meshSize */
    mesh_ptrtype createMesh( double meshSize );

    /** run the application */
    void run();

private:

    /**
     * solve the system
     */
    template<typename elem_type>
    void solve( sparse_matrix_ptrtype const& D, elem_type& u, vector_ptrtype const& F );


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u,
                        element_type_cont& v,
                        element_type_cont& e );


	void writeResults( value_type e1,
                       value_type e2,
                       value_type e3,
                       value_type e4,
                       value_type e5,
					   value_type& theta,
					   value_type& delta,
					   value_type& gamma,
					   int& anisomesh
					   );

private:

    backend_ptrtype M_backend;
    double meshSize;
    mesh_ptrtype mesh;

    space_ptrtype Xh;
    ptrtype_cont  Xch;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;

}; // Laplacian

template<int Order>
Laplacian<Order>::Laplacian( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers()
{
/*	int anisomesh = this->vm()["anisomesh"].template as<int>();
	// This is a convention to separate folders for data of the computations
	if (anisomesh) { meshSize = 2; }
*/
	if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }

    int example = this->vm()["example"].template as<int>();

    if( example==3 ) {
		this->changeRepository( boost::format( "%1%/Ex%2%/alpha_%3%/%4%/P%5%/%6%/" )
                                % this->about().appName()
                                % example
                                % this->vm()["alpha"].template as<double>()
                                % entity_type::name()
                                % Order
                                % this->vm()["hsize"].template as<double>()
                                );
	} else {
		this->changeRepository( boost::format( "%1%/Ex%2%/%3%/P%4%/%5%/" )
                                % this->about().appName()
                                % example
                                % entity_type::name()
                                % Order
                                % this->vm()["hsize"].template as<double>()
                                );
	}

	Log() << "[Laplacian] hsize = " << meshSize << "\n";
    Log() << "[Laplacian] export = " << this->vm().count("export") << "\n";

    mesh = createMesh( meshSize );
    Xh = space_type::New( mesh );
    Log() << "Xh information\n";
    Xh->printInfo();
    Log() << "Xch information\n";
    Xch = type_cont::New( mesh );
    Xch->printInfo();


}

template<int Order>
typename Laplacian<Order>::mesh_ptrtype
Laplacian<Order>::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptrtype _mesh( new mesh_type );
	std::string mesh_name;

	int anisomesh = this->vm()["anisomesh"].template as<int>();

	if(anisomesh==0) {
		GmshHypercubeDomain td(entity_type::nDim,entity_type::nOrder,entity_type::nRealDim,entity_type::is_hypercube);
		td.setCharacteristicLength( meshSize );
		td.setX( std::make_pair( -1, 1 ) );
		td.setY( std::make_pair( -1, 1 ) );
		mesh_name = td.generate( entity_type::name().c_str() );
	} else {
		// Todo: load anisotropic mesh with two elements

		// Todo: copy "Hypercube_2_1_aniso.msh"-file from ../../../../.. (where this file should be copied once) to the actual repository if file is not yet existing
		mesh_name = "Hypercube_2_1_aniso.msh";
	}

	ImporterGmsh<mesh_type> import( mesh_name );
	_mesh->accept( import );

	timers["mesh"].second = timers["mesh"].first.elapsed();
    Log() << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    Log() << "[Laplacian] meshtype = " << anisomesh << "\n";
    return _mesh;
} // Laplacian::createMesh


template<int Order>
void
Laplacian<Order>::run()
{


    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;



    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();

    element_type u( Xh, "U" );
    element_type v( Xh, "V" );
    /*    element_0_type u( U.template element<0>() );
          element_1_type tau( U.template element<1>() );
          element_0_type v( V.template element<0>() );
          element_1_type phi( V.template element<1>() );*/

    timers["init"].second = timers["init"].first.elapsed();

    /*
     * a quadrature rule for numerical integration
     */
    im_type im;
	im_type_norm im_norm;

    value_type gamma = this->vm()["gamma"].template as<value_type>();
    value_type delta = this->vm()["delta"].template as<value_type>();
    value_type theta = this->vm()["theta"].template as<value_type>();
    value_type diff = this->vm()["diff"].template as<double>();
	//value_type hsize = this->vm()["hsize"].template as<value_type>();
    value_type alpha = this->vm()["alpha"].template as<value_type>();
	int example = this->vm()["example"].template as<int>();
	int anisomesh = this->vm()["anisomesh"].template as<int>();

	value_type pi = 4.0*math::atan(1.0);

	// ONLY 2D
	AUTO( g, sin(pi*Px())*sin(pi*Py())* chi(example==1) + chi(example==3)*( (1-Px()*Px())*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0)) ) );
	AUTO( grad_g, ((pi*cos(pi*Px())*sin(pi*Py()))*oneX() + (pi*sin(pi*Px())*cos(pi*Py()))*oneY())* chi(example==1)
          + chi(example==3)*(
                             ( alpha*Px()*(1-Px()*Px())*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0-1.0))
                               -2*Px()*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0))
                               )*oneX()
                             +( alpha*Py()*(1-Px()*Px())*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0-1.0))
                                -2*Py()*(1-Px()*Px())*pow((Px()*Px()+Py()*Py()),(alpha/2.0))
                                )*oneY()
                             )
          );
	AUTO( f,        2*pi*pi*sin(pi*Px())*sin(pi*Py())* chi(example==1) + chi(example==3)*(
                                                                                          2*((1-Px()*Px())+(1-Py()*Py()))*pow((Px()*Px()+Py()*Py()),(alpha/2.0))
                                                                                          -2*alpha*(1-Px()*Px())*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0-1.0))
                                                                                          +4*alpha*pow((Px()*Px()+Py()*Py()),(alpha/2.0))
                                                                                          -8*alpha*(Px()*Px()*Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0-1.0))
                                                                                          - alpha*(alpha-2)*(1-Px()*Px())*(1-Py()*Py())*pow((Px()*Px()+Py()*Py()),(alpha/2.0-1.0))
                                                                                          )
          );


    /*    AUTO( g, val( chi(example==1)*sin(pi*Px())*sin(pi*Py())
          +chi(example==2)*(1-cos()*cos*(vf::pow(abs(Px()),alpha)))
          ));
          AUTO( grad_g, val( chi(example==1)*(
          pi*cos(pi*Px())*sin(pi*Py())*oneX() + pi*sin(pi*Px())*cos(pi*Py())*oneY() )
          +chi(example==2)*(
          alpha*(vf::pow(abs(Px()),alpha-1)))*oneX()

          ));
          AUTO( f, 		val( 2*pi*pi*sin(pi*Px())*sin(pi*Py()) - chi(example==2)*(alpha*(alpha-1)*(vf::pow(abs(Px()),alpha-2)))
          ));*/

    vector_ptrtype F( M_backend->newVector( Xh ) );
    timers["assembly"].first.restart();


    form1( Xh, F, _init=true )  = integrate( _range=elements(mesh), _quad=im_norm, _expr=f*id(v) );
	// ATTENTION: MESH MIGHT HAVE NEGATIVE X-COORDINATES (1E-18)
    form1( Xh, F ) += integrate( _range=boundaryfaces(mesh), _quad=im, _expr=g*( - theta*grad(v)*N() + gamma*id(v)*Order*Order/hFace() ) );
    //    form( Xh, F, false ) += integrate( boundaryfaces(mesh), im, -g*trans(N())*(id(phi)) );

    if ( this->vm().count( "export-matlab" ) )
        F->printMatlab( "F.m" );
    timers["assembly"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    sparse_matrix_ptrtype D( M_backend->newMatrix( Xh, Xh ) );

    timers["assembly"].first.restart();
    size_type pattern = (Pattern::COUPLED|Pattern::EXTENDED );
    form2( Xh, Xh, D, _init=true, _pattern=pattern ) = integrate( _range=elements(mesh), _quad=im, _expr=( diff*gradt(u)*trans(grad(v)))
                                                                  //                                                             +trans(idt(tau))*id(phi)
                                                                  //                                                             +0.*idt(u)*div(phi)
                                                                  //                                                             +0.*divt(tau)*id(v)

                                                                  );
    timers["assembly_D_elements"].second += timers["assembly"].first.elapsed();
    timers["assembly"].first.restart();
	// -----------------------------------------
	// TODO: Write correct lifting stabilization
	// -----------------------------------------
    form2( Xh, Xh, D ) +=integrate( internalfaces(mesh),
                                    // - {grad(u)} . [v]
                                    - averaget(gradt(u))*jump(id(v))
                                    // - theta*[u] . {grad(v)}
                                    - theta * average(grad(v))*jumpt(idt(u))
                                    // gamma*[u] . [v]*p/h_face
                                    + gamma * (trans(jumpt(idt(u)))*jump(id(v)) )*Order*Order/hFace()
                                    // delta*{tau} . [v]
                                    //                                        + delta * (trans(averaget(idt(tau)))*jump(id(v)) )
                                    // - [u] . {phi}
                                    //                                        - trans(jumpt(idt(u)))*average(id(phi))
                                    );
    timers["assembly_D_internalfaces"].second += timers["assembly"].first.elapsed();
    timers["assembly"].first.restart();
    form2( Xh, Xh, D ) += integrate( boundaryfaces(mesh),
                                     - trans(id(v))*(gradt(u)*N())
                                     - theta * trans(idt(u))*(grad(v)*N())
                                     + gamma * (trans(idt(u))*id(v))*Order*Order/hFace()
                                     //										  + delta * (trans(idt(tau))*id(v)*N() )
                                     //										  - trans(idt(u)*N())*id(phi)
                                     );
    timers["assembly_D_boundaryfaces"].second += timers["assembly"].first.elapsed();
    timers["assembly"].first.restart();

    sparse_matrix_ptrtype Mdelta;
    if ( delta > 0 )
        {
            // TODO: pass graph of D
            Mdelta = sparse_matrix_ptrtype( M_backend->newMatrix( Xh, Xh ) );
            //size_type pattern = (Pattern::COUPLED|Pattern::EXTENDED );
            //form2( Xh, Xh, Mdelta, _init=true, _pattern=pattern );

            vectorial_space_ptrtype Wh( vectorial_space_type::New( mesh ) );

            std::vector<bool> face_done( mesh->numFaces() );
            std::fill( face_done.begin(), face_done.end(), false );
            typename mesh_type::face_iterator it, en;
            it = mesh->beginFace();
            en = mesh->endFace();

            uint16_type nLocalDof = Xh->fe()->nLocalDof;
            std::vector<typename mesh_type::element_type> localelts( 2 );
            for ( ; it!=en; ++it )
                {
                    size_type face_id = it->id();

                    boost::shared_ptr<mesh_type> localmesh ( new mesh_type );
                    std::vector<int> lface;
                    std::vector<size_type> lelt;
                    localelts.clear();

                    Log() << "\t Dealing with face " << face_id << " is on boundary? " << it->isOnBoundary() << "\n";
                    // create submesh : the two elements
                    // associated to the face
                    localelts.push_back( it->element0() );
                    lelt.push_back( it->ad_first() );
                    lface.push_back( it->pos_first() );
                    if ( it->isConnectedTo1() )
                        {
                            localelts.push_back( it->element1() );
                            lelt.push_back( it->ad_second() );
                            lface.push_back( it->pos_second() );
                        }

                    // be careful do not renumber the mesh we want to keep the correspondance between the local and global
                    // mesh below
                    if ( it->isOnBoundary() )
                        Log() << "\t\t creating submesh from elements (" << localelts[0].id()  << ")\n";
                    else
                        Log() << "\t\t creating submesh from elements (" << localelts[0].id() << "," << localelts[1].id() << ")\n";
                    mesh->createSubmesh( *localmesh, localelts.begin(), localelts.end() );
                    Log() << "\t\t creating submesh done. Num elements: " << localmesh->numElements() << "\n";


                    vectorial_space_ptrtype Wvh = vectorial_space_type::New( localmesh );
                    space_ptrtype Wsh = space_type::New( localmesh );

                    FEEL_ASSERT
                        ( ublas::norm_frobenius( localmesh->element( 0 ).G() -
                                                 it->element0().G() ) < 1e-10 )

                        ( it->id() )( localmesh->element( 0 ).id() ).error( "invalid local element");
                    if ( !it->isOnBoundary() )
                        {
                            FEEL_ASSERT
                                ( ublas::norm_frobenius( localmesh->element( 1 ).G() -
                                                         it->element1().G() ) < 1e-10 )

                                ( it->id() )( localmesh->element( 1 ).id() ).error( "invalid local element for internal facexs");
                        }
                    std::map<size_type,std::vector<vectorial_element_ptrtype> > wlocal;
                    std::map<size_type,std::vector<vectorial_element_ptrtype> > wglobal;
                    element_type wg( Xh, "w" );
                    element_type w( Wsh, "w" );
                    element_type v( Wsh, "v" );
                    vectorial_element_type wv( Wvh, "w" );
                    vectorial_element_type vv( Wvh, "v" );

                    sparse_matrix_ptrtype Mloc( M_backend->newMatrix( Wvh, Wvh ) );
                    // mass matrix for L2 Projection
                    form2( Wvh, Wvh, Mloc, _init=true ) = integrate( elements( localmesh ), trans(idt( wv )) * id( vv ) );
                    Mloc->close();
                    if ( this->vm().count( "export-matlab" ) )
                        {
                            std::ostringstream ostr;
                            ostr << "Mloc-" << face_id << ".m";
                            Mloc->printMatlab( ostr.str() );
                        }
                    vector_ptrtype Floc( M_backend->newVector( Wvh ) );

                    mesh_type::element_iterator itl, enl;
                    //typename std::vector<typename mesh_type::element_type>::iterator itl, enl;
                    boost::tie( itl, enl ) = localmesh->elementsRange();
                    FEEL_ASSERT( std::distance( itl, enl ) ==  localmesh->numElements() )( localmesh->numElements() )( std::distance( itl, enl ) ).error( "invalid element range");

                    // for each basis function in each element
                    // associated to the face, construct the
                    // associated L2 projection which is first stored
                    // in wlocal[elt_id][i] then in
                    // wglobal[elt_id][i],
                    //
                    // BE CAREFULL wlocal and wglobal must encode the
                    // global id of the element and the basis function
                    // id hence the usage of std::map
                    for( size_type n = 0; itl != enl; ++ itl, ++n )
                        {
                            wlocal[lelt[n]].resize( nLocalDof );
                            wglobal[lelt[n]].resize( nLocalDof );
                            for( uint16_type i = 0; i < nLocalDof; ++i )
                                {
                                    uint16_type k = boost::get<0>( Wsh->dof()->localToGlobal( itl->id(), i , 0 ) );
                                    wlocal[lelt[n]][i] = vectorial_element_ptrtype( new vectorial_element_type( Wvh, "wlocal" ) );
                                    wglobal[lelt[n]][i] = vectorial_element_ptrtype( new vectorial_element_type( Wh, "wglobal" ) );
                                    w.zero();
                                    w( k ) = 1;
                                    FEEL_ASSERT( w.l2Norm() == 1 )( lelt[n] )( i )( w.l2Norm() )( k ).error( "invalid local basis function" );

                                    FEEL_ASSERT
                                        ( ( ublas::norm_frobenius( itl->face( lface[n] ).G() - it->G() ) < 1e-10 ) &&
                                          ( ublas::norm_frobenius( localmesh->face( itl->face( lface[n] ).id() ).G() - it->G() ) < 1e-10 ) )
                                        ( itl->face( lface[n] ).G() )
                                        ( localmesh->face( itl->face( lface[n] ).id() ).G() )
                                        ( it->G() ).error( "invalid face" );


                                    if ( it->isOnBoundary() )
                                        {
                                            FEEL_ASSERT( localmesh->face( itl->face( lface[n] ).id() ).isOnBoundary() ).error( "invalid local boundary face" );
                                            form1( Wvh, Floc, _init=true ) = integrate( _range=idedfaces( localmesh, itl->face( lface[n] ).id() ), _quad=im,
                                                                                        _expr=trans(idv( w )*N())*id( vv ) );
                                        }
                                    else
                                        {
                                            FEEL_ASSERT( ( std::distance( localmesh->beginInternalFace(),
                                                                          localmesh->endInternalFace() ) == 1 ) &&
                                                         ( ublas::norm_frobenius( localmesh->beginInternalFace()->G() -
                                                                                  it->G() ) < 1e-10 ) ).error( "invalid local internal face");
                                            //form1( Wvh, Floc, _init=true ) = integrate( idedfaces( localmesh, itl->face( lface[n] ).id() ), im,
                                            form1( Wvh, Floc, _init=true ) = integrate( _range=internalfaces( localmesh ), _quad=im,
                                                                                        _expr=trans(jumpv( idv( w ) ))*average( id( vv ) ) );
                                        }
                                    Floc->close();
                                    Log() << "Floc element " << lelt[n] << " w_" << i << " flocnorm= " << Floc->l2Norm() << "\n";
                                    if ( this->vm().count( "export-matlab" ) )
                                        {
                                            std::ostringstream ostr;
                                            ostr << "Floc-" << lelt[n] << "-" << i << ".m";
                                            Floc->printMatlab( ostr.str() );
                                        }
                                    // construct the L2 projection, we should use a direct solver here and actually
                                    // construct the LU factorisation prior to the solve
                                    solve( Mloc, *wlocal[lelt[n]][i], Floc );
                                    Log() << "element " << lelt[n] << " basis_" << i << " norm= " << wlocal[lelt[n]][i]->l2Norm() << "\n";
                                    if ( this->vm().count( "export-matlab" ) )
                                        {
                                            std::ostringstream ostr;
                                            ostr << "wlocal-" << lelt[n] << "-" << i << ".m";
                                            wlocal[lelt[n]][i]->printMatlab( ostr.str() );
                                        }

                                    wglobal[lelt[n]][i]->zero();
                                    for ( size_type e = 0; e < lelt.size(); ++e )
                                        for( uint16_type j = 0; j < nLocalDof; ++j )
                                            {
                                                for( uint16_type c = 0; c < 2; ++c )
                                                    {
                                                        size_type klocal = boost::get<0>( Wvh->dof()->localToGlobal( e, j , c ) );
                                                        size_type kglobal = boost::get<0>( Wh->dof()->localToGlobal( lelt[e], j , c ) );
                                                        (*wglobal[lelt[n]][i])( kglobal ) = (*wlocal[lelt[n]][i])( klocal );

                                                    }
                                            }
                                    Log() << "element " << lelt[n] << " global basis_" << i << " norm= " << wglobal[lelt[n]][i]->l2Norm() << "\n";
#if 0
                                    FEEL_ASSERT( wlocal[lelt[n]][i]->l2Norm() > 1e-10 &&
                                                 wglobal[lelt[n]][i]->l2Norm() > 1e-10 )
                                        ( lelt[n] )( i )
                                        ( wlocal[lelt[n]][i]->l2Norm()  )
                                        ( wglobal[lelt[n]][i]->l2Norm()  ).error( "invalid basis function" );
#endif

                                }
                        }
                    // global assembly: we pass the vector of local/global basis functions
                    form2( Xh, Xh, D ) += integrate( _range=idedelements(mesh,lelt[0]), _quad=im,
                                                     _expr=delta*trans(basist(wglobal))*basis(wglobal) );

                    if ( !it->isOnBoundary() )
                        form2( Xh, Xh, D ) += integrate( _range=idedelements(mesh,lelt[1]), _quad=im,
                                                         _expr=delta*trans(basist(wglobal))*basis(wglobal) );
                    face_done[ face_id ] = true;
                    Log() << "\t Done with face " << face_id << "\n";
                }
            //Mdelta->close();
            if ( this->vm().count( "export-matlab" ) )
                Mdelta->printMatlab( "Mdelta.m" );
        }
    D->close();
//if ( delta > 0 )
//D->addMatrix( 1.0, Mdelta );
    timers["assembly_D_close"].second += timers["assembly"].first.elapsed();

    timers["assembly"].second += ( timers["assembly_D_elements"].second +
                                   timers["assembly_D_internalfaces"].second +
                                   timers["assembly_D_boundaryfaces"].second +
                                   timers["assembly_D_close"].second );
    if ( this->vm().count( "export-matlab" ) )
        D->printMatlab( "D" );

    this->solve( D, u, F );

    element_type_cont uEx( Xch, "uEx" );

    sparse_matrix_ptrtype M( M_backend->newMatrix( Xch, Xch ) );
    form2( Xch, Xch, M, _init=true ) = integrate( elements( mesh ), trans(idt(uEx))*id(uEx) );
    M->close();
    vector_ptrtype L( M_backend->newVector( Xch ) );
    form1( Xch, L ) = integrate( elements( mesh ), im_norm, trans(g)*id(uEx) );
    this->solve( M, uEx, L );

    //    Log() << "||error||_0 = " << math::sqrt(integrate( elements(mesh), im, val( trans(idv(u)-g)*(idv(u)-g) ) ).evaluate()(0,0)) << "\n";
    //    Log() << "||error||_0 = " << math::sqrt(integrate( elements(mesh), im, trans(idv(u)-idv(uEx))*(idv(u)-idv(uEx)) ).evaluate()(0,0)) << "\n";

	// Norm computations
	value_type n1 = integrate( elements(mesh), im_norm, val( (gradv(u)-trans(grad_g))*(trans(gradv(u))-grad_g) ) ).evaluate()(0,0);
    //value_type n2 = delta*integrate( elements(mesh), im_norm, val( trans(idv(tau))*idv(tau) ) ).evaluate()(0,0);

    value_type n2 = 0;
    if ( delta > 0 )
        {
            vector_ptrtype U( M_backend->newVector( Xh ) );
            *U = u;
            //n2 = Mdelta->energy( U, U );
        }

	value_type n3 = gamma*Order*Order*integrate( internalfaces(mesh), im_norm,  (trans(jumpv(idv(u)))*jumpv(idv(u)))/hFace()  ).evaluate()(0,0);
	value_type n4 = gamma*Order*Order*integrate( boundaryfaces(mesh), im_norm,  (trans(idv(u)-g)*(idv(u)-g))/hFace()  ).evaluate()(0,0);
	value_type n5 = integrate( elements(mesh), im_norm, val( (idv(u)-g)^2 ) ).evaluate()(0,0);


    // 	n1 = integrate( elements(mesh), im, val( (gradv(u)-trans(grad_g))*(trans(gradv(u))-grad_g) ) ).evaluate()(0,0);
    // 	n2 = delta*integrate( elements(mesh), im, val( trans(idv(tau))*idv(tau) ) ).evaluate()(0,0);
    // 	n3 = gamma*/hsize*integrate( internalfaces(mesh), im, (trans(jumpv(idv(u)))*jumpv(idv(u))) ).evaluate()(0,0);
    // 	n4 = gamma*Order/hsize*integrate( boundaryfaces(mesh), im, (trans(idv(u))*idv(u)) ).evaluate()(0,0);

    //	Log() << "||dg-error||_dg = " << math::sqrt(n1+n2+n3+n4) << "\n";
	Log() << "||dg-error||_dg = " << math::sqrt(n1+n3+n4) << "\n";
    Log() << "||u-error||_0 = " << math::sqrt(n5) << "\n";
	Log() << "||u-error||_1 = " << math::sqrt(n1) << "\n";
    Log() << "||L-error||_0 = " << math::sqrt(n2) << "\n";
	Log() << "||jump-error||_0 = " << math::sqrt(n3+n4) << "\n";

    //this->writeResults(math::sqrt(n1+n2+n3+n4), math::sqrt(n5), math::sqrt(n1), -1, math::sqrt(n3+n4), theta, delta, gamma );
	this->writeResults(math::sqrt(n1+n2+n3+n4), math::sqrt(n5), math::sqrt(n1), (n2>0)?math::sqrt(n2):-1, math::sqrt(n3+n4), theta, delta, gamma, anisomesh );

    form1( Xch, L, _init=true ) = integrate( elements( mesh ), trans(idv(u))*id(uEx) );
    element_type_cont uc( Xch, "uc" );
    this->solve( M, uc, L );

    this->exportResults( u, uc, uEx );

    Log() << "        run():numElements: " << mesh->numElements() << "\n";

    Log() << "[im]    run():  nPoints: " << im.nPoints() << "\n";
    Log() << "[im]    run():    Order: " << im.nOrder << "\n";
    Log() << "[im_norm]run(): nPoints: " << im_norm.nPoints() << "\n";
    Log() << "[im_norm]run():   Order: " << im_norm.nOrder << "\n";
    Log() << "[timer] run():     init: " << timers["init"].second << "\n";
    Log() << "[timer] run(): assembly: " << timers["assembly"].second << "\n";
    Log() << "[timer] run():     o D elements : " << timers["assembly_D_elements"].second << "\n";
    Log() << "[timer] run():     o D internalfaces : " << timers["assembly_D_internalfaces"].second << "\n";
    Log() << "[timer] run():     o D boundaryfaces : " << timers["assembly_D_boundaryfaces"].second << "\n";
    Log() << "[timer] run():     o D close : " << timers["assembly_D_close"].second << "\n";
    //Log() << "[timer] run():     o F : " << timers["assembly_F"].second << "\n";
    //Log() << "[timer] run():     o M : " << timers["assembly_M"].second << "\n";
    //Log() << "[timer] run():     o L : " << timers["assembly_L"].second << "\n";
    //Log() << "[timer] run():     o i : " << timers["assembly_evaluate"].second << "\n";
    Log() << "[timer] run():   solver: " << timers["solver"].second << "\n";
    Log() << "[timer] run():   solver: " << timers["export"].second << "\n";

} // Laplacian::run

template<int Order>
template<typename elem_type>
void
Laplacian<Order>::solve( sparse_matrix_ptrtype const& D, elem_type& u, vector_ptrtype const& F  )
{
    timers["solver"].first.restart();


    backend_ptrtype b( backend_type::build( this->vm() ) );
    vector_ptrtype U( b->newVector( u.functionSpace() ) );
    Log() << "[solve] D sizes: " << D->size1() << ", " << D->size2() << "\n";
    Log() << "[solve] u  size: " << u.size() << "\n";
    Log() << "[solve] U  size: " << U->size() << "\n";
    Log() << "[solve] F  size: " << F->size() << "\n";


    b->solve( D, D, U, F );
    u = *U;

    timers["solver"].second = timers["solver"].first.elapsed();
    Log() << "[timer] solve: " << timers["solver"].second << "\n";
} // Laplacian::solve


template<int Order>
void
Laplacian<Order>::exportResults( element_type& U, element_type_cont& V, element_type_cont& E )
{
    timers["export"].first.restart();

    exporter->step(1.)->setMesh( createMesh( this->vm()["hvisu"].template as<value_type>() ) );
    exporter->step(1.)->add( "u", U );
    exporter->step(1.)->add( "v", V );
    exporter->step(1.)->add( "e", E );
    exporter->save();



    timers["export"].second = timers["export"].first.elapsed();
    Log() << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Laplacian::export

template<int Order>
void
Laplacian<Order>::writeResults( value_type e1,
                                value_type e2,
                                value_type e3,
                                value_type e4,
                                value_type e5,
                                value_type& theta,
                                value_type& delta,
                                value_type& gamma,
								int& anisomesh
                                )
{
    timers["write"].first.restart();

//     exporter->step(1.)->setMesh( U.functionSpace()->mesh() );
//     exporter->step(1.)->add( "u", U );
//     exporter->step(1.)->add( "v", V );
//     exporter->step(1.)->add( "e", E );
//     exporter->save();


	std::ostringstream fname_errors;
	fname_errors << "errors_" << theta << "_" << delta << "_" << gamma << "_" << anisomesh << ".dat";
	std::ofstream ofs( fname_errors.str().c_str() );

	int wspace=15;
	ofs << std::setw( wspace ) << "theta" << " "
		<< std::setw( wspace ) << "delta" << " "
		<< std::setw( wspace ) << "gamma" << " "
		<< std::setw( wspace ) << "anisomesh" << " \n";
	ofs << std::setw( wspace ) << theta << " "
		<< std::setw( wspace ) << delta << " "
		<< std::setw( wspace ) << gamma << " "
		<< std::setw( wspace ) << anisomesh << " \n\n";
	ofs << std::setw( wspace ) << "||u||_DG" << " "
		<< std::setw( wspace ) << "||u||_0" << " "
		<< std::setw( wspace ) << "||u||_1" << " "
		<< std::setw( wspace ) << "||L(u)||_0" << " "
		<< std::setw( wspace ) << "||[u]||_0" << " \n";
	ofs << std::setw( wspace ) << e1 << " "
		<< std::setw( wspace ) << e2 << " "
		<< std::setw( wspace ) << e3 << " "
		<< std::setw( wspace ) << e4 << " "
		<< std::setw( wspace ) << e5 << " ";
	ofs << "\n\n";
	ofs << e1 << " "
		<< e2 << " "
		<< e3 << " "
		<< e4 << " "
		<< e5 << " ";
	ofs << "\n";
	ofs.close();

    timers["write"].second = timers["write"].first.elapsed();
    Log() << "[timer] writeResults(): " << timers["write"].second << "\n";
} // Laplacian::writeResults

} // Feel




int
main( int argc, char** argv )
{
    using namespace Feel;

    /* change parameters below */
    const int nOrder = 1;

    //typedef Continuous MyContinuity;  no continuous version
    typedef Discontinuous MyContinuity;
    typedef Feel::Laplacian<nOrder> laplacian_type;
    //typedef Feel::Laplacian<nDim, nOrder, MyContinuity, Simplex, Scalar> laplacian_type;

    /* assertions handling */
    Feel::Assert::setLog( "laplacian.assert");

    /* define and run application */
    laplacian_type laplacian( argc, argv, makeAbout(), makeOptions() );
    laplacian.run();
}




