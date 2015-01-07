/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-02-05

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file dar.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-02-05
 */
#include <feel/feel.hpp>

namespace Feel
{
/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description DARoptions( "DAR options" );
    DARoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ( "bx", Feel::po::value<double>()->default_value( 1.0 ), "convection X component" )
    ( "by", Feel::po::value<double>()->default_value( 0.0 ), "convection Y component" )
    ( "bz", Feel::po::value<double>()->default_value( 0.0 ), "convection Z component" )
    ( "epsilon", po::value<double>()->default_value( 1 ), "diffusion coefficient" )
    ( "mu", po::value<double>()->default_value( 1 ), "reaction coefficient" )
    ( "weakdir", po::value<int>()->default_value( 1 ), "use weak Dirichlet condition" )
    ( "penaldir", Feel::po::value<double>()->default_value( 10 ),
      "penalisation parameter for the weak boundary Dirichlet formulation" )
    ( "stab", po::value<int>()->default_value( 1 ), "use stabilisation=1, no stabilisation=0" )
    ( "stabcoeff", Feel::po::value<double>()->default_value( 2.5e-2 ),
      "stabilisation coefficient" )
    ;
    return DARoptions;
}

/**
 * \class DAR
 *
 * DAR Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
template<int Dim>
class DAR
    :
public Simget
{
    typedef Simget super;
public:

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 2;

    //! numerical type is double
    typedef double value_type;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<Dim> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    DAR()
        :
        super(),
        meshSize( doption("hsize") ),
        shape(    soption("shape") )
    {
    }

    void run();

private:

    //! mesh characteristic size
    double meshSize;

    //! shape of the domain
    std::string shape;
}; // DAR

template<int Dim> const uint16_type DAR<Dim>::Order;

template<int Dim>
void
DAR<Dim>::run()
{
    LOG(INFO) << "------------------------------------------------------------\n";
    LOG(INFO) << "Execute DAR<" << Dim << ">\n";

    Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % Order
                                   % meshSize );

    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                      _usenames=true,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=meshSize ) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    space_ptrtype Xh = space_type::New( mesh );
    element_type u( Xh, "u" );
    element_type v( Xh, "v" );
    element_type gproj( Xh, "v" );
    /** \endcode */

<<<<<<< Updated upstream
    bool weak_dirichlet =  ioption("weakdir"  );
    value_type penaldir =  doption("penaldir" );
    bool stab =            ioption("stab"     );
    value_type stabcoeff = doption("stabcoeff");
    value_type epsilon =   doption("epsilon"  );
    value_type mu =        doption("mu"       );
    value_type bx =        doption("bx"       );
    value_type by =        doption("by"       );
=======
    bool weak_dirichlet = ioption("weakdir");
    value_type penaldir = doption("penaldir");
    bool stab = ioption("stab");
    value_type stabcoeff = doption("stabcoeff");
    value_type epsilon = doption("epsilon");
    value_type mu = doption("mu");
    value_type bx = doption("bx");
    value_type by = doption("by");
>>>>>>> Stashed changes

    LOG(INFO) << "[DAR] hsize = " << meshSize << "\n";
    LOG(INFO) << "[DAR] bx = " << bx << "\n";
    LOG(INFO) << "[DAR] by = " << by << "\n";
    LOG(INFO) << "[DAR] mu = " << mu << "\n";
    LOG(INFO) << "[DAR] epsilon = " << epsilon << "\n";
    LOG(INFO) << "[DAR] bccoeff = " << penaldir << "\n";
    LOG(INFO) << "[DAR] bctype = " << weak_dirichlet << "\n";
    LOG(INFO) << "[DAR] stab = " << stab << "\n";
    LOG(INFO) << "[DAR] stabcoeff = " << stabcoeff << "\n";



    /** define \f$g\f$ the expression of the exact solution and
     * \f$f\f$ the expression of the right hand side such that \f$g\f$
     * is the exact solution
     */
    /** \code */
    //# marker1 #
    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto g = sin( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() );
    gproj = vf::project( Xh, elements( mesh ), g );
    auto grad_g = vec( +pi*cos( pi*Px() )*cos( pi*Py() )*cos( pi*Pz() ),
                       -pi*sin( pi*Px() )*sin( pi*Py() )*cos( pi*Pz() ) );
    auto beta = vec( cst( bx ),cst( by ) );
    //! deduce from expression the type of f (thanks to keyword 'auto')
    auto f = ( pi*pi*Dim*epsilon*g+trans( grad_g )*beta + mu*g );
    //# endmarker1 #
    /** \endcode */





    /**
     * Construction of the right hand side. F is the vector that holds
     * the algebraic representation of the right habd side of the
     * problem
     */
    /** \code */
    //# marker2 #
    auto F = backend( _vm=this->vm() )->newVector( Xh );
    form1( _test=Xh, _vector=F, _init=true ) =
        integrate( _range=elements( mesh ), _expr=f*id( v ) );

    //# endmarker2 #
    if ( weak_dirichlet )
    {
        //# marker41 #
#if 1
        form1( _test=Xh, _vector=F ) +=
            integrate( _range=boundaryfaces( mesh ),
                       _expr=g*( -epsilon*grad( v )*vf::N()+penaldir*id( v )/hFace() ) );
#endif
        //# endmarker41 #
    }

    /** \endcode */

    /**
     * create the matrix that will hold the algebraic representation
     * of the left hand side
     */
    //# marker3 #
    /** \code */
    size_type pattern=Pattern::COUPLED;

    if ( stab )
        pattern = Pattern::COUPLED|Pattern::EXTENDED;

    auto D = backend()->newMatrix( _test=Xh, _trial=Xh, _pattern=pattern  );
    /** \endcode */

    //! assemble $\int_\Omega \epsilon \nabla u \cdot \nabla v$
    /** \code */
    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements( mesh ),
                   _expr=( epsilon*gradt( u )*trans( grad( v ) )+( gradt( u )*beta )*id( v )+mu*idt( u )*id( v ) ) );
    /** \endcode */
    //# endmarker3 #

    if ( stab )
    {
        // define the stabilisation coefficient expression
        auto stab_coeff = ( stabcoeff*abs( bx*Nx()+by*Ny() )* vf::pow( hFace(),2.0 ) );
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=internalfaces( mesh ),
                       _expr=stab_coeff*( trans( jumpt( gradt( u ) ) )*jump( grad( v ) ) ) );
    }

    if ( weak_dirichlet )
    {
        /** weak dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# assemble \f$\int_{\partial \Omega} -\nabla u \cdot \mathbf{n} v\f$
         * -# assemble \f$\int_{\partial \Omega} -\nabla v \cdot \mathbf{n} u\f$
         * -# assemble \f$\int_{\partial \Omega} \frac{\gamma}{h} u v\f$
         */
        /** \code */
        //# marker10 #
#if 1
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=boundaryfaces( mesh ),
                       _expr= ( -( epsilon*gradt( u )*vf::N() )*id( v )
                                -( epsilon*grad( v )*vf::N() )*idt( u )
                                +penaldir*id( v )*idt( u )/hFace() ) );
#else
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            integrate( _range=boundaryfaces( mesh ),
                       _expr= ( -( epsilon*gradt( u )*vf::N() )*id( v ) ) );


#endif
        //# endmarker10 #
        /** \endcode */
    }

    else
    {
        /** strong(algebraic) dirichlet conditions treatment for the boundaries marked 1 and 3
         * -# first close the matrix (the matrix must be closed first before any manipulation )
         * -# modify the matrix by cancelling out the rows and columns of D that are associated with the Dirichlet dof
         */
        /** \code */
        //# marker5 #
#if 1
        form2( _test=Xh, _trial=Xh, _matrix=D ) +=
            on( _range=boundaryfaces( mesh ),_element=u, _rhs=F, _expr=g );
#endif
        //# endmarker5 #
        /** \endcode */

    }

    /** \endcode */

    //! solve the system
    /** \code */
    //# marker6 #
    backend( _rebuild=true,_vm=this->vm() )->solve( _matrix=D, _solution=u, _rhs=F );
    //# endmarker6 #
    /** \endcode */

    //! compute the \f$L_2$ norm of the error
    /** \code */
    //# marker7 #

    LOG(INFO) << "||error||_L2=" << normL2( _range=elements(mesh), _expr=idv(u)-g ) << "\n";

    //# endmarker7 #
    /** \endcode */

    //! save the results
    /** \code */
    //! project the exact solution
    element_type e( Xh, "e" );
    e = vf::project( Xh, elements( mesh ), g );

    export_ptrtype exporter( export_type::New( this->vm(),
                             ( boost::format( "%1%-%2%-%3%" )
                               % this->about().appName()
                               % shape
                               % Dim ).str() ) );

    if ( exporter->doExport() )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( mesh );

        exporter->step( 0 )->add( "u", u );
        exporter->step( 0 )->add( "g", e );

        exporter->save();
        LOG(INFO) << "exportResults done\n";
    }

    /** \endcode */
} // DAR::run

} // Feel

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;
    /**
     * Initialize Feel++ Environment
     */
    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="dar",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    /**
     * create an application
     */
    /** \code */
    Application app;

    /** \endcode */

    /**
     * register the simgets
     */
    /** \code */
    app.add( new DAR<2>() );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}






