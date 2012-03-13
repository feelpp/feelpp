/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-11-13

  Copyright (C) 2009 Université Joseph Fourier (Grenoble I)

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
   \file RBheat.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-11-13
 */
#ifndef __RBHEAT_H
#define __RBHEAT_H 1

#include <boost/timer.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelcrb/parameterspace.hpp>



namespace Feel
{

po::options_description
makeRbHeatOptions()
{
    po::options_description rbheatoptions("RbHeat options");
    rbheatoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.01 ), "mesh size")
        ("mu1", po::value<double>()->default_value( 0.2 ), "mu1")
        ("mu2", po::value<double>()->default_value( 0.2 ), "mu2")
        ("mu3", po::value<double>()->default_value( -1 ), "mu3")
        ("mu4", po::value<double>()->default_value( 0.1 ), "mu4")
        ("no-export", "don't export results")
        ;
    return rbheatoptions.add( Feel::feel_options() );
}
AboutData
makeRbHeatAbout( std::string const& str = "rbheat" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "RbHeat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//auto mesh = GeoTool::createMeshFromGeoFile<mesh_type>("/u/c/calandrj/Feel.src/examples/heat/MCS_heat.geo", "nameExport",meshSize);
/////////////////////////////////////////////////////////////////////////////////////////////////////////

gmsh_ptrtype
createGeo( double meshSize  )
{
std::ostringstream ostr;


 ostr <<"Point(1) = {-0.5,0,0,"<<0.01<<"};\n"
      <<"Point(2) = {0.5,0,0,"<<0.01<<"};\n"
      <<"Point(3) = {0.5,0.5,0,"<<0.01<<"};\n"
      <<"Point(4) = {-0.5,0.5,0,"<<0.01<<"};\n"
      <<"\n"
      <<"Line(1) = {1,2};\n"
      <<"Line(2) = {2,3};\n"
      <<"Line(3) = {3,4};\n"
      <<"Line(4) = {4,1};\n"
   
      <<"Line Loop(10) = {1,2,3,4};\n"
      <<"Plane Surface(20) = {10};\n"

      <<"Physical Line(\"left\") = {4};\n"
      <<"Physical Line(\"right\")={2};\n"
      <<"Physical Line(\"hautbas\")={1,4};\n"
      <<"Physical Surface(\"maille1\") = {20};\n"
      <<"\n";

#if 0
    ostr << "p=" << meshSize << ";\n"
	 <<"Point(1) = {-0.5,0,0,p};\n"
	 <<"Point(2) = {0,0,0,p};\n"
	 <<"Point(3) = {0.5,0,0,p};\n"
	 <<"Point(4) = {0.5,0.5,0,p};\n"
	 <<"Point(5) = {0,0.5,0,p};\n"
	 <<"Point(6) = {-0.5,0.5,0,p};\n"
	 <<"\n"
	 <<"Line(1) = {1,2};\n"
	 <<"Line(2) = {2,3};\n"
	 <<"Line(3) = {3,4};\n"
	 <<"Line(4) = {4,5};\n"
	 <<"Line(5) = {5,6};\n"
	 <<"Line(6) = {6,1};\n"
	 <<"Line(7) = {2,5};\n"

	 <<"Line Loop(1) = {1,7,5,6};\n"
	 <<"Line Loop(2) = {2,3,4,-7};\n"
	 <<"Plane Surface(3) = {1};\n"
	 <<"Plane Surface(4) = {2};\n"

	 <<"Physical Line(\"left\") = {6};\n"
	 <<"Physical Line(\"right\")={3};\n"
	 <<"Physical Line(\"hautbas\")={1,2,4,5};\n"
	 <<"Physical Line(\"interface\") = {7};\n"
	 <<"Physical Surface(\"maille1\") = {3};\n"
	 <<"Physical Surface(\"maille2\") = {4};\n"
	 <<"\n";
#endif
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( "rbheat" );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


class RbHeat
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 5;
    static const uint16_type ParameterSpaceDimension = 4;
    static const bool is_time_dependent = false;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef boost::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    /*mesh*/
    typedef Simplex<1,1> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;


    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef Eigen::VectorXd theta_vector_type;


    typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    RbHeat();

    //! constructor from command line
    RbHeat( po::variables_map const& vm );


    //! copy constructor
    //Heat1D( Heat1D const & );
    //! destructor
    ~RbHeat(){}

    //! initialisation of the model
    void init();
    //@}

    /** @name Operator overloads
     */
    //@{

    //@}

    /** @name Accessors
     */
    //@{

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const { return 3; }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const { return 2; }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const { if ( l == 0 ) return 2; return 1; }

    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace() { return Xh; }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const { return M_Dmu;}

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu, double time=0 )
        {
	  std::cout<<" computeThetaQ ======= "<<std::endl;
            M_thetaAq.resize( Qa() );
            M_thetaAq( 0 ) = 1;
            M_thetaAq( 1 ) = mu( 0 ); // k_1
            M_thetaAq( 2 ) = mu( 1 ); // k_2

            M_thetaFq.resize( Nl() );
            M_thetaFq[0].resize( Ql(0) );
            M_thetaFq[0]( 0 ) = mu(2); // delta
            M_thetaFq[0]( 1 ) = mu(3); // phi

            M_thetaFq[1].resize( Ql(1) );
            M_thetaFq[1]( 0 ) = 1;

            return boost::make_tuple( M_thetaAq, M_thetaFq );
        }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const { return M_thetaAq; }

    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const { return M_thetaFq; }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaAq( int q ) const
        {
            return M_thetaAq( q );
        }

    /**
     * \return the \p q -th term of the \p l -th output
     */
    value_type thetaL( int l, int q ) const
        {
            return M_thetaFq[l]( q );
        }

    //@}

    /** @name  Mutators
     */
    //@{

    /**
     * set the mesh characteristic length to \p s
     */
    void setMeshSize( double s ) { meshSize = s; }


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * run the convergence test
     */

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& T );

    /**
     * solve for a given parameter \p mu
     */
    void solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );


    /**
     * update the PDE system with respect to \param mu
     */
    void update( parameter_type const& mu );
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( element_type& u );

    void solve( sparse_matrix_ptrtype& ,element_type& ,vector_ptrtype&  );

    /**
     * returns the scalar product of the boost::shared_ptr vector x and
     * boost::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );

    /**
     * specific interface for OpenTURNS
     *
     * \param X input vector of size N
     * \param N size of input vector X
     * \param Y input vector of size P
     * \param P size of input vector Y
     */
    void run( const double * X, unsigned long N, double * Y, unsigned long P );

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu );


private:

    po::variables_map M_vm;
    backend_ptrtype backend;

    double meshSize;

    bool M_use_weak_dirichlet;
    double M_gammabc;

    bool M_do_export;
    export_ptrtype exporter;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;
    element_ptrtype pT;

    std::vector<sparse_matrix_ptrtype> M_Aq;
    std::vector<std::vector<vector_ptrtype> > M_Fq;

    parameterspace_ptrtype M_Dmu;
    theta_vector_type M_thetaAq;
    std::vector<theta_vector_type> M_thetaFq;
};

RbHeat::RbHeat()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_do_export( false ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}


RbHeat::RbHeat( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "rbHeat" ) ),
    M_Dmu( new parameterspace_type )
{
  this->init();
}
void
RbHeat::init()
{
    /*
     * First we create the mesh
     */
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=createGeo(meshSize) ,
			   _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    M_Aq.resize( 3 );
    M_Aq[0] = backend->newMatrix( Xh, Xh );
    M_Aq[1] = backend->newMatrix( Xh, Xh );
    M_Aq[2] = backend->newMatrix( Xh, Xh );

    M_Fq.resize( 2 );
    M_Fq[0].resize( 2 );
    M_Fq[0][0] = backend->newVector( Xh );
    M_Fq[0][1] = backend->newVector( Xh );

    M_Fq[1].resize( 1 );
    M_Fq[1][0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    using namespace Feel::vf;
    static const int N = 2;
    Feel::ParameterSpace<4>::Element mu_min( M_Dmu );
    mu_min << 0.2, 0.2, 0.01, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<4>::Element mu_max( M_Dmu );
    mu_max << 50, 50, 5, 5;
    M_Dmu->setMax( mu_max );


    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    Log() << "Number of dof " << Xh->nLocalDof() << "\n";
    std::cout<<"Number of dof " << Xh->nLocalDof() << "\n";




    double surface_left = integrate( _range= markedfaces(mesh,"left"), _expr=cst(1.) ).evaluate()(0,0);
    std::cout<<"left = "<<surface_left<<std::endl;

 double surface_right = integrate( _range= markedfaces(mesh,"right"), _expr=cst(1.) ).evaluate()(0,0);
    std::cout<<"right = "<<surface_right<<std::endl;

 double surface__haut_bas = integrate( _range= markedfaces(mesh,"hautbas"), _expr=cst(1.) ).evaluate()(0,0);
    std::cout<<"hautbas = "<<surface__haut_bas<<std::endl;

 double surface__maille1 = integrate( _range= markedelements(mesh,"maille1"), _expr=cst(1.) ).evaluate()(0,0);
    std::cout<<"maille1 = "<<surface__maille1<<std::endl;
    exit(0);

 double surface__maille2 = integrate( _range= markedelements(mesh,"maille2"), _expr=cst(1.) ).evaluate()(0,0);
    std::cout<<"maille2 = "<<surface__maille2<<std::endl;

    // right hand side
    form1( Xh, M_Fq[0][0], _init=true ) = integrate( markedfaces(mesh,"left"), id(v) );

    form1( _test=Xh, _vector=M_Fq[0][1], _init=true ) = integrate( elements(mesh), id(v) );
    M_Fq[0][0]->close();
    M_Fq[0][1]->close();

    // output
    form1( Xh, M_Fq[1][0], _init=true ) = integrate( markedelements(mesh,"maille1"), id(v)/0.2 );
    form1( Xh, M_Fq[1][0] ) += integrate( markedelements(mesh,"maille2"), id(v)/0.2 );
    M_Fq[1][0]->close();

    form2( Xh, Xh, M_Aq[0], _init=true ) = integrate( elements(mesh), 0.1*(gradt(u)*trans(grad(v)) ) );
    form2( Xh, Xh, M_Aq[0] ) += integrate( markedfaces(mesh,"right"), id(u)*idt(v) );
    M_Aq[0]->close();

    form2( Xh, Xh, M_Aq[1], _init=true ) = integrate( markedelements(mesh,"maille1"), (gradt(u)*trans(grad(v)) ) );
    M_Aq[1]->close();

    form2( Xh, Xh, M_Aq[2], _init=true ) = integrate( markedelements(mesh,"maille2"), (gradt(u)*trans(grad(v)) ) );
    M_Aq[2]->close();
    M = backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, M, _init=true ) =
        integrate( elements(mesh), id(u)*idt(v) + grad(u)*trans(gradt(u)) );
    M->close();

////////////////////////////////////////////////////////////////////////////////


} // RbHeat::init

RbHeat::sparse_matrix_ptrtype
RbHeat::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

RbHeat::affine_decomposition_type
RbHeat::computeAffineDecomposition()
{
  std::cout<<"computeAffineDecomposition du model"<<std::endl;
    return boost::make_tuple( M_Aq, M_Fq );
}


void
RbHeat::solve( sparse_matrix_ptrtype& D,
               element_type& t,
               vector_ptrtype& F )
{

    vector_ptrtype T( backend->newVector( t.functionSpace() ) );
    backend->solve( D, D, T, F );
    t = *T;
} // Heat1d::solve


void
RbHeat::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        Log() << "exportResults starts\n";

        exporter->step(0)->setMesh( U.functionSpace()->mesh() );

        exporter->step(0)->add( "u", U );

        exporter->save();
    }
} // Heat1d::export

void
RbHeat::update( parameter_type const& mu )
{

  std::cout<<"================================================================= "<<std::endl;

    *D = *M_Aq[0];
    //std::cout<<"*M_Aq[0] = "<<std::endl;
    //std::cout<<*M_Aq[0]<<std::endl;
    for( size_type q = 1;q < M_Aq.size(); ++q )
    {
      //std::cout<<"\n\n"<<std::endl;
      //std::cout<<"*M_Aq["<<q<<"] = "<<std::endl;
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }

    F->close();
    F->zero();
    for( size_type q = 0;q < M_Fq[0].size(); ++q )
    {
      std::cout<<"M_Fq[0]["<<q<<"] = "<<std::endl;
      std::cout<<*M_Fq[0][q]<<std::endl;
        F->add( M_thetaFq[0][q], M_Fq[0][q] );
    }
}
void
RbHeat::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    //this->exportResults( *T );

}

void
RbHeat::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeThetaq( mu );
    this->update( mu );
    backend->solve( _matrix=D,  _solution=T, _rhs=F );
}

void
RbHeat::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f);
    //std::cout << "l2solve(u,f) done\n";
}

double
RbHeat::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
RbHeat::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
RbHeat::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
  std::cout<<" model RUN"<<std::endl;
    using namespace vf;
    Feel::ParameterSpace<4>::Element mu( M_Dmu );
    std::cout<<"model RUN -1- "<<std::endl;
    mu << X[0], X[1], X[2], X[3];
    std::cout<<"model RUN -2- "<<std::endl;
    static int do_init = true;
    if ( do_init )
    {
    std::cout<<"model RUN -3- "<<std::endl;
        meshSize = X[4];
        this->init();
	std::cout<<"model RUN -4- "<<std::endl;
        do_init = false;
    }
    std::cout<<"model RUN call solve"<<std::endl;
    this->solve( mu, pT );
    std::cout<<"model RUN END "<<std::endl;
}


double
RbHeat::output( int output_index, parameter_type const& mu )
{

    std::cout<<"model output"<<std::endl;

    using namespace vf;
    this->solve( mu, pT );
    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    // right hand side (compliant)
    if( output_index == 0 )
    {
        double s1 = M_thetaFq[0](0)*dot( M_Fq[0][0], U )+M_thetaFq[0](1)*dot( M_Fq[0][1], U );

        return s1;
    }
    // output
    if ( output_index == 1 )
    {
	    double s1 = M_thetaFq[1](0)*dot( M_Fq[1][0], U )+M_thetaFq[1](1)*dot( M_Fq[1][1], U );
        double mean = integrate( elements(mesh),
                                 chi( (Px() >= -0.1) && (Px() <= 0.1) )*idv(*pT) ).evaluate()(0,0)/0.2;
        //std::cout<<"output1 c1 = "<<mean<<std::endl;

        double meanT = ( integrate( markedelements(mesh,"maille1"),idv(*pT) ).evaluate()(0,0)+
                         integrate( markedelements(mesh,"maille2"),idv(*pT) ).evaluate()(0,0) )/0.2;
        //std::cout<<"output1 c2 = "<<meanT<<std::endl;
        //std::cout<<"output1 c3= "<< dot( M_Fq[1][0], U ) <<std::endl;
	std ::cout <<"sortie"<<endl;
        return s1;

    }

}

}

#endif /* __RBHEAT_H */



