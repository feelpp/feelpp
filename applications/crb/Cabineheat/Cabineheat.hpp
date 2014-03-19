/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \file Cabineheat.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __CABINEHEAT_H
#define __CABINEHEAT_H 1

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
makeCabineHeatOptions()
{
    po::options_description rbheatoptions( "CabineHeat options" );
    rbheatoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.01 ), "mesh size" )
    ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    ( "mu2", po::value<double>()->default_value( -1 ), "mu2" )
    ( "mu3", po::value<double>()->default_value( 0.1 ), "mu3" )
    ( "no-export", "don't export results" )
    ;
    return rbheatoptions;
}
AboutData
makeCabineHeatAbout( std::string const& str = "CabineHeat" )
{
    Feel::AboutData about( str.c_str(),
                           str.c_str(),
                           "0.1",
                           "CabineHeat Benchmark",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2010,2011 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//auto mesh = GeoTool::createMeshFromGeoFile<mesh_type>("/u/c/calandrj/Feel.src/examples/heat/MCS_heat.geo", "nameExport",meshSize);
/////////////////////////////////////////////////////////////////////////////////////////////////////////

gmsh_ptrtype
createGeo( double meshSize  )
{
    std::ostringstream ostr;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////// MAILLAGE CABINE TOP ////////////////////////////////////////////////////////

    ostr << "lc1=" <<meshSize <<";\n"

         <<"Point(1) = {-1.9, -0.05, 0, lc1};\n"
         <<"Point(2) = {1.9, -0.05,  0, lc1};\n"
         <<"Point(3) = {1.9, 0.05, 0, lc1} ;\n"
         <<"Point(4) = {1.9993749, 0.05,0, lc1} ;\n"
         <<"Point(5) = {1.812616, 0.845236,  0, lc1} ;\n"
         <<"Point(6) = {1.3, 1.519868, 0, lc1} ;\n"
         <<"Point(7) = {0.46915, 2, 0, lc1} ;\n"	      // point on major axis of right side ellipse
         <<"Point(8) = {0.125, 1.996089, 0, lc1};\n"
         <<"Point(9) = {0.05, 1.9993749, 0, lc1};\n"
         <<"Point(10) = {0, 2, 0, lc1};\n"
         <<"Point(11) = {-1.9, 0.05, 0, lc1};\n"
         <<"Point(12) = {-1.9993749,0.05,0, lc1} ;\n"
         <<"Point(13) = {-1.812616, 0.845236,  0, lc1} ;\n"
         <<"Point(14) = {-1.3, 1.519868, 0, lc1} ;\n"
         <<"Point(15) = {-0.46915, 2, 0, lc1} ;\n"	// point on major axis of left side ellipse
         <<"Point(16) = {-0.125, 1.996089, 0, lc1} ;\n"
         <<"Point(17) = {-0.05, 1.9993749, 0, lc1} ;\n"
         <<"Point(18) = {-1.3, 2, 0, lc1};\n"		// left side ellipse center
         <<"Point(19) = {1.3, 2, 0, lc1};	\n"	// right side ellipse center
         <<"Point(20) = {0, 0, 0, lc1};	\n"	// center of circle
         <<"Point(21) = {0.05, 1.875, 0, lc1};\n"
         <<"Point(22) = {-0.05, 1.875, 0, lc1};\n"

         <<"Point (23) = {-1.55, 0.50, 0, lc1};\n"	//center of P1
         <<"Point (24) = {-1.55, 0.80, 0, lc1};\n"	// point on major axis of P1
         <<"Point (25) = {-1.55, 1.00, 0, lc1};\n"
         <<"Point (26) = {-1.55, 0.00, 0, lc1};\n"
         <<"Point (27) = {-1.35, 0.50, 0, lc1};\n"
         <<"Point (28) = {-1.75, 0.50, 0, lc1};\n"

         <<"Point (29) = {-1.05, 0.50, 0, lc1};\n"	//center of P2
         <<"Point (30) = {-1.05, 0.80, 0, lc1};\n"	// point on major axis of P2
         <<"Point (31) = {-1.05, 1.00, 0, lc1};\n"
         <<"Point (32) = {-1.05, 0.00, 0, lc1};\n"
         <<"Point (33) = {-0.85, 0.50, 0, lc1};\n"
         <<"Point (34) = {-1.25, 0.50, 0, lc1};\n"

         <<"Point (35) = {-0.55, 0.50, 0, lc1};\n"	//center of P3
         <<"Point (36) = {-0.55, 0.80, 0, lc1};\n"	// point on major axis of P3
         <<"Point (37) = {-0.55, 1.00, 0, lc1};\n"
         <<"Point (38) = {-0.55, 0.00, 0, lc1};\n"
         <<"Point (39) = {-0.35, 0.50, 0, lc1};\n"
         <<"Point (40) = {-0.75, 0.50, 0, lc1};\n"

         <<"Point (41) = { 0.55, 0.50, 0, lc1};\n"	//center of P4
         <<"Point (42) = { 0.55, 0.80, 0, lc1};\n"	// point on major axis of P4
         <<"Point (43) = { 0.55, 1.00, 0, lc1};\n"
         <<"Point (44) = { 0.55, 0.00, 0, lc1};\n"
         <<"Point (45) = { 0.75, 0.50, 0, lc1};\n"
         <<"Point (46) = { 0.35, 0.50, 0, lc1};\n"

         <<"Point (47) = { 1.05, 0.50, 0, lc1};\n"	//center of P5
         <<"Point (48) = { 1.05, 0.80, 0, lc1};\n"	// point on major axis of P5
         <<"Point (49) = { 1.05, 1.00, 0, lc1};\n"
         <<"Point (50) = { 1.05, 0.00, 0, lc1};\n"
         <<"Point (51) = { 1.25, 0.50, 0, lc1};\n"
         <<"Point (52) = { 0.85, 0.50, 0, lc1};\n"

         <<"Point (53) = { 1.55, 0.50, 0, lc1};\n"	//center of P6
         <<"Point (54) = { 1.55, 0.80, 0, lc1};\n"// point on major axis of P6
         <<"Point (55) = { 1.55, 1.00, 0, lc1};\n"
         <<"Point (56) = { 1.55, 0.00, 0, lc1};\n"
         <<"Point (57) = { 1.75, 0.50, 0, lc1};\n"
         <<"Point (58) = { 1.35, 0.50, 0, lc1};\n"

         <<"Point(59)={-1.9993749,-0.05, 0, lc1};\n"		//Point(59)={-1.9993749,-0.05,0, lc1};
         <<"Point(60)={-1.130265, -1.65, 0, lc1};\n" //Point(60)={-1.130265, -1.65, 0, lc2};
         <<"Point(61)={ 1.130265, -1.65, 0, lc1};\n" //Point(61)={ 1.130265, -1.65, 0, lc2};
         <<"Point(62)={ 1.9993749,-0.05, 0, lc1};\n"		//Point(62)={ 1.9993749,-0.05,0, lc1};

         <<"Point(63)={-0.15,-1.65,0, lc1}; \n"         // bay inlet
         <<"Point(64)={ 0.20,-1.65,0, lc1};\n"
         <<"Point(65)={ 0.20,-1.55,0, lc1};\n"
         <<"Point(66)={-0.15,-1.55,0, lc1};\n"

         <<"Line(1) = {1,2} ;\n"
         <<"Line(2) = {2,3} ;\n"
         <<"Line(3) = {3,4} ;\n"
         <<"Line(4) = {11,1} ;\n"
         <<"Line(5) = {12,11} ;\n"
         <<"Circle (6) = {4, 20, 6};\n"
         <<"Ellipse (7) = {6, 19, 7, 8};\n"	     // Right HL ellipse
         <<"Circle (8) = {8, 20, 9};\n"
         <<"Circle (9) = {14, 20, 12};\n"
         <<"Ellipse (10) = {16, 18, 15, 14};\n"	 // Left HL ellipse
         <<"Circle (11) = {17, 20, 16};\n"
         <<"Line(12) = {9,21} ;\n"
         <<"Line(13) = {21,22} ;\n"
         <<"Line(14) = {22,17} ;\n"
         <<"Circle (15) = {6, 20, 8};\n"
         <<"Circle (16) = {9, 20, 17};\n"
         <<"Circle (17) = {16, 20, 14};\n"
         //creation passagers
         <<"Ellipse (18) = {26, 23, 24, 27};\n"     // P1
         <<"Ellipse (19) = {27, 23, 24, 25};\n"
         <<"Ellipse (20) = {25, 23, 24, 28};\n"
         <<"Ellipse (21) = {28, 23, 24, 26};\n"
         //
         <<"Ellipse (22) = {32, 29, 30, 33}; \n"    // p2
         <<"Ellipse (23) = {33, 29, 30, 31};\n"
         <<"Ellipse (24) = {31, 29, 30, 34};\n"
         <<"Ellipse (25) = {34, 29, 30, 32};\n"
         //
         <<"Ellipse (26) = {38, 35, 36, 39}; \n"    // p3
         <<"Ellipse (27) = {39, 35, 36, 37};\n"
         <<"Ellipse (28) = {37, 35, 36, 40};\n"
         <<"Ellipse (29) = {40, 35, 36, 38};\n"
         //
         <<"Ellipse (30) = {44, 41, 42, 45};\n"     // p4
         <<"Ellipse (31) = {45, 41, 42, 43};\n"
         <<"Ellipse (32) = {43, 41, 42, 46};\n"
         <<"Ellipse (33) = {46, 41, 42, 44};\n"
         //
         <<"Ellipse (34) = {50, 47, 48, 51};  \n"   // p5
         <<"Ellipse (35) = {51, 47, 48, 49};\n"
         <<"Ellipse (36) = {49, 47, 48, 52};\n"
         <<"Ellipse (37) = {52, 47, 48, 50};\n"
         //
         <<"Ellipse (38) = {56, 53, 54, 57}; \n"     // p6
         <<"Ellipse (39) = {57, 53, 54, 55};\n"
         <<"Ellipse (40) = {55, 53, 54, 58};\n"
         <<"Ellipse (41) = {58, 53, 54, 56};\n"

         //creation contour bay+entrees
         <<"Line(42)={59,1};\n"
         <<"Line(43)={2,62};\n"

         <<"Circle(44)={12,20,59};\n"
         <<"Circle(45)={59,20,60};\n"	          //Circle(45)={59,20,60};

         <<"Line(46)={60,63};\n"
         <<"Line(47)={63,64};\n"
         <<"Line(48)={64,61};\n"

         <<"Circle(49)={61,20,62};  \n" 	          //Circle(49)={61,20,62};
         <<"Circle(50)={62,20,4}; \n"

         <<"Line(51)={64,65}; \n"                    // bay inlet
         <<"Line(52)={65,66};\n"
         <<"Line(53)={66,63};\n"

         <<"Line Loop(1) = {1,2,3,6,7,8,12,13,14,11,10,9,5,4} ; \n"    // Cabin Fluid
         <<"Line Loop(2) = {18,19,20,21} ;    \n"                      // P1
         <<"Line Loop(3) = {22,23,24,25} ; \n"                         // P2
         <<"Line Loop(4) = {26,27,28,29}; \n"                          // P3
         <<"Line Loop(5) = {30,31,32,33}; \n"                          // P4
         <<"Line Loop(6) = {34,35,36,37}; \n"                          // P5
         <<"Line Loop(7) = {38,39,40,41}; \n"                          // P6

         <<"Plane Surface(8) = {1,2,3,4,5,6,7};\n"		                // cabin fluid

         <<"Line Loop(15) = {45,46,-53,-52,-51,48,49,-43,-1,-42};\n"   // eBay

         <<"Physical Surface(\"cabin\") = {8};\n"
         <<"Physical Surface(\"bay_Surface\") = {26};\n"

         <<"Physical Line(\"passenger1\") = {20, 19, 21, 18};\n"
         <<"Physical Line(\"passenger2\") = {24, 23, 25, 22};\n"
         <<"Physical Line(\"passenger3\") = {28, 27, 29, 26};\n"
         <<"Physical Line(\"passenger4\") = {32, 31, 33, 30};\n"
         <<"Physical Line(\"passenger5\") = {36, 35, 37, 34};\n"
         <<"Physical Line(\"passenger6\") = {40, 39, 41, 38};\n"
         <<"Physical Line(\"cabin_wall\") = {9, 6, 10, 7, 8, 11};\n"
         <<"Physical Line(\"plancher\") = {1};\n"
         <<"Physical Line(\"bay_wall\") = {45, 49, 48, 46};\n"
         <<"Physical Line(\"bay_inlet\") = {53, 52, 51};\n"
         <<"Physical Line(\"left_outlet\") = {5, 4};\n"
         <<"Physical Line(\"right_outlet\") = {3, 2};\n"
         <<"Physical Line(\"cabine_inlet\") = {14, 13, 12};\n"
         <<"\n";

    ///////////////////FIN MAILLAGE CABINE /////////////////////////////////////////////
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( "Cabineheat" );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


class CabineHeat
{
public:


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 5;
    static const uint16_type ParameterSpaceDimension = 5;
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
    typedef Simplex<2,1> entity_type;
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
    CabineHeat();

    //! constructor from command line
    CabineHeat( po::variables_map const& vm );


    //! copy constructor
    //Heat1D( Heat1D const & );
    //! destructor
    ~CabineHeat() {}

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
    int Qa() const
    {
        return 1;
    }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     */
    int Nl() const
    {
        return 2;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     */
    int Ql( int l ) const
    {
        if ( l == 0 ) return 4;

        return 1;
    }

    /**
     * \brief Returns the function space
     */
    space_ptrtype functionSpace()
    {
        return Xh;
    }

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<theta_vector_type, std::vector<theta_vector_type> >
    computeThetaq( parameter_type const& mu, double time=0 )
    {
        M_thetaAq.resize( Qa() );
        M_thetaAq( 0 ) = mu( 0 ); //k
        M_thetaFq.resize( Nl() );
        M_thetaFq[0].resize( Ql( 0 ) );
        M_thetaFq[0]( 0 ) = mu( 1 ); // delta 1Flux_passenger(0->175)
        M_thetaFq[0]( 1 ) = mu( 2 ); // delta2
        M_thetaFq[0]( 2 ) = mu( 3 ); // delta3
        M_thetaFq[0]( 3 ) = mu( 4 ); // phi
        M_thetaFq[1].resize( Ql( 1 ) );
        M_thetaFq[1]( 0 ) = 1; //output 1
        return boost::make_tuple( M_thetaAq, M_thetaFq );
    }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaAq() const
    {
        return M_thetaAq;
    }

    /**
     * \brief return the coefficient vector
     */
    std::vector<theta_vector_type> const& thetaFq() const
    {
        return M_thetaFq;
    }

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
    void setMeshSize( double s )
    {
        meshSize = s;
    }


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

CabineHeat::CabineHeat()
    :
    backend( backend_type::build( BACKEND_PETSC ) ),
    meshSize( 0.01 ),
    M_do_export( false ),
    exporter( Exporter<mesh_type>::New( "ensight" ) ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}


CabineHeat::CabineHeat( po::variables_map const& vm )
    :
    M_vm( vm ),
    backend( backend_type::build( vm ) ),
    meshSize( vm["hsize"].as<double>() ),
    M_do_export( !vm.count( "no-export" ) ),
    exporter( Exporter<mesh_type>::New( vm, "cabineheat" ) ),
    M_Dmu( new parameterspace_type )
{
    this->init();
}
void
CabineHeat::init()
{
    /*
     * First we create the mesh
     */
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _desc=createGeo( meshSize ) ,
                           _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );

    /*
     * The function space and some associate elements are then defined
     */
    Xh = space_type::New( mesh );
    // allocate an element of Xh
    pT = element_ptrtype( new element_type( Xh ) );

    //  initialisation de A1 et A2
    M_Aq.resize( 1 );
    M_Aq[0] = backend->newMatrix( Xh, Xh );

    M_Fq.resize( 4 );
    M_Fq[0].resize( 4 );
    M_Fq[0][0] = backend->newVector( Xh );
    M_Fq[0][1] = backend->newVector( Xh );
    M_Fq[0][2] = backend->newVector( Xh );
    M_Fq[0][3] = backend->newVector( Xh );

    M_Fq[1].resize( 1 );
    M_Fq[1][0] = backend->newVector( Xh );

    D = backend->newMatrix( Xh, Xh );
    F = backend->newVector( Xh );

    using namespace Feel::vf;
    static const int N = 2;
    Feel::ParameterSpace<5>::Element mu_min( M_Dmu );
    mu_min << 0.2,100, 0.1, 0.1, 0.1;
    M_Dmu->setMin( mu_min );
    Feel::ParameterSpace<5>::Element mu_max( M_Dmu );
    mu_max << 50, 175, 5, 5, 5;
    M_Dmu->setMax( mu_max );


    element_type u( Xh, "u" );
    element_type v( Xh, "v" );

    LOG(INFO) << "Number of dof " << Xh->nLocalDof() << "\n";
    //std::cout<<"Number of dof " << Xh->nLocalDof() << "\n";


    // right hand side
    form1( Xh, M_Fq[0][0], _init=true ) = integrate( markedfaces( mesh,"passenger1" ), id( v ) );
    form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mesh,"passenger2" ), id( v ) );
    form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mesh,"passenger3" ), id( v ) );
    form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mesh,"passenger4" ), id( v ) );
    form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mesh,"passenger5" ), id( v ) );
    form1( Xh, M_Fq[0][0] ) += integrate( markedfaces( mesh,"passenger6" ), id( v ) );

    form1( Xh, M_Fq[0][1], _init=true ) = integrate( markedfaces( mesh,"cabin_wall" ), id( v ) );

    form1( Xh, M_Fq[0][2], _init=true ) = integrate( markedfaces( mesh,"right_outlet" ), id( v ) );
    form1( Xh, M_Fq[0][2] ) += integrate( markedfaces( mesh,"left_outlet" ), id( v ) );

    form1( Xh, M_Fq[0][3], _init=true ) = integrate( elements( mesh ), id( v ) );
    M_Fq[0][0]->close();
    M_Fq[0][1]->close();
    M_Fq[0][2]->close();
    M_Fq[0][3]->close();

    // output non compliant : mean temperature on plancher
    double onplancher = integrate( markedfaces( mesh,"plancher" ), cst( 1. ) ).evaluate()( 0,0 );
    std::cout<<"onplancher"<<onplancher<<std::endl;
    form1( Xh, M_Fq[1][0], _init=true ) = integrate( markedfaces( mesh,"plancher" ), id( v )/onplancher );
    M_Fq[1][0]->close();

    form2( Xh, Xh, M_Aq[0], _init=true ) = integrate( elements( mesh ),( gradt( u )*trans( grad( v ) ) ) );
    form2( Xh, Xh, M_Aq[0] ) += on( markedfaces( mesh,"cabine_inlet" ), u, F, cst( 20. ) );
    M_Aq[0]->close();


    M = backend->newMatrix( Xh, Xh );

    form2( Xh, Xh, M, _init=true ) =
        integrate( elements( mesh ), id( u )*idt( v ) + grad( u )*trans( gradt( u ) ) );
    M->close();



    ////////////////////////////////////////////////////////////////////////////////


} // Cabineheat::init

CabineHeat::sparse_matrix_ptrtype
CabineHeat::newMatrix() const
{
    return backend->newMatrix( Xh, Xh );
}

CabineHeat::affine_decomposition_type
CabineHeat::computeAffineDecomposition()
{
    return boost::make_tuple( M_Aq, M_Fq );
}


void
CabineHeat::solve( sparse_matrix_ptrtype& D,
                   element_type& t,
                   vector_ptrtype& F )
{
    vector_ptrtype T( backend->newVector( t.functionSpace() ) );
    backend->solve( D, D, T, F );
    t = *T;
} // Heat1d::solve


void
CabineHeat::exportResults( element_type& U )
{
    if ( M_do_export )
    {
        LOG(INFO) << "exportResults starts\n";

        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );

        exporter->step( 0 )->add( "u", U );

        exporter->save();
    }
} // Heat1d::export

void
CabineHeat::update( parameter_type const& mu )
{

    *D = *M_Aq[0];

    for ( size_type q = 0; q < M_Aq.size(); ++q )
    {
        D->addMatrix( M_thetaAq[q], M_Aq[q] );
    }

    F->close();
    F->zero();

    for ( size_type q = 0; q < M_Fq[0].size(); ++q )
    {
        F->add( M_thetaFq[0][q], M_Fq[0][q] );
    }
}
void
CabineHeat::solve( parameter_type const& mu )
{
    element_ptrtype T( new element_type( Xh ) );
    this->solve( mu, T );
    //this->exportResults( *T );

}

void
CabineHeat::solve( parameter_type const& mu, element_ptrtype& T )
{
    this->computeThetaq( mu );
    this->update( mu );
    backend->solve( _matrix=D,  _solution=T, _rhs=F );
}

void
CabineHeat::l2solve( vector_ptrtype& u, vector_ptrtype const& f )
{
    //std::cout << "l2solve(u,f)\n";
    backend->solve( _matrix=M,  _solution=u, _rhs=f );
    //std::cout << "l2solve(u,f) done\n";
}

double
CabineHeat::scalarProduct( vector_ptrtype const& x, vector_ptrtype const& y )
{
    return M->energy( x, y );
}
double
CabineHeat::scalarProduct( vector_type const& x, vector_type const& y )
{
    return M->energy( x, y );
}

void
CabineHeat::run( const double * X, unsigned long N, double * Y, unsigned long P )
{
    using namespace vf;
    Feel::ParameterSpace<5>::Element mu( M_Dmu );
    mu << X[0], X[1], X[2], X[3], X[4];
    static int do_init = true;

    if ( do_init )
    {
        meshSize = X[5];
        this->init();
        do_init = false;
    }

    this->solve( mu, pT );
}


double
CabineHeat::output( int output_index, parameter_type const& mu )
{
    std::cout<<"model output"<<std::endl;

    using namespace vf;

    this->solve( mu, pT );

    vector_ptrtype U( backend->newVector( Xh ) );
    *U = *pT;

    // (compliant) and (non compliant: mean temperature on plancher)
    double s=0;

    if ( output_index<4 )
    {
        for ( int i=0; i<Ql( output_index ); i++ )  s += M_thetaFq[output_index]( i )*dot( M_Fq[output_index][i], U );
    }

    else
    {
        throw std::logic_error( "[Cabineheat::output] error with output_index : only 0 or 1 " );
    }

    return s;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
}


}

#endif /* __CABINEHEAT_H */


