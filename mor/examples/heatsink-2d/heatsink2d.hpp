/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Universite Joseph Fourier (Grenoble I)

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
   \file heatSink.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __HeatSink2D_H
#define __HeatSink2D_H 1

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
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelalg/solvereigen.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelmor/parameterspace.hpp>

#include <feel/feelts/bdf.hpp>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>


#include <feel/feelmor/modelcrbbase.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description makeHeatSink2DOptions();
FEELPP_EXPORT AboutData makeHeatSink2DAbout( std::string const& str = "heatSink" );

class FEELPP_EXPORT FunctionSpaceDefinition
{
public :
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef fusion::vector<Lagrange<1, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    /*element*/
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};

/**
 * \class HeatSink2D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class FEELPP_EXPORT HeatSink2D : public ModelCrbBase< ParameterSpaceX, FunctionSpaceDefinition, 1 >
{
public:

    typedef ModelCrbBase<ParameterSpaceX,FunctionSpaceDefinition, 1> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;


    /** @name Constants
     */
    //@{

    static const uint16_type Order = 1;
    //@}

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_type vector_type;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix_type;
    typedef eigen_matrix_type ematrix_type;
    typedef std::shared_ptr<eigen_matrix_type> eigen_matrix_ptrtype;

    /*mesh*/
    typedef Simplex<2,Order> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef FunctionSpace<mesh_type, fusion::vector<Lagrange<0, Scalar> >, Discontinuous> p0_space_type;
    typedef typename p0_space_type::element_type p0_element_type;

    /*basis*/
    typedef fusion::vector<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef space_type functionspace_type;
    typedef space_ptrtype functionspace_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef std::shared_ptr<element_type> element_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef std::shared_ptr<export_type> export_ptrtype;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;


    //typedef std::vector< std::vector< double > > beta_vector_type;

    /*typedef boost::tuple<
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<sparse_matrix_ptrtype> >,
        std::vector< std::vector<std::vector<vector_ptrtype> > >
     > affine_decomposition_type;*/

    //typedef boost::tuple<std::vector<sparse_matrix_ptrtype>, std::vector<std::vector<vector_ptrtype>  > > affine_decomposition_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    HeatSink2D();

    //! copy constructor
    HeatSink2D( HeatSink2D const & );
    //! destructor
    ~HeatSink2D() {}

    //! initialization of the model
    void initModel();
    //@}
    /** @name Accessors
     */
    //@{

    // \return the number of terms in affine decomposition of left hand
    // side bilinear form
    int Qa() const
    {
        return 4;
    }

    // \return the number of terms in affine decomposition of bilinear form
    // associated to mass matrix
    int Qm() const
    {
        return 2;
    }

    /**
     * there is at least one output which is the right hand side of the
     * primal problem
     *
     * \return number of outputs associated to the model
     * in our case we have a compliant output and 2 others outputs : average temperature on boundaries
     */
    int Nl() const
    {
        return 2;
    }

    /**
     * \param l the index of output
     * \return number of terms  in affine decomposition of the \p q th output term
     * in our case no outputs depend on parameters
     */
    int Ql( int l ) const
    {
        return 1*( l==0 ) + 1*( l>0 );
    }

    int mMaxA( int q )
    {
        if( q < this->Qa() )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxM( int q )
    {
        if ( q < this->Qm() )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxM(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q < this->Ql(output_index) )
            return 1;
        else
            throw std::logic_error( "[Model heatsink] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }

    /**
     * \brief compute the beta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */

    betaqm_type
    computeBetaQm( element_type const& T, parameter_type const& mu, double time , bool only_terms_time_dependent=false )
    {
        return computeBetaQm( mu, time, only_terms_time_dependent );
    }

    betaqm_type
    computeBetaQm( vectorN_type const& urb, parameter_type const& mu, double time , bool only_terms_time_dependent=false  )
    {
        return computeBetaQm( mu, time, only_terms_time_dependent );
    }

    betaqm_type
    computeBetaQm( parameter_type const& mu , double time , bool only_terms_time_dependent=false )
    {
        double biot      = mu( 0 );
        double L         = mu( 1 );
        double k         = mu( 2 );
        double detJ = L/Lref;
        //coefficient from J^{-T}J^{-1}[4]
        double JJ4 = ( Lref*Lref )/( L*L );
        double t = M_bdf->time();

        M_betaMqm.resize( Qm() );
        M_betaMqm[0].resize(1);
        M_betaMqm[1].resize(1);
        M_betaMqm[0][0]=rho*C/k_fin;
        M_betaMqm[1][0]=rho*C/k_fin * detJ;

        M_betaAqm.resize( Qa() );
        for(int i=0; i<Qa(); i++)
            M_betaAqm[i].resize(1);
        M_betaAqm[0][0] = k ;
        M_betaAqm[1][0] = 1 ;
        M_betaAqm[2][0] = JJ4*detJ;
        M_betaAqm[3][0] = biot * detJ;

        M_betaFqm.resize( Nl() );
        for(int i=0; i<Nl(); i++)
        {
            M_betaFqm[i].resize( Ql(i) );
            for(int j=0;j<Ql(i);j++)
                M_betaFqm[i][j].resize(1);
        }

        M_betaFqm[0][0][0] = 1-math::exp( -t );
        M_betaFqm[1][0][0] = 2;

        return boost::make_tuple( M_betaMqm, M_betaAqm, M_betaFqm );
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * create a new matrix
     * \return the newly created matrix
     */
    sparse_matrix_ptrtype newMatrix() const;

    vector_ptrtype newVector() const;

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition();

    std::vector< std::vector< element_ptrtype > > computeInitialGuessAffineDecomposition()
    {
        std::vector< std::vector<element_ptrtype> > q;
        q.resize(1);
        q[0].resize(1);
        element_ptrtype elt ( new element_type ( Xh ) );
        q[0][0] = elt;
        return q;
    }

    /**
     * \brief solve the model for parameter \p mu
     * \param mu the model parameter
     * \param T the temperature field
     */
    void solve( parameter_type const& mu, element_ptrtype& T, int output_index=0 );


    int computeNumberOfSnapshots();

    void assemble();
    double timeFinal()
    {
        return M_bdf->timeFinal();
    }
    double timeStep()
    {
        return  M_bdf->timeStep();
    }
    double timeInitial()
    {
        return M_bdf->timeInitial();
    }
    int timeOrder()
    {
        return M_bdf->timeOrder();
    }
    void initializationField( element_ptrtype& initial_field, parameter_type const& mu );
    bool isSteady()
    {
        return M_is_steady;
    }


    /**
     * solve for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu );

    /**
     * solve \f$ M u = f \f$
     */
    void l2solve( vector_ptrtype& u, vector_ptrtype const& f );

    /**
     * H1 scalar product
     */
    sparse_matrix_ptrtype energyMatrix ( void )
    {
        return M;
    }


    /**
     * update the PDE system with respect to \param mu
     */
    void update( parameter_type const& mu,double bdf_coeff, element_type const& bdf_poly, int output_index=0 ) ;
    //@}

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( double time, element_type& T , parameter_type const& mu );

    void exportOutputs( double time, double output1, double output2 );

    void solve( sparse_matrix_ptrtype& ,element_type& ,vector_ptrtype&  );

    /**
     * returns the scalar product of the std::shared_ptr vector x and
     * std::shared_ptr vector y
     */
    double scalarProduct( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product of the vector x and vector y
     */
    double scalarProduct( vector_type const& x, vector_type const& y );

    /**
     * returns the scalar product used in POD of the std::shared_ptr vector x and
     * std::shared_ptr vector y
     */
    double scalarProductForPod( vector_ptrtype const& X, vector_ptrtype const& Y );

    /**
     * returns the scalar product used in POD of the vector x and vector y
     */
    double scalarProductForPod( vector_type const& x, vector_type const& y );

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
    value_type output( int output_index, parameter_type const& mu, element_type& u, bool need_to_solve=false );

    gmsh_ptrtype createGeo( double mu2 );

    parameter_type refParameter()
    {
        return Dmu->min();
    }

    bdf_ptrtype bdfModel(){ return M_bdf; }


private:
    bool M_is_steady ;

    double Lref;//reference value for geometric parameter

    int export_number;

    bool do_export;

    double rho;
    double C;
    double k_fin;

    /* surfaces*/
    double surface_gamma4, surface_gamma1;

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;

    sparse_matrix_ptrtype D,M,Mpod;
    vector_ptrtype F;

    element_ptrtype pT;

    bdf_ptrtype M_bdf;
    int M_Nsnap;

};

}

#endif /* __HeatSink2D_H */
