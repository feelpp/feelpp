/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2014-01-19

  Copyright (C) 2011-2014 Feel++ Consortium

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
   \file benchmarkgrepl.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \author Cecile Daversin <daversin@math.unistra.fr>

   date 2014-01-19
 */
#ifndef FEELPP_BENCHMARKGREPLNONLINEARELLIPTIC_HPP
#define FEELPP_BENCHMARKGREPLNONLINEARELLIPTIC_HPP 1

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelcrb/modelcrbbase.hpp>
#include <benchmarkgrepl-options.hpp>

namespace Feel
{

namespace BenchmarkGreplNonlinearElliptic_Definition
{
class ParameterDefinition
{
public :
    //static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace</*ParameterSpaceDimension*/> parameterspace_type;
};


template<int Order, int Dim>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<Dim,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;
    typedef bases<Lagrange<Order + 2, Scalar> > basis_type_eimg;
    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef FunctionSpace<mesh_type, basis_type_eimg, value_type> space_type_eimg;

    typedef typename space_type::element_type element_type;

};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef typename FunctionSpaceDefinition::space_type_eimg space_type_eimg;

    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type_eimg, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type_eimg, space_type , parameterspace_type> fund_type;

};

} // namespace BenchmarkGreplNonlinearElliptic_Definition

/**
 * \class BenchmarkGreplNonlinearElliptic
 * \brief brief description
 *
 * This is from the paper
 * EFFICIENT REDUCED-BASIS TREATMENT OF NONAFFINE
 * AND NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS
 *  authors :
 * Martin A. Grepl, Yvon Maday, Ngoc C. Nguyen and Anthony T. Patera
 * ESAIM: Mathematical Modelling and Numerical Analysis
 * --
 * 5. Nonaffine monotonic elliptic equations
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
template<int Order, int Dim>
class FEELPP_EXPORT BenchmarkGreplNonlinearElliptic :
        public ModelCrbBase< BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                             BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>,
                             NonLinear,
                             BenchmarkGreplNonlinearElliptic_Definition::EimDefinition<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim> > >
{
    typedef ModelCrbBase<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,
                         BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>,
                         NonLinear,
                         BenchmarkGreplNonlinearElliptic_Definition::EimDefinition<BenchmarkGreplNonlinearElliptic_Definition::ParameterDefinition,BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim> > > super_type;
    typedef BenchmarkGreplNonlinearElliptic<Order,Dim> self_type;

 public :
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    typedef double value_type;

    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>::space_type space_type;
    typedef typename std::shared_ptr<space_type> space_ptrtype;
    typedef typename BenchmarkGreplNonlinearElliptic_Definition::FunctionSpaceDefinition<Order,Dim>::space_type_eimg space_type_eimg;
    typedef typename std::shared_ptr<space_type_eimg> space_ptrtype_eimg;

    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::monolithic_type monolithic_type;
    typedef typename super_type::eim_interpolation_error_type eim_interpolation_error_type;

    typedef std::vector< std::vector<sparse_matrix_ptrtype> > vector_sparse_matrix;

    using super_type::computeBetaQm;
    typedef boost::tuple<beta_vector_type,  std::vector<beta_vector_type> > beta_type;



    BenchmarkGreplNonlinearElliptic(std::string const& prefix = "");

    //! initialization of the model
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    //@}

    vector_ptrtype assembleForDEIMnl( parameter_type const& mu, element_type const& u, int const& tag ) override
    {
        auto Xh = this->Xh;
        vector_ptrtype V = this->M_backend->newVector(Xh);
        auto temp = Xh->element(V,0);
        mesh = Xh->mesh();
        temp.on( _range=elements(mesh), _expr=cst(mu(0))/cst(mu(1))*( exp( cst(mu(1))*idv(u) ) - 1  ) );
        return V;
    }


    virtual beta_vector_type computeBetaInitialGuess( parameter_type const& mu ) override
    {

        this->M_betaInitialGuess.resize( 1 );
        this->M_betaInitialGuess[0].resize( 1 );

        this->M_betaInitialGuess[0][0] = 1;

        return this->M_betaInitialGuess;
    }


    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    beta_type
    computeBetaQm( element_type const& T,parameter_type const& mu ) override
    {
        std::vector<vectorN_type*> betas;
        if ( M_use_deim )
        {
            auto beta = this->deim()->beta(mu,T);
            betas.push_back( &beta );
        }

        else
        {
            auto beta = this->scalarContinuousEim()[0]->beta( mu , T );
            betas.push_back( &beta );
        }

        fillBetaQm(betas, mu);

        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computeBetaQm( vectorN_type const& urb, parameter_type const& mu ) override
    {
        std::vector<vectorN_type*> betas;
        if ( M_use_deim )
        {
            auto beta = this->deim()->beta(mu,urb);
            betas.push_back( &beta );
        }

        else
        {
            auto beta = this->scalarContinuousEim()[0]->beta( mu ,urb );
            betas.push_back( &beta );
        }

        fillBetaQm(betas, mu);

        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computeBetaQm( parameter_type const& mu ) override
    {
        std::vector<vectorN_type*> betas;
        if ( M_use_deim )
        {
            auto beta = this->deim()->beta(mu);
            betas.push_back( &beta );
        }

        else
        {
            auto beta = this->scalarContinuousEim()[0]->beta(mu);
            betas.push_back( &beta );
        }

        fillBetaQm(betas, mu);

        if( M_use_newton )
            return boost::make_tuple( this->M_betaJqm, this->M_betaRqm);
        else
            return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);
    }

    beta_type
    computePicardBetaQm( element_type const& T,parameter_type const& mu ) override
    {
        std::vector<vectorN_type*> betas;
        if ( M_use_deim )
        {
            auto beta = this->deim()->beta(mu,T);
            betas.push_back( &beta );
        }

        else
        {
            auto beta = this->scalarContinuousEim()[0]->beta(mu,T);
            betas.push_back( &beta );
        }

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    beta_type
    computePicardBetaQm( parameter_type const& mu ) override
    {
        std::vector<vectorN_type*> betas;
        if ( M_use_deim )
        {
            auto beta = this->deim()->beta(mu);
            betas.push_back( &beta );
        }

        else
        {
            auto beta = this->scalarContinuousEim()[0]->beta(mu);
            betas.push_back( &beta );
        }

        fillBetaQm(betas, mu);
        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm);

    }

    void initBetaQm( int M )
    {
        if( M_use_newton )
        {
            this->M_betaJqm.resize( 3 );
            this->M_betaJqm[0].resize( 1 );
            this->M_betaJqm[1].resize( M );
            this->M_betaJqm[2].resize( 1 );
            this->M_betaRqm.resize( 2 );
            this->M_betaRqm[0].resize( 3 );
            this->M_betaRqm[0][0].resize( 1 );
            this->M_betaRqm[0][1].resize( M );
            this->M_betaRqm[0][2].resize( 1 );
            this->M_betaRqm[1].resize( 1 );
            this->M_betaRqm[1][0].resize( 1 );
        }
        if( !M_use_newton || M_useSerErrorEstimation )
        {
            this->M_betaAqm.resize( 1 );
            this->M_betaAqm[0].resize( 1 );
            this->M_betaFqm.resize(2);
            this->M_betaFqm[0].resize( 2 );
            this->M_betaFqm[0][0].resize( M );
            this->M_betaFqm[0][1].resize( 1 );
            this->M_betaFqm[1].resize( 1 );
            this->M_betaFqm[1][0].resize( 1 );
        }
    }

    void fillBetaQm(std::vector<vectorN_type*> betas, parameter_type const& mu)
    {
        int M;
        if ( M_use_deim )
            M = this->deim()->size();
        else
            M = this->scalarContinuousEim()[0]->mMax();

        auto beta_g=*betas[0];

        if( M_use_newton )
        {
            if ( this->M_betaJqm.empty() )
                this->initBetaQm( M );
            else
            {
                //needed if cobuild
                this->M_betaJqm[1].resize( M );
                this->M_betaRqm[0][1].resize( M );
            }

            this->M_betaJqm[0][0] = 1;
            for(int m=0; m<M; m++)
            {
                this->M_betaJqm[1][m] = mu(1)*beta_g(m);
            }
            this->M_betaJqm[2][0] = mu(0);

            this->M_betaRqm[0][0][0] = this->computeBetaInitialGuess( mu )[0][0];
            for(int m=0; m<M; m++)
            {
                this->M_betaRqm[0][1][m] = beta_g(m);
            }
            this->M_betaRqm[0][2][0] = -100;
            //output
            this->M_betaRqm[1][0][0] = 1;
        }
        //else
        if( !M_use_newton || M_useSerErrorEstimation )
        {
            if ( this->M_betaAq.empty() )
                this->initBetaQm( M );
            else
                this->M_betaFqm[0][0].resize( M ); //needed if cobuild

            this->M_betaAqm[0][0]=1;
            for(int m=0; m<M; m++)
            {
                this->M_betaFqm[0][0][m]=-beta_g(m) ;
            }
            this->M_betaFqm[0][1][0]=100;
            this->M_betaFqm[1][0][0]=1;
        }

    }

    beta_vector_type computeBetaLinearDecompositionA( parameter_type const& mu , double time=1e30 ) override
    {
        beta_vector_type beta;
        beta.resize(1);
        beta[0].resize(1);
        beta[0][0]=1;
        return beta;
    }

    /**
     * \brief Returns the affine decomposition
     */
    affine_decomposition_type computeAffineDecomposition() override;
    affine_decomposition_type computePicardAffineDecomposition() override;
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;

    void assemble() override;

    std::vector< std::vector<element_ptrtype> > computeInitialGuessAffineDecomposition() override;


    //@}


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false ) override;

    gmsh_ptrtype createGeo( double hsize );

    void assembleResidualWithAffineDecomposition(std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm);
    void assembleJacobianWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & Jqm);
    void assembleFunctionalWithAffineDecomposition(std::vector<std::vector<sparse_matrix_ptrtype> > & RF_Aqm,
                                                   std::vector< std::vector<std::vector<vector_ptrtype> > >& RF_Fqm);
    bool updateResidual(element_type const& X, std::vector< std::vector<std::vector<vector_ptrtype> > >& Rqm) override;
    void updateResidualMonolithic(vector_ptrtype const& X, vector_ptrtype & R, parameter_type const& mu);
    void updateJacobianMonolithic(vector_ptrtype const& X, sparse_matrix_ptrtype & J, parameter_type const& mu);
    monolithic_type computeMonolithicFormulationU( parameter_type const& mu , element_type const& solution ) override;

    element_type solve( parameter_type const& mu ) override;

private:

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;
    element_ptrtype pT;
    parameter_type M_mu;

    sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    bool M_use_newton;
    bool M_useSerErrorEstimation;

    std::vector< std::vector< element_ptrtype > > M_InitialGuess;

    bool M_use_deim;

};


}

#endif /* __BenchmarkGreplNonlinearElliptic_H */
