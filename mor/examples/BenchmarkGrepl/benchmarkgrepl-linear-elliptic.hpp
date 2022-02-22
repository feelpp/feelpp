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

   date 2014-01-19
 */
#ifndef FEELPP_BENCHMARKGREPLLINEARELLIPTIC_HPP
#define FEELPP_BENCHMARKGREPLLINEARELLIPTIC_HPP 1

#include <feel/feelmor/modelcrbbase.hpp>

#include <benchmarkgrepl-options.hpp>

namespace Feel
{

namespace BenchmarkGreplLinearElliptic_Definition
{

class ParameterDefinition
{
public :
    //static const uint16_type ParameterSpaceDimension = 2;
    typedef ParameterSpace</*ParameterSpaceDimension*/> parameterspace_type;
};

template<int Order>
class FunctionSpaceDefinition
{
public :
    //static const uint16_type Order = 1;
    typedef double value_type;

    /*mesh*/
    typedef Simplex<2,1> entity_type; /*dim,order*/
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    typedef typename space_type::element_type element_type;

};

template <typename ParameterDefinition, typename FunctionSpaceDefinition >
class EimDefinition
{
public :
    typedef typename ParameterDefinition::parameterspace_type parameterspace_type;
    typedef typename FunctionSpaceDefinition::space_type space_type;


    /* EIM */
    // Scalar continuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fun_type;
    // Scalar Discontinuous //
    typedef EIMFunctionBase<space_type, space_type , parameterspace_type> fund_type;

};

} // namespace BenchmarkGreplLinearElliptic_Definition

/**
 * \class BenchmarkGreplLinearElliptic
 * \brief brief description
 *
 * This is from the paper
 * EFFICIENT REDUCED-BASIS TREATMENT OF NONAFFINE
 * AND NONLINEAR PARTIAL DIFFERENTIAL EQUATIONS
 *  authors :
 * Martin A. Grepl, Yvon Maday, Ngoc C. Nguyen and Anthony T. Patera
 * ESAIM: Mathematical Modelling and Numerical Analysis
 * --
 * 3. Nonaffine linear coercive elliptic equations
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 * @see
 */
template<int Order>
class FEELPP_EXPORT BenchmarkGreplLinearElliptic : public ModelCrbBase< BenchmarkGreplLinearElliptic_Definition::ParameterDefinition,
                                                                        BenchmarkGreplLinearElliptic_Definition::FunctionSpaceDefinition<Order>,
                                                                        TimeIndependent,
                                                                        BenchmarkGreplLinearElliptic_Definition::EimDefinition<BenchmarkGreplLinearElliptic_Definition::ParameterDefinition,
                                                                                                                               BenchmarkGreplLinearElliptic_Definition::FunctionSpaceDefinition<Order> > >
{
public:

    typedef ModelCrbBase<BenchmarkGreplLinearElliptic_Definition::ParameterDefinition,
                         BenchmarkGreplLinearElliptic_Definition::FunctionSpaceDefinition<Order>,
                         TimeIndependent,
                         BenchmarkGreplLinearElliptic_Definition::EimDefinition<BenchmarkGreplLinearElliptic_Definition::ParameterDefinition,
                                                                                BenchmarkGreplLinearElliptic_Definition::FunctionSpaceDefinition<Order> > > super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;

    typedef double value_type;

    typedef typename super_type::element_ptrtype element_ptrtype;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::parameter_type parameter_type;
    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    typedef typename BenchmarkGreplLinearElliptic_Definition::FunctionSpaceDefinition<Order>::space_type space_type;

    typedef typename super_type::beta_vector_type beta_vector_type;
    typedef typename super_type::affine_decomposition_type affine_decomposition_type;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vectorN_type vectorN_type;
    typedef typename super_type::monolithic_type monolithic_type;
    typedef typename super_type::eim_interpolation_error_type eim_interpolation_error_type;

    typedef std::vector< std::vector<sparse_matrix_ptrtype> > vector_sparse_matrix;

    using super_type::computeBetaQm;

    BenchmarkGreplLinearElliptic();

    //! initialization of the model
    void initModel() override;
    void setupSpecificityModel( boost::property_tree::ptree const& ptree, std::string const& dbDir ) override;

    //@}


    int mMaxA( int q )
    {
        if ( q==0 )
            return 1;
        if( q==1 )
        {
            auto eim_g = this->scalarContinuousEim()[0];
            return eim_g->mMax();
        }
        else
            throw std::logic_error( "[Model Benchmark Grepl] ERROR : try to acces to mMaxA(q) with a bad value of q");
    }

    int mMaxF( int output_index, int q)
    {
        if ( q == 0 )
            return 1;
        if( q == 1 )
        {
            auto eim_g = this->scalarContinuousEim()[0];
            return eim_g->mMax();
        }
        else
            throw std::logic_error( "[Model Benchmark Grepl] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }


    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    boost::tuple<beta_vector_type,  std::vector<beta_vector_type> >
    computeBetaQm( element_type const& T,parameter_type const& mu ) override
    {
        return computeBetaQm( mu );
    }

    boost::tuple<beta_vector_type,  std::vector<beta_vector_type>  >
    computeBetaQm( parameter_type const& mu ) override
    {
        this->M_betaAqm[0][0]=1;
        auto eim_g = this->scalarContinuousEim()[0];
        int eim_g_size = eim_g->mMax();
        vectorN_type beta_g = eim_g->beta( mu );
        this->M_betaAqm[1].resize( eim_g_size );
        for(int m=0; m<eim_g_size; m++)
        {
            this->M_betaAqm[1][m] = beta_g(m);
        }

        this->M_betaFqm[0][0].resize( eim_g_size );
        for(int m=0; m<eim_g_size; m++)
        {
            this->M_betaFqm[0][0][m] = beta_g(m);
        }
        this->M_betaFqm[1][0][0] = 1;

        return boost::make_tuple( this->M_betaAqm, this->M_betaFqm );
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
    std::vector< std::vector<sparse_matrix_ptrtype> > computeLinearDecompositionA() override;
    monolithic_type computeMonolithicFormulation( parameter_type const& mu ) override;

    void assemble() override;


    //@}


    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu, element_type &T, bool need_to_solve=false ) override;

    bool referenceParametersGivenByUser() override { return true; }
    parameter_type refParameter() override
    {
        return this->Dmu->min();
    }

    gmsh_ptrtype createGeo( double hsize );

    void checkEimExpansion();

    eim_interpolation_error_type eimInterpolationErrorEstimation() override
    {
        std::map<int,double> eim_error_aq;
        std::vector< std::map<int,double> > eim_error_fq;

        //in that case we make an error on the eim approximation
        //M_Aqm[1] contains the eim approximation
        eim_error_aq.insert( std::pair<int,double>(1 , 0 ) );
        eim_error_fq.resize(2);
        //M_Fqm[0][0] contains the eim approximation
        eim_error_fq[0].insert( std::pair<int,double>(0 , 0 ) );
        return boost::make_tuple( eim_error_aq, eim_error_fq );
    }

    eim_interpolation_error_type eimInterpolationErrorEstimation( parameter_type const& mu , vectorN_type const& uN ) override
    {

        std::map<int,double> eim_error_aq;
        std::vector< std::map<int,double> > eim_error_fq;

        auto eim_g = this->scalarContinuousEim()[0];
        bool error;
        int max = eim_g->mMax(error);
        if( error )
        {
            //in that case we make an error on the eim approximation
            int size=uN.size();
            auto solution = Feel::expansion( this->XN->primalRB(), uN , size);
            double eim_error = eim_g->interpolationErrorEstimation(mu,solution,max).template get<0>();
            //M_Aqm[1] contains the eim approximation
            eim_error_aq.insert( std::pair<int,double>(1 , eim_error) );
            eim_error_fq.resize(2);
            //M_Fqm[0][0] contains the eim approximation
            eim_error_fq[0].insert( std::pair<int,double>(0 , eim_error) );
        }
        return boost::make_tuple( eim_error_aq, eim_error_fq );
    }


private:

    /* mesh, pointers and spaces */
    mesh_ptrtype mesh;

    sparse_matrix_ptrtype D,M;
    vector_ptrtype F;

    element_ptrtype pT;

    sparse_matrix_ptrtype M_monoA;
    std::vector<vector_ptrtype> M_monoF;

    parameter_type M_mu;

};



template<int Order>
gmsh_ptrtype
BenchmarkGreplLinearElliptic<Order>::createGeo( double hsize )
{
    gmsh_ptrtype gmshp( new Gmsh );
    std::ostringstream ostr;
    double H = hsize;
    ostr <<"Point (1) = {0, 0, 0,"<< H <<"};\n"
         <<"Point (2) = {1, 0, 0,"<< H <<"};\n"
         <<"Point (3) = {1, 1, 0,"<< H <<"};\n"
         <<"Point (4) = {0, 1, 0,"<< H <<"};\n"
         <<"Line (11) = {1,2};\n"
         <<"Line (12) = {2,3};\n"
         <<"Line (13) = {3,4};\n"
         <<"Line (14) = {4,1};\n"
         <<"Line Loop (21) = {11, 12, 13, 14};\n"
         <<"Plane Surface (30) = {21};\n"
         <<"Physical Line (\"boundaries\") = {11,12,13,14};\n"
         <<"Physical Surface (\"Omega\") = {30};\n"
         ;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << "benchmarkgrepl_geo";
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
}


}

#endif /* __BenchmarkGreplLinearElliptic_H */


