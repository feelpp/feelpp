/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-08-09

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
 * @file   crbmodelsaddlepoint.hpp
 * @author wahl
 * @date   Tue Nov 25 16:24:55 2014
 */


#ifndef __CRBModelSaddlePoint_H
#define __CRBModelSaddlePoint_H 1

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbmodelbase.hpp>



namespace Feel
{
class FunctionSpaceDefinitionSaddlePoint
{
public :
    /*mesh*/
    typedef Simplex<1> entity_type ;
    typedef Mesh<entity_type > mesh_type ;

    /*basis*/
    typedef Lagrange< 2, Vectorial, Continuous >  basis_u_type ;
    typedef Lagrange< 1, Scalar, Continuous > basis_p_type ;
    typedef bases< basis_u_type, basis_p_type > basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type , basis_type > space_type ;
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent=false;
    static const bool is_linear=true;
};


template < typename ParameterDefinition=ParameterSpace<1>,
           typename FunctionSpaceDefinition=FunctionSpaceDefinitionSaddlePoint,
           int _Options=0,
           typename EimDefinition=EimDefinitionBase<ParameterDefinition,FunctionSpaceDefinition>
           >
class CRBModelSaddlePoint :
        public CRBModelBase< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition >
{
    typedef CRBModelBase< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition > super_type;
    typedef CRBModelSaddlePoint< ParameterDefinition, FunctionSpaceDefinition, _Options, EimDefinition > self_type;

public :
    typedef typename super_type::value_type value_type;

    typedef typename super_type::mesh_type mesh_type;
    typedef typename super_type::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename super_type::space_type space_type;
    typedef typename super_type::element_type element_type;
    typedef typename super_type::element_ptrtype element_ptrtype;

    //! reduced basis function space type
    typedef typename super_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename super_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;

    typedef typename super_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename super_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename super_type::vector_ptrtype vector_ptrtype;
    typedef typename super_type::vector_type vector_type;

    typedef typename super_type::parameterspace_type parameterspace_type;
    typedef typename super_type::parameter_type parameter_type;

    typedef Eigen::VectorXd vectorN_type;
    typedef std::vector< std::vector< double > > beta_vector_type;

    typedef typename boost::tuple<sparse_matrix_ptrtype,
                                  sparse_matrix_ptrtype,
                                  std::vector<vector_ptrtype>
                                  > offline_merge_type;

    typedef typename boost::tuple<std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector<sparse_matrix_ptrtype> >,
                                  std::vector< std::vector< std::vector<vector_ptrtype> > >
                                  > affine_decomposition_type;

    typedef typename boost::tuple< beta_vector_type,
                                   beta_vector_type,
                                   std::vector<beta_vector_type>
                                   > betaqm_type;




    //@{ /// Constructors
    /// Default
    CRBModelSaddlePoint() :
        super_type()
        {}


protected :


}; // class CRBModelSaddlepoint

} // namespace Feel

#endif // __CRBModelSaddlePoint_H
