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
#include <feel/feelcrb/crbmodel.hpp>

namespace Feel
{
template<typename ModelType>
class CRBModelSaddlePoint :
        public CRBModel<ModelType>
{
    typedef CRBModel<ModelType> super;
public :
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    //! reduced basis function space type
    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename model_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::vector_type vector_type;

    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename model_type::parameter_type parameter_type;

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


    static const int nb_spaces = space_type::nSpaces;

    //@{ /// Constructors
    /// Default
    CRBModelSaddlePoint() :
        super()
        {
            this->init();
        }

    CRBModelSaddlePoint( po::variables_map const& vm,
                         CRBModelMode mode=CRBModelMode::PFEM ) :
        super( vm, mode )
        {
            this->init();
        }

    CRBModelSaddlePoint( model_ptrtype & model ) :
        super( model )
        {
            this->init();
        }
    CRBModelSaddlePoint( CRBModelSaddlePoint const & o ) :
        super( o )
        {
            this->init();
        }
    //@}

    /// Destructor
    virtual ~CRBModelSaddlePoint() {}

    virtual size_type Qb() const
        {
            return M_Qb;
        }

    virtual size_type Qg( int output_index) const
        {
            return M_Qg[output_index];
        }

protected :
    int M_Qb;
    std::vector<int> M_Qg;


}; // class CRBModelSaddlepoint

template<typename ModelType>
void
CRBModelSaddlePoint<ModelType>::countAffineDecompositionTerms()
{
    if( this->M_alreadyCountAffineDecompositionTerms )
        return;
    else
        this->M_alreadyCountAffineDecompositionTerms=true;

    if ( M_Aqm.size() > 0)
    {

    }
} // countAffineDecompositionterms

} // namespace Feel

#endif // __CRBModelSaddlePoint_H
