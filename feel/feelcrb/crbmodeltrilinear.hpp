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
   \file crbmodeltrilinear.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-14
 */
#ifndef __CRBModelTrilinear_H
#define __CRBModelTrilinear_H 1

#include <boost/shared_ptr.hpp>

#include <vector>


#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbmodel.hpp>


namespace Feel
{

/**
 * \class CRBModelTrilinear
 * \brief Certified Reduced Basis Model trilinear class
 *
 * This class implements the requirements over a model to be usable by the
 * certified reduced basis method.
 *
 * \tparam ModelType the type of the finite element method model
 *
 * The FEM model type should derive from this class and fill the vector of
 * matrices M_Aq and vector of vectors M_Fq
 *
 * @author Christophe Prud'homme
 * @see crb
 */
template<typename ModelType>
class CRBModelTrilinear : public CRBModel<ModelType>
{
    typedef  CRBModel<ModelType> super;

public:


    /** @name Constants
     */
    //@{

    //static const uint16_type ParameterSpaceDimension = ModelType::ParameterSpaceDimension;
    //static const bool is_time_dependent = ModelType::is_time_dependent;

    //@}

    /** @name Typedefs
     */
    //@{

    //! model type
    typedef ModelType model_type;
    typedef boost::shared_ptr<ModelType> model_ptrtype;

    //! value_type
    typedef typename model_type::value_type value_type;
    //! mesh type
    typedef typename ModelType::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    //! space_type
    typedef typename ModelType::space_type space_type;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    //! reduced basis function space type
    typedef typename model_type::rbfunctionspace_type rbfunctionspace_type;
    typedef typename model_type::rbfunctionspace_ptrtype rbfunctionspace_ptrtype;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename model_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename model_type::vector_ptrtype vector_ptrtype;
    typedef typename model_type::vector_type vector_type;

    typedef typename model_type::eigen_matrix_type eigen_matrix_type;

    typedef typename model_type::parameterspace_type parameterspace_type;
    typedef typename model_type::parameterspace_ptrtype parameterspace_ptrtype;
    typedef typename model_type::parameter_type parameter_type;
    typedef typename model_type::parameter_ptrtype parameter_ptrtype;
    typedef typename model_type::sampling_type sampling_type;
    typedef typename model_type::sampling_ptrtype sampling_ptrtype;


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


    static const int nb_spaces = functionspace_type::nSpaces;
    typedef typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<2> > , fusion::vector< mpl::int_<0>, mpl::int_<1> >  ,
                       typename mpl::if_ < boost::is_same< mpl::int_<nb_spaces> , mpl::int_<3> > , fusion::vector < mpl::int_<0> , mpl::int_<1> , mpl::int_<2> > ,
                                  typename mpl::if_< boost::is_same< mpl::int_<nb_spaces> , mpl::int_<4> >, fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3> >,
                                                     fusion::vector< mpl::int_<0>, mpl::int_<1>, mpl::int_<2>, mpl::int_<3>, mpl::int_<4> >
                                                     >::type >::type >::type index_vector_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    CRBModelTrilinear()
        :
        super()
    {
        this->init();
    }

    CRBModelTrilinear( po::variables_map const& vm, CRBModelMode mode = CRBModelMode::PFEM  )
        :
        super( vm, mode )
    {
        this->init();
    }

    /**
     * \param model the model to be used
     */
    CRBModelTrilinear( model_ptrtype & model )
        :
        super( model )
    {
        this->init();
    }

    /**
     * copy constructor
     */
    CRBModelTrilinear( CRBModelTrilinear const & o )
        :
        super( o )
    {
        this->init();
    }

    //! destructor
    virtual ~CRBModelTrilinear()
    {}


    //@}

    /** @name Operator overloads
     */
    //@{


    /** @name Accessors
     */
    //@{


    //! return the number of \f$\mu\f$ independent terms for the bilinear form
    size_type Qa() const
    {
        return this->M_model->Qa();
    }


    //! return the number of \f$\mu\f$ independent terms for the trilinear form
    size_type QaTri() const
    {
        return this->M_model->QaTri();
    }

    size_type Ql( int l ) const
    {
        return this->M_model->Ql( l );
    }


    sparse_matrix_ptrtype computeTrilinearForm( element_type const& xi )
    {
        return this->M_model->computeTrilinearForm( xi );
    }

    sparse_matrix_ptrtype jacobian( element_type const& xi )
    {
        return this->M_model->jacobian( xi );
    }

    vector_ptrtype residual( element_type const& xi )
    {
        return this->M_model->residual( xi );
    }



    /**
     * solve the model for a given parameter \p mu
     */
    element_type solve( parameter_type const& mu )
    {
        return this->M_model->solve( mu );
    }



protected:

    //! affine decomposition terms for the left hand side ( bilinear )
    //std::vector< std::vector<sparse_matrix_ptrtype> > M_Aqm;

    //! affine decomposition terms for the right hand side
    //std::vector< std::vector<std::vector<vector_ptrtype> > > M_Fqm;


};


}
#endif /* __CRBModelTrilinear_H */
