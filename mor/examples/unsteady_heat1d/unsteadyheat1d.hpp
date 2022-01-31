/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-13

  Copyright (C) 2009 Universit� Joseph Fourier (Grenoble I)

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
   \file unsteadyHeat1d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-11-13
 */
#ifndef __UnsteadyHeat1D_H
#define __UnsteadyHeat1D_H 1

#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>

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
#include <feel/feelmor/reducedbasisspace.hpp>

namespace Feel
{

FEELPP_EXPORT po::options_description
makeUnsteadyHeat1DOptions();

FEELPP_EXPORT AboutData
makeUnsteadyHeat1DAbout( std::string const& str = "unsteadyHeat1d" );


class FEELPP_EXPORT FunctionSpaceDefinition
{
public :
    static const uint16_type Order = 1;

    typedef double value_type;

    /*mesh*/
    typedef Simplex<1,1> entity_type;
    typedef Mesh<entity_type> mesh_type;

    /*basis*/
    typedef bases<Lagrange<Order, Scalar> > basis_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef typename space_type::element_type element_type;

    static const bool is_time_dependent = true;
    static const bool is_linear = true;

};


/**
 * \class UnsteadyHeat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class FEELPP_EXPORT UnsteadyHeat1D : public ModelCrbBase<ParameterSpaceX,FunctionSpaceDefinition>
{
public:
    typedef ModelCrbBase<ParameterSpaceX,FunctionSpaceDefinition> super_type;
    typedef typename super_type::funs_type funs_type;
    typedef typename super_type::funsd_type funsd_type;
    typedef typename super_type::beta_vector_light_type beta_vector_light_type;
    typedef typename super_type::affine_decomposition_light_type affine_decomposition_light_type;
    using super_type::computeBetaQm;

    /** @name Typedefs
     */
    //@{

    typedef double value_type;

    typedef typename FunctionSpaceDefinition::mesh_type mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename FunctionSpaceDefinition::basis_type basis_type;

    typedef Backend<value_type> backend_type;
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*space*/
    typedef typename FunctionSpaceDefinition::space_type space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    /* parameter space */
    using parameterspace_type = ParameterSpaceX;
    typedef std::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;

    /* time discretization */
    typedef Bdf<space_type>  bdf_type;
    typedef std::shared_ptr<bdf_type> bdf_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    UnsteadyHeat1D();

    //! constructor from command line
    UnsteadyHeat1D( po::variables_map const& vm );


    //! copy constructor
    UnsteadyHeat1D( UnsteadyHeat1D const & );
    //! destructor
    ~UnsteadyHeat1D() {}

    //! initialisation of the model
    void initModel();
    //@}


    /** @name  Methods
     */
    //@{


    void assemble();

    value_type output( int output_index, parameter_type const& mu, element_type &u, bool need_to_solve=false );

    bdf_ptrtype bdfModel(){ return M_bdf; }

private:

    double alpha;

    double meshSize;
    mesh_ptrtype mesh;
    bdf_ptrtype M_bdf;
};


}

#endif /* __Heat1D_H */
