/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-12-10

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file opusmodel.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-12-10
 */
#ifndef __OpusModel_H
#define __OpusModel_H 1

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeldiscr/operatorlagrangep1.hpp>
#include <feel/feelpoly/im.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelfilters/exporter.hpp>


#include <opusmodelbase.hpp>
#include <opusdefs.hpp>

namespace Feel
{
/**
 * \addtogroup Models
 * \\@{
 */
template<typename SpaceType> class OpusModelThermal;
template<typename SpaceType> class OpusModelFluidPoiseuille;
template<typename SpaceType> class OpusModelFluidOseen;


/**
 * \class OpusModel
 * \brief The Opus model coupling heat transfer and fluid
 *
 * @author Christophe Prud'homme
 * @see
 */

template<int OrderU=2, int OrderP=OrderU-1, int OrderT=OrderP>
class OpusModel : public OpusModelBase
{
    typedef OpusModelBase super;
public:

    /** @name Typedefs
     */
    //@{
    static const uint16_type Dim = 2;

    //@}
    /** @name Typedefs
     */
    //@{


    typedef double value_type;


    /*mesh*/
    typedef Simplex<Dim> entity_type;
    typedef Mesh<entity_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Simplex<1,1,2> entity12_type;
    typedef Mesh<entity12_type> mesh12_type;
    typedef boost::shared_ptr<mesh12_type> mesh12_ptrtype;


    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    /* temperature */
    typedef DiscontinuousInterfaces<fusion::vector<mpl::vector<mpl::int_<3>, mpl::int_<11>, mpl::int_<13> >,
            mpl::vector<mpl::int_<4>, mpl::int_<11>, mpl::int_<14> >
            > >  discontinuity_type;
#if defined( OPUS_WITH_THERMAL_DISCONTINUITY )
    typedef bases<Lagrange<OrderT, Scalar, discontinuity_type> > temp_basis_type;
#else
    typedef bases<Lagrange<OrderT, Scalar> > temp_basis_type;
#endif
    typedef bases<Lagrange<OrderT, Vectorial> > grad_temp_basis_type;
    //typedef Periodic<1,2,value_type> periodic_type;
    typedef Periodic<> periodic_type;


    typedef FunctionSpace<mesh_type, temp_basis_type, Periodicity<Periodic<> > > temp_functionspace_type;
    typedef boost::shared_ptr<temp_functionspace_type> temp_functionspace_ptrtype;
    typedef typename temp_functionspace_type::element_type temp_element_type;
    typedef boost::shared_ptr<temp_element_type> temp_element_ptrtype;

    typedef FunctionSpace<mesh_type, grad_temp_basis_type> grad_temp_functionspace_type;
    typedef boost::shared_ptr<grad_temp_functionspace_type> grad_temp_functionspace_ptrtype;
    typedef typename grad_temp_functionspace_type::element_type grad_temp_element_type;
    typedef boost::shared_ptr<grad_temp_element_type> grad_temp_element_ptrtype;

    /* 1D profil */
    typedef bases<Lagrange<1, Scalar> > p1_basis_type;
    typedef FunctionSpace<mesh12_type, p1_basis_type, value_type> p1_functionspace_type;
    typedef boost::shared_ptr<p1_functionspace_type> p1_functionspace_ptrtype;
    typedef typename p1_functionspace_type::element_type p1_element_type;
    typedef boost::shared_ptr<p1_element_type> p1_element_ptrtype;

    /* P0 */
    typedef FunctionSpace<mesh_type, bases<Lagrange<0, Scalar, Discontinuous> > > p0_space_type;
    typedef boost::shared_ptr<p0_space_type> p0_space_ptrtype;
    typedef typename p0_space_type::element_type p0_element_type;
    typedef boost::shared_ptr<p0_element_type> p0_element_ptrtype;

    typedef Feel::OpusModelThermal<temp_functionspace_type> thermal_operator_type;
    typedef boost::shared_ptr<thermal_operator_type> thermal_operator_ptrtype;

    /* velocity */
    typedef bases<Lagrange<OrderU, Vectorial> > velocity_basis_type;
    typedef FunctionSpace<mesh_type, velocity_basis_type, value_type> velocity_functionspace_type;
    typedef boost::shared_ptr<velocity_functionspace_type> velocity_functionspace_ptrtype;
    typedef typename velocity_functionspace_type::element_type velocity_element_type;
    typedef boost::shared_ptr<velocity_element_type> velocity_element_ptrtype;

    /* pressure */
    typedef bases<Lagrange<OrderP, Scalar> > pressure_basis_type;
    typedef FunctionSpace<mesh_type, pressure_basis_type, value_type> pressure_functionspace_type;
    typedef boost::shared_ptr<pressure_functionspace_type> pressure_functionspace_ptrtype;
    typedef typename pressure_functionspace_type::element_type pressure_element_type;
    typedef boost::shared_ptr<pressure_element_type> pressure_element_ptrtype;
    /* fluid */
    typedef bases<Lagrange<OrderU, Vectorial>,Lagrange<OrderP, Scalar> > fluid_basis_type;

    typedef FunctionSpace<mesh_type, fluid_basis_type> fluid_functionspace_type;
    typedef boost::shared_ptr<fluid_functionspace_type> fluid_functionspace_ptrtype;
    typedef typename fluid_functionspace_type::element_type fluid_element_type;
    typedef boost::shared_ptr<fluid_element_type> fluid_element_ptrtype;
    typedef typename fluid_element_type::template sub_element<0>::type fluid_element_0_type;
    typedef typename fluid_element_type::template sub_element<1>::type fluid_element_1_type;

#if 1
    typedef Feel::OpusModelFluidPoiseuille<fluid_functionspace_type> fluid_operator_type;
#else
    typedef Feel::OpusModelFluidOseen<fluid_functionspace_type> fluid_operator_type;
#endif
    typedef boost::shared_ptr<fluid_operator_type> fluid_operator_ptrtype;

    /* time */
    typedef Bdf<temp_functionspace_type>  temp_bdf_type;
    typedef boost::shared_ptr<temp_bdf_type> temp_bdf_ptrtype;
    typedef Bdf<fluid_functionspace_type>  fluid_bdf_type;
    typedef boost::shared_ptr<fluid_bdf_type> fluid_bdf_ptrtype;

    /* export */
    typedef Exporter<mesh_type> export_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    OpusModel();
    OpusModel( po::variables_map const& vm );
    OpusModel( OpusModel const & );
    ~OpusModel();

    // initialize the opus model
    void init();

    //@}

    /** @name Operator overloads
    fa     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P );
    void run();

    void exportResults( double time, temp_element_type& T, fluid_element_type& U, bool force_export = false );



    //@}



protected:

private:

private:
    std::string mu_file;
    double M_dt;
    double s1,s2;

    double M_meshsize;
    bool M_force_rebuild;

    double M_thermal_conductance;

    mesh_ptrtype M_mesh;
    mesh_ptrtype M_mesh_air;
    mesh12_ptrtype M_mesh_line;
    mesh12_ptrtype M_mesh_cross_section_2;



    p0_space_ptrtype M_P0h;
    p0_element_ptrtype domains;
    p0_element_ptrtype rhoC;
    p0_element_ptrtype k;
    p0_element_ptrtype Q;

    p1_functionspace_ptrtype M_P1h;

    temp_bdf_ptrtype M_temp_bdf;
    temp_functionspace_ptrtype M_Th;
    grad_temp_functionspace_ptrtype M_grad_Th;
    thermal_operator_ptrtype M_thermal;

    fluid_bdf_ptrtype M_fluid_bdf;
    fluid_functionspace_ptrtype M_Xh;
    fluid_operator_ptrtype M_fluid;

    boost::shared_ptr<export_type> M_exporter,M_exporter_fluid;
};

/** \\@} */
}
#endif /* __OpusModel_H */
