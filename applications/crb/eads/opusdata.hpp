/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-11-15

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
   \file opusdata.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-11-15
 */
#ifndef __OpusData_H
#define __OpusData_H 1

//#include <boost/property_tree/ptree.hpp>
//#include <boost/property_tree/xml_parser.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typetraits.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelalg/glas.hpp>

#include "opuscomponent.hpp"

namespace Feel
{
/**
 * \class OpusData
 * \brief Data for the OPUS EADS Benchmark
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class OpusData
{
public:


    /** @name Typedefs
     */
    //@{

    enum { NO_EXPORT = 0, EXPORT_SAME_MESH, EXPORT_LAGP1_MESH };

    enum { FLOW_POISEUILLE = 0, FLOW_NAVIERSTOKES };
    enum { INIT_FLOW_WITH_ZERO = 0, INIT_FLOW_WITH_STOKES };

    typedef Feel::node<double>::type node_type;

    typedef OpusData opusdata_type;
    typedef boost::shared_ptr<opusdata_type> opusdata_ptrtype;

    //@}

    /** @name Structures
     */
    //@{

    /**
     * \class Inflow
     * \brief Inflow functor
     */
    struct Inflow
    {
        typedef double value_type;
        typedef Feel::node<value_type>::type node_type;
        typedef Feel::uint16_type uint16_type;
        static const uint16_type rank = 1;

        Inflow( OpusData const& opusdata, double time );
        Inflow( Inflow const& i ): M_opusdata( i.M_opusdata ), M_time( i.M_time ) {}

        double operator()( uint16_type, uint16_type, node_type const& p, node_type const& n ) const;

    private:
        OpusData const& M_opusdata;
        double M_time;
    };

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * constructor
     * \param d dimension (2 or 3)
     */
    OpusData( int d );

    /**
     * constructor
     * \param d dimension (2 or 3)
     * \param vm variables map from the command line options
     */
    OpusData( int d, Feel::po::variables_map const& vm );

    //! copy constructor
    OpusData( OpusData const & );

    //! destructor
    ~OpusData();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    Feel::po::variables_map vm() const
    {
        return M_vm;
    }

    bool isSteady() const
    {
        return M_is_steady;
    }

    double d() const
    {
        return M_dimension;
    }
    double h() const
    {
        return M_h;
    }

    double orderTime() const
    {
        return M_order_time;
    }
    double orderTemp() const
    {
        return M_order_temp;
    }
    double orderU() const
    {
        return M_order_u;
    }
    double orderP() const
    {
        return M_order_p;
    }

    double gammaBc() const
    {
        return M_gamma_bc;
    }
    double gammaU() const
    {
        return M_gamma_u;
    }
    double gammaP() const
    {
        return M_gamma_p;
    }
    double gammaTemp() const
    {
        return M_gamma_conv_T;
    }
    double gammaDivDiv() const
    {
        return M_gamma_divdiv;
    }
    double deltaDivDiv() const
    {
        return M_delta_divdiv;
    }
    double epsPseudoCompressibility() const
    {
        return M_eps_compress;
    }

    OpusComponent const& component( std::string const& comp ) const
    {
        return M_components.find( comp )->second;
    }
    OpusComponent& component( std::string const& comp )
    {
        return M_components.find( comp )->second;
    }

    std::vector<std::string> const& dirichletTemperatureMarkers() const
    {
        return M_dirichlet_temp;
    }
    std::vector<std::string> const& dirichletVelocityMarkers() const
    {
        return M_dirichlet_velocity;
    }
    std::vector<std::string> const& dirichletPressureMarkers() const
    {
        return M_dirichlet_pressure;
    }

    /**
     * \return whether to use the same preconditioner
     */
    int useSamePreconditioner() const
    {
        return M_use_same_prec;
    }


    /**
     * \return how to initialize the Navier-Stokes solver
     */
    int init() const
    {
        return M_init;
    }

    /**
     * \return the export strategy
     * 0 no export, 1 export with same mesh, 2 export on finer mesh
     */
    int doExport() const
    {
        return M_export;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //@}

    /** @name  Methods
     */
    //@{

    gmsh_ptrtype createMesh( double h, bool ref = false );
    gmsh_ptrtype createMeshLine( double h );
    gmsh_ptrtype createMeshCrossSection2( double h );
    gmsh_ptrtype createMeshAir( double h );

    void createMatlabScript();

    Inflow inflow( double& time ) const
    {
        return Inflow( *this, time );
    }

    void print() const;

    static Feel::po::options_description makeOptions();

    //! load the property tree
    void load( const std::string &filename );

    //! save the property tree
    void save( const std::string &filename );

    //@}




protected:

    //! variables map from the command line
    Feel::po::variables_map M_vm;

    //! steady or unsteady case
    bool M_is_steady;

    //! geometric dimension of the problem
    int M_dimension;

    //! geometric order
    int M_order_g;

    //! time order
    int M_order_time;

    //! temperature order
    int M_order_temp;

    //! velocity order
    int M_order_u;

    //! pressure order
    int M_order_p;

    //! true if use equal order (velocity/pressure), false otherwise
    bool M_equal_order;

    //! h size
    double M_h;

    //! penalisation coefficient for weakly imposed BC
    double M_gamma_bc;

    //! parameter for equal order stabilization
    double M_gamma_p;

    //! parameter for convection stabilization
    double M_gamma_u;

    //! parameter for convection stabilization
    double M_gamma_conv_T;

    //! parameter for divergence jump stabilization
    double M_gamma_divdiv;

    //! parameter for divergence stabilization
    double M_delta_divdiv;

    // pseudo compressibility
    double M_eps_compress;

    //! preconditioner strategy
    int M_use_same_prec;

    //! initialisation strategy, 0: with zero, 1: with Stokes
    int M_init;

    //! export strategy : 0 no export, 1 export with (P1) pressure mesh, 2 export with (P1) velocity mesh
    int M_export;

    //! component dictionary
    std::map<std::string,OpusComponent> M_components;

    //! vector of Dirichlet markers for temperature
    std::vector<std::string> M_dirichlet_temp;

    //! vector of Dirichlet markers for velocity
    std::vector<std::string> M_dirichlet_velocity;

    //! vector of Dirichlet markers for pressure
    std::vector<std::string> M_dirichlet_pressure;


};
}
#endif /* __OpusData_H */

