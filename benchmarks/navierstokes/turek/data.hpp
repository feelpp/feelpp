/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-05-02

  Copyright (C) 2008,2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file data.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-05-02
 */
#ifndef __Data_H
#define __Data_H 1


#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/typetraits.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/application.hpp>


/**
 * \class Data
 * \brief Data for the Turek benchmark
 *
 *  @author Christophe Prud'homme
 *  @see
 */
class Data
{
public:


    /** @name Typedefs
     */
    //@{

    enum { CROSS_SECTION_CIRCULAR = 0, CROSS_SECTION_SQUARE };

    enum { NO_EXPORT = 0, EXPORT_SAME_MESH, EXPORT_LAGP1_MESH };

    enum { INIT_WITH_ZERO = 0, INIT_WITH_STOKES };

    enum { INFLOW_STEADY = 0, INFLOW_UNSTEADY };

    typedef Feel::node<double>::type node_type;

    typedef Data data_type;
    typedef boost::shared_ptr<data_type> data_ptrtype;

    //@}

    /** @name Structures
     */
    //@{

    static data_ptrtype New( Feel::po::variables_map const& vm );

#if 1
    struct Inflow
    {
        typedef double value_type;
        typedef Feel::node<value_type>::type node_type;
        typedef Feel::uint16_type uint16_type;
        static const uint16_type rank = 1;

        Inflow( Data const& data, double time );
        Inflow( Inflow const& i ): M_data( i.M_data ), M_time( i.M_time ) {}

        double operator()( uint16_type, uint16_type, node_type const& p, node_type const& n ) const;

    private:
        Data const& M_data;
        double M_time;
    };
#else
    template<int D>
    struct Inflow;

    template<>
    struct Inflow<2>
    {
        typedef double value_type;
        typedef Feel::node<value_type>::type node_type;
        typedef Feel::uint16_type uint16_type;
        static const uint16_type rank = 1;

        double operator()( uint16_type, uint16_type, node_type const& p, node_type const& n ) const;
    };

    template<>
    struct Inflow<3>
    {
        typedef double value_type;
        typedef Feel::node<value_type>::type node_type;
        typedef Feel::uint16_type uint16_type;
        static const uint16_type rank = 1;

        double operator()( uint16_type, uint16_type, node_type const& p, node_type const& n ) const;
    };
#endif

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * constructor
     * \param d dimension (2 or 3)
     */
    Data( int d );

    /**
     * constructor
     * \param d dimension (2 or 3)
     * \param vm variables map from the command line options
     */
    Data( int d, Feel::po::variables_map const& vm );

    //! copy constructor
    Data( Data const & );

    //! destructor
    ~Data();

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

    double d() const
    {
        return M_dimension;
    }
    double h() const
    {
        return M_h;
    }

    //! returns the scale by which h() is divided on the cylinder
    double hCylinderScale() const
    {
        return M_hcyl_scale;
    }

    double orderG() const
    {
        return M_order_g;
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

    double H() const
    {
        return M_H;
    }
    double D() const
    {
        return M_D;
    }
    double Re() const
    {
        return M_Re;
    }
    double rho() const
    {
        return M_rho;
    }
    double nu() const
    {
        return M_nu;
    }
    double mu() const
    {
        return M_rho*M_nu;
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
     * coordinates of the front point where the pressure difference is
     * calculated
     *
     * @return coordinates of the front point
     */
    node_type xa() const;

    /**
     * coordinates of the back point where the pressure difference is
     * calculated
     *
     * @return coordinates of the back point
     */
    node_type xe() const;

    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    double Ubar() const
    {
        return M_nu*M_Re/M_D;
    }

    /**
     * get the magnitude of the profile velocity for the 2D cylinder
     * 3*Ubar/2 in 2D and 9*Ubar/4 in 3D
     *
     * @return the magnitude of the profile velocity
     */
    double Um() const
    {
        return Feel::math::pow( 3./2., M_dimension-1 )*Ubar();
    }

    /**
     * \return the scalar coefficient for the dimensionalise the
     * adminensional force
     */
    double scalingForce() const;

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

    //! return the type of inflow : steady or unsteady
    int inflowType() const
    {
        return M_inflow_type;
    }

    //! \return the cross section type : circular or square
    int crossSectionType() const
    {
        return M_cross_section_type;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    std::pair<std::string,std::string> createCylinder();

    void createMatlabScript();

    Inflow inflow( double& time ) const
    {
        return Inflow( *this, time );
    }

    void print() const;

    static Feel::AboutData makeAbout();
    static Feel::po::options_description makeOptions();

    virtual void run() = 0;
    //@}




protected:

    //! variables map from the command line
    Feel::po::variables_map M_vm;

    //! geometric dimension of the problem
    int M_dimension;

    //! geometric order
    int M_order_g;

    //! velocity order
    int M_order_u;

    //! pressure order
    int M_order_p;

    //! true if use equal order (velocity/pressure), false otherwise
    bool M_equal_order;

    //! h size
    double M_h;

    //! h size scale around the obstable
    double M_hcyl_scale;

    //! penalisation coefficient for weakly imposed BC
    double M_gamma_bc;

    //! parameter for equal order stabilization
    double M_gamma_p;

    //! parameter for convection stabilization
    double M_gamma_u;

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

    //! height and width of the domain (in m)
    double M_H;

    //! diameter of the cylinder (in m)
    double M_D;

    //! Reynolds number
    double M_Re;

    //! viscosity is constant (in m^2/s)
    double M_nu;

    //! density is constant (in kg/m^3)
    double M_rho;

    //! inflow type
    int M_inflow_type;

    //! cross section type
    int M_cross_section_type;

    //! vector of Dirichlet markers for velocity
    std::vector<std::string> M_dirichlet_velocity;

    //! vector of Dirichlet markers for pressure
    std::vector<std::string> M_dirichlet_pressure;
};
#endif /* __Data_H */
