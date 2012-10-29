/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-07-20

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
   \file opuscomponent.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-07-20
 */
#ifndef __OpusComponent_H
#define __OpusComponent_H 1

//#include <boost/property_tree/ptree.hpp>

namespace Feel
{
/**
 * \class OpusComponent
 * \brief describes a component of the Opus EADS Benchmark
 *
 *  @author Christophe Prud'homme
 */
class OpusComponent
{
public:
    OpusComponent()
        :
        M_name(),
        M_k( ),
        M_rhoC( ),
        M_Q( ),
        M_h ( ),
        M_e( ),
        M_flow_rate( 0 )
    {
    }

    OpusComponent( std::string name, Feel::po::variables_map const& vm )
        :
        M_name( name ),
        M_k( ),
        M_rhoC( ),
        M_Q( ),
        M_h ( ),
        M_e( ),
        M_flow_rate( 0 )
    {
        std::ostringstream ostr;
        ostr << name;
        M_k = vm[ostr.str()+".k"].as<double>();
        M_rhoC = vm[ostr.str()+".rhoC"].as<double>();
        M_Q = vm[ostr.str()+".Q"].as<double>();
        M_h = vm[ostr.str()+".h"].as<double>();
        M_e = vm[ostr.str()+".e"].as<double>();

        if ( vm.find( ostr.str()+".flow-rate" ) != vm.end() )
            M_flow_rate = vm[ostr.str()+".flow-rate"].as<double>();
    }
    OpusComponent( std::string name, double k, double rhoC, double Q, double h, double e, double fr = 0. )
        :
        M_name( name ),
        M_k( k ),
        M_rhoC( rhoC ),
        M_Q( Q ),
        M_h ( h ),
        M_e( e ),
        M_flow_rate( fr )
    {}

    OpusComponent( OpusComponent const& oc )
        :
        M_name( oc.M_name ),
        M_k( oc.M_k ),
        M_rhoC( oc.M_rhoC ),
        M_Q( oc.M_Q ),
        M_h( oc.M_h ),
        M_e( oc.M_e ),
        M_flow_rate( oc.M_flow_rate )
    {

    }
    ~OpusComponent()
    {}

    OpusComponent& operator=( OpusComponent const& oc )
    {
        if ( this != &oc )
        {
            M_name = oc.M_name;
            M_k = oc.M_k;
            M_rhoC = oc.M_rhoC;
            M_Q = oc.M_Q;
            M_h = oc.M_h;
            M_e = oc.M_e;
            M_flow_rate = oc.M_flow_rate;
        }

        return *this;
    }
    //! \return name of the component
    std::string name() const
    {
        return M_name;
    }

    //! \return conductivity of the component
    double k() const
    {
        return M_k;
    }
    //! set the conductivity of the component
    void setK( double k )
    {
        M_k = k;
    }

    //! \return heat capacity of the component
    double rhoC() const
    {
        return M_rhoC;
    }
    //! set heat capacity of the component
    void setRhoC( double r )
    {
        M_rhoC = r;
    }

    //! \return heat dissipated by the component
    double Q() const
    {
        return M_Q;
    }
    //! set heat dissipated by the component
    void setQ( double q )
    {
        M_Q = q;
    }

    //! \return height of the component
    double h() const
    {
        return M_h;
    }
    //! set height of the component
    void setH( double h )
    {
        M_h=h;
    }

    //! \return width of the component
    double e() const
    {
        return M_e;
    }
    //! set width of the component
    void setE( double e )
    {
        M_e = e;
    }

    //! \return the flow rate
    double flowRate() const
    {
        return M_flow_rate;
    }

    //! set the flow rate
    void setFlowRate( double f )
    {
        M_flow_rate = f;
    }

    // load the property tree in the structure
    //void load( const boost::property_tree::ptree& pt );

    //! save the property tree
    //void save( boost::property_tree::ptree& pt );

private:

    std::string M_name;
    double M_k;
    double M_rhoC;
    double M_Q;
    double M_h;
    double M_e;
    double M_flow_rate;
};
Feel::po::options_description makeComponentOptions();

std::ostream& operator<<( std::ostream& os, OpusComponent const& );
}
#endif /* __OpusComponent_H */
