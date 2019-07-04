//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @author Idrissa Niakh <>
//! @date 14 May 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef FEELPP_CRB_SENSORS_HPP
#define FEELPP_CRB_SENSORS_HPP 1

#include <limits>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>
#if BOOST_VERSION >= 104700
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#endif
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/base_object.hpp>

#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/parameterspace.hpp>

#include <feel/feelvf/vf.hpp>

#include <feel/feelcrb/crb.hpp>
#include <feel/feelcrb/crbmodel.hpp>
#include <feel/feeldiscr/reducedbasisspace.hpp>
#include <feel/feeldiscr/geometricspace.hpp>

#include <Eigen/Core>
#include <feel/feelcrb/sensordesc.hpp>
namespace Feel
{

//!
//! base class for sensors
//!
template<typename Space, typename ValueT = double>
class SensorBase: public FsFunctionalLinear<Space>
{
public:

    // -- TYPEDEFS --
    typedef SensorBase<Space> this_type;
    typedef FsFunctionalLinear<Space> super_type;

    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
    using node_t = typename space_type::mesh_type::node_type;

    SensorBase( space_ptrtype const& space, std::string const& n = "sensor" ):
        super_type(space),
        M_name( n )
    {}

    virtual ~SensorBase() = default;

    void setName( std::string const& n ) { M_name = n; }
    std::string const& name() const { return M_name; }

    //!
    //! other interface may include:
    //!  - DB storage for past measurements
    //!  - DB storage for future measurements
    //!
private:
    std::string M_name;

};

//!
//! pointwise type sensor
//!
template<typename Space>
class SensorPointwise: public SensorBase<Space>
{
public:

    // -- TYPEDEFS --
    typedef SensorPointwise<Space> this_type;
    typedef SensorBase<Space> super_type;

    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
    using node_t = typename space_type::mesh_type::node_type;
    
    SensorPointwise( space_ptrtype const& space, node_t const& p, std::string const& n = "pointwise"):
        super_type(space,n),
        M_point(p)
    {
        init();
    }

    virtual ~SensorPointwise(){}

    void setPoint( node_t const& p )
    {
        M_point = p;
    }

    node_t const& point() const
    {
        return M_point;
    }

    void init()
    {
        auto v=this->space()->element();
        //auto expr=integrate(_range=elements(M_space->mesh()), _expr=id(v)*phi);
        //super_type::operator=( expr );
        this->close();
    }
private:
    node_t M_point;
};


//!
//! gaussian type sensor
//!
template<typename Space>
class SensorGaussian: public SensorBase<Space>
{
public:

    // -- TYPEDEFS --
    typedef SensorGaussian<Space> this_type;
    typedef SensorBase<Space> super_type;

    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
    using node_t = typename space_type::mesh_type::node_type;

    SensorGaussian( space_ptrtype const& space, node_t const& center, double radius = 0., std::string const& n = "gaussian"):
        super_type(space,n),
        M_center(center),
        M_radius(radius)
    {
        init();
    }

    virtual ~SensorGaussian(){}

    void setCenter( node_t const& center )
    {
        M_center = center;
    }

    void setRadius( double radius )
    {
        M_radius = radius;
    }

    std::vector<double> const& center()
    {
        return M_center;
    }

    double radius()
    {
        return M_radius;
    }

    void init()
    {
        auto v = this->space()->element();
        auto phi = this->phiExpr( mpl::int_< space_type::nDim >() );
        if constexpr ( space_type::nDim == 1 )
            form1( _test=this->space(), _vector=this->containerPtr() ) = integrate(_range=elements(this->space()->mesh()), _expr=id(v)*exp( pow(Px()-M_center[0],2)/(2*std::pow(M_radius,2))) );
        if constexpr ( space_type::nDim == 2 )
            form1( _test=this->space(), _vector=this->containerPtr() ) = integrate(_range=elements(this->space()->mesh()), _expr=id(v)*exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2))/(2*std::pow(M_radius,2))) );
        if constexpr ( space_type::nDim == 3 )
            form1( _test=this->space(), _vector=this->containerPtr() ) = integrate(_range=elements(this->space()->mesh()), _expr=id(v)*exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2)+pow(Pz()-M_center[2],2))/(2*std::pow(M_radius,2))) );
        //super_type::operator= ( expr );
        this->close();
    }

private:
    auto phiExpr( mpl::int_<1> /**/ )
        {
         return exp( pow(Px()-M_center[0],2)/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<2> /**/ )
    {
         return exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2))/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<3> /**/ )
    {
         return exp(( pow(Px()-M_center[0],2)+pow(Py()-M_center[1],2)+pow(Pz()-M_center[2],2))/(2*std::pow(M_radius,2)));
    }
    node_t M_center;
    double M_radius;
};

template<typename Space>
class SensorMap : public std::map<std::string,std::shared_ptr<SensorBase<Space>>>
{

public:

    // -- TYPEDEFS --
    typedef SensorMap<Space> this_type;
    typedef SensorBase<Space> element_type;

    typedef Space space_type;

    typedef std::shared_ptr<space_type> space_ptrtype;

    static const int nDim = space_type::mesh_type::nDim;
    
    using sensor_type = SensorBase<Space>;
    using node_t = typename sensor_type::node_t;
    using sensor_ptrtype = std::shared_ptr<sensor_type>;

    SensorMap( space_ptrtype const& space, SensorDescriptionMap<nDim>  const& sensor_desc):
        M_sensor_desc(sensor_desc),
        M_space(space)
    {
        this->init();
    }

    space_ptrtype const&  space()
    {
        return M_space;
    }

    SensorDescriptionMap<nDim> const& sensor_desc()
    {
        return M_sensor_desc;
    }

    void setSpace( space_ptrtype const& space )
    {
        M_space = space;
    }

    void setSensor_desc( SensorDescriptionMap<nDim> const& sensor_desc)
    {
        M_sensor_desc = sensor_desc;
    }


    virtual ~SensorMap() {}

    int size() const
    {
        return M_sensor_desc.size();
    }

    void init()
    {
        Feel::cout << "init " << M_sensor_desc.size() << " sensors" << std::endl;
        //this->reserve( this->size() );
        for( auto const& [sensor_name, sensor_desc]  : M_sensor_desc )
        {
            Feel::cout << "sensor " << sensor_name << std::endl;
            sensor_ptrtype newElement;
            node_t n( sensor_desc.position().size() );
            for ( int i = 0; i < sensor_desc.position().size(); ++i )
                n( i ) = sensor_desc.position()[i];
            if (sensor_desc.type() == "gaussian" )
            {
                newElement = std::make_shared<SensorGaussian<space_type>> (M_space, n, sensor_desc.radius(),sensor_desc.name());
            }
            else
            {
               if (sensor_desc.type() == "pointwise" )
               {
                   newElement = std::make_shared<SensorPointwise<space_type>> (M_space, n,sensor_desc.name());
               }
            }
            this->insert( std::pair( sensor_name, newElement ) );
            Feel::cout << "Sensor " << sensor_name << " added to map" << std::endl;
        }
    }

private:
    
    SensorDescriptionMap<nDim> M_sensor_desc;
    space_ptrtype M_space;

};


}
#endif /* _FEELPP_CRB_SENSORS_HPP */
