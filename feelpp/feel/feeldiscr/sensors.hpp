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

#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/fsfunctionallinear.hpp>

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
    using mesh_type = typename space_type::mesh_type;


    SensorBase( space_ptrtype const& space, node_t const& p, std::string const& name = "sensor" , std::string const& t = ""):
        super_type( space ),
        M_name( name ),
        M_position( p ),
        M_type( t )
    {}

    virtual ~SensorBase() = default;

    virtual void init() = 0;

    void setName( std::string const& n ) { M_name = n; }
    std::string const& name() const { return M_name; }

    void setPosition( node_t const& p ) { M_position = p; this->init(); }
    node_t const& position() const { return M_position; }

    std::string const& type() const { return M_type; }

    //!
    //! other interface may include:
    //!  - DB storage for past measurements
    //!  - DB storage for future measurements
    //!
protected:
    std::string M_name;
    node_t M_position;
    std::string M_type;
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
    using node_t = typename super_type::node_t;

    SensorPointwise( space_ptrtype const& space, node_t const& p, std::string const& n = "pointwise"):
        super_type( space, p, n, "pointwise" )
    {
        this->init();
    }

    virtual ~SensorPointwise(){}

    void init() override
    {
        // auto v=this->space()->element();
        //auto expr=integrate(_range=elements(M_space->mesh()), _expr=id(v)*phi);
        //super_type::operator=( expr );
        // this->close();
    }
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
    using mesh_type = typename super_type::mesh_type;
    using node_t = typename super_type::node_t;
    static const int nDim = space_type::mesh_type::nDim;

    SensorGaussian( space_ptrtype const& space, node_t const& center,
                    double radius = 1., std::string const& n = "gaussian"):
        super_type( space, center, n, "gaussian" ),
        M_radius( radius )
    {
        this->init();
    }

    virtual ~SensorGaussian(){}

    void setRadius( double radius ) { M_radius = radius; this->init(); }
    double radius() { return M_radius; }

    void init() override
    {
        auto v = this->space()->element();
        auto phi = this->phiExpr( mpl::int_< space_type::nDim >() );
        auto n = integrate(_range=elements(this->space()->mesh()), _expr=phi).evaluate()(0,0);
        form1( _test=this->space(), _vector=this->containerPtr() ) =
            integrate(_range=elements(this->space()->mesh()), _expr=id(v)*phi/n );
        this->close();
    }

private:
    auto phiExpr( mpl::int_<1> /**/ )
    {
        return exp( -inner(P()-vec(cst(this->M_position[0])))/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<2> /**/ )
    {
        return exp( -inner(P()-vec(cst(this->M_position[0]),cst(this->M_position[1])))/(2*std::pow(M_radius,2)));
    }
    auto phiExpr( mpl::int_<3> /**/ )
    {
        return exp( -inner(P()-vec(cst(this->M_position[0]),cst(this->M_position[1]),cst(this->M_position[2])))/(2*std::pow(M_radius,2)));
    }

    double M_radius;
};

#if 0
template<typename Space>
class SensorMap : public std::map<std::string,std::shared_ptr<SensorBase<Space>>>
{

public:

    // -- TYPEDEFS --
    typedef SensorMap<Space> this_type;
    typedef SensorBase<Space> element_type;

    typedef Space space_type;

    typedef std::shared_ptr<space_type> space_ptrtype;
    using mesh_type = typename space_type::mesh_type;
    using p0_space_ptrtype = Pdh_ptrtype<mesh_type,0>;

    static const int nDim = space_type::mesh_type::nDim;
    
    using sensor_type = SensorBase<Space>;
    using node_t = typename sensor_type::node_t;
    using sensor_ptrtype = std::shared_ptr<sensor_type>;

    SensorMap() = default;
    SensorMap( space_ptrtype const& space, SensorDescriptionMap<nDim>  const& sensor_desc):
        M_sensor_desc(sensor_desc),
        M_space(space),
        M_space_p0( Pdh<0>( space->mesh() ) )
    {
        this->init();
    }

    space_ptrtype const&  space()
    {
        return M_space;
    }
    p0_space_ptrtype const&  spaceP0()
    {
        return M_space_p0;
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
                newElement = std::make_shared<SensorGaussian<space_type>> (M_space, M_space_p0, n, sensor_desc.radius(),sensor_desc.name());
            }
            else
            {
               if (sensor_desc.type() == "pointwise" )
               {
                   newElement = std::make_shared<SensorPointwise<space_type>> (M_space, M_space_p0, n,sensor_desc.name());
               }
            }
            this->insert( std::pair( sensor_name, newElement ) );
            Feel::cout << "Sensor " << sensor_name << " added to map" << std::endl;
        }
    }

private:
    
    SensorDescriptionMap<nDim> M_sensor_desc;
    space_ptrtype M_space;
    p0_space_ptrtype M_space_p0;
};
#endif

}
#endif /* _FEELPP_CRB_SENSORS_HPP */
