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
#include <feel/feeldiscr/geometricspace.hpp>
#include <feel/feelcore/json.hpp>

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
    using mesh_type = typename space_type::mesh_type;


    SensorBase() = default;
    SensorBase( space_ptrtype const& space, std::string const& name = "sensor" , std::string const& t = ""):
        super_type( space ),
        M_name( name ),
        M_type( t )
    {}

    virtual ~SensorBase() = default;

    void setSpace( space_ptrtype const& space ) override
    {
        super_type::setSpace(space);
        this->init();
    }
    virtual void init() = 0;

    void setName( std::string const& n ) { M_name = n; }
    std::string const& name() const { return M_name; }

    std::string const& type() const { return M_type; }
    virtual json to_json() const
    {
        json j;
        j["type"] = M_type;
        return j;
    }

    //!
    //! other interface may include:
    //!  - DB storage for past measurements
    //!  - DB storage for future measurements
    //!
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & M_name;
        ar & M_type;
    }

protected:
    std::string M_name;
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
    using mesh_type = typename space_type::mesh_type;
    using node_t = typename mesh_type::node_type;
    using geometricspace_type = GeometricSpace<mesh_type>;
    using geometricspacectx_type = typename geometricspace_type::Context;
    using gmc_type = typename geometricspace_type::gmc_type;
    using gmc_ptrtype = std::shared_ptr<gmc_type>;
    using gm_type = typename gmc_type::gm_type;
    using geoelement_type = typename gmc_type::element_type;
    using fe_type = typename space_type::fe_type;

    SensorPointwise() = default;
    SensorPointwise( space_ptrtype const& space, node_t const& p, std::string const& n = "pointwise"):
        super_type( space, n, "pointwise" ),
        M_position(p)
    {
        this->init();
    }

    virtual ~SensorPointwise(){}

    void setPosition( node_t const& p ) { M_position = p; this->init(); }
    node_t const& position() const { return M_position; }

    void init() override
    {
        auto space = this->space();
        auto mesh = space->mesh();
        auto geospace = std::make_shared<geometricspace_type>(mesh);
        auto geospacectx = std::make_shared<geometricspacectx_type>( geospace );
        geospacectx->add(M_position, false);
        geospacectx->updateForUse();

        auto u = space->element();
        auto expr = id(u);
        geospacectx->template updateGmcContext<std::decay_t<decltype(expr)>::context>();

        using expr_type = std::decay_t<decltype(expr)>;
        using expr_basis_t = typename expr_type::test_basis;
        static constexpr size_type expr_context = expr_type::context;
        // fe context
        using fecontext_type = typename space_type::fe_type::template Context< expr_context, fe_type, gm_type, geoelement_type>;
        using fecontext_ptrtype = std::shared_ptr<fecontext_type>;

        for ( auto& [geoCtxId,geoCtx] : *geospacectx )
        {
            auto gmc = std::get<0>(geoCtx)->gmContext();
            auto const& curCtxIdToPointIds = std::get<1>(geoCtx);
            auto const& elt = gmc->element();
            auto fepc = space->fe()->preCompute( space->fe(), gmc->xRefs() );
            auto fec = std::make_shared<fecontext_type>( space->fe(), gmc, fepc );
            auto mapgmc = Feel::vf::mapgmc( gmc );
            auto mapfec = Feel::vf::mapfec( fec );
            auto tExpr = expr.evaluator(mapgmc, mapfec);

            // TODO check if next lines are really useful?
            fec->update( gmc );
            tExpr.update( mapgmc, mapfec );

            //Eigen::MatrixXd M_IhLoc = Vh->fe()->localInterpolants(1);
            Eigen::MatrixXd IhLoc = Eigen::MatrixXd::Zero(expr_basis_t::nLocalDof,gmc->nPoints());
            space->fe()->interpolateBasisFunction( tExpr, IhLoc );


            for( auto const& ldof : space->dof()->localDof( elt.id() ) )
            {
                index_type index = ldof.second.index();

                for ( uint16_type q=0;q<curCtxIdToPointIds.size();++q )
                {
                    this->containerPtr()->set( index, IhLoc( ldof.first.localDof(), q ) ); // TODO use datamap mapping
                }
            }
        }

        this->containerPtr()->close();
    }

    json to_json() const override
    {
        json j = super_type::to_json();
        j["coord"] = json::array();
        for(int i = 0; i < space_type::nDim; ++i)
            j["coord"].push_back(M_position(i));
        return j;
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<super_type>(*this);
        ar & M_position;
    }


    node_t M_position;
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
    using node_t = typename space_type::mesh_type::node_type;
    static const int nDim = space_type::mesh_type::nDim;

    SensorGaussian() = default;
    SensorGaussian( space_ptrtype const& space, node_t const& center,
                    double radius = 1., std::string const& n = "gaussian"):
        super_type( space, n, "gaussian" ),
        M_radius( radius ),
        M_position( center )
    {
        this->init();
    }

    virtual ~SensorGaussian(){}

    void setPosition( node_t const& p ) { M_position = p; this->init(); }
    node_t const& position() const { return M_position; }
    void setRadius( double radius ) { M_radius = radius; this->init(); }
    double radius() { return M_radius; }

    void init() override
    {
        auto v = this->space()->element();
        auto phi = this->phiExpr( mpl::int_< space_type::nDim >() );
        auto n = integrate(_range=elements(support(this->space())), _expr=phi).evaluate()(0,0);
        form1( _test=this->space(), _vector=this->containerPtr() ) =
            integrate(_range=elements(support(this->space())), _expr=id(v)*phi/n );
        this->close();
    }

    json to_json() const override
    {
        json j = super_type::to_json();
        j["radius"] = M_radius;
        j["coord"] = json::array();
        for(int i = 0; i < space_type::nDim; ++i)
            j["coord"].push_back(M_position(i));
        return j;
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

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<super_type>(*this);
        ar & M_position;
        ar & M_radius;
    }

    double M_radius;
    node_t M_position;
};

//!
//! surface type sensor
//!
template<typename Space>
class SensorSurface: public SensorBase<Space>
{
public:

    // -- TYPEDEFS --
    typedef SensorSurface<Space> this_type;
    typedef SensorBase<Space> super_type;

    using space_type = typename super_type::space_type;
    using space_ptrtype = typename super_type::space_ptrtype;
    using mesh_type = typename super_type::mesh_type;
    static const int nDim = space_type::mesh_type::nDim;

    SensorSurface() = default;
    SensorSurface( space_ptrtype const& space, std::set<std::string> const& markers, std::string const& n = "surface"):
        super_type( space, n, "surface" ),
        M_markers( markers )
    {
        this->init();
    }

    virtual ~SensorSurface(){}

    void setMarkers( std::set<std::string> const& markers ) { M_markers = markers; this->init(); }
    std::set<std::string> markers() { return M_markers; }

    void init() override
    {
        auto v = this->space()->element();
        auto n = integrate(_range=markedfaces(support(this->space()), M_markers), _expr=cst(1.)).evaluate()(0,0);
        form1( _test=this->space(), _vector=this->containerPtr() ) =
            integrate(_range=markedfaces(support(this->space()), M_markers), _expr=id(v)/n );
        this->close();
    }

    json to_json() const override
    {
        json j = super_type::to_json();
        j["markers"] = json::array();
        for(auto const& m : M_markers)
            j["markers"].push_back(m);
        return j;
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<super_type>(*this);
        ar & M_markers;
    }

    std::set<std::string> M_markers;
};

template<typename Space>
class SensorMap : public std::map<std::string, std::shared_ptr<SensorBase<Space>>>
{
public:
    using this_type = SensorMap<Space>;
    using super_type = std::map<std::string, std::shared_ptr<SensorBase<Space>>>;
    using value_type = typename super_type::value_type;
    using element_type = SensorBase<Space>;

    using space_type = Space;
    using space_ptrtype = std::shared_ptr<space_type>;
    using fselement_type = typename space_type::element_type;

    SensorMap() = default;
    explicit SensorMap( space_ptrtype const& Xh ): M_Xh(Xh) {}
    explicit SensorMap( space_ptrtype const& Xh, json j)
        : M_Xh(Xh)
    {
        this->loadJson(j);
    }

    void loadJson(json const& j)
    {
        M_j = j;
        for( auto const& /*it*/[name, sensor] : M_j.items() )
        {
            // auto name = it.first;
            // auto sensor = it.second;
            std::string type = sensor.value("type", "");

            if( type == "pointwise" )
            {
                using node_t = typename SensorPointwise<space_type>::node_t;
                node_t n(space_type::nDim);
                auto coords = sensor["coord"].template get<std::vector<double>>();
                for( int i = 0; i < space_type::nDim; ++i )
                    n(i) = coords[i];
                auto s = std::make_shared<SensorPointwise<space_type>>(M_Xh, n, name);
                this->insert(value_type(name, s));
            }
            else if( type == "gaussian" )
            {
                double radius = sensor.value("radius", 0.1);
                using node_t = typename SensorGaussian<space_type>::node_t;
                node_t n(space_type::nDim);
                auto coords = sensor["coord"].template get<std::vector<double>>();
                for( int i = 0; i < space_type::nDim; ++i )
                    n(i) = coords[i];
                auto s = std::make_shared<SensorGaussian<space_type>>(M_Xh, n, radius, name);
                this->insert(value_type(name, s));
            }
            else if( type == "surface" )
            {
                auto markers = sensor["markers"].template get<std::set<std::string>>();
                auto s = std::make_shared<SensorSurface<space_type>>(M_Xh, markers, name);
                this->insert(value_type(name, s));
            }
        }
    }

    Eigen::VectorXd apply(fselement_type const& u) const
    {
        Eigen::VectorXd v(this->size());
        int i = 0;
        for( auto const& it /*[name, sensor]*/ : *this )
            v(i++) = (*it.second)(u);
        return v;
    }

    json to_json()
    {
        if( M_j.empty() )
        {
            for( auto const& it /*[name, sensor]*/ : *this )
            {
                auto name = it.first;
                auto sensor = it.second;
                M_j[name] = sensor->to_json();
            }
        }
        return M_j;
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // add any subtype of SensorBase here
        ar.template register_type<SensorPointwise<space_type>>();
        ar.template register_type<SensorGaussian<space_type>>();
        ar.template register_type<SensorSurface<space_type>>();
        ar & boost::serialization::base_object<super_type>(*this);
        if( Archive::is_loading::value )
        {
            for( auto& it /*[name, sensor]*/ : *this)
                it.second->setSpace(M_Xh);
        }
    }

    space_ptrtype M_Xh;
    json M_j;
};

}
#endif /* _FEELPP_CRB_SENSORS_HPP */
