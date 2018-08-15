/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 20 Aug 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_CONVOLVE_HPP
#define FEELPP_CONVOLVE_HPP 1

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>

namespace Feel {

namespace detail
{
template<typename Args>
struct convolve_type
{
    using space_type = decay_type<typename parameter::binding<Args,tag::space>::type>;
    using element_type =  typename space_type::element_type;
};
} // detail

BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::convolve_type<Args>::element_type ), // return type
    convolve,    // 2. function name

    tag,           // 3. namespace of tag types

    ( required
      ( range, *  )
      ( expr,   * )
      ( space, *( boost::is_convertible<mpl::_,std::shared_ptr<FunctionSpaceBase> > ) )
    ) // 4. one required parameter, and

    ( optional
      ( quad,   *, typename vf::detail::integrate_type<Args>::_quad_type() )
      ( geomap, *, (vf::detail::integrate_type<Args>::geoOrder > 1 )?GeomapStrategyType::GEOMAP_OPT:GeomapStrategyType::GEOMAP_HO )
      ( quad1,   *, typename vf::detail::integrate_type<Args>::_quad1_type() )
      ( use_tbb,   ( bool ), false )
      ( use_harts,   ( bool ), false )
      ( grainsize,   ( int ), 100 )
      ( partitioner,   *, "auto" )
      ( verbose,   ( bool ), false )
      ( quadptloc, *, typename vf::detail::integrate_type<Args>::_quadptloc_ptrtype() )
    )
)
{
    using T = double;
    constexpr int nDim = decay_type<space_type>::nDim;
    std::vector<Eigen::Matrix<T, nDim,1>> X,Z;
    Eigen::Matrix<T, nDim,1> Y;

    LOG(INFO) << "convolve::dof:" << space->dof()->dofPoints().size();
    for( auto const& d :  space->dof()->dofPoints() )
    {
        auto const& y = d.second.template get<0>();
        for ( int i = 0; i < nDim; ++i ) Y(i)=y(i);
        X.push_back( Y );
    }
    LOG(INFO) << "convolve::X(" << X.size() << ")=" << X;
    auto const& wc = worldComm(range);
    mpi::all_gather( wc.localComm(), X.data(), X.size(), Z  );
    LOG(INFO) << "convolve::Z(" << Z.size() << ")=" << Z;
    auto ret = integrate( _range=range, _expr=expr, _quad=quad, _geomap=geomap,
                          _quad1=quad1,_use_tbb=use_tbb,_grainsize=grainsize,
                          _partitioner=partitioner,_verbose=verbose,_quadptloc=quadptloc).evaluate( Z );
    LOG(INFO) << "convolve::ret(" << ret.size() <<")=" << ret << std::endl;
    auto v = space->element();
    for( int i = 0; i < space->nLocalDof(); ++i )
        v(i) = ret[space->dof()->firstDof()+i]( 0, 0);
    v.close();
    return v;
}

}

#endif
