/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 This file is part of the Feel library

 Copyright (C) 2010 University of Coimbra

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
 \file ale.hpp
 \author Goncalo Pena <gpena@mat.uc.pt>
 \date 2010-10-12
 */

#ifndef FEELPP_MODELS_ALE_H
#define FEELPP_MODELS_ALE_H 1

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelmodels/modelcore/modelbase.hpp>

namespace Feel
{
namespace FeelModels
{

template <class Convex, int Order = 1>
class ALE : public ModelBase
{
  public:
    typedef ModelBase super_type;

    typedef ALE<Convex, Order> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

    typedef Convex convex_type;
    static const uint16_type Dim = convex_type::nDim;
    static const uint16_type Order_low = convex_type::nOrder;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef std::map<std::string, std::vector<flag_type>> flagSet_type;

  protected:
    template <int i>
    struct MyReferenceFunctionSpace
    {
        typedef bases<Lagrange<i, Vectorial>> basis_type;
        typedef FunctionSpace<mesh_type, basis_type> type;
        typedef boost::shared_ptr<type> ptrtype;
        typedef typename type::element_type elt_type;
        typedef boost::shared_ptr<elt_type> elt_ptrtype;
    };

    static const int orderUse = mpl::if_<mpl::greater<mpl::int_<Order>, mpl::int_<Order_low>>,
                                         mpl::int_<Order>,
                                         mpl::int_<Order_low>>::type::value;

  public:
    typedef typename MyReferenceFunctionSpace<orderUse>::basis_type ale_map_basis_type;
    typedef typename MyReferenceFunctionSpace<orderUse>::type ale_map_functionspace_type;
    typedef typename MyReferenceFunctionSpace<orderUse>::ptrtype ale_map_functionspace_ptrtype;
    typedef typename MyReferenceFunctionSpace<orderUse>::elt_type ale_map_element_type;

    /**
     * constructor,copy,desctructor
     */
    ALE( mesh_ptrtype mesh, std::string prefix = "",
         WorldComm const& worldcomm = WorldComm(),
         bool moveGhostEltFromExtendedStencil = false );
    ALE( ALE const& tc );
    //~ALE();

    /**
     * static builder
     */
    static self_ptrtype build( mesh_ptrtype mesh, std::string prefix = "",
                               WorldComm const& worldcomm = WorldComm(),
                               bool moveGhostEltFromExtendedStencil = false );

    /**
     * Add the set of flags that mark the boundary
     */
    void addBoundaryFlags( flagSet_type flags );
    void addBoundaryFlags( std::string str, flag_type flag );
    void clearFlagSet();
    void setFlagSet( flagSet_type const& fl );
    flagSet_type const& flagSet() const;
    std::vector<flag_type> const& flagSet( std::string key ) const;
    flag_type flagSet( std::string key, int k ) const;

    virtual ale_map_functionspace_ptrtype const& functionSpace() const = 0;
    virtual ale_map_element_type const& displacement() const = 0;
    virtual ale_map_element_type const& map() const = 0;

    virtual void generateMap( ale_map_element_type const& dispOnBoundary,
                              ale_map_element_type const& oldDisp ) = 0;

  private:
    flagSet_type M_flagSet;
};
} // namespace FeelModels
} // namespace Feel
#endif // FEELPP_MODELS_ALE_H
