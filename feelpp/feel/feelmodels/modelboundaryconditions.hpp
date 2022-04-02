/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Mar 2015

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
#ifndef FEELPP_MODELBOUNDARYCONDITIONS_HPP
#define FEELPP_MODELBOUNDARYCONDITIONS_HPP 1

#include <boost/property_tree/ptree.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelmarkers.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel
{

namespace pt = boost::property_tree;

class FEELPP_EXPORT ModelBoundaryCondition : public CommObject
{
  public:
    using super = CommObject;
    explicit ModelBoundaryCondition( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ModelBoundaryCondition( ModelBoundaryCondition const& ) = default;
    ModelBoundaryCondition( ModelBoundaryCondition&& ) = default;
    ModelBoundaryCondition( pt::ptree const& p, std::string const& name,
                            std::string const& material, std::string const& e1,
                            std::string const& e2,
                            worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ModelBoundaryCondition& operator=( ModelBoundaryCondition const& ) = default;
    ModelBoundaryCondition& operator=( ModelBoundaryCondition && ) = default;

    std::string const& name() const { return M_name; }
    std::string const& marker() const { return *M_markers.begin(); }
    std::set<std::string> const& markers() const { return M_markers; }
    bool emptyMarkers() const { return M_markers.empty(); }
    std::string const& material() const { return M_material; }
    std::string const& expression() const { return M_expr1; }
    std::string const& expression1() const { return M_expr1; }
    std::string const& expression2() const { return M_expr2; }
    ModelExpression const& modelExpr() const { return M_modelExpr1; }
    ModelExpression const& modelExpr1() const { return M_modelExpr1; }
    ModelExpression const& modelExpr2() const { return M_modelExpr2; }

    void setName(std::string const& name) { M_name = name; }
    void setMarker(std::string const& marker) { M_markers = ModelMarkers(marker); }
    void setMarkers(std::set<std::string> const& markers) { M_markers = ModelMarkers(markers); }
    void setMaterial(std::string const& material) { M_material = material; }
    void setExpr(std::string const& expr, WorldComm const& worldComm = Environment::worldComm(), std::string const& directoryLibExpr = "" ) { M_expr1 = expr; M_modelExpr1.setExpr(expr, worldComm, directoryLibExpr); }
    void setExpr1(std::string const& expr, WorldComm const& worldComm = Environment::worldComm(), std::string const& directoryLibExpr = "" ) { M_expr1 = expr; M_modelExpr1.setExpr(expr, worldComm, directoryLibExpr); }
    void setExpr2(std::string const& expr, WorldComm const& worldComm = Environment::worldComm(), std::string const& directoryLibExpr = "" ) { M_expr2 = expr; M_modelExpr2.setExpr(expr, worldComm, directoryLibExpr); }

    bool hasExpression() const { return !M_expr1.empty(); }
    bool hasExpression1() const { return !M_expr1.empty(); }
    bool hasExpression2() const { return !M_expr2.empty(); }
    bool isExpressionScalar() const { return M_modelExpr1.isScalar(); }
    bool isExpressionVector() const { return M_modelExpr1.isVector(); }
    bool isExpressionMatrix() const { return M_modelExpr1.isMatrix(); }
    bool isExpression1Scalar() const { return M_modelExpr1.isScalar(); }
    bool isExpression1Vector() const { return M_modelExpr1.isVector(); }
    bool isExpression1Matrix() const { return M_modelExpr1.isMatrix(); }
    bool isExpression2Scalar() const { return M_modelExpr2.isScalar(); }
    bool isExpression2Vector() const { return M_modelExpr2.isVector(); }
    bool isExpression2Matrix() const { return M_modelExpr2.isMatrix(); }

    template<int M=1, int N=1>
    auto expr() const { return M_modelExpr1.expr<M,N>(); }
    template<int M=1, int N=1>
    auto expr1() const { return M_modelExpr1.expr<M,N>(); }
    template<int M=1, int N=1>
    auto expr2() const { return M_modelExpr2.expr<M,N>(); }

    void setParameterValues( std::map<std::string,double> const& mp );

protected:
    pt::ptree M_pt;
    std::string M_name;
    std::string M_material;
    ModelMarkers M_markers;
    std::string M_expr1;
    std::string M_expr2;
    ModelExpression M_modelExpr1;
    ModelExpression M_modelExpr2;
};

struct FEELPP_EXPORT ModelBoundaryId : public std::tuple<std::string,std::string,std::string>
{
  using parent = std::tuple<std::string,std::string,std::string>;
  ModelBoundaryId( parent const& t ) : parent( t ) {}
  std::string const& type() const { return std::get<1>( *this ); }
  // 
  std::string const& field() const { return std::get<0>( *this ); }
  // the name of the bc, can be also the name of the marker
  std::string const& name() const { return std::get<2>( *this ); }

};
inline std::ostream&
operator<<(std::ostream& os, ModelBoundaryId const& bcid )
{
  os << "(" << bcid.field() << "," << bcid.type() << "," << bcid.name() << ")";
  return os;
}
/**
 * a map of boundary conditions
 */
class FEELPP_EXPORT ModelBoundaryConditions : public std::map<std::string,std::map<std::string,std::map<std::string,ModelBoundaryCondition> > >, public CommObject
{
  public:
    using super = CommObject;
    explicit ModelBoundaryConditions( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    virtual ~ModelBoundaryConditions() = default;
    void setPTree( pt::ptree const& p );
    void setParameterValues( std::map<std::string,double> const& mp );
    // flatten the boundary condition map so that we can easily iterate in a one level loop
    std::map<ModelBoundaryId,ModelBoundaryCondition> flatten() const;
    std::map<std::string,std::map<std::string,ModelBoundaryCondition> > const& byField( std::string const& field ) const;
    std::map<std::string,ModelBoundaryCondition> const& byFieldType( std::string const& field, std::string const& type ) const;

  private:
    void setup();

    pt::ptree M_pt;
    std::map<std::string,std::map<std::string,ModelBoundaryCondition> > M_emptyField;
    std::map<std::string,ModelBoundaryCondition> M_emptyFieldType;
};


class FEELPP_EXPORT ModelBoundaryConditionsNEW// : public std::map<std::string,nl::json>
{
  public:
  ModelBoundaryConditionsNEW() = default;
  ModelBoundaryConditionsNEW( ModelBoundaryConditionsNEW && ) = default;
  ModelBoundaryConditionsNEW( ModelBoundaryConditionsNEW const& ) = default;

  bool hasSection( std::string const& sn ) const { return M_sections.find( sn ) != M_sections.end(); }
  nl::json const& section( std::string const& sn ) const { return M_sections.at( sn ); }

  void setup( nl::json const& jarg )
  {
      if ( !jarg.is_object() )
          return;
      for ( auto const& [jargkey,jargval] : jarg.items() )
      {
          M_sections.emplace( jargkey,jargval );
      }
  }
  private:
  std::map<std::string,nl::json> M_sections;
};

}

#endif
