/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 11 Apr 2015
 
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
#ifndef FEELPP_MODELPOSTPROCESS_HPP
#define FEELPP_MODELPOSTPROCESS_HPP 1


#include <map>

#include <boost/property_tree/ptree.hpp>

#include <feel/feelmodels/modelexpression.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

class FEELPP_EXPORT ModelPostprocessExports
{
  public :
    ModelPostprocessExports()
        :
        M_format( "ensightgold" )
    {}
    ModelPostprocessExports( ModelPostprocessExports const& ) = default;
    ModelPostprocessExports( ModelPostprocessExports&& ) = default;
    ModelPostprocessExports& operator=( ModelPostprocessExports const& ) = default;
    ModelPostprocessExports& operator=( ModelPostprocessExports && ) = default;

    std::set<std::string> const& fields() const { return M_fields; }
    std::string const& format() const { return M_format; }

    void setup( pt::ptree const& p );

  private :
    std::set<std::string> M_fields;
    std::string M_format;
};

class FEELPP_EXPORT ModelPointPosition
{
public :
    typedef Eigen::MatrixXd coord_value_type;
    typedef vector_field_expression<3> coord_expr_type;

    ModelPointPosition()
        :
        M_value( coord_value_type::Zero(3,1) )
        {}

    ModelPointPosition( ModelPointPosition const& ) = default;
    ModelPointPosition( ModelPointPosition&& ) = default;
    ModelPointPosition& operator=( ModelPointPosition const& ) = default;
    ModelPointPosition& operator=( ModelPointPosition && ) = default;

    std::string const& name() const { return M_name; }
    coord_value_type const& value() const { return M_value; }
    std::string const& meshMarker() const { return M_meshMarker; }
    coord_expr_type const& expression() const { CHECK( this->hasExpression() ) << "no expression defined"; return *M_expr; }
    bool hasExpression() const { return M_expr.get_ptr() != 0; }

    void setName( std::string const& s ) { M_name = s; }
    void setValue( coord_value_type const& v ) { M_value = v; }
    void setMeshMarker( std::string const& s ) { M_meshMarker = s; }
    void setExpression( std::string const& expression, std::string const& dirLibExpr = "",
                        WorldComm const& world = Environment::worldComm() )
        {
            M_expr = expr<3,1>( expression,"",world,dirLibExpr );
            M_value = M_expr->evaluate();
        }

    void setParameterValues( std::map<std::string,double> const& mp )
        {
            if ( !this->hasExpression() )
                return;
            M_expr->setParameterValues( mp );
            M_value = M_expr->evaluate();
        }

private:
    std::string M_name;
    coord_value_type M_value;
    boost::optional<coord_expr_type> M_expr;
    std::string M_meshMarker;

};

class FEELPP_EXPORT ModelPostprocessPointPosition : public std::pair< ModelPointPosition, std::set<std::string> >
{
    typedef std::pair< ModelPointPosition, std::set<std::string> > super_type;
public :
    ModelPostprocessPointPosition( WorldComm const& world = Environment::worldComm() )
        :
        super_type(),
        M_worldComm( world )
        {}
    ModelPostprocessPointPosition( ModelPointPosition const& ptPos, WorldComm const& world = Environment::worldComm() )
        :
        super_type( ptPos,std::set<std::string>() ),
        M_worldComm( world )
        {}
    ModelPostprocessPointPosition( ModelPointPosition const& ptPos, std::set<std::string> const& fields,
                                   WorldComm const& world = Environment::worldComm() )
        :
        super_type( ptPos,fields ),
        M_worldComm( world )
        {}
    ModelPostprocessPointPosition( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition( ModelPostprocessPointPosition&& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition && ) = default;

    ModelPointPosition const& pointPosition() const { return this->first; }
    ModelPointPosition & pointPosition() { return this->first; }
    std::set<std::string> const& fields() const { return this->second; }
    std::set<std::string> & fields() { return this->second; }

    void setPTree( pt::ptree const& _p, std::string const& name ) { M_p = _p; this->setup( name ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { this->second = fields; }
    void addFields( std::string const& field ) { this->second.insert( field ); }
    void setParameterValues( std::map<std::string,double> const& mp ) { this->pointPosition().setParameterValues( mp ); }
private:
    void setup( std::string const& name );
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};

class FEELPP_EXPORT ModelExtremum
{
public :
    ModelExtremum() = default;
    ModelExtremum( ModelExtremum const& ) = default;
    ModelExtremum( ModelExtremum&& ) = default;
    ModelExtremum& operator=( ModelExtremum const& ) = default;
    ModelExtremum& operator=( ModelExtremum && ) = default;

    std::string const& name() const { return M_name; }
    std::string const& type() const { return M_type; }
    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }

    void setName( std::string const& s ) { M_name = s; }
    void setType( std::string const& s ) { CHECK(s == "max" || s=="min" || s=="mean" ) << "invalid type " << s; M_type = s; }
    void addMarker( std::string const& mark )
        {
            if ( std::find( M_meshMarkers.begin(),M_meshMarkers.end(), mark ) == M_meshMarkers.end() )
                M_meshMarkers.push_back( mark );
        }
private :
    std::string M_name;
    std::string M_type;
    std::list<std::string> M_meshMarkers;
};
class FEELPP_EXPORT ModelPostprocessExtremum : public std::pair< ModelExtremum, std::set<std::string> >
{
public :
    ModelPostprocessExtremum( WorldComm const& world = Environment::worldComm() )
        :
        M_worldComm( world )
        {}
    ModelPostprocessExtremum( ModelPostprocessExtremum const& ) = default;
    ModelPostprocessExtremum( ModelPostprocessExtremum&& ) = default;
    ModelPostprocessExtremum& operator=( ModelPostprocessExtremum const& ) = default;
    ModelPostprocessExtremum& operator=( ModelPostprocessExtremum && ) = default;

    ModelExtremum const& extremum() const { return this->first; }
    ModelExtremum & extremum() { return this->first; }
    std::set<std::string> const& fields() const { return this->second; }
    std::set<std::string> & fields() { return this->second; }

    void setPTree( pt::ptree const& _p, std::string const& name ) { M_p = _p; this->setup( name ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { this->second = fields; }
    void addFields( std::string const& field ) { this->second.insert( field ); }
    //void setParameterValues( std::map<std::string,double> const& mp ) { this->extremum().setParameterValues( mp ); }

private:
    void setup( std::string const& name );
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
};

//! store informations require by the postprocessing Norm
class FEELPP_EXPORT ModelPostprocessNorm
{
public :
    ModelPostprocessNorm( WorldComm const& world = Environment::worldComm() )
        :
        M_worldComm( world ),
        M_quadOrder( 5 ),
        M_quad1Order( 5 )
        {}
    ModelPostprocessNorm( ModelPostprocessNorm const& ) = default;
    ModelPostprocessNorm( ModelPostprocessNorm&& ) = default;
    ModelPostprocessNorm& operator=( ModelPostprocessNorm const& ) = default;
    ModelPostprocessNorm& operator=( ModelPostprocessNorm && ) = default;

    //! name given to a Norm postprocessing
    std::string const& name() const { return M_name; }
    //! set of norm types (L2,H1,...) that the postprocessing compute
    std::set<std::string> const& types() const { return M_types; }
    //! name of fe field for which norms is computed (not require if expr is given)
    std::string const& field() const { return M_field; }
    //! mesh markers where norms is computed
    std::set<std::string> const& markers() const { return M_markers; }
    //! an expression represented an analytical solution
    ModelExpression const& solution() const { return M_solution; }
    //! expression of the gradient of the solution
    ModelExpression const& gradSolution() const { return M_gradSolution; }
    //! an expression for which norms is computed (not require if field is given)
    ModelExpression const& expr() const { return M_expr; }
    //! expression of the gradient of the expression
    ModelExpression const& gradExpr() const { return M_gradExpr; }
    //! quad order used in norms computations (not use if the quadrature can be exact)
    uint16_type quadOrder() const { return M_quadOrder; }
    //! quad1 order used with ho geometry and the optimized geomap
    uint16_type quad1Order() const { return M_quad1Order; }

    //! return true if an expression has been given
    bool hasExpr() const { return M_expr.hasAtLeastOneExpr(); }
    //! return true if a field has been given
    bool hasField() const { return !M_field.empty(); }

    void setPTree( pt::ptree const& _p, std::string const& name ) { M_p = _p; this->setup( name ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup( std::string const& name );
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types, M_markers;
    ModelExpression M_solution, M_gradSolution, M_expr, M_gradExpr;
    uint16_type M_quadOrder, M_quad1Order;
};

//! store informations require by the postprocessing Statistics
class FEELPP_EXPORT ModelPostprocessStatistics
{
public :
    ModelPostprocessStatistics( WorldComm const& world = Environment::worldComm() )
        :
        M_worldComm( world ),
        M_quadOrder( 5 ),
        M_quad1Order( 5 )
        {}
    ModelPostprocessStatistics( ModelPostprocessStatistics const& ) = default;
    ModelPostprocessStatistics( ModelPostprocessStatistics&& ) = default;
    ModelPostprocessStatistics& operator=( ModelPostprocessStatistics const& ) = default;
    ModelPostprocessStatistics& operator=( ModelPostprocessStatistics && ) = default;

    //! name given to a Statistics postprocessing
    std::string const& name() const { return M_name; }
    //! set of Statistics types (min,max,...) that the postprocessing compute
    std::set<std::string> const& types() const { return M_types; }
    //! name of fe field for which Statistics is computed (not require if expr is given)
    std::string const& field() const { return M_field; }
    //! mesh markers where Statistics is computed
    std::set<std::string> const& markers() const { return M_markers; }
    //! an expression for which Statistics is computed (not require if field is given)
    ModelExpression const& expr() const { return M_expr; }
    //! quad order used in Statistics computations (not use if the quadrature can be exact)
    uint16_type quadOrder() const { return M_quadOrder; }
    //! quad1 order used with ho geometry and the optimized geomap
    uint16_type quad1Order() const { return M_quad1Order; }

    //! return true if an expression has been given
    bool hasExpr() const { return M_expr.hasAtLeastOneExpr(); }
    //! return true if a field has been given
    bool hasField() const { return !M_field.empty(); }

    void setPTree( pt::ptree const& _p, std::string const& name ) { M_p = _p; this->setup( name ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup( std::string const& name );
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types, M_markers;
    ModelExpression M_expr;
    uint16_type M_quadOrder, M_quad1Order;
};

class FEELPP_EXPORT ModelPostprocess
{
public:

    ModelPostprocess( WorldComm const& world = Environment::worldComm() );
    ModelPostprocess( pt::ptree const& p, WorldComm const& world = Environment::worldComm() );
    virtual ~ModelPostprocess();
    pt::ptree const& pTree() const { return M_p; }
    pt::ptree & pTree() { return M_p; }
    pt::ptree pTree( std::string const& name ) const;
    bool useModelName() const { return M_useModelName; }
    std::map<std::string,ModelPostprocessExports> allExports() const { return M_exports; }
    std::map<std::string,std::vector<ModelPostprocessPointPosition> > const& allMeasuresPoint() const { return M_measuresPoint; }
    std::map<std::string,std::vector<ModelPostprocessNorm> > const& allMeasuresNorm() const { return M_measuresNorm; }
    std::map<std::string,std::vector<ModelPostprocessStatistics> > const& allMeasuresStatistics() const { return M_measuresStatistics; }
    bool hasExports( std::string const& name = "" ) const;
    bool hasMeasuresPoint( std::string const& name = "" ) const;
    bool hasMeasuresExtremum( std::string const& name = "" ) const;
    bool hasMeasuresNorm( std::string const& name = "" ) const;
    bool hasMeasuresStatistics( std::string const& name = "" ) const;
    ModelPostprocessExports const& exports( std::string const& name = "" ) const;
    std::vector<ModelPostprocessPointPosition> const& measuresPoint( std::string const& name = "" ) const;
    std::vector<ModelPostprocessExtremum> const& measuresExtremum( std::string const& name = "" ) const;
    std::vector<ModelPostprocessNorm> const& measuresNorm( std::string const& name = "" ) const;
    std::vector<ModelPostprocessStatistics> const& measuresStatistics( std::string const& name = "" ) const;

    void setPTree( pt::ptree const& _p );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setParameterValues( std::map<std::string,double> const& mp );

    std::map<std::string,double> toPostprocessValues() const
        {
            std::map<std::string,double> pv;
            #if 0
            for( auto const& p : *this )
                pv[p.first]=p.second.value();
            #endif
            return pv;
        }
    void saveMD(std::ostream &os);
private:
    void setup();
    void setup( std::string const& name, pt::ptree const& p );
private:
    WorldComm const& M_worldComm;
    pt::ptree M_p;
    bool M_useModelName;
    std::map<std::string,ModelPostprocessExports> M_exports;
    std::map<std::string,std::vector< ModelPostprocessPointPosition > > M_measuresPoint;
    std::map<std::string,std::vector< ModelPostprocessExtremum > > M_measuresExtremum;
    std::map<std::string,std::vector< ModelPostprocessNorm > > M_measuresNorm;
    std::map<std::string,std::vector< ModelPostprocessStatistics > > M_measuresStatistics;
    ModelPostprocessExports M_emptyExports;
    std::vector< ModelPostprocessPointPosition > M_emptyMeasuresPoint;
    std::vector< ModelPostprocessExtremum > M_emptyMeasuresExtremum;
    std::vector< ModelPostprocessNorm > M_emptyMeasuresNorm;
    std::vector< ModelPostprocessStatistics > M_emptyMeasuresStatistics;
    std::string M_directoryLibExpr;
};


}
#endif

