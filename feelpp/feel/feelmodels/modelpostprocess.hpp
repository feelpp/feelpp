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

#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

class FEELPP_EXPORT ModelPostprocessExports : public CommObject
{
  public :
    using super = CommObject;
    // (name,expr,markers,representations,tags)
    using vector_export_expr_type = std::vector<std::tuple<std::string,ModelExpression,ModelMarkers,std::set<std::string>,std::set<std::string>>>;

    ModelPostprocessExports( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super( world ),
        M_format( "ensightgold" )
    {}
    ModelPostprocessExports( ModelPostprocessExports const& ) = default;
    ModelPostprocessExports( ModelPostprocessExports&& ) = default;
    ModelPostprocessExports& operator=( ModelPostprocessExports const& ) = default;
    ModelPostprocessExports& operator=( ModelPostprocessExports && ) = default;

    std::set<std::string> const& fields() const { return M_fields; }
    vector_export_expr_type const& expressions() const { return M_exprs; }
    std::string const& format() const { return M_format; }

    void setup( pt::ptree const& p );

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );
  private :
    std::string M_directoryLibExpr;
    std::set<std::string> M_fields;
    std::string M_format;
    vector_export_expr_type M_exprs;
};

class FEELPP_EXPORT ModelPostprocessSave
{
  public :
    ModelPostprocessSave() = default;
    ModelPostprocessSave( ModelPostprocessSave const& ) = default;
    ModelPostprocessSave( ModelPostprocessSave&& ) = default;
    ModelPostprocessSave& operator=( ModelPostprocessSave const& ) = default;
    ModelPostprocessSave& operator=( ModelPostprocessSave && ) = default;

    std::set<std::string> const& fieldsNames() const { return M_fieldsNames; }
    std::string const& fieldsFormat() const { return M_fieldsFormat; }

    void setup( pt::ptree const& p );
  private :
    std::set<std::string> M_fieldsNames;
    std::string M_fieldsFormat;
};

class FEELPP_EXPORT ModelPostprocessQuantities : public CommObject
{
    using super = CommObject;
  public :
    ModelPostprocessQuantities( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super( world )
    {}
    ModelPostprocessQuantities( ModelPostprocessQuantities const& ) = default;
    ModelPostprocessQuantities( ModelPostprocessQuantities&& ) = default;
    ModelPostprocessQuantities& operator=( ModelPostprocessQuantities const& ) = default;
    ModelPostprocessQuantities& operator=( ModelPostprocessQuantities && ) = default;

    std::set<std::string> const& quantities() const { return M_quantities; }
    std::map<std::string,ModelExpression> const& expressions() const { return M_exprs; }

    void setup( pt::ptree const& p );

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

  private :
    std::string M_directoryLibExpr;
    std::set<std::string> M_quantities;
    std::map<std::string,ModelExpression> M_exprs;
};

class FEELPP_EXPORT ModelPostprocessPointPosition : public CommObject
{
    using super_type = CommObject;
public :

    struct PointPosition
    {
        using coord_value_type = Eigen::MatrixXd;

        PointPosition() = default;
        PointPosition( PointPosition const& ) = default;
        PointPosition( PointPosition && ) = default;
        PointPosition& operator=( PointPosition const& ) = default;
        PointPosition& operator=( PointPosition && ) = default;

        coord_value_type const& coordinatesEvaluated() const { return M_coordinatesEvaluated; }

        void setExpression( ModelExpression && mexpr ) { M_coordinatesExpr = mexpr; this->updateForUse(); }
        void setMeshMarker( std::string const& s ) { M_meshMarker = s; }

        void setParameterValues( std::map<std::string,double> const& mp )
            {
                if ( !M_coordinatesExpr.hasAtLeastOneExpr() )
                    return;
                M_coordinatesExpr.setParameterValues( mp );
                this->updateForUse();
            }

        void updateForUse()
            {
                if ( !M_coordinatesExpr.hasAtLeastOneExpr() )
                    return;
                M_coordinatesEvaluated = M_coordinatesExpr.evaluate();
            }
    private:
        coord_value_type M_coordinatesEvaluated;
        ModelExpression M_coordinatesExpr;
        std::string M_meshMarker;
    };

    ModelPostprocessPointPosition( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super_type( world )
        {}
    ModelPostprocessPointPosition( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition( ModelPostprocessPointPosition&& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition && ) = default;

    std::string const& name() const { return M_name; }
    std::vector<PointPosition> const& pointsSampling() const { return M_pointsSampling; }

    std::set<std::string> const& fields() const { return M_fields; }
    std::set<std::string> & fields() { return M_fields; }

    std::map<std::string,std::tuple<ModelExpression,std::string> > const& expressions() const { return M_exprs; }

    void setPTree( pt::ptree const& _p, std::string const& name, ModelIndexes const& indexes ) { M_p = _p; this->setup( name, indexes ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { M_fields = fields; }
    void addFields( std::string const& field ) { M_fields.insert( field ); }
    void setParameterValues( std::map<std::string,double> const& mp );
private:
    void setup( std::string const& name, ModelIndexes const& indexes );
private:
    std::string M_name;
    pt::ptree M_p;
    std::string M_directoryLibExpr;
    std::vector<PointPosition> M_pointsSampling;
    std::set<std::string> M_fields;
    std::map<std::string,std::tuple<ModelExpression,std::string> > M_exprs; // name -> ( expr, tag )
};

//! store informations require by the postprocessing Norm
class FEELPP_EXPORT ModelPostprocessNorm : public CommObject
{
public :
    using super = CommObject;
    ModelPostprocessNorm( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super( world ),
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

    void setPTree( pt::ptree const& _p, std::string const& name, ModelIndexes const& indexes ) { M_p = _p; this->setup( name, indexes ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup( std::string const& name, ModelIndexes const& indexes );
private:
    pt::ptree M_p;
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types;
    ModelMarkers M_markers;
    ModelExpression M_solution, M_gradSolution, M_expr, M_gradExpr;
    uint16_type M_quadOrder, M_quad1Order;
};

//! store informations require by the postprocessing Statistics
class FEELPP_EXPORT ModelPostprocessStatistics : public CommObject
{
public :
    using super = CommObject;
    ModelPostprocessStatistics( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super( world ),
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

    void setPTree( pt::ptree const& _p, std::string const& name, ModelIndexes const& indexes ) { M_p = _p; this->setup( name, indexes ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    void setup( std::string const& name, ModelIndexes const& indexes );
private:
    pt::ptree M_p;
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types;
    ModelMarkers M_markers;
    ModelExpression M_expr;
    uint16_type M_quadOrder, M_quad1Order;
};

class FEELPP_EXPORT ModelPostprocessCheckerMeasure : public CommObject
{
    using super = CommObject;
  public :
    ModelPostprocessCheckerMeasure( worldcomm_ptr_t const& world = Environment::worldCommPtr(), double value=0., double tol=0.01 )
        :
        super( world ),
        M_value( value ),
        M_tolerance( tol )
        {}
    ModelPostprocessCheckerMeasure( ModelPostprocessCheckerMeasure const& ) = default;
    ModelPostprocessCheckerMeasure( ModelPostprocessCheckerMeasure&& ) = default;
    ModelPostprocessCheckerMeasure& operator=( ModelPostprocessCheckerMeasure const& ) = default;
    ModelPostprocessCheckerMeasure& operator=( ModelPostprocessCheckerMeasure && ) = default;

    //! name of a measure to check
    std::string const& name() const { return M_name; }
    //! reference value
    double value() const { return M_value; }
    //! tolerance used in comparison
    double tolerance() const { return M_tolerance; }
    //! check if the arg val is close at tolerance to the reference value
    //! \return tuple ( check is true, diff value )
    std::tuple<bool,double> run( double val ) const;
    //! check if the arg val is close at tolerance to the reference value (the target value can be an expression which depending on symbols expr given in arg)
    //! \return tuple ( check is true, diff value )
    template <typename SymbolsExprType>
    std::tuple<bool,double> run( double val, SymbolsExprType const& se ) const
    {
        M_value = expr( M_valueExpr.exprScalar(), se ).evaluate()(0,0);
        return this->runImpl( val );
    }

    void setPTree( pt::ptree const& _p, std::string const& name, ModelIndexes const& indexes ) { M_p = _p; this->setup( name, indexes ); }
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );
  private:
    std::tuple<bool,double> runImpl( double val ) const;
    void setup( std::string const& name, ModelIndexes const& indexes );
  private :
    pt::ptree M_p;
    std::string M_name;
    std::string M_directoryLibExpr;
    ModelExpression M_valueExpr;
    mutable double M_value;
    double M_tolerance;
};

class FEELPP_EXPORT ModelPostprocess : public CommObject
{
public:
    using super = CommObject;
    ModelPostprocess( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    ModelPostprocess( pt::ptree const& p, worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    virtual ~ModelPostprocess();
    pt::ptree const& pTree() const { return M_p; }
    pt::ptree & pTree() { return M_p; }
    pt::ptree pTree( std::string const& name ) const;
    bool useModelName() const { return M_useModelName; }
    std::map<std::string,ModelPostprocessExports> allExports() const { return M_exports; }
    std::map<std::string,ModelPostprocessSave> allSave() const { return M_save; }
    std::map<std::string,ModelPostprocessQuantities> allMeasuresQuantities() const { return M_measuresQuantities; }
    std::map<std::string,std::vector<ModelPostprocessPointPosition> > const& allMeasuresPoint() const { return M_measuresPoint; }
    std::map<std::string,std::vector<ModelPostprocessNorm> > const& allMeasuresNorm() const { return M_measuresNorm; }
    std::map<std::string,std::vector<ModelPostprocessStatistics> > const& allMeasuresStatistics() const { return M_measuresStatistics; }
    std::map<std::string,std::vector<ModelPostprocessCheckerMeasure> > const& allCheckersMeasure() const { return M_checkersMeasure; }
    bool hasExports( std::string const& name = "" ) const;
    bool hasSave( std::string const& name = "" ) const;
    bool hasMeasuresQuantities( std::string const& name = "" ) const;
    bool hasMeasuresPoint( std::string const& name = "" ) const;
    bool hasMeasuresNorm( std::string const& name = "" ) const;
    bool hasMeasuresStatistics( std::string const& name = "" ) const;
    bool hasCheckersMeasure( std::string const& name = "" ) const;
    ModelPostprocessExports const& exports( std::string const& name = "" ) const;
    ModelPostprocessSave const& save( std::string const& name = "" ) const;
    ModelPostprocessQuantities const& measuresQuantities( std::string const& name = "" ) const;
    std::vector<ModelPostprocessPointPosition> const& measuresPoint( std::string const& name = "" ) const;
    std::vector<ModelPostprocessNorm> const& measuresNorm( std::string const& name = "" ) const;
    std::vector<ModelPostprocessStatistics> const& measuresStatistics( std::string const& name = "" ) const;
    std::vector<ModelPostprocessCheckerMeasure> const& checkersMeasure( std::string const& name = "" ) const;

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
    pt::ptree M_p;
    bool M_useModelName;
    std::map<std::string,ModelPostprocessExports> M_exports;
    std::map<std::string,ModelPostprocessSave> M_save;
    std::map<std::string,ModelPostprocessQuantities> M_measuresQuantities;
    std::map<std::string,std::vector< ModelPostprocessPointPosition > > M_measuresPoint;
    std::map<std::string,std::vector< ModelPostprocessNorm > > M_measuresNorm;
    std::map<std::string,std::vector< ModelPostprocessStatistics > > M_measuresStatistics;
    std::map<std::string,std::vector< ModelPostprocessCheckerMeasure > > M_checkersMeasure;
    ModelPostprocessExports M_emptyExports;
    ModelPostprocessSave M_emptySave;
    ModelPostprocessQuantities M_emptyMeasuresQuantities;
    std::vector< ModelPostprocessPointPosition > M_emptyMeasuresPoint;
    std::vector< ModelPostprocessNorm > M_emptyMeasuresNorm;
    std::vector< ModelPostprocessStatistics > M_emptyMeasuresStatistics;
    std::vector< ModelPostprocessCheckerMeasure > M_emptyCheckersMeasure;
    std::string M_directoryLibExpr;
};


}
#endif

