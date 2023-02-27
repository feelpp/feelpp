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

    void setup( nl::json const& jarg );

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

    void setup( nl::json const& jarg );
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

    void setup( nl::json const& jarg );

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

    struct PointsOverGeometry
    {
        using coord_value_type = Eigen::MatrixXd;

        virtual ~PointsOverGeometry() {}

        virtual void setParameterValues( std::map<std::string,double> const& mp ) = 0;

        //! return the coordinates of all points (as a vector of coord)
        std::vector<coord_value_type> const& coordinates() const { return M_coordinates; }

    protected :
        std::vector<coord_value_type> M_coordinates;
    };

    struct PointPosition
    {
        using coord_value_type = Eigen::MatrixXd;

        PointPosition() = default;
        PointPosition( PointPosition const& ) = default;
        PointPosition( PointPosition && ) = default;
        PointPosition& operator=( PointPosition const& ) = default;
        PointPosition& operator=( PointPosition && ) = default;

        coord_value_type const& coordinatesEvaluated() const { return M_coordinatesEvaluated; }

        void setup( std::string const& coordExprStr, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes );
        void setup( nl::json const& jData, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes );

        void setExpression( ModelExpression && mexpr ) { M_coordinatesExpr = mexpr/*; this->updateForUse()*/; }
        void setMeshMarker( std::string const& s ) { M_meshMarker = s; }

        void setParameterValues( std::map<std::string,double> const& mp )
            {
                if ( !M_coordinatesExpr.hasAtLeastOneExpr() )
                    return;
                M_coordinatesExpr.setParameterValues( mp );
                //this->updateForUse();
            }

        template <typename SymbolsExprType>
        void updateForUse( SymbolsExprType const& se )
            {
                if ( !M_coordinatesExpr.hasAtLeastOneExpr() )
                    return;
                M_coordinatesEvaluated = M_coordinatesExpr.evaluate( se );
            }
    private:
        coord_value_type M_coordinatesEvaluated;
        ModelExpression M_coordinatesExpr;
        std::string M_meshMarker;
    };

    struct PointsOverCoordinates : PointsOverGeometry
    {
        PointsOverCoordinates() = default;
        PointsOverCoordinates( PointsOverCoordinates const& ) = default;
        PointsOverCoordinates( PointsOverCoordinates && ) = default;

        void setup( nl::json const& jarg, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes );

        void setParameterValues( std::map<std::string,double> const& mp ) override
            {
                for (auto & ptPos : M_pointPositions )
                    ptPos.setParameterValues( mp );
            }

        template <typename SymbolsExprType>
        void updateForUse( SymbolsExprType const& se )
            {
                this->M_coordinates.resize( M_pointPositions.size() );
                int k=0;
                for (auto & ptPos : M_pointPositions )
                {
                    ptPos.updateForUse( se );
                    this->M_coordinates[k++] = ptPos.coordinatesEvaluated();
                }
            }
    private :
        std::vector<PointPosition> M_pointPositions;
    };


    struct PointsOverSegment : PointsOverGeometry
    {
        PointsOverSegment() : M_nPoints( 10 ) {}
        PointsOverSegment( PointsOverSegment const& ) = default;
        PointsOverSegment( PointsOverSegment && ) = default;

        void setup( nl::json const& jarg, worldcomm_t const& world, std::string const& directoryLibExpr, ModelIndexes const& indexes );

        void setParameterValues( std::map<std::string,double> const& mp ) override
            {
                M_point1.setParameterValues( mp );
                M_point2.setParameterValues( mp );
            }

        template <typename SymbolsExprType>
        void updateForUse( SymbolsExprType const& se )
            {
                M_point1.updateForUse( se );
                M_point2.updateForUse( se );
                auto const& pt1 = M_point1.coordinatesEvaluated();
                auto const& pt2 = M_point2.coordinatesEvaluated();

                this->updatePointsSampling( pt1, pt2 );
            }
    private:
        void updatePointsSampling( coord_value_type const& pt1, coord_value_type const& pt2 );
    private:
        PointPosition M_point1, M_point2;
        int M_nPoints;
    };

    struct MeasuresOutput
    {
        MeasuresOutput()
            :
            M_type( "values" )
            {}
        MeasuresOutput( MeasuresOutput const& ) = default;
        MeasuresOutput( MeasuresOutput && ) = default;
        MeasuresOutput& operator=( MeasuresOutput const& ) = default;
        MeasuresOutput& operator=( MeasuresOutput && ) = default;
        void setup( nl::json const& jarg );

        std::string const& name() const { return M_name; }
        std::string const& type() const { return M_type; }

        void setName( std::string const& name ) { M_name = name; }
    private :
        std::string M_name, M_type;
    };

    ModelPostprocessPointPosition( worldcomm_ptr_t const& world = Environment::worldCommPtr() )
        :
        super_type( world ),
        M_includeCoordinates( false )
        {}
    ModelPostprocessPointPosition( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition( ModelPostprocessPointPosition&& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition const& ) = default;
    ModelPostprocessPointPosition& operator=( ModelPostprocessPointPosition && ) = default;

    //! name given to the postprocessing
    std::string const& name() const { return M_name; }
    //! points (over each geo) where postprocessing is applied
    std::vector<std::shared_ptr<PointsOverGeometry>> const& pointsOverAllGeometry() const { return M_pointsOverAllGeometry; }
    //! return fields names used in postprocessing
    std::set<std::string> const& fields() const { return M_fields; }
    //! expressions used in postprocessing
    std::map<std::string,ModelExpression> const& expressions() const { return M_exprs; }
    //! include coordinates measure in post processing
    bool includeCoordinates() const { return M_includeCoordinates; }
    //! return info about measures output
    MeasuresOutput const& measuresOutput() const { return M_measuresOutput; }
    //! return a related tag (for example : allow to identify the used mesh)
    std::string const& tag() const { return M_tag; }


    void setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setFields( std::set<std::string> const& fields ) { M_fields = fields; }
    void addFields( std::string const& field ) { M_fields.insert( field ); }

    void setParameterValues( std::map<std::string,double> const& mp );

    template <typename SymbolsExprType>
    void updateForUse( SymbolsExprType const& se )
    {
         for ( auto & pointsOverGeometry : M_pointsOverAllGeometry )
         {
             if ( !pointsOverGeometry )
                 continue;
             if ( auto pointsOverCoordinates = std::dynamic_pointer_cast<PointsOverCoordinates>( pointsOverGeometry ) )
                 pointsOverCoordinates->updateForUse( se );
             if ( auto pointsOverSegment = std::dynamic_pointer_cast<PointsOverSegment>( pointsOverGeometry ) )
                 pointsOverSegment->updateForUse( se );
         }
    }

  private:
    std::string M_name;
    std::string M_directoryLibExpr;
    std::vector<std::shared_ptr<PointsOverGeometry>> M_pointsOverAllGeometry;
    std::set<std::string> M_fields;
    std::map<std::string,ModelExpression> M_exprs; // name -> expr
    bool M_includeCoordinates;
    MeasuresOutput M_measuresOutput;
    std::string M_tag;
};

//! store information require by the postprocessing Norm
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

    void setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types;
    ModelMarkers M_markers;
    ModelExpression M_solution, M_gradSolution, M_expr, M_gradExpr;
    uint16_type M_quadOrder, M_quad1Order;
};

//! store information require by the postprocessing Statistics
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
    //! return requires markers connection
    ModelMarkers const& requiresMarkersConnection() const { return M_requiresMarkersConnection; }
    //! return internalfaces evalutation type (i.e. mean,sum,max,...)
    std::string const& internalFacesEvalutationType() const { return M_internalFacesEvalutationType; }

    //! return true if an expression has been given
    bool hasExpr() const { return M_expr.hasAtLeastOneExpr(); }
    //! return true if a field has been given
    bool hasField() const { return !M_field.empty(); }

    void setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );

private:
    std::string M_directoryLibExpr;
    std::string M_name, M_field;
    std::set<std::string> M_types;
    ModelMarkers M_markers;
    ModelExpression M_expr;
    uint16_type M_quadOrder, M_quad1Order;
    ModelMarkers M_requiresMarkersConnection;
    std::string M_internalFacesEvalutationType;
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
    //std::tuple<bool,double> run( double val ) const;
    //! check if the arg val is close at tolerance to the reference value (the target value can be an expression which depending on symbols expr given in arg)
    //! \return tuple ( check is true, diff value )
    template <typename SymbolsExprType>
    std::tuple<bool,double> run( double val, SymbolsExprType const& se ) const
    {
        M_value = expr( M_valueExpr.exprScalar(), se ).evaluate()(0,0);
        M_tolerance = expr( M_toleranceExpr.exprScalar(), se ).evaluate()(0,0);
        return this->runImpl( val );
    }

    void setup( nl::json const& jarg, std::string const& name, ModelIndexes const& indexes );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }
    void setParameterValues( std::map<std::string,double> const& mp );
  private:
    std::tuple<bool,double> runImpl( double val ) const;
  private :
    std::string M_name;
    std::string M_directoryLibExpr;
    ModelExpression M_valueExpr;
    mutable double M_value;
    ModelExpression M_toleranceExpr;
    mutable double M_tolerance;
};

class FEELPP_EXPORT ModelPostprocess : public CommObject
{
public:
    using super = CommObject;
    ModelPostprocess( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    virtual ~ModelPostprocess();
    bool hasJsonProperties( std::string const& name = "" ) const;
    nl::json const& jsonProperties( std::string const& name = "" ) const;
    bool useModelName() const { return M_useModelName; }
    std::map<std::string,ModelPostprocessExports> const& allExports() const { return M_exports; }
    std::map<std::string,ModelPostprocessSave> const& allSave() const { return M_save; }
    std::map<std::string,ModelPostprocessQuantities> const& allMeasuresQuantities() const { return M_measuresQuantities; }
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
    std::vector<ModelPostprocessPointPosition> & measuresPoint( std::string const& name = "" );
    std::vector<ModelPostprocessNorm> const& measuresNorm( std::string const& name = "" ) const;
    std::vector<ModelPostprocessStatistics> const& measuresStatistics( std::string const& name = "" ) const;
    std::vector<ModelPostprocessCheckerMeasure> const& checkersMeasure( std::string const& name = "" ) const;

    void setPTree( nl::json const& jarg );
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
    void setup( std::string const& name, nl::json const& jarg );
private:
    nl::json M_p;
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

