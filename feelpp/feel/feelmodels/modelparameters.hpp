/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015
 
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
#ifndef FEELPP_MODELPARAMETERS_HPP
#define FEELPP_MODELPARAMETERS_HPP 1

#include <map>

#include <feel/feelcore/commobject.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelfit/interpolator.hpp>
#include <feel/feelfit/fit.hpp>
//#include <feel/feelvf/symbolsexpr.hpp>
#include <feel/feelmodels/modelexpression.hpp>

namespace Feel {

//namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelParameter
{
    ModelParameter() = default;
    ModelParameter( ModelParameter const& ) = default;
    ModelParameter( ModelParameter&& ) = default;
    ModelParameter& operator=( ModelParameter const& ) = default;
    ModelParameter& operator=( ModelParameter && ) = default;

    ModelParameter( std::string const& name, double value, double min = 0., double max = 0., std::string const& desc = "" )
        :
        M_name( name ),
        M_type( "value" ),
        M_min( min ),
        M_max( max ),
        M_expr( value ),
        M_desc( desc )
        {}

    ModelParameter( std::string const& name, ModelExpression const& mexpr,
                    double min = 0., double max = 0., std::string const& desc = "" )
        :
        M_name( name ),
        M_type( "expression" ),
        M_min( min ),
        M_max( max ),
        M_expr( mexpr ),
        M_desc( desc )
        {}
    ModelParameter( std::string const& name, std::shared_ptr<Interpolator> interpolator, ModelExpression const& mexpr, std::string const& desc = "" )
        :
        M_name( name ),
        M_type( "fit" ),
        M_min( 0. ),
        M_max( 0. ),
        M_expr( mexpr ),
        M_interpolator( interpolator ),
        M_desc( desc )
        {}

    std::string const& name() const { return M_name; }
    void setName( std::string const& name ) { M_name = name; }

    std::string const& type() const { return M_type; }

    double value() const { return M_expr.value(); }
    void setValue( double v ) { M_expr.setValue( v ); }
    double min() const { return M_min; }
    void setMin( double v ) { M_min = v; }
    double max() const { return M_max; }
    void setMax( double v ) { M_max = v; }
    bool hasMinMax() const { return M_min != 0 || M_max != 0; }

    std::string const& description() const { return M_desc; }
    void setDescription( std::string const& desc ) { M_desc = desc; }

    //bool hasExpression() const { return M_expr.hasAtLeastOneExpr(); }
    template <int M=1,int N=1>
    bool hasExpression() const { return M_expr.hasExpr<M,N>(); }

    template <int M=1,int N=1>
        auto const& expression() const
    {
        if ( !this->hasExpression<M,N>() )
            CHECK( false ) << "no expression defined";
        return M_expr.template expr<M,N>();
    }

    bool isEvaluable() const
    {
        if ( this->type() == "expression" || this->type() == "value" )
            return M_expr.isEvaluable();
        else
            return false;
    }

    auto evaluate() const
    {
        return M_expr.evaluate();
    }

    bool hasFitInterpolator() const { return ( M_interpolator )? true : false; }
    std::shared_ptr<Interpolator> fitInterpolator() const { return M_interpolator; }

    void setParameterValues( std::map<std::string,double> const& mp )
    {
        M_expr.setParameterValues( mp );
    }
    void eraseParameterValues( std::set<std::string> const& symbNames )
    {
         M_expr.eraseParameterValues( symbNames );
    }
    void updateParameterValues( std::map<std::string,double> & pv ) const
    {
        M_expr.updateParameterValues( this->name(), pv );
    }
    void updateSymbolNames( std::set<std::string> & outputSymbNames ) const
    {
        M_expr.updateSymbolNames( this->name(), outputSymbNames );
    }
    void updateInformationObject(  std::string const& symbol, nl::json::array_t & ja ) const;
private:
    std::string M_name, M_type;
    double M_min, M_max;
    ModelExpression M_expr;
    std::shared_ptr<Interpolator> M_interpolator;
    std::string M_desc;

};

//!
//! class for Model Parameters
//!
class ModelParameters: public std::map<std::string,ModelParameter>, public CommObject
{
public:
    using super=CommObject;
    ModelParameters( worldcomm_ptr_t const& world = Environment::worldCommPtr() );
    ModelParameters( ModelParameters const& ) = default;
    virtual ~ModelParameters();
    void setPTree( nl::json const& jarg );
    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void updateParameterValues();
    void setParameterValues( std::map<std::string,double> const& mp );
    std::map<std::string,double> toParameterValues() const;


    auto symbolsExpr() const
        {
            auto tupleSymbolExprs = hana::transform( ModelExpression::expr_shapes, [this](auto const& e_ij)
                                                     {
                                                            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
                                                            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;

                                                            using _expr_type = std::decay_t< decltype( ModelExpression{}.template expr<ni,nj>() ) >;
                                                            symbol_expression_t<_expr_type> seParamValue;
                                                            for( auto const& [symbName, mparam] : *this )
                                                            {
                                                                if ( mparam.isEvaluable() )
                                                                    continue;

                                                                if ( mparam.type() == "expression" )
                                                                {
                                                                    if ( !mparam.template hasExpression<ni,nj>() )
                                                                        continue;
                                                                    VLOG(1) << "add parameter symbolsexpr " << symbName;
                                                                    auto const& theexpr = mparam.template expression<ni,nj>();
                                                                    seParamValue.add( symbName, theexpr, SymbolExprComponentSuffix( ni, nj ) );
                                                                }
                                                            }
                                                            return seParamValue;
                                                        } );

            using _fit_expr_type = Expr< Fit<typename ModelExpression::expr_scalar_type,0> >;
            symbol_expression_t<_fit_expr_type> seFit;
            for( auto const& p : *this )
            {
                auto const& mparam = p.second;
                if ( mparam.type() != "fit" )
                    continue;
                auto const& theexpr = mparam.template expression<1,1>();
                std::string symbName = p.first;
                seFit.add( symbName, fit( theexpr, mparam.fitInterpolator() ) );
            }

            return Feel::vf::symbolsExpr( SymbolsExpr( tupleSymbolExprs ), seFit );
        }

    void updateInformationObject( nl::json & p ) const;

   void saveMD(std::ostream &os);
private:
    void setup();
private:
    nl::json M_p;
    std::string M_directoryLibExpr;
};


}
#endif
