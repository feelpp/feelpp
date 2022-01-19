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
#ifndef FEELPP_MODELMATERIALS_HPP
#define FEELPP_MODELMATERIALS_HPP 1


#include <vector>
//#include <boost/property_tree/ptree.hpp>
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelMaterial : public CommObject
{
    using super = CommObject;
    typedef ModelExpression mat_property_expr_type;
    static const uint16_type expr_order = mat_property_expr_type::expr_order;
    typedef mat_property_expr_type::expr_scalar_type expr_scalar_type;
    typedef mat_property_expr_type::expr_vectorial2_type expr_vectorial2_type;
    typedef mat_property_expr_type::expr_vectorial3_type expr_vectorial3_type;
    typedef mat_property_expr_type::expr_matrix22_type expr_matrix22_type;
    typedef mat_property_expr_type::expr_matrix33_type expr_matrix33_type;

    ModelMaterial( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    ModelMaterial( ModelMaterial const& ) = default;
    ModelMaterial( ModelMaterial&& ) = default;
    ModelMaterial& operator=( ModelMaterial const& ) = default;
    ModelMaterial& operator=( ModelMaterial && ) = default;
    ModelMaterial( std::string const& name, nl::json const& jarg,
                   worldcomm_ptr_t const& worldComm = Environment::worldCommPtr(),
                   std::string const& directoryLibExpr = "",
                   ModelIndexes const& indexes = ModelIndexes() );

    std::string const& name() const { return M_name; }
    std::set<std::string> const& meshMarkers() const { return M_meshMarkers; }
    std::set<std::string> const& physics() const { return M_physics; }
    std::string const physic() const { return M_physics.empty() ? "" : *(M_physics.begin()); }
    //! Set Name
    void setName( std::string const& name ) { M_name = name; }

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setPhysics( std::set<std::string> const& s) { M_physics = s; }
    void addPhysics( std::string const& s) { M_physics.insert( s ); }

    void setProperty( std::string const& property, nl::json const& jarg, ModelIndexes const& indexes = ModelIndexes() );
    void setProperty( std::string const& property, std::string const& e );

    bool hasProperty( std::string const& prop ) const;
    bool hasPropertyConstant( std::string const& prop ) const;
    bool hasPropertyExprScalar( std::string const& prop ) const;
    bool hasPropertyExprVectorial2( std::string const& prop ) const;
    bool hasPropertyExprVectorial3( std::string const& prop ) const;
    template <int M,int N>
    bool hasPropertyExprMatrix( std::string const& prop ) const
    {
        auto itFindProp = M_materialProperties.find( prop );
        if ( itFindProp == M_materialProperties.end() )
            return false;
        auto const& matProp = itFindProp->second;
        return matProp.template hasExprMatrix<M,N>();
    }

    std::map<std::string,mat_property_expr_type>& properties() { return M_materialProperties; }
    std::map<std::string,mat_property_expr_type> const& properties() const { return M_materialProperties; }
    mat_property_expr_type const& property( std::string const& prop ) const;
    double propertyConstant( std::string const& prop ) const;
    expr_scalar_type const& propertyExprScalar( std::string const& prop ) const;
    expr_vectorial2_type const& propertyExprVectorial2( std::string const& prop ) const;
    expr_vectorial3_type const& propertyExprVectorial3( std::string const& prop ) const;
    template <int M,int N>
    auto const& propertyExprMatrix( std::string const& prop ) const
    {
        bool hasProp = hasPropertyExprMatrix<M,N>( prop );
        CHECK( hasProp ) << "no matrix expr";
        return M_materialProperties.find( prop )->second.template exprMatrix<M,N>();
    }

    bool hasPhysics() const { return !M_physics.empty(); }
    bool hasPhysics( std::string const& physic ) const { return M_physics.find(physic) != M_physics.end(); }
    bool hasPhysics( std::initializer_list<std::string> const& physics ) const { return this->hasPhysics( std::set<std::string>( physics ) ); }
    bool hasPhysics( std::set<std::string> const& physics ) const
    {
        for ( std::string const& s : physics )
            if ( this->hasPhysics( s ) )
                return true;
        return false;
    }
    /*! Material mass density
     */
    double rho() const { return this->propertyConstant( "rho" ); }

    /*! Molecular(dynamic) viscosity
     */
    double mu() const { return this->propertyConstant( "mu" ); }

    /*! Specify the constant-pressure specific heat Cp.
     */
    double Cp() const { return this->propertyConstant( "Cp" ); }

    /*! Specify the constant-volume specific heat Cv.
     */
    double Cv() const { return this->propertyConstant( "Cv" ); }

    /*! heat diffusion coefficients
     */
    double k11() const {  return this->propertyConstant( "k11" ); }
    double k12() const { return this->propertyConstant( "k12" ); }
    double k13() const { return this->propertyConstant( "k13" ); }
    double k22() const { return this->propertyConstant( "k22" ); }
    double k23() const { return this->propertyConstant( "k23" ); }
    double k33() const { return this->propertyConstant( "k33" ); }

    /*! Material Reference temperature
     */
    double Tref() const { return this->propertyConstant( "Tref" ); }

    /*! Material coefficient for thermal expansion
     */
    double beta() const { return this->propertyConstant( "beta" ); }

    /*! heat capacity
     */
    double C() const { return this->propertyConstant( "C" ); }

    double Cs() const { return this->propertyConstant( "Cs" ); }
    double Cl() const { return this->propertyConstant( "Cl" ); }
    double L() const { return this->propertyConstant( "L" ); }
    double Ks() const { return this->propertyConstant( "Ks" ); }
    double Kl() const { return this->propertyConstant( "Kl" ); }
    double Tsol() const { return this->propertyConstant( "Tsol" ); }
    double Tliq() const { return this->propertyConstant( "Tliq" ); }

    // Mechanical properties
    /*! Young's Modulus
     */
    double E() const { return this->propertyConstant( "E" ); }

    /*! Poisson's ratio
     */
    double nu() const { return this->propertyConstant( "nu" ); }

    /*! Electrical conductivity
     */
    double sigma() const { return this->propertyConstant( "sigma" ); }

    void load( std::string const& );

    void setParameterValues( std::map<std::string,double> const& mp );
private:

    std::string M_name; /*!< Material name*/
    std::string M_directoryLibExpr;

    //! mat propeteries
    std::map<std::string, mat_property_expr_type > M_materialProperties;
    //! material physics
    std::set<std::string> M_physics;
    //! mesh markers
    ModelMarkers M_meshMarkers;

};

std::ostream& operator<<( std::ostream& os, ModelMaterial const& m );

/**
 * @brief a set of materials
 * key: mesh marker
 * name -> name of the materials - can be different
 */
class FEELPP_EXPORT ModelMaterials: public std::map<std::string,ModelMaterial>, public CommObject
{
public:
    using super = CommObject;
    using value_type = std::map<std::string,ModelMaterial>::value_type;
    ModelMaterials( worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    //ModelMaterials( pt::ptree const& p, worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() );
    virtual ~ModelMaterials() = default;
    void setPTree( nl::json const& jarg ) { M_p = jarg; setup(); }

    ModelMaterial const& material( std::string const& m ) const;

    /** return all the materials which physic is physic
     *
     */
    std::map<std::string,ModelMaterial> materialWithPhysic(std::string const& physic) const;

    /** return all the materials which physic is one of physics
     *
     */
    std::map<std::string,ModelMaterial> materialWithPhysic(std::vector<std::string> const& physics) const;
    std::set<std::string> markersWithPhysic(std::string const& physic) const;
    std::set<std::string> markersWithPhysic(std::vector<std::string> const& physic) const;

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setParameterValues( std::map<std::string,double> const& mp );

    void saveMD(std::ostream &os);

    /**
     * get the list of symbols associated to expressions
     */
    auto symbolsExpr( bool add_evaluable = false ) const
    {
        auto extract = [this,&add_evaluable](auto const& e_ij) {
            constexpr int ni = std::decay_t<decltype(hana::at_c<0>(e_ij))>::value;
            constexpr int nj = std::decay_t<decltype(hana::at_c<1>(e_ij))>::value;
            using _expr_type = std::decay_t< decltype( ModelExpression{}.expr<ni,nj>() ) >;
            symbol_expression_t<_expr_type> seParamValue;
            for ( auto const& [cname, mat] : *this )
            {
                for( auto const& [symbName,mparam] : mat.properties() )
                {
                    if ( mparam.isEvaluable() && !add_evaluable )
                        continue;
                    if ( !mparam.template hasExpr<ni,nj>() )
                        continue;
                    auto const& theexpr = mparam.template expr<ni,nj>();
                    VLOG(1) << "material " << cname << " has property " << symbName;
                    seParamValue.add( cname+"_"+symbName, theexpr, SymbolExprComponentSuffix( ni, nj ) );
                }
            }
            return seParamValue;
        };
        auto tupleSymbolExprs = hana::transform( ModelExpression::expr_shapes, extract );

        return Feel::vf::symbolsExpr( SymbolsExpr( tupleSymbolExprs ) );
    }
private:
    void setup();
private:
    nl::json M_p;
    std::string M_directoryLibExpr;

};



FEELPP_EXPORT inline ModelMaterial
material( ModelMaterials::value_type const& m )
{
    return m.second;
}

FEELPP_EXPORT inline std::string
name( ModelMaterials::value_type const& m )
{
    return m.first;
}

}
#endif
