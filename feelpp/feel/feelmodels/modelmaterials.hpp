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
#include <boost/property_tree/ptree.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelMaterial
{
    typedef ModelExpression mat_property_expr_type;
    static const uint16_type expr_order = mat_property_expr_type::expr_order;
    typedef mat_property_expr_type::expr_scalar_type expr_scalar_type;
    typedef mat_property_expr_type::expr_vectorial2_type expr_vectorial2_type;
    typedef mat_property_expr_type::expr_vectorial3_type expr_vectorial3_type;
    typedef mat_property_expr_type::expr_matrix22_type expr_matrix22_type;
    typedef mat_property_expr_type::expr_matrix33_type expr_matrix33_type;

    ModelMaterial( WorldComm const& worldComm = Environment::worldComm() );
    ModelMaterial( ModelMaterial const& ) = default;
    ModelMaterial( ModelMaterial&& ) = default;
    ModelMaterial& operator=( ModelMaterial const& ) = default;
    ModelMaterial& operator=( ModelMaterial && ) = default;
    ModelMaterial( std::string const& name, pt::ptree const& p,
                   WorldComm const& worldComm = Environment::worldComm(),
                   std::string const& directoryLibExpr = "" );

    std::string const& name() const { return M_name; }
    std::set<std::string> const& meshMarkers() const { return M_meshMarkers; }
    std::set<std::string> const& physics() const { return M_physics; }
    std::string const physic() const { return M_physics.empty() ? "" : *(M_physics.begin()); }
    //! return the property tree
    pt::ptree const& pTree() const { return M_p; }
    /*! Set Name
     */
    void setName( std::string const& name ) { M_name = name; }

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setPhysics( std::set<std::string> const& s) { M_physics = s; }
    void addPhysics( std::string const& s) { M_physics.insert( s ); }

    void setProperty( std::string const& property, pt::ptree const& p );
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
    bool hasPhysics( std::initializer_list<std::string> const& physics ) const
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

    /**
     *
     */
    pt::ptree getNode( std::string const& key ) const {
        try { return M_p.get_child( key ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    std::string getString( std::string const& key ) const {
        try { return M_p.get<std::string>( key ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    int getInt( std::string const& key ) const {
        try { return M_p.get<int>( key ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    double getDouble( std::string const& key ) const {
        try { return M_p.get<double>( key ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    std::vector<std::string> getVecString( std::string const& key ) const {
        try {
            std::vector<std::string> res;
            for( auto const& item : M_p.get_child( key ) )
                res.push_back(item.second.template get_value<std::string>());
            return res;
        }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    std::vector<double> getVecDouble( std::string const& key ) const {
        try {
            std::vector<double> res;
            for( auto const& item : M_p.get_child( key ) )
                res.push_back(item.second.template get_value<double>());
            return res;
        }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    Expr<GinacEx<2> > getScalar( std::string const& key ) const {
        try { return expr( M_p.get<std::string>( key ) ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    // Expr<GinacEx<2> > getScalar( std::string const& key, std::pair<std::string,double> const& params ) { return expr( M_p.get( key, "0" ), params ); }
    /**
     *
     */
    Expr<GinacEx<2> > getScalar( std::string const& key, std::map<std::string,double> const& params ) const {
        try { return expr( M_p.get<std::string>( key ), params ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    template<typename ExprT> auto/*Expr<GinacExVF<ExprT> >*/ getScalar( std::string const& key,
                                                                std::string const& sym, ExprT e ) const {
        try { return expr( M_p.get<std::string>( key ), sym, e ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    template<typename ExprT> auto/*Expr<GinacExVF<ExprT> >*/ getScalar( std::string const& key,
                                                                std::initializer_list<std::string> const& sym,
                                                                std::initializer_list<ExprT> e ) const {
        try { return expr( M_p.get<std::string>( key ), sym, e ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    template<typename ExprT> auto/*Expr<GinacExVF<ExprT> >*/ getScalar( std::string const& key,
                                                                std::vector<std::string> const& sym,
                                                                std::vector<ExprT> e ) const {
        try { return expr( M_p.get<std::string>( key ), sym, e ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }

    /**
     *
     */
    template<typename ExprT> auto/*Expr<GinacExVF<ExprT> >*/ getScalar( std::string const& key,
                                                                std::initializer_list<std::string> const& sym,
                                                                std::initializer_list<ExprT> e,
                                                                std::map<std::string, double> params ) const {
        try {
            auto ex = expr( M_p.get<std::string>( key ), sym, e );
            ex.setParameterValues( params );
            return ex;
        }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    template<int T> Expr<GinacMatrix<T,1,2> > getVector( std::string const& key ) const {
        try { return expr<T,1>( M_p.get<std::string>( key ) ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    template<int T> Expr<GinacMatrix<T,1,2> > getVector( std::string const& key,
                                                         std::pair<std::string,double> const& params ) const {
        try { return expr<T,1>( M_p.get<std::string>( key ), params ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    template<int T> Expr<GinacMatrix<T,1,2> > getVector( std::string const& key,
                                                         std::map<std::string,double> const& params ) const {
        try { return expr<T,1>( M_p.get<std::string>( key ), params ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    // template<int T, typename ExprT> Expr<GinacMatrixVF<ExprT> > getVector( std::string const& key, std::string const& sym, ExprT e );
    // template<int T, typename ExprT> Expr<GinacMatrix<T,1,2> > getVector( std::string const& key, std::initializer_list<std::string> const& sym, std::initializer_list<ExprT> e );
    // template<int T, typename ExprT> Expr<GinacMatrix<T,1,2> > getVector( std::string const& key, std::vector<std::string> const& sym, std::vector<ExprT> e );
    template<int T1, int T2=T1> Expr<GinacMatrix<T1,T2,2> > getMatrix( std::string const& key ) const {
        try { return expr<T1,T2>( M_p.get<std::string>( key ) ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    template<int T1, int T2=T1> Expr<GinacMatrix<T1,T2,2> > getMatrix( std::string const& key,
                                                                       std::pair<std::string,double> const& params ) const {
        try { return expr<T1,T2>( M_p.get<std::string>( key ), params ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    template<int T1, int T2=T1> Expr<GinacMatrix<T1,T2,2> > getMatrix( std::string const& key,
                                                                       std::map<std::string,double> const& params ) const {
        try { return expr<T1,T2>( M_p.get<std::string>( key ), params ); }
        catch( pt::ptree_error const& e ) {
            cerr << "key " << key << ": " << e.what() << std::endl;
            exit(1);
        }
    }
    // template<int T1, int T2, typename ExprT> Expr<GinacExVF<ExprT> > getMatrix( std::string const& key, std::string const& sym, ExprT e );
    // template<int T1, int T2, typename ExprT> Expr<GinacExVF<ExprT> > getMatrix( std::string const& key, std::initializer_list<std::string> const& sym, std::initializer_list<ExprT> e );
    // template<int T1, int T2, typename ExprT> Expr<GinacExVF<ExprT> > getMatrix( std::string const& key, std::vector<std::string> const& sym, std::vector<ExprT> e );

private:

    WorldComm const * M_worldComm;
    std::string M_name; /*!< Material name*/
    pt::ptree M_p;
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
class FEELPP_EXPORT ModelMaterials: public std::map<std::string,ModelMaterial>
{
public:
    using value_type = std::map<std::string,ModelMaterial>::value_type;
    ModelMaterials( WorldComm const& worldComm = Environment::worldComm() );
    ModelMaterials( pt::ptree const& p, WorldComm const& worldComm = Environment::worldComm() );
    virtual ~ModelMaterials() = default;
    void setPTree( pt::ptree const& _p ) { M_p = _p; setup(); }

    ModelMaterial const&
    material( std::string const& m ) const
        {
            auto it = this->find( m );
            if ( it == this->end() )
                throw std::invalid_argument( std::string("ModelMaterial: Invalid material name ") + m );
            return it->second;
        }

    /** return all the materials which physic is physic
     *
     */
    std::map<std::string,ModelMaterial> materialWithPhysic(std::string const& physic) const
    {
        std::map<std::string,ModelMaterial> mat;
        std::copy_if(this->begin(),this->end(),std::inserter(mat,mat.begin()),
                     [physic](std::pair<std::string,ModelMaterial> const& mp)
                     { return mp.second.hasPhysics(physic); } );
        return mat;
    }

    /** return all the materials which physic is one of physics
     *
     */
    std::map<std::string,ModelMaterial> materialWithPhysic(std::vector<std::string> const& physics) const
    {
        std::map<std::string,ModelMaterial> mat;
        std::copy_if(this->begin(),this->end(),std::inserter(mat,mat.begin()),
                     [physics](std::pair<std::string,ModelMaterial> const& mp)
                     {
                         bool b = false;
                         for( auto const& p : physics )
                             b = b || mp.second.hasPhysics(p);
                         return b;
                     });
        return mat;
    }

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    void setParameterValues( std::map<std::string,double> const& mp );

    void saveMD(std::ostream &os);
private:
    void setup();
private:
    WorldComm const* M_worldComm;
    pt::ptree M_p;
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
