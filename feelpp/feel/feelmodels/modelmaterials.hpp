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
#include <feel/feelcore/commobject.hpp>
#include <feel/feelmodels/modelexpression.hpp>
#include <feel/feelmodels/modelmarkers.hpp>

namespace Feel {

namespace pt =  boost::property_tree;

struct FEELPP_EXPORT ModelMaterialProperty
{
    ModelMaterialProperty() = default;
    ModelMaterialProperty( ModelMaterialProperty const& ) = default;
    ModelMaterialProperty( ModelMaterialProperty && ) = default;

    void setup( nl::json const& jarg, worldcomm_t const& worldComm, std::string const& directoryLibExpr, ModelIndexes const& indexes = ModelIndexes() );
    ModelExpression const& mexpr() const { return M_mexpr; }

    template <int M=1,int N=1>
    auto const& expr() const { return M_mexpr.template expr<M,N>(); }

    nl::json const& data() const { return M_data; }

    void setParameterValues( std::map<std::string,double> const& mp ) { M_mexpr.setParameterValues( mp ); }
  private:
    ModelExpression M_mexpr;
    nl::json M_data;
};

struct FEELPP_EXPORT ModelMaterial : public CommObject
{
    using super = CommObject;
    using material_property_type = ModelMaterialProperty;
    using material_properties_type = std::map<std::string, material_property_type >;

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

    material_properties_type & properties() { return M_materialProperties; }
    material_properties_type const& properties() const { return M_materialProperties; }
    material_property_type const& property( std::string const& prop ) const;

    bool hasSubMaterial( std::string const& subMat ) const;
    std::map<std::string, material_properties_type> & subMaterialProperties() { return M_subMaterialProperties; }
    std::map<std::string, material_properties_type> const& subMaterialProperties() const { return M_subMaterialProperties; }
    material_properties_type const& subMaterialProperties( std::string const& subMat ) const;
    void setSubMaterialProperty( std::string const& subMat, std::string const& property, nl::json const& jarg, ModelIndexes const& indexes = ModelIndexes() );
    void setSubMaterialProperty( std::string const& subMat, std::string const& property, std::string const& e );
    bool hasSubMaterialProperty( std::string const& submat, std::string const& prop ) const;

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

    void load( std::string const& );

    void setParameterValues( std::map<std::string,double> const& mp );
private:

    std::string M_name; /*!< Material name*/
    std::string M_directoryLibExpr;

    //! mat properties
    material_properties_type M_materialProperties;
    std::map<std::string, material_properties_type> M_subMaterialProperties;
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
                    auto const& mexpr = mparam.mexpr();
                    if ( mexpr.isEvaluable() && !add_evaluable )
                        continue;
                    if ( !mexpr.template hasExpr<ni,nj>() )
                        continue;
                    auto const& theexpr = mexpr.template expr<ni,nj>();
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
