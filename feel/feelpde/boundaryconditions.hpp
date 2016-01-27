/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 25 Jan 2015
 
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
#ifndef FEELPP_BOUNDARYCONDITIONS_HPP
#define FEELPP_BOUNDARYCONDITIONS_HPP 1

#include <map>
#include <set>
#include <string>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelvf/ginac.hpp>

#include <boost/property_tree/ptree.hpp>

namespace Feel
{
namespace pt =  boost::property_tree;
struct ExpressionStringAtMarker : public std::tuple<std::string,std::string,std::string>
{
    typedef std::tuple<std::string,std::string,std::string> super;

    ExpressionStringAtMarker( super && s ) : super( s ) { M_meshMarkers.push_back( this->marker() ); }
    
    /**
     * @return the marker
     */
    std::string const& marker() const { return std::get<0>( *this ); }
    
    /**
     * @return the expression
     */
    std::string const& expression() const { return std::get<1>( *this ); }
    
    /**
     * @return the expression
     */
    std::string const& expression1() const { return std::get<1>( *this ); }

    /**
     * @return the expression
     */
    std::string const& expression2() const { return std::get<2>( *this ); }

    bool hasExpression() const { return !std::get<1>( *this ).empty(); } 
    bool hasExpression1() const { return !std::get<1>( *this ).empty(); } 
    bool hasExpression2() const { return !std::get<2>( *this ).empty(); }

    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }

    void setMeshMarkers( std::list<std::string> const& s ) { M_meshMarkers=s; }

private :
    std::list<std::string> M_meshMarkers;
};
/**
 * Defines boundary conditions dictionary
 */
class BoundaryConditions
    :
        public std::map<std::string,std::map<std::string,std::vector<ExpressionStringAtMarker>>>
{
    typedef std::map<std::string,std::map<std::string,std::vector<ExpressionStringAtMarker>>> super;
    
  public:
    using value_type = typename super::value_type;

    BoundaryConditions( WorldComm const& world = Environment::worldComm() );
    
    /**
     * constructor from an \c initializer_list<>
     */
    BoundaryConditions( std::initializer_list<value_type> l, WorldComm const& world = Environment::worldComm() )
        :
        super(l),
        M_worldComm( world ),
        M_prefix()
        {}
    
    BoundaryConditions( std::string const& prefix, WorldComm const& world = Environment::worldComm() );
    BoundaryConditions( BoundaryConditions const& b ) = default;
    BoundaryConditions( BoundaryConditions && b ) = default;
    BoundaryConditions& operator=( BoundaryConditions const& bc ) = default;
    BoundaryConditions& operator=( BoundaryConditions && bc ) = default;
    ~BoundaryConditions() = default;
    
    /**
     * @return prefix 
     */
    std::string const& prefix() const { return M_prefix; }
    
    /**
     * \param p prefix to be set
     */
    void setPrefix( std::string const& p ) { M_prefix = p; }

    void setPTree( pt::ptree const& p );

    void setDirectoryLibExpr( std::string const& directoryLibExpr ) { M_directoryLibExpr = directoryLibExpr; }

    /**
     * load property tree from file \p filename 
     */
    void load(const std::string &filename);

    void saveMD(std::ostream &os);

    std::pair<bool,int> iparam( std::string const& field,std::string const& bc, std::string const& marker, std::string const& param ) const;
    std::pair<bool,double> dparam( std::string const& field,std::string const& bc, std::string const& marker, std::string const& param ) const;
    std::pair<bool,std::string> sparam( std::string const& field,std::string const& bc, std::string const& marker, std::string const& param ) const;

    std::list<std::string> markers( std::string const& field, std::string const& type ) const;
    std::list<std::string> markers( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const;

    /**
     * retrieve scalar field \p field with boundary conditions of type \p type
     */
    template<int Order=2> map_scalar_field<Order> getScalarFields( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
    {
        using namespace Feel::vf;
        map_scalar_field<Order> m_f;
        for ( auto const& it : listKeys )
        {
            auto curFields = getScalarFields<Order>( std::string(it.first), std::string(it.second) );
            for ( auto const& curField : curFields )
                m_f[curField.first] = std::move(curField.second);
        }
        return std::move(m_f);
    }
    template<int Order=2> map_scalar_field<Order> getScalarFields( std::string && field, std::string && type ) const
    {
        using namespace Feel::vf;
        map_scalar_field<Order> m_f;
        auto const& itFindField = this->find(field);
        if ( itFindField == this->end() ) return std::move(m_f);
        auto const& itFindType = itFindField->second.find(type);
        if ( itFindType == itFindField->second.end() ) return std::move(m_f);
        for ( auto f : itFindType->second )
        {
            LOG(INFO) << "Building expr " << f.expression() << " for " << f.marker();
            m_f[std::get<0>(f)] = expr<Order>( f.expression(), "", M_worldComm, M_directoryLibExpr );
        }
        return std::move(m_f);
    }
    /**
     * retrieve scalar field pair \p field with boundary conditions of type \p type
     */
    template<int Order=2> map_scalar_fields<Order> getScalarFieldsList( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
    {
        using namespace Feel::vf;
        map_scalar_fields<Order> m_f;
        for ( auto const& it : listKeys )
        {
            auto curFieldsList = getScalarFieldsList<Order>( std::string(it.first), std::string(it.second) );
            for ( auto const& curFieldList : curFieldsList )
                for ( auto const& curField : curFieldList.second )
                    m_f[curFieldList.first].push_back( std::move(curField) );
        }
        return std::move(m_f);
    }
    template<int Order=2> map_scalar_fields<Order> getScalarFieldsList( std::string && field, std::string && type ) const
        {
            using namespace Feel::vf;
            map_scalar_fields<Order> m_f;
            auto const& itFindField = this->find(field);
            if ( itFindField == this->end() ) return std::move(m_f);
            auto const& itFindType = itFindField->second.find(type);
            if ( itFindType == itFindField->second.end() ) return std::move(m_f);
            for ( auto f : itFindType->second )
            {
                CHECK( f.hasExpression1() && f.hasExpression2() ) << "Invalid call";
                LOG(INFO) << "Building expr1 " << f.expression1() << " for " << f.marker();
                m_f[f.marker()].push_back( expr<Order>( f.expression1(), "", M_worldComm, M_directoryLibExpr ) );
                LOG(INFO) << "Building expr2 " << f.expression2() << " for " << f.marker();
                m_f[f.marker()].push_back( expr<Order>( f.expression2(), "", M_worldComm, M_directoryLibExpr ) );
            }
            return std::move(m_f);
        }
    template<int d> map_vector_field<d> getVectorFields( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
    {
        using namespace Feel::vf;
        map_vector_field<d> m_f;
        for ( auto const& it : listKeys )
        {
            auto curFields = getVectorFields<d>( std::string(it.first), std::string(it.second) );
            for ( auto const& curField : curFields )
                m_f[curField.first] = std::move(curField.second);
        }
        return std::move(m_f);
    }
    template<int d> map_vector_field<d> getVectorFields( std::string && field, std::string && type )  const
    {
        using namespace Feel::vf;
        map_vector_field<d> m_f;
        auto const& itFindField = this->find(field);
        if ( itFindField == this->end() ) return std::move(m_f);
        auto const& itFindType = itFindField->second.find(type);
        if ( itFindType == itFindField->second.end() ) return std::move(m_f);
        for ( auto f : itFindType->second )
        {
            LOG(INFO) << "Building expr " << f.expression() << " for " << std::get<0>(f);
            m_f[std::get<0>(f)] = expr<d,1,2>( f.expression(), "", M_worldComm, M_directoryLibExpr );
        }
        return std::move(m_f);
    }
    template<int d> map_vector_fields<d> getVectorFieldsList( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
    {
        using namespace Feel::vf;
        map_vector_fields<d> m_f;
        for ( auto const& it : listKeys )
        {
            auto curFieldsList = getVectorFieldsList<d>( std::string(it.first), std::string(it.second) );
            for ( auto const& curFieldList : curFieldsList )
                for ( auto const& curField : curFieldList.second )
                    m_f[curFieldList.first].push_back( std::move(curField) );
        }
        return std::move(m_f);
    }
    template<int d> map_vector_fields<d> getVectorFieldsList( std::string && field, std::string && type )  const
    {
        using namespace Feel::vf;
        map_vector_fields<d> m_f;
        auto const& itFindField = this->find(field);
        if ( itFindField == this->end() ) return std::move(m_f);
        auto const& itFindType = itFindField->second.find(type);
        if ( itFindType == itFindField->second.end() ) return std::move(m_f);
        for ( auto f : itFindType->second )
        {
            CHECK( f.hasExpression1() && f.hasExpression2() ) << "Invalid call";
            LOG(INFO) << "Building expr " << f.expression() << " for " << std::get<0>(f);
            m_f[std::get<0>(f)].push_back( expr<d,1,2>( f.expression1(), "", M_worldComm, M_directoryLibExpr ) );
            m_f[std::get<0>(f)].push_back( expr<d,1,2>( f.expression2(), "", M_worldComm, M_directoryLibExpr ) );
        }
        return std::move(m_f);
    }
    template<int d> map_matrix_field<d,d> getMatrixFields( std::initializer_list< std::pair<std::string,std::string > > const& listKeys ) const
    {
        using namespace Feel::vf;
        map_matrix_field<d,d> m_f;
        for ( auto const& it : listKeys )
        {
            auto curFields = getMatrixFields<d>( std::string(it.first), std::string(it.second) );
            for ( auto const& curField : curFields )
                m_f[curField.first] = std::move(curField.second);
        }
        return std::move(m_f);
    }
    template<int d> map_matrix_field<d,d> getMatrixFields( std::string && field, std::string && type )  const
    {
        using namespace Feel::vf;
        map_matrix_field<d,d> m_f;
        auto const& itFindField = this->find(field);
        if ( itFindField == this->end() ) return std::move(m_f);
        auto const& itFindType = itFindField->second.find(type);
        if ( itFindType == itFindField->second.end() ) return std::move(m_f);
        for ( auto f : itFindType->second )
        {
            LOG(INFO) << "Building expr " << f.expression() << " for " << std::get<0>(f);
            m_f[std::get<0>(f)] = expr<d,d,2>( f.expression(), "", M_worldComm, M_directoryLibExpr );
        }
        return std::move(m_f);
    }
  private:
    void setup();

    template <typename CastType>
    std::pair<bool,CastType> param( std::string const& field,std::string const& bc, std::string const& marker, std::string const& param, CastType const& defaultValue ) const;

  private:
    WorldComm const& M_worldComm;
    std::string M_prefix;
    pt::ptree M_pt;
    std::string M_directoryLibExpr;

};

using BoundaryConditionFactory = Singleton<BoundaryConditions>;
}
#endif
