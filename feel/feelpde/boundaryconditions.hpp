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

    ExpressionStringAtMarker( super && s ) : super( s ) {}
    
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

    BoundaryConditions();
    
    /**
     * constructor from an \c initializer_list<>
     */
    BoundaryConditions( std::initializer_list<value_type> l )
        : super(l),M_prefix() {}
    
    BoundaryConditions( std::string const& prefix );
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
    void setPrefix( std::string p ) { M_prefix = p; }

    void setPTree( pt::ptree const& p );
        
    /**
     * load property tree from file \p filename 
     */
    void load(const std::string &filename);

    /**
     * retrieve scalar field \p field with boundary conditions of type \p type
     */
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
            m_f[std::get<0>(f)] = expr<Order>( f.expression() );
        }
        return std::move(m_f);
    }
    /**
     * retrieve scalar field pair \p field with boundary conditions of type \p type
     */
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
                LOG(INFO) << "Building expr " << f.expression() << " for " << f.marker();
                m_f[f.marker()].push_back( expr<Order>( f.expression1() ) );
                m_f[f.marker()].push_back( expr<Order>( f.expression2() ) );
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
            m_f[std::get<0>(f)] = expr<d,1,2>( f.expression() );
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
            m_f[std::get<0>(f)] = expr<d,d,2>( f.expression() );
        }
        return std::move(m_f);
    }
  private:
    void setup();
  private:

    std::string M_prefix;
    pt::ptree M_pt;
};

using BoundaryConditionFactory = Singleton<BoundaryConditions>;
}
#endif
