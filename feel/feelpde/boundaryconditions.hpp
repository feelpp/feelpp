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
#include <feel/feelfilters/loadcsv.hpp>

#include <boost/property_tree/ptree.hpp>

namespace Feel
{
namespace pt =  boost::property_tree;

//! Boundary condition at marker is an e≈pression string
class ExpressionStringAtMarker : public std::tuple<std::string,std::string,std::string,std::string,std::string>
{
public:
    typedef std::tuple<std::string,std::string,std::string,std::string,std::string> super;
    enum boundarycondition_t { EXPRESSION = 0, FILE = 1 };
    
    ExpressionStringAtMarker( super && s )
        :
        super( s ),
        M_type( EXPRESSION )
        
        {
            M_meshMarkers.push_back( this->marker() );
            if ( typeStr() == "file" )
            {
                M_type = FILE;
                M_filename = Environment::expand( std::get<2>( *this ) );
                M_data = loadXYFromCSV( M_filename, std::get<3>( *this ), std::get<4>( *this ) );
            }
        }

    //! type of boundary condition : expression or data
    std::string typeStr() const { return std::get<0>( *this ); }

    //! return type
    boundarycondition_t type() const { return M_type; }

    //! @return true if boundary condition is an expression, false otherwise
    bool isExpression() const { return M_type == EXPRESSION; }

    //! @return true
    bool isFile() const { return M_type == FILE; }
    
    /**
     * @return the marker
     */
    std::string const& marker() const { return std::get<1>( *this ); }
    
    /**
     * @return the expression
     */
    std::string const& expression() const { return std::get<2>( *this ); }
    
    /**
     * @return the expression
     */
    std::string const& expression1() const { return std::get<2>( *this ); }

    /**
     * @return the expression
     */
    std::string const& expression2() const { return std::get<3>( *this ); }

    bool hasExpression() const { return isExpression() && !std::get<2>( *this ).empty(); } 
    bool hasExpression1() const { return isExpression() && !std::get<2>( *this ).empty(); } 
    bool hasExpression2() const { return isExpression() && !std::get<3>( *this ).empty(); }

    std::string filename() const { LOG_IF( ERROR, !isFile() ) << "boundary condition is not given by a file"; return M_filename; }
    bool hasFilename() const { return isFile() && !M_filename.empty(); }
    
    std::list<std::string> const& meshMarkers() const { return M_meshMarkers; }

    void setMeshMarkers( std::list<std::string> const& s ) { M_meshMarkers=s; }

    double data( double time, double epsilon = 1e-7 ) const;
private :
    
    std::list<std::string> M_meshMarkers;
    boundarycondition_t M_type;
    std::string M_filename;
    std::map<double,double> M_data;
};
// works on both `const` and non-`const` associative containers:
template<class Container>
inline
auto floating_equal_range( Container&& container, double target, double epsilon = 0.00001 )
    -> decltype( container.equal_range(target) )
{
    auto lower = container.lower_bound( target-epsilon );
    auto upper = container.upper_bound( target+epsilon );
    return std::make_pair(lower, upper);
}
template<typename Container>
inline
auto floating_key_exists( Container const& container, double target, double epsilon = 0.00001 )
    -> decltype( std::make_pair( true, floating_equal_range( container, target, epsilon ) ) )
{
    auto range = floating_equal_range(container, target, epsilon);
    return std::make_pair( range.first != range.second, range );
}
inline double
ExpressionStringAtMarker::data( double time, double epsilon ) const
{
    auto r = floating_key_exists( M_data, time, epsilon );
    LOG(INFO) << "Look for " << time << " in data file, found: " << r.first;
    if ( r.first )
        return r.second.first->second;
    throw std::logic_error( "invalid time "+std::to_string(time)+" not found in data file " + this->filename() );
}
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
            m_f[f.marker()] = expr<Order>( f.expression(), "", M_worldComm, M_directoryLibExpr );
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
            LOG(INFO) << "Building expr " << f.expression() << " for " << f.marker();
            m_f[f.marker()] = expr<d,1,2>( f.expression(), "", M_worldComm, M_directoryLibExpr );
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
            LOG(INFO) << "Building expr " << f.expression() << " for " << f.marker();
            m_f[f.marker()].push_back( expr<d,1,2>( f.expression1(), "", M_worldComm, M_directoryLibExpr ) );
            m_f[f.marker()].push_back( expr<d,1,2>( f.expression2(), "", M_worldComm, M_directoryLibExpr ) );
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
            LOG(INFO) << "Building expr " << f.expression() << " for " << f.marker();
            m_f[f.marker()] = expr<d,d,2>( f.expression(), "", M_worldComm, M_directoryLibExpr );
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
