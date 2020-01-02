/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 Dec 2019

 Copyright (C) 2019 Feel++ Consortium

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
#ifndef FEELPP_DOFRELATION_HPP
#define FEELPP_DOFRELATION_HPP 1

#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <feel/feeldiscr/dof.hpp>
namespace Feel {

template<typename FeT, typename IndexT=uint32_type>
class BimapDofRelation
{
public:
    using fe_t = FeT;
    static constexpr uint16_type nDofComponents() { return fe_t::is_product?fe_t::nComponents:1; }
    typedef Dof globaldof_type;
    typedef LocalDof<nDofComponents()> localdof_type;
    typedef boost::bimap<bimaps::set_of<localdof_type>, bimaps::multiset_of<Dof> > dof_table;
    typedef typename dof_table::value_type dof_relation;
    typedef typename dof_table::left_iterator local_dof_iterator;
    typedef typename dof_table::left_const_iterator local_dof_const_iterator;
    typedef typename dof_table::right_iterator global_dof_iterator;
    typedef typename dof_table::right_const_iterator global_dof_const_iterator;

    BimapDofRelation() = default;

    bool hasLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.left.find( ldof  );
            if ( eit == M_el_l2g.left.end() )
                return false;
            else
                return true;
        }
    auto findLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.left.find( ldof  );
            return std::pair(eit,eit!=M_el_l2g.left.end());
        }
    std::optional<globaldof_type> findGlobalDofFromLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.left.find( ldof  );
            if ( eit!=M_el_l2g.left.end())
                return eit->second;
            else
                return std::nullopt;
        }
    auto insertDofRelation( localdof_type const& ldof, globaldof_type const& gdof )
        {
            auto res = M_el_l2g.insert( dof_relation( ldof, gdof ) );
            return res;
        }
    template<typename Iterator>
    bool modifyDofRelation( Iterator it, globaldof_type const& gdof )
        {
            return M_el_l2g.left.modify_data( it, bimaps::_data = gdof );
        }

    void renumberDofRelation( std::vector<size_type> const& previousGlobalIdToNewGlobalId )
        {
            for( auto it = M_el_l2g.left.begin(), en = M_el_l2g.left.end(); it != en; ++it )
            {
                auto const& previousGDof=it->second;
                Dof newGDof( previousGDof );
                newGDof.setIndex( previousGlobalIdToNewGlobalId[previousGDof.index()] );
                bool successfulModify = M_el_l2g.left.modify_data( it, boost::bimaps::_data = newGDof );
                CHECK( successfulModify ) << "modify global dof id fails";
            }
        } 


    std::pair<global_dof_const_iterator,global_dof_const_iterator>  globalDof()  const
        {
            return std::make_pair( M_el_l2g.right.begin(), M_el_l2g.right.end() );
        }
    std::pair<global_dof_const_iterator,global_dof_const_iterator> globalDof( size_type GlobalDofId ) const
        {
            auto lower = M_el_l2g.right.lower_bound( globaldof_type(GlobalDofId,-1) );
            auto upper = M_el_l2g.right.upper_bound( globaldof_type(GlobalDofId,2) );
            return std::make_pair( lower, upper );
        }
    /**
     * \return the specified entries of the globalToLocal table
     *
     * \param DofId the Dof ID
     *
     * \return the element id and local dof id
     */
    localdof_type const& globalToLocal( size_type dof )  const
        {
            auto it = M_el_l2g.right.find( Dof( dof ) );
            DCHECK( it != M_el_l2g.right.end() ) << "Invalid global dof entry ( " << dof << ")";
            return it->second;
        }

    std::pair<local_dof_const_iterator,local_dof_const_iterator> localDof() const
        {
            return std::make_pair( M_el_l2g.left.begin(), M_el_l2g.left.end() );
        }

    /**
     * @code
     * for( auto d:localDof( elid )
     * {
     *   // do something with d.first (local dof) or d.second (global dof)
     * )
     * @endcode
     * @return the pair of iterators [beg,end) of Dof associated to element \p ElId
     */
    std::pair<local_dof_const_iterator,local_dof_const_iterator> localDof( size_type ElId ) const
        {
            auto lower = M_el_l2g.left.lower_bound( localdof_type(ElId) );
            auto upper = M_el_l2g.left.upper_bound( localdof_type(ElId,invalid_uint16_type_value) );
            //DCHECK( it.first != M_el_l2g.left.end() ) << "Invalid element dof entry " << ElId;
            return std::make_pair( lower, upper );
        }
    /**
     * @return the number of Dof associated to element \p ElId
     */
    int localDofSize( size_type ElId ) const
        {
            auto const& [ beg,end] = localDof( ElId );
            return std::distance( beg,end );
        }
    globaldof_type const& localToGlobalDof( const size_type ElId,
                                            const uint16_type localNode,
                                            const uint16_type c = 0 ) const 
        {
            auto it = M_el_l2g.left.find( localdof_type(ElId,fe_t::nLocalDof * c  + localNode ) );
            DCHECK( it != M_el_l2g.left.end() ) << "Invalid dof entry ( " << ElId << ", " << fe_t::nLocalDof * c  + localNode << ")";
            //DCHECK( it->second.index() < nDof() && nDof() > 0 ) << "Invalid Dof Entry: " << it->second.index() << " > " << this->nDof();
            return it->second;
        }


private:
    dof_table M_el_l2g;
    
};


template<typename FeT, typename IndexT=uint32_type>
class MapDofRelation
{
public:
    using fe_t = FeT;
    static constexpr uint16_type nDofComponents() { return fe_t::is_product?fe_t::nComponents:1; }
    typedef Dof globaldof_type;
    typedef LocalDof<nDofComponents()> localdof_type;
    typedef std::unordered_map<IndexT,std::vector<Dof>> dof_table;
    typedef typename dof_table::value_type dof_relation;
    typedef typename dof_table::iterator local_dof_iterator;
    typedef typename dof_table::const_iterator local_dof_const_iterator;
    typedef typename dof_table::iterator global_dof_iterator;
    typedef typename dof_table::const_iterator global_dof_const_iterator;

    MapDofRelation() = default;
    ~MapDofRelation() = default;
    
    bool hasLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.find( ldof.elementId()  );
            if ( eit == M_el_l2g.end() || eit->second.at( ldof.localDof() ) == globaldof_type(-1) )
                return false;
            else
                return true;
        }
    auto findLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.find( ldof  );
            if ( eit == M_el_l2g.end() || eit->second.at( ldof.localDof() ) == globaldof_type(-1) )
                return std::pair(eit,false);
            else
                return std::pair(std::begin(eit->second)+ldof.localDof(),true);
        }
    std::optional<globaldof_type> findGlobalDofFromLocalDof( localdof_type const& ldof ) const
        {
            auto eit = M_el_l2g.find( ldof.elementId()  );
            if ( eit!=M_el_l2g.end())
                return eit->second[ldof.localDof()];
            else
                return std::nullopt;
        }
    auto insertDofRelation( localdof_type const& ldof, globaldof_type const& gdof )
        {
            auto eit = M_el_l2g.find(ldof.elementId()); 
            if ( eit == M_el_l2g.end() )
            {
                M_el_l2g[ldof.elementId()].resize( fe_t::nLocalDof*nDofComponents(), Dof(-1) );
                eit = M_el_l2g.find(ldof.elementId());
            }
            eit->second.at(ldof.localDof()) = gdof;
            return std::pair(eit,true);
        }
    template<typename Iterator>
    bool modifyDofRelation( Iterator it, globaldof_type const& gdof )
        {
            *it = gdof;
            return true;
        }

    void renumberDofRelation( std::vector<size_type> const& previousGlobalIdToNewGlobalId )
        {
            for( auto it = M_el_l2g.begin(), en = M_el_l2g.end(); it != en; ++it )
            {
                for( auto& gdof : it->second )
                {
                    gdof.setIndex( previousGlobalIdToNewGlobalId[gdof.index()] );
                }
            }
        } 


    std::pair<global_dof_const_iterator,global_dof_const_iterator>  globalDof()  const
        {
            return std::make_pair( M_el_l2g.begin(), M_el_l2g.end() );
        }
    std::multimap<globaldof_type,localdof_type> globalDof( size_type GlobalDofId ) const
        {
            std::multimap<globaldof_type,localdof_type> m;
            
            for( auto const& el: M_el_l2g )
            {
                int ldof = 0;
                for ( auto const& gdof: el.second )
                {
                    if ( gdof.index() == GlobalDofId )
                        m.insert( std::pair{gdof,localdof_type(el.first,ldof)} );
                    ldof++;
                }
            } 
            return m;
        }
    /**
     * \return the specified entries of the globalToLocal table
     *
     * \param DofId the Dof ID
     *
     * \return the element id and local dof id
     */
    localdof_type const& globalToLocal( size_type the_gdof )  const
        {
            CHECK(false) << "call globalToLocal error " << the_gdof;
        }

    std::map<localdof_type,globaldof_type> localDof() const
        {
            std::map<localdof_type,globaldof_type> m;
            
            for( auto const& el: M_el_l2g )
            {
                int ldof = 0;
                for ( auto const& gdof: el.second )
                    m.insert( std::pair{localdof_type(el.first,ldof++), gdof} );

            } 
            return m;
        }

    std::map<localdof_type,globaldof_type> localDof( size_type ElId ) const
        {
            auto const&vdof = M_el_l2g.at( ElId );
            std::map<localdof_type,globaldof_type> m;
            //m.reserve( vdof.size() );
            int i = 0;
            std::for_each( vdof.begin(), vdof.end(), [&m,ElId,&i]( auto l ) { m.insert( std::pair{localdof_type(ElId,i++), Dof(l) } ); } );
            return m;
        }
    /**
     * @return the number of Dof associated to element \p ElId
     */
    int localDofSize( size_type ElId ) const
        {
            return M_el_l2g.at( ElId ).size();
        }
    globaldof_type const& localToGlobalDof( const size_type ElId,
                                            const uint16_type localNode,
                                            const uint16_type c = 0 ) const 
        {
            return M_el_l2g.at( ElId ).at( fe_t::nLocalDof * c  + localNode );
        }


private:
    dof_table M_el_l2g;
    
};


}
#endif
