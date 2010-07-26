/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2006-11-26

  Copyright (C) 2006 Université Joseph Fourier

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file gmshtensorizeddomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2006-11-26
 */
#ifndef __GmshTensorizedDomain_H
#define __GmshTensorizedDomain_H 1

#include <life/lifefilters/gmsh.hpp>

namespace Life
{
/**
 * \class GmshTensorizedDomain
 * \brief Tensorized Domain description for gmsh mesh generation
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
class GmshTensorizedDomain : public Gmsh
{
    typedef Gmsh super;
public:


    /** @name Constants and Typedefs
     */
    //@{
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RDim;

    typedef Entity<Dim,Order, nRealDim> entity_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GmshTensorizedDomain()
        :
        super(),
        _M_I( nRealDim ),
        _M_h( 0.1 ),
        _M_descr()
    {
        if ( nRealDim >= 1 )
            _M_I[0] = std::make_pair( 0, 1 );
        if ( nRealDim >= 2 )
            _M_I[1] = std::make_pair( 0, 1 );
        if ( nRealDim >= 3 )
            _M_I[2] = std::make_pair( 0, 1 );
        this->setOrder( (GMSH_ORDER) nOrder );
    }

    GmshTensorizedDomain( GmshTensorizedDomain const & td )
        :
        super( td ),
        _M_I( td._M_I ),
        _M_h( td._M_h ),
        _M_descr( td._M_descr )
    {
        this->setOrder( (GMSH_ORDER) nOrder );
    }
    ~GmshTensorizedDomain()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return bounding box
     */
    std::vector<std::pair<double,double> > const& boundingBox() const {return _M_I;}
    /**
     * \return characteristic length
     */
    double const& h() const {return _M_h; }

    /**
     * \return the geometry description
     */
    std::string description() const { return this->getDescription(); }

    //@}

    /** @name  Mutators
     */
    //@{

    void setX( std::pair<double,double> const& x )
    {
        LIFE_ASSERT( nRealDim >= 1 )( nDim )( nRealDim ) .error( "invalid dimension" );
        _M_I[0] = x;
    }
    void setY( std::pair<double,double> const& y )
    {
        LIFE_ASSERT( nRealDim >= 2 )( nDim )( nRealDim ).error( "invalid dimension" );
        _M_I[1] = y;
    }
    void setZ( std::pair<double,double> const& z )
    {
        LIFE_ASSERT( nRealDim >= 3 )( nDim )( nRealDim ).error( "invalid dimension" );
        _M_I[2] = z;
    }
    void setReferenceDomain()
    {
        if ( nDim >= 1 )
            _M_I[0] = std::make_pair( -1, 1 );
        if ( nDim >= 2 )
            _M_I[1] = std::make_pair( -1, 1 );
        if ( nDim >= 3 )
            _M_I[2] = std::make_pair( -1, 1 );
    }

    void setCharacteristicLength( double h ) { _M_h = h; }

    //@}

    /** @name  Methods
     */
    //@{

    std::string generate( std::string const& name ) const
    {
        std::string descr = getDescription();
        return super::generate( name, descr );
    }

    //@}



private:
    std::string getDescription() const
    {return getDescription( mpl::int_<nDim>(), mpl::bool_<entity_type::is_simplex_product>() );}
    // 1D
    std::string getDescription( mpl::int_<1>, mpl::bool_<false> ) const;
    std::string getDescription( mpl::int_<1>, mpl::bool_<true> ) const
    { return getDescription( mpl::int_<1>(), mpl::bool_<false>() ); }
    // 2D
    std::string getDescription( mpl::int_<2>, mpl::bool_<false> ) const;
    std::string getDescription( mpl::int_<2>, mpl::bool_<true> ) const;
    // 3D
    std::string getDescription( mpl::int_<3>, mpl::bool_<false>, bool do_recombine = false ) const;
    std::string getDescription( mpl::int_<3>, mpl::bool_<true> ) const;

private:

    std::vector<std::pair<double,double> > _M_I;
    double _M_h;
    std::string _M_descr;
};

} // Life

#if !defined( LIFE_INSTANTIATION_MODE )
# include <life/lifefilters/gmshtensorizeddomain.cpp>
#endif // LIFE_INSTANTIATION_MODE

#endif /* __GmshTensorizedDomain_H */
