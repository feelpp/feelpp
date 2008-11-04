/* -*- mode: c++ -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-10-11

  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file gmshsimplexdomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-11
 */
#ifndef __GmshSimplexDomain_H
#define __GmshSimplexDomain_H 1

#include <life/lifefilters/gmsh.hpp>
#include <life/lifemesh/simplex.hpp>

namespace Life
{
/**
 * \class GmshSimplexDomain
 * \brief Simplex Domain description for gmsh mesh generation
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
template<int Dim, int Order>
class GmshSimplexDomain : public Gmsh
{
    typedef Gmsh super;
public:


    /** @name Constants and Typedefs
     */
    //@{
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef Simplex<Dim,Order, Dim> entity_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GmshSimplexDomain()
        :
        super(),
        _M_I( nDim ),
        _M_h( 0.1 ),
        _M_descr()
    {
        if ( nDim >= 1 )
            _M_I[0] = std::make_pair( 0, 1 );
        if ( nDim >= 2 )
            _M_I[1] = std::make_pair( 0, 1 );
        if ( nDim >= 3 )
            _M_I[2] = std::make_pair( 0, 1 );
        this->setOrder( (GMSH_ORDER) nOrder );
    }

    GmshSimplexDomain( GmshSimplexDomain const & td )
        :
        super( td ),
        _M_I( td._M_I ),
        _M_h( td._M_h ),
        _M_descr( td._M_descr )
    {
        this->setOrder( (GMSH_ORDER) nOrder );
    }
    ~GmshSimplexDomain()
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


    //@}

    /** @name  Mutators
     */
    //@{

    void setX( std::pair<double,double> const& x )
    {
        LIFE_ASSERT( nDim >= 1 )( nDim ).error( "invalid dimension" );
        _M_I[0] = x;
    }
    void setY( std::pair<double,double> const& y )
    {
        LIFE_ASSERT( nDim >= 2 )( nDim ).error( "invalid dimension" );
        _M_I[1] = y;
    }
    void setZ( std::pair<double,double> const& z )
    {
        LIFE_ASSERT( nDim >= 3 )( nDim ).error( "invalid dimension" );
        _M_I[2] = z;
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
    {return getDescription( mpl::int_<nDim>() );}
    // 1D
    std::string getDescription( mpl::int_<1> ) const;
    // 2D
    std::string getDescription( mpl::int_<2> ) const;
    // 3D
    std::string getDescription( mpl::int_<3> ) const;

private:

    std::vector<std::pair<double,double> > _M_I;
    double _M_h;
    std::string _M_descr;


};

} // Life
#endif /* __GmshSimplexDomain_H */
