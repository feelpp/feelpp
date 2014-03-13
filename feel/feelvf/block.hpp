/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-01-26

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2006-2012 Universite Joseph Fourier (Grenoble I)

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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file block.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-01-26
 */
#ifndef __Block_H
#define __Block_H 1

#include <sstream>
#include <list>

#include <feel/feelcore/feel.hpp>

namespace Feel
{

    //template<typename T> class MatrixSparse;

namespace vf
{
/// \cond detail
/*!
  \class Block
  \brief class that describes a matrix or vector block

  @author Christophe Prud'homme
  @see
*/
class Block
{
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    Block( uint16_type __ic = 0, uint16_type __jc = 0, size_type __gic = 0, size_type __gjc = 0 )
        :
        M_lr( __ic ),
        M_lc( __jc ),
        M_gr( __gic ),
        M_gc( __gjc )
    {}
    Block( Block const & __b )
        :
        M_lr( __b.M_lr ),
        M_lc( __b.M_lc ),
        M_gr( __b.M_gr ),
        M_gc( __b.M_gc )
    {}
    ~Block()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    Block& operator=( Block const& __b )
    {
        if ( FEELPP_ISLIKELY( this != &__b ) )
        {
            M_lr = __b.M_lr;
            M_lc = __b.M_lc;
            M_gr = __b.M_gr;
            M_gc = __b.M_gc;
        }

        return *this;
    }


    //@}

    /** @name Accessors
     */
    //@{

    uint16_type localRow() const
    {
        return M_lr;
    }
    uint16_type localColumn() const
    {
        return M_lc;
    }

    size_type globalRowStart() const
    {
        return M_gr;
    }
    size_type globalColumnStart() const
    {
        return M_gc;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    void setLocalRow( uint16_type __lr )
    {
        M_lr = __lr;
    }
    void setLocalColumn( uint16_type __lc )
    {
        M_lc = __lc;
    }

    void setGlobalRowStart( size_type __r )
    {
        M_gr = __r;
    }
    void setGlobalColumnStart( size_type __c )
    {
        M_gc = __c;
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}

protected:

private:

    uint16_type M_lr;
    uint16_type M_lc;
    size_type M_gr;
    size_type M_gc;
};
typedef std::list<Block> list_block_type;

inline
std::ostream&
operator<<( std::ostream& __os, Block const& __b )
{
    __os << "Block [ "
         << __b.localRow() << ","
         << __b.localColumn() << ","
         << __b.globalRowStart() << ","
         << __b.globalColumnStart() << "]";

    return __os;
}
/// \endcond


template <typename T>
struct BlocksBase
{
    typedef T block_type;
    typedef uint16_type index_type;
    BlocksBase()
        :
        M_nRow( 0 ),
        M_nCol( 0 ),
        M_vec( 1 ),
        M_cptToBuild( 0 )
    {}

    BlocksBase( index_type nr,index_type nc )
        :
        M_nRow( nr ),
        M_nCol( nc ),
        M_vec( nr*nc ),
        M_cptToBuild( 0 )
    {}

    BlocksBase( index_type nr,index_type nc,block_type const& a )
        :
        M_nRow( nr ),
        M_nCol( nc ),
        M_vec( nr*nc, a ),
        M_cptToBuild( 0 )
    {}

    BlocksBase( BlocksBase<T> const& b )
        :
        M_nRow( b.M_nRow ),
        M_nCol( b.M_nCol ),
        M_vec( b.M_vec ),
        M_cptToBuild( b.M_cptToBuild )
    {}

    BlocksBase<T>
    operator<<( block_type const& m )  const
    {
        BlocksBase<T> newBlock( *this );
        newBlock.M_vec[M_cptToBuild]=m;
        ++( newBlock.M_cptToBuild );
        return newBlock;
    }

    void
    push_back( block_type const& m )
    {
        M_vec[M_cptToBuild]=m;
        ++M_cptToBuild;
    }
#if 0
    void
    operator<<( block_type const& m )
    {
        M_vec[M_cptToBuild]=m;
        ++M_cptToBuild;
        //return *this;
    }
#endif

    block_type &
    operator()( index_type c1,index_type c2=0 )
    {
        return M_vec[c1*M_nCol+c2];
    }

    block_type
    operator()( index_type c1,index_type c2=0 ) const
    {
        return M_vec[c1*M_nCol+c2];
    }

    std::vector<block_type> const&
    getSetOfBlocks() const
    {
        return M_vec;
    }

    index_type nRow() const
    {
        return M_nRow;
    }
    index_type nCol() const
    {
        return M_nCol;
    }

    void reset()
    {
        M_cptToBuild=0;
        M_vec.clear();
        M_vec.resize( this->nRow()*this->nCol() );
    }

    void resize( index_type nr,index_type nc=1)
    {
        M_nRow=nr;
        M_nCol=nc;
        M_vec.resize( nr*nc );
    }

    void merge( index_type c1,index_type c2,BlocksBase<T> const& b )
    {
        const index_type nRb = b.nRow(), nCb = b.nCol();
        for ( index_type i=0; i<nRb; ++i )
            for ( index_type j=0; j<nCb; ++j )
                this->operator()( c1+i,c2+j ) = b( i,j );
    }

    BlocksBase<T> subBlock(index_type x1r,index_type x1c,index_type x2r,index_type x2c) const
    {
        const index_type nrow = x2r-x1r+1;
        const index_type ncol = x2c-x1c+1;
        BlocksBase<T> res(nrow,ncol);
        for (index_type kr=0;kr<nrow;++kr)
            for (index_type kc=0;kc<ncol;++kc)
                res(kr,kc)=this->operator()(x1r+kr,x1c+kc);
        return res;
    }

    BlocksBase<T> transpose() const
    {
        BlocksBase<T> res(this->nCol(),this->nRow());
        for (index_type i=0;i<res.nRow();++i)
            for (index_type j=0;j<res.nCol();++j)
                res(i,j) = this->operator()(j,i);
        return res;
    }

private :
    index_type M_nRow,M_nCol;
    std::vector<block_type> M_vec;
    index_type M_cptToBuild;
};



template <int NR, int NC, typename T>
struct Blocks : public BlocksBase<T>
{
    typedef BlocksBase<T> super_type;
    typedef typename super_type::index_type index_type;
    static const index_type NBLOCKROWS = NR;
    static const index_type NBLOCKCOLS = NC;

    typedef T block_type;

    Blocks()
        :
        super_type( NR,NC )
    {}

    Blocks( block_type const& a )
        :
        super_type( NR,NC, a )
    {}

    Blocks( super_type const& a )
        :
        super_type( a )
    {}

    Blocks<NR,NC,T>
    operator<<( block_type const& m ) const
    {
        Blocks<NR,NC,T> newBlock = super_type::operator<<( m ) ;
        return newBlock;
    }

    void
    push_back( block_type const& m )
    {
        super_type::push_back( m );
    }

#if 0
    void
    operator<<( block_type const& m )
    {
        super_type::operator<<( m );
    }
#endif

};

} // vf

class BlocksStencilPattern : public vf::BlocksBase<size_type>
{
    typedef vf::BlocksBase<size_type> super_type;
    typedef super_type::index_type index_type;
    typedef BlocksStencilPattern self_type;
public :

    BlocksStencilPattern(index_type nr,index_type nc)
        :
        super_type(nr,nc)
    {}

    BlocksStencilPattern(index_type nr,index_type nc, size_type pat)
        :
        super_type(nr,nc,pat)
    {}

    BlocksStencilPattern(super_type const & b)
        :
        super_type(b)
    {}

    self_type
    operator<<( size_type const& m ) const
    {
        return super_type::operator<<( m );
    }

};


} // feel
#endif /* __Block_H */
