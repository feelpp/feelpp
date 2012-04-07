/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2005-01-26
 */
#ifndef __Block_H
#define __Block_H 1

#include <sstream>
#include <list>

#include <feel/feelcore/feel.hpp>
#include <feel/feelalg/matrixsparse.hpp>

namespace Feel
{
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
        _M_lr( __ic ),
        _M_lc( __jc ),
        _M_gr( __gic ),
        _M_gc( __gjc )
        {}
    Block( Block const & __b )
        :
        _M_lr( __b._M_lr ),
        _M_lc( __b._M_lc ),
        _M_gr( __b._M_gr ),
        _M_gc( __b._M_gc )
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
                _M_lr = __b._M_lr;
                _M_lc = __b._M_lc;
                _M_gr = __b._M_gr;
                _M_gc = __b._M_gc;
            }
            return *this;
        }


    //@}

    /** @name Accessors
     */
    //@{

    uint16_type localRow() const { return _M_lr; }
    uint16_type localColumn() const { return _M_lc; }

    size_type globalRowStart() const { return _M_gr; }
    size_type globalColumnStart() const { return _M_gc; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setLocalRow( uint16_type __lr ) { _M_lr = __lr; }
    void setLocalColumn( uint16_type __lc ) { _M_lc = __lc; }

    void setGlobalRowStart( size_type __r ) { _M_gr = __r; }
    void setGlobalColumnStart( size_type __c ) { _M_gc = __c; }

    //@}

    /** @name  Methods
     */
    //@{


    //@}

protected:

private:

    uint16_type _M_lr;
    uint16_type _M_lc;
    size_type _M_gr;
    size_type _M_gc;
};
typedef std::list<Block> list_block_type;

inline
DebugStream&
operator<<( DebugStream& __os, Block const& __b )
{
    std::ostringstream __str;
    __str << "Block [ "
          << __b.localRow() << ","
          << __b.localColumn() << ","
          << __b.globalRowStart() << ","
          << __b.globalColumnStart() << "]";

    __os << __str.str() << "\n";
    return __os;
}

inline
NdebugStream&
operator<<( NdebugStream& __os, Block const& /*__b*/ )
{
    return __os;
}

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


template <typename T= boost::shared_ptr< MatrixSparse<double> > >
struct BlocksBase
{
    typedef T block_type;
    //typedef boost::shared_ptr<block_type> block_ptrtype;
    //typedef MatrixSparse<T> matrix_type;
    //typedef boost::shared_ptr<matrix_type> matrix_ptrtype;
    BlocksBase()
        :
        M_nRow(0),
        M_nCol(0),
        M_vec(1),
        M_cptToBuild(0)
    {}

    BlocksBase(uint16_type nr,uint16_type nc)
        :
        M_nRow(nr),
        M_nCol(nc),
        M_vec(nr*nc),
        M_cptToBuild(0)
    {}

    BlocksBase(uint16_type nr,uint16_type nc,block_type /*const&*/ a)
        :
        M_nRow(nr),
        M_nCol(nc),
        M_vec(nr*nc, a),
        M_cptToBuild(0)
    {}

    BlocksBase(BlocksBase<T> const& b)
        :
        M_nRow(b.M_nRow),
        M_nCol(b.M_nCol),
        M_vec(b.M_vec),
        M_cptToBuild(b.M_cptToBuild)
    {}

    BlocksBase<T>
    operator<<(block_type const& m)  const
    {
        BlocksBase<T> newBlock(*this);
        newBlock.M_vec[M_cptToBuild]=m;
        ++(newBlock.M_cptToBuild);
        return newBlock;
    }

    void
    push_back(block_type const& m)
    {
        M_vec[M_cptToBuild]=m;
        ++M_cptToBuild;
    }
#if 0
    void
    operator<<(block_type const& m)
    {
        M_vec[M_cptToBuild]=m;
        ++M_cptToBuild;
        //return *this;
    }
#endif

    block_type &
    operator()(uint16_type c1,uint16_type c2)
    {
        return M_vec[c1*M_nCol+c2];
    }

    block_type
    operator()(uint16_type c1,uint16_type c2) const
    {
        return M_vec[c1*M_nCol+c2];
    }

    std::vector<block_type>
    getSetOfBlocks() const { return M_vec; }

    uint16_type nRow() const { return M_nRow; }
    uint16_type nCol() const { return M_nCol; }

    void reset()
    {
        M_cptToBuild=0;
        M_vec.clear();
        M_vec.resize(this->nRow()*this->nCol());
    }

    void merge(uint16_type c1,uint16_type c2,BlocksBase<T> /*const&*/ b)
    {
        uint16_type nRb = b.nRow(), nCb = b.nCol();
        for (uint16_type i=0;i<nRb;++i)
            for (uint16_type j=0;j<nCb;++j)
                this->operator()(c1+i,c2+j) = b(i,j);
    }

private :
    uint16_type M_nRow,M_nCol;
    std::vector<block_type> M_vec;
    uint16_type M_cptToBuild;
};



template <int NR, int NC, typename T= boost::shared_ptr< MatrixSparse<double> > >
struct Blocks : public BlocksBase<T>
{
    static const uint16_type NBLOCKROWS = NR;
    static const uint16_type NBLOCKCOLS = NC;

    typedef BlocksBase<T> super_type;
    typedef T block_type;

    Blocks()
        :
        super_type(NR,NC)
    {}

    Blocks(block_type const& a)
        :
        super_type(NR,NC, a)
    {}

    Blocks(super_type const& a)
        :
        super_type(a)
    {}

    Blocks<NR,NC,T>
    operator<<(block_type const& m) const
    {
        Blocks<NR,NC,T> newBlock = super_type::operator<<(m) ;
        return newBlock;
    }

    void
    push_back(block_type const& m)
    {
        super_type::push_back(m);
    }

#if 0
    void
    operator<<(block_type const& m)
    {
        super_type::operator<<(m);
    }
#endif

};








} // vf
} // feel
#endif /* __Block_H */
