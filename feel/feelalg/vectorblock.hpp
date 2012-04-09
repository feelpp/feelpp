/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@imag.fr>
       Date: 2012-01-18

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file matrixblock.hpp
   \author Vincent Chabannes <vincent.chabannes@imag.fr>
   \date 2012-01-18
 */

#ifndef __VectorBlock_H
#define __VectorBlock_H 1

#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelvf/block.hpp>


namespace Feel
{


template<typename T> class Backend;

template <int NR, typename T=double>
class BlocksVector : public vf::Blocks<NR,1,boost::shared_ptr<Vector<T> > >
{
    typedef vf::Blocks<NR,1,boost::shared_ptr<Vector<T> > > super_type;
public :
    typedef BlocksVector<NR,T> self_type;

    typedef Vector<T> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    template <typename VectorType>
    self_type
    operator<<( VectorType const& m ) const
    {
        //typedef typename mpl::or_<is_shared_ptr<VectorType>, boost::is_pointer<VectorType> >::type is_ptr_or_shared_ptr;
        self_type newBlock( *this );
        typedef is_shared_ptr<VectorType> is_shared_ptr;
        newBlock.push_back( m,is_shared_ptr() );
        return newBlock;
    }

    template <typename VectorType>
    void
    push_back( VectorType const& m, mpl::bool_<true> )
    {
        super_type::push_back( m );
    }

#if 0 //Don't work
    template <typename VectorType>
    void push_back( VectorType const& m, mpl::bool_<false> )
    {
        //self_type::operator<<(m.shared_from_this());
        auto m2=m.shared_from_this();
        self_type::operator<<( m2 );
    }
#endif

    void localize( vector_ptrtype const& vb )
    {
        size_type _start_i=0;

        for ( uint16_type i=0; i<this->nRow(); ++i )
        {
            size_type nBlockRow = this->operator()( i,0 )->size();

            for ( size_type k=0; k<nBlockRow; ++k )
            {
                this->operator()( i,0 )->set( k, vb->operator()( _start_i+k ) );
            }

            _start_i += nBlockRow;
        }
    }



};

/**
 * \class VectorBlock
 * \brief block of vector
 *
 * <code>
 * auto myBlocks = BlocksVector<2,2>() << F1 << F2
 *
 * auto A = backend->newBlockVector(myBlocks);
 * </code>
 *
 * @author Vincent Chabannes
 */

template< typename T>
class VectorBlockBase
{
    //typedef Vector<T> super;
public:

    /** @name Typedefs
     */
    //@{

    typedef VectorBlockBase<T> self_type;
    typedef T value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef Vector<T> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    VectorBlockBase( vf::BlocksBase<vector_ptrtype > const & blockVec,
                     backend_type &backend,
                     bool copy_values=true );


    VectorBlockBase( VectorBlockBase const & vb )
        :
        M_vec( vb.M_vec )
    {}

    ~VectorBlockBase()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    VectorBlockBase operator=( VectorBlockBase const& vb )
    {
        if ( this != &vb )
        {
            M_vec = vb.M_vec;
        }

        return *this;
    }

    //@}

    /** @name Accessors
     */
    //@{
    vector_ptrtype getVector()
    {
        return M_vec;
    }
    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    void updateBlockVec( vector_ptrtype const & m, size_type start_i );

    //@}



protected:

private:

    vector_ptrtype M_vec;
};


template<int NR, typename T>
class VectorBlock : public VectorBlockBase<T>
{
    typedef VectorBlockBase<T> super_type;

public:

    static const uint16_type NBLOCKROWS = NR;

    typedef typename super_type::value_type value_type;
    typedef typename super_type::vector_ptrtype vector_ptrtype;
    typedef typename super_type::backend_type backend_type;
    typedef vf::Blocks<NBLOCKROWS,1,vector_ptrtype > blocks_type;

    VectorBlock(  blocks_type const & blockVec,
                  backend_type &backend,
                  bool copy_values=true )
        :
        super_type( blockVec,backend,copy_values )
    {}

    VectorBlock( VectorBlock const & vb )
        :
        super_type( vb )
    {}

    VectorBlock operator=( VectorBlock const& vb )
    {
        super_type::operator=( vb );
        return *this;
    }

    VectorBlock & operator = ( vector_ptrtype const& F )
    {
        super_type::operator=( F );
        return *this;
    }

}; // VectorBlock


} // Feel

#endif /* __VectorBlock_H */
