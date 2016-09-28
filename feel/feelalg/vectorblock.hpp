/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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
#include <feel/feelalg/products.hpp>

namespace Feel
{


template<typename T> class Backend;



template <typename T=double>
class BlocksBaseVector : public vf::BlocksBase<boost::shared_ptr<Vector<T> > >
{
public :
    typedef vf::BlocksBase<boost::shared_ptr<Vector<T> > > super_type;
    typedef BlocksBaseVector<T> self_type;
    typedef Vector<T> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    typedef boost::shared_ptr<Backend<T> > backend_ptrtype;
    //using local_vector_type = Eigen::Matrix<value_type,Eigen::Dynamic,Eigen::Dynamic>;

    BlocksBaseVector(uint16_type nr = 0,
                     backend_ptr_t<T> b = backend() )
        :
        super_type(nr,1),
        M_backend( b )
    {}

    BlocksBaseVector(self_type const & b) = default;

    BlocksBaseVector(super_type const & b)
        :
        super_type(b),
        M_backend( backend() )
    {}

    template<typename PS>
    BlocksBaseVector( PS&& ps, backend_ptr_t<T> b = backend(),
                      std::enable_if_t<std::is_base_of<ProductSpacesBase,std::remove_reference_t<PS>>::value>* = nullptr )
        :
        super_type( ps.numberOfSpaces(), 1 ),
        M_backend( b )
        {
            int n = 0;
            hana::for_each( ps, [&]( auto const& e )
                            {

                                hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(e)>>{},
                                          [&]( auto&& x )
                                          {
                                              for( int i = 0; i < x->numberOfSpaces(); ++i, ++n )
                                              {
                                                  cout << "creating dyn vector block (" << n << "," << i  << ")\n";
                                                  (*this)(n,0) = (*x)[i]->elementPtr();
                                              }

                                          },
                                          [&]( auto&& x )
                                          {
                                              cout << "creating vector block (" << n  << ")\n";
                                              (*this)(n,0) = x->elementPtr();
                                              ++n;
                                          } )(e);

                            });


        }

    template<typename PS>
    BlocksBaseVector( PS&& ps, backend_ptr_t<T> b = backend(),
                      std::enable_if_t<std::is_base_of<ProductSpaceBase,std::remove_reference_t<PS>>::value>* = nullptr )
        :
        super_type( ps.numberOfSpaces(), 1 ),
        M_backend( b )
        {
            for( int i = 0; i < ps.numberOfSpaces(); ++i )
            {
                cout << "creating dyn vector block (" << i  << ")\n";
                (*this)(i,0) = ps[i]->elementPtr();
            }

        }

    /**
     * push_back method
     */
    self_type
    operator<<( vector_ptrtype const& m ) const
    {
        return super_type::operator<<( m );
    }

    /**
     * copy subblock and M_vector if init from global vector
     */
    void localize( vector_ptrtype const& vb, size_type _start_i=0 );

    /**
     * copy subblock from M_vector ( contained in this object )
     */
    void localize();

    /**
     * build vector representing all blocks
     */
    void buildVector( backend_ptrtype _backend = Feel::backend(_rebuild=false));

    /**
     * set values of VectorBlock from subvector
     */
    void setVector( vector_type & vec, vector_type const& subvec , int blockId ) const;
    /**
     * set values of subvector from global vector build with VectorBlock
     */
    void setSubVector( vector_type & subvec, vector_type const& vec , int blockId ) const;

    vector_ptrtype& vector() { return M_vector; }
    vector_ptrtype const& vector() const { return M_vector; }

    /**
     * termination function to fill
     */
    void fill( int r  ) {}

    /**
     * termination function to fillWithNewVector
     */
    void fillWithNewVector( int r  ) {}

    /**
     * recurse through entries of subvector
     */
    template<typename Arg1, typename ...Args>
    void
    fillWithNewVector( int n, const boost::shared_ptr<Arg1>& arg1, const boost::shared_ptr<Args>&... args )
        {
            this->operator()( n ) = M_backend->newBlockVector( arg1 );
            // do submatrix (n+1,n+1)
            fillWithNewVector( ++n, args... );
        }

    /**
     * recurse through entries of subvector
     */
    template<typename Arg1, typename ...Args>
    void
    fill( int n, const boost::shared_ptr<Arg1>& arg1, const boost::shared_ptr<Args>&... args )
        {
            this->operator()( n ) = arg1;
            // do submatrix (n+1,n+1)
            fill( ++n, args... );
        }


private :
    vector_ptrtype M_vector;
    backend_ptr_t<T> M_backend;
    //std::unordered_map<size_type, local_vector_type> M_local_vec;
};

template <int NR, typename T=double>
class BlocksVector : public BlocksBaseVector<T>
{
public :
    static const uint16_type NBLOCKROWS = NR;
    static const uint16_type NBLOCKCOLS = 1;

    typedef BlocksBaseVector<T> super_type;

    BlocksVector()
        :
        super_type(NBLOCKROWS)
    {}
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
    typedef vf::BlocksBase<vector_ptrtype> blocksbase_type;

    VectorBlock(  blocksbase_type const & blockVec,
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


/**
 * Build blocks of CSR graphs with function spaces \p
 * (args1,args2,argn). Variadic template is used to handle an arbitrary number of
 * function spaces.
 * The blocks are organized then matrix wise with the stencil associated of pairs of function spaces in \p (arg1,...,argn)
 *
 */
template<typename Arg1, typename ...Args>
BlocksBaseVector<typename decay_type<Arg1>::value_type>
newVectorBlocks( const Arg1& arg1, const Args&... args )
{
    const int size = sizeof...(Args)+1;
    BlocksBaseVector<typename decay_type<Arg1>::value_type> g( size, backend() );
    int n = 0;
    g.fillWithNewVector( n, arg1, args... );
    return g;
}


/**
 * Build blocks of CSR graphs with function spaces \p
 * (args1,args2,argn). Variadic template is used to handle an arbitrary number of
 * function spaces.
 * The blocks are organized then matrix wise with the stencil associated of pairs of function spaces in \p (arg1,...,argn)
 *
 */
template<typename Arg1, typename ...Args>
BlocksBaseVector<typename decay_type<Arg1>::value_type>
vectorBlocks( const Arg1& arg1, const Args&... args,
              std::enable_if_t<std::is_base_of<ProductSpacesBase,std::remove_reference_t<Arg1>>::value>* = nullptr )
{
    const int size = sizeof...(Args)+1;
    BlocksBaseVector<typename decay_type<Arg1>::value_type> g( size, backend() );
    int n = 0;
    g.fill( n, arg1, args... );
    return g;
}

template<typename Arg>
BlocksBaseVector<typename decay_type<Arg>::value_type>
vectorBlocks( const Arg& arg,
              std::enable_if_t<std::is_base_of<ProductSpaceBase,std::remove_reference_t<Arg>>::value>* = nullptr )
{
    BlocksBaseVector<typename decay_type<Arg>::value_type> g( arg.numberOfSpaces(), backend() );
    for( int i = 0; i < arg.numberOfSpaces(); ++i )
        g(i) = arg.elementPtr();
    return g;
}

/**
 * Build blocks of CSR graphs with function spaces \p
 * (args1,args2,argn). Variadic template is used to handle an arbitrary number of
 * function spaces.
 * The blocks are organized then matrix wise with the stencil associated of pairs of function spaces in \p (arg1,...,argn)
 *
 */
template<typename PS>
//BlocksBaseVector<typename decay_type<PS>::value_type>
BlocksBaseVector<double>
blockVector( PS && ps, backend_ptrtype b = backend(),
             std::enable_if_t<std::is_base_of<ProductSpacesBase,std::remove_reference_t<PS>>::value>* = nullptr )
{
    const int size = ps.numberOfSpaces();
    //BlocksBaseVector<typename decay_type<PS>::value_type> g( size, backend() );
    BlocksBaseVector<double> g( size, b );

    int n = 0;
    hana::for_each( ps, [&]( auto const& e )
                    {

                        hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(e)>>{},
                                  [&] (auto&& x) {
                                      for( int i = 0; i < x->numberOfSpaces(); ++i, ++n  )
                                      {
                                          cout << "creating dyn vector block (" << n  << ")\n";
                                          g(n,0) = b->newVector( (*x)[i] );
                                      }
                                  },
                                  [&] (auto&& x){
                                      cout << "creating vector block (" << n  << ")\n";
                                      g(n++,0) = b->newVector( x );
                                  })(e);
                    });
    return g;
}

template<typename PS>
BlocksBaseVector<double>
blockVector( PS && ps, backend_ptrtype b = backend(),
             std::enable_if_t<std::is_base_of<ProductSpaceBase,std::remove_reference_t<PS>>::value>* = nullptr )
{
    BlocksBaseVector<double> g( ps.numberOfSpaces(), b );

    for( int i = 0; i < ps.numberOfSpaces(); ++i )
        g(i,0) = b->newVector( ps[i] );
    return g;
}

template<typename PS>
//BlocksBaseVector<typename decay_type<PS>::value_type>
BlocksBaseVector<double>
blockElement( PS && ps, backend_ptrtype b = backend(),
              std::enable_if_t<std::is_base_of<ProductSpacesBase,std::remove_reference_t<PS>>::value>* = nullptr )
{
    const int size = hana::size(ps);
    //BlocksBaseVector<typename decay_type<PS>::value_type> g( size, backend() );
    BlocksBaseVector<double> g( size, b );

    int n = 0;
    hana::for_each( ps, [&]( auto const& e )
                    {
                        cout << "creating vector element (" << n  << ")\n";
                        g(n,0) = e->elementPtr();
                        ++n;
                    });
    return g;
}

template<typename PS>
//BlocksBaseVector<typename decay_type<PS>::value_type>
BlocksBaseVector<double>
blockElement( PS && ps, backend_ptrtype b = backend(),
              std::enable_if_t<std::is_base_of<ProductSpaceBase,std::remove_reference_t<PS>>::value>* = nullptr )
{
    BlocksBaseVector<double> g( ps.numberOfSpaces(), b );

    int n = 0;
    for( int i = 0;i < ps.numberOfSpaces(); ++i )
        g(i,0) = ps[i]->elementPtr();
    return g;
}

} // Feel

#endif /* __VectorBlock_H */
