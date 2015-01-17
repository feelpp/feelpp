/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-13

  Copyright (C) 2005,2006 EPFL

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
   \file matrixublas.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-13
 */
#ifndef __MatrixUBlas_H
#define __MatrixUBlas_H 1

#include <set>
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/bindings/traits/traits.hpp>
//#include <boost/numeric/bindings/traits/ublas_sparse.hpp>


namespace Feel
{
namespace ublas = boost::numeric::ublas;

/*!
 * \class MatrixUBlas
 * \brief interface to ublas sparse matrix
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename T, typename LayoutType>
class MatrixUBlas
{
public:


    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef ublas::compressed_matrix<value_type, LayoutType> matrix_type;

    typedef typename boost::numeric::bindings::traits::sparse_matrix_traits<matrix_type>::ordering_type ordering_type;
    typedef typename boost::numeric::bindings::traits::sparse_matrix_traits<matrix_type>::layout_type layout_type;

    static const bool is_row_major = boost::is_same<ordering_type,
                      boost::numeric::bindings::traits::row_major_t>::value;


    typedef std::vector<std::set<size_type> > pattern_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    MatrixUBlas()
        :
        M_is_initialized( false ),
        M_mat()
    {}
    MatrixUBlas( MatrixUBlas const & m )
        :
        M_is_initialized( m.M_is_initialized ),
        M_mat( m.M_mat )
    {}

    ~MatrixUBlas()
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    value_type& operator()( size_type i, size_type j )
    {
        return M_mat( i, j );
    }

    value_type const& operator()( size_type i, size_type j ) const
    {
        return M_mat( i, j );
    }

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    unsigned int size1 () const
    {
        return M_mat.size1();
    }

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    unsigned int size2 () const
    {
        return M_mat.size2();
    }

    /**
     * \return the number of non-zeros entries in the matrix
     */
    size_type nnz() const
    {
        return M_mat.nnz();
    }

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    unsigned int rowStart () const
    {
        return 0;
    }

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    unsigned int rowStop () const
    {
        return 0;
    }

    /**
     * \return true if matrix is initialized/usable, false otherwise
     */
    bool isInitialized() const
    {
        return M_is_initialized;
    }

    /**
     * Returns the raw ublas matrix.
     */
    matrix_type const& mat () const
    {
        return M_mat;
    }

    /**
     * Returns the raw ublas matrix.
     */
    matrix_type & mat ()
    {
        return M_mat;
    }



    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Initialize a Ublas matrix that is of global
     * dimension \f$ m \times  n \f$ with local dimensions
     * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
     * nonzeros per row (defaults to 30).
     * \p noz is the number of on-processor
     * nonzeros per row (defaults to 30).
     */
    void init ( const unsigned int m,
                const unsigned int n,
                const unsigned int m_l,
                const unsigned int n_l,
                const unsigned int nnz=30,
                const unsigned int noz=10 );

    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear ()
    {
        M_mat.clear();
    }

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero ()
    {
        M_mat = ublas::zero_matrix<value_type>( M_mat.size1(), M_mat.size2() );
    }

    void zero ( size_type start1, size_type stop1, size_type start2, size_type stop2 )
    {
        ublas::subrange( M_mat, start1, stop1, start2, stop2 )  = ublas::zero_matrix<value_type>( stop1-start1, stop2-start2 );
    }

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void add ( const unsigned int i,
               const unsigned int j,
               const value_type value )
    {
        M_mat( i, j ) += value;
    }

    /**
     * Set \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void set ( const unsigned int i,
               const unsigned int j,
               const value_type value )
    {
        M_mat( i, j ) = value;
    }


    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL" ) const;


    /**
     * fill sparse matrix with non zero entries
     */
    void fill( pattern_type const& );

    void resize( size_type nr, size_type nc, bool preserve = false )
    {
        M_mat.resize( nr, nc, preserve );
    }

    /**
     * \return \f$ v^T M u \f$
     */
    template<typename VE1, typename VE2>
    value_type
    energy( ublas::vector_expression<VE1> const& __v,
            ublas::vector_expression<VE2> const& __u ) const
    {
        return ublas::inner_prod( __v, ublas::prod( M_mat, __u ) );
    }

    /**
     *
     */
    void diagonalize( size_type );
    //@}



protected:

private:

    bool M_is_initialized;

    /**
     * the ublas sparse matrix data structure
     */
    matrix_type M_mat;

};
template <typename T, typename LayoutType>
void MatrixUBlas<T,LayoutType>::init ( const unsigned int m,
                                       const unsigned int n,
                                       const unsigned int /*m_l*/,
                                       const unsigned int /*n_l*/,
                                       const unsigned int /*nnz*/,
                                       const unsigned int /*noz*/ )
{
    if ( ( m==0 ) || ( n==0 ) )
        return;

    M_mat.resize( m,n,false );
    this->zero ();
}

template<typename T, typename LayoutType>
void
MatrixUBlas<T, LayoutType>::diagonalize( size_type __dof_index )
{
    // eliminating row
    ublas::matrix_row<matrix_type> mr ( M_mat, __dof_index );
    typedef typename ublas::matrix_row<matrix_type>::iterator r_it;

    for ( r_it __r = mr.begin(); __r != mr.end(); ++__r )
    {
        *__r = 0.0;
    }

    // eliminating column
    ublas::matrix_column<matrix_type> mc ( M_mat, __dof_index );

    for ( typename ublas::matrix_column<matrix_type>::iterator therow = mc.begin();
            therow != mc.end(); ++therow )
    {
        *therow = 0.0;
    }

    // 1 on the diagonal
    M_mat( __dof_index, __dof_index ) = 1.0;
}
template<typename T, typename LayoutType>
void
MatrixUBlas<T, LayoutType>::fill( pattern_type const& __pattern )
{
    namespace bindings = boost::numeric::bindings;

    boost::timer chrono;
    typename pattern_type::const_iterator __itl = __pattern.begin();
    size_type __nnz = 0;
#if 0
    std::for_each( __pattern.begin(),
                   __pattern.end(),
                   ( boost::lambda::var( __nnz ) += boost::lambda::bind( &std::set<size_type>::size,
                           boost::lambda::_1 ) ) );
#else

    for ( size_type __i = 0; __i < __pattern.size(); ++__i )
    {
        __nnz += __pattern[__i].size();
        //         std::cout << "line " << __i << "\n";
        //         std::for_each( __pattern[__i].begin(), __pattern[__i].end(), std::cout << lambda::_1 << " " );
        //         std::cout << "\n";
    }

#endif

    //     ublas::unbounded_array<value_type> __val( __nnz );
    //     std::for_each( __val.begin(), __val.end(), boost::lambda::_1 = 0.0 );

    FEELPP_ASSERT( __nnz >= M_mat.nnz() )( __nnz )( M_mat.nnz() ).error( "incompatible sizes" );

    DVLOG(2) << "number of nnz in old M : " << M_mat.nnz() << ", " << M_mat.nnz_capacity() <<"\n";
    DVLOG(2) << "size M.value_data() :  " << M_mat.value_data().size() << "\n";
    //     DVLOG(2) << "           size val :  " << __val.size() << "\n";
    // save current nonzero entrie of M in the vector val
    // std::copy( M_mat.value_data().begin(), M_mat.value_data().begin()+M_mat.nnz(), __val.begin() );

    //std::cout << "values=";
    // std::for_each( __val.begin(), __val.end(), std::cout << boost::lambda::_1 << "\n" );
    //std::cout << "\n";

    matrix_type M_mat_backup( M_mat );

    DVLOG(2) << "resizing M old : " << M_mat_backup.size1() << "," << M_mat_backup.size2() << " nnz = " << M_mat_backup.nnz() <<"\n";

    M_mat.reserve( __nnz, false );

    DVLOG(2) << "resizing M new : " << M_mat.size1() << "," << M_mat.size2() << " nnz = " << M_mat.nnz() <<"\n";

    std::set<size_type>::const_iterator __it;
    std::set<size_type>::const_iterator __en;

    DVLOG(2) << "*** counting Nnz  : " << chrono.elapsed() <<"\n";
    //DVLOG(2) << "*** l size  : " << __pattern.size() <<"\n";
    DVLOG(2) << "*** nnz size  : " << __nnz <<"\n";
    chrono.restart();

    uint32_type __max_nnz_per_line = 0;
    uint __row = 0;
    uint __nnz_entry = 0;

    size_type __filled1 = 1;
    size_type __filled2 = 0;

    size_type thesize;

    if ( is_row_major )
        thesize = M_mat.size1();

    else
        thesize = M_mat.size2();

    while ( __row < thesize )
    {
        __it = __pattern[__row].begin();
        __en = __pattern[__row].end();

        M_mat.index1_data()[__row] = __nnz_entry;

        ++__filled1;
        uint32_type __nnz_line = 0;

        while ( __it != __en )
        {
            M_mat.index2_data()[__nnz_entry] = *__it;

            ++__filled2;

            namespace bindings = boost::numeric::bindings;

            typename matrix_type::value_type* __pv = 0;

            if ( boost::is_same<typename bindings::traits::sparse_matrix_traits<matrix_type>::ordering_type,
                    bindings::traits::row_major_t>::value )
            {
                if ( __row < M_mat_backup.size1() && *__it < M_mat_backup.size2() )
                    __pv = M_mat_backup.find_element( __row, *__it );
            }

            else if ( boost::is_same<typename bindings::traits::sparse_matrix_traits<matrix_type>::ordering_type,
                      bindings::traits::column_major_t>::value )
            {
                if ( __row < M_mat_backup.size2() && *__it < M_mat_backup.size1() )
                    __pv = M_mat_backup.find_element( *__it, __row );
            }

            else
            {
                std::cout << "ERROR " << __FILE__ << ": " << __LINE__ << "\n";
            }

            if ( __pv )
                M_mat.value_data()[__nnz_entry] = *__pv;

            else
                M_mat.value_data()[__nnz_entry] =  value_type( 0 );

            ++__nnz_entry;
            ++__it;
            ++__nnz_line;
        }

        __max_nnz_per_line = std::max( __nnz_line, __max_nnz_per_line );
        ++__row;
    }

    M_mat.index1_data()[thesize] = __filled2;
    FEELPP_ASSERT( thesize+1 == __filled1 )( thesize )( __filled1 ).error( "invalid matrix storage" );
    M_mat.set_filled( __filled1, __filled2 );
    FEELPP_ASSERT( M_mat.nnz() == __filled2 )( M_mat.nnz() )( __filled2 ).error( "inconsistent matrix storage" );

    DVLOG(2) << "***  value data size  : " << M_mat.value_data().size() << "\n";
    DVLOG(2) << "***              nnz  : " << M_mat.nnz() << "\n";
    DVLOG(2) << "*** max nnz per line  : " << __max_nnz_per_line << "\n";
    DVLOG(2) << "*** fillMatrixFromPattern() done in " << chrono.elapsed() <<"s\n";
}

template<typename T, typename LayoutType>
void
MatrixUBlas<T, LayoutType>::printMatlab( const std::string filename ) const
{
    std::string name = filename;
    std::string separator = " , ";

    // check on the file name
    int i = filename.find( "." );

    if ( i <= 0 )
        name = filename + ".m";

    else
    {
        if ( ( unsigned int ) i != filename.size() - 2 ||
                filename[ i + 1 ] != 'm' )
        {
            std::cerr << "Wrong file name ";
            name = filename + ".m";
        }
    }

    std::ofstream file_out( name.c_str() );

    FEELPP_ASSERT( file_out )( filename ).error( "[Feel::spy] ERROR: File cannot be opened for writing." );

    file_out << "var_"<<filename<<" = [ ";

    for ( typename matrix_type::const_iterator1 i1=M_mat.begin1();
            i1!=M_mat.end1(); ++i1 )
    {
        for ( typename matrix_type::const_iterator2 i2=i1.begin();
                i2!=i1.end(); ++i2 )
            file_out << i2.index1() + 1 << separator
                     << i2.index2() + 1 << separator
                     << *i2  << std::endl;
    }

    file_out << "];" << std::endl;
    file_out << "I=var_"<<filename<<"(:,1); J=var_"<<filename<<"(:,2); var_"<<filename<<"=var_"<<filename<<"(:,3);" << std::endl;
    file_out << "A=sparse(I,J,var_"<<filename<<"); spy(A);" << std::endl;
}

} // Feel
#endif /* __MatrixUBlas_H */
