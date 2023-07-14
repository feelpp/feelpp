// $Id: numeric_vector.h,v 1.11 2005/02/22 22:17:34 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __numeric_vector_h__
#define __numeric_vector_h__

#include <vector>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <feel/feelcore/disablewarnings.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/export.hpp>
#include <feel/feelcore/reenablewarnings.hpp>

#include <feel/feelcore/traits.hpp>

#include <feel/feelalg/datamap.hpp>
#include <feel/feelalg/glas.hpp>
#include <range/v3/view/take.hpp>
namespace Feel
{
namespace ublas = boost::numeric::ublas;

// forward declarations
template <typename T, typename SizeT> class Vector;
template <typename T> class MatrixSparse;
template <typename T> class MatrixShell;

/**
 * Numeric vector. Provides a uniform interface
 * to vector storage schemes for different linear
 * algebra libraries.
 *
 * @author Benjamin S. Kirk, 2003
 * @author Christophe Prud'homme 2005
 */
template <typename T,typename SizeT = uint32_type>
class FEELPP_EXPORT Vector : public std::enable_shared_from_this<Vector<T,SizeT> >
{
public:

    typedef T value_type;
    using size_type = SizeT;
    typedef typename type_traits<T>::real_type real_type;

    typedef Vector<T> self_type;
    typedef std::shared_ptr<Vector<T> > self_ptrtype;
    typedef std::shared_ptr<Vector<T> > clone_ptrtype;

    typedef DataMap<size_type> datamap_type;
    typedef std::shared_ptr<datamap_type> datamap_ptrtype;
    using vector_ptrtype = std::shared_ptr<Vector<T>>;

    /**
     *  Dummy-Constructor. Dimension=0
     */
    Vector ( worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() );

    Vector ( datamap_ptrtype const& n );

    /**
     * Constructor. Set dimension to \p n and initialize all elements with zero.
     */
    Vector ( const size_type n, worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() );

    /**
     * Constructor. Set local dimension to \p n_local, the global dimension
     * to \p n, and initialize all elements with zero.
     */
    Vector ( const size_type n,
             const size_type n_local,
             worldcomm_ptr_t const& _worldComm = Environment::worldCommPtr() );

    Vector ( Vector const& v );

    /**
     * Destructor, deallocates memory. Made virtual to allow
     * for derived classes to behave properly.
     */
    virtual ~Vector ();

    datamap_type const& map() const
    {
        return *M_map;
    }

    datamap_ptrtype const& mapPtr() const
    {
        return M_map;
    }

    void setMap( datamap_ptrtype const& d )
    {
        M_map=d;
    }

    /**
     * @returns true if the vector has been initialized,
     * false otherwise.
     */
    virtual bool isInitialized() const
    {
        return M_is_initialized;
    }

    /**
     * @returns true if the vector is closed and ready for
     * computation, false otherwise.
     */
    virtual bool closed() const
    {
        return M_is_closed;
    }

    /**
     * Call the assemble functions
     */
    virtual void close () = 0;

    /**
     * @returns the \p Vector<T> to a pristine state.
     */
    virtual void clear ();

    /**
     *
     */
    void localize(const Vector<T>& V);

    /**
     * set initialized only for subclasses
     */
    void setInitialized( bool b )
    {
        M_is_initialized = b;
    }

    /**
     * @ set false if the vector is in assembly state and need to be closed
     * for some next used (global operation) , false otherwise.
     */
    void setIsClosed( bool b )
    {
        M_is_closed = b;
    }

    /**
     * Set all entries to zero. Equivalent to \p v = 0, but more obvious and
     * faster.
     */
    virtual void zero () = 0;

    /**
     * Set entries to zero between \p start and \p stop
     */
    virtual void zero ( size_type /*start*/,  size_type /*stop*/ ) = 0;

    /**
     * set the entries to the constant \p v
     */
    virtual void setConstant( value_type v ) = 0;

    /**
     * set the entries to 0
     */
    virtual void setZero()
    {
        this->zero();
    }

    /**
     * set the entries to 1
     */
    virtual void setOnes()
    {
        this->setConstant( 1 );
    }

    /**
     * Replaces each component of a vector by its reciprocal.
     */
    virtual int reciprocal();

    /**
     * Creates a copy of this vector and returns it in an \p shared_ptr<>.
     * This must be overloaded in the derived classes.
     */
    virtual clone_ptrtype clone () const = 0;

    /**
     * Change the dimension of the vector to \p N. The reserved memory for
     * this vector remains unchanged if possible, to make things faster, but
     * this may waste some memory, so take this in the back of your head.
     * However, if \p N==0 all memory is freed, i.e. if you want to resize
     * the vector and release the memory not needed, you have to first call
     * \p init(0) and then \p init(N). This cited behaviour is analogous
     * to that of the STL containers.
     *
     * On \p fast==false, the vector is filled by
     * zeros.
     */

    virtual void init ( const size_type,
                        const size_type,
                        const bool = false );

    /**
     * call init with n_local = N,
     */
    virtual void init ( const size_type,
                        const bool = false );


    /**
     * call init with datamap,
     */
    virtual void init ( datamap_ptrtype const& dm )
    {
        M_is_closed = false;
        M_is_initialized = false;
        M_map = dm;
    }



    //   /**
    //    * Change the dimension to that of the
    //    * vector \p V. The same applies as for
    //    * the other \p init function.
    //    *
    //    * The elements of \p V are not copied, i.e.
    //    * this function is the same as calling
    //    * \p init(V.size(),fast).
    //    */
    //   virtual void init (const Vector<T>&,
    // 		     const bool = false) {}

    /**
     * \f$U(0-N) = s\f$: fill all components.
     */
    Vector<T> & operator= ( const T s );

    /**
     *  \f$U = V\f$: copy all components.
     */
    virtual Vector<T> & operator= ( const Vector<T> &V );

    /**
     *  \f$U = V\f$: copy all components.
     */
    Vector<T> & operator= ( const std::vector<T> &v );

    /**
     *  \f$v = x*y\f$: coefficient-wise multiplication
     */
    virtual void pointwiseMult ( Vector<T> const& x, Vector<T> const& y ) {}

    /**
     *  \f$v = x/y\f$: coefficient-wise divide
     */
    virtual void pointwiseDivide ( Vector<T> const& x, Vector<T> const& y ) {}

    /**
     * \return the sum of the components of the vector
     */
    virtual value_type sum() const = 0;

    /**
     * @returns the minimum element in the vector.
     * In case of complex numbers, this returns the minimum
     * Real part.
     */
    virtual real_type min () const = 0;

    /**
     * retrieve the min component as well as the index of the min component
     */
    //virtual real_type min( size_type& index ) const = 0;

    /**
     * @returns the maximum element in the vector.
     * In case of complex numbers, this returns the maximum
     * Real part.
     */
    virtual real_type max () const = 0;
    virtual real_type maxWithIndex( int* index=nullptr ) const;

    //! Replaces every element in a vector with its absolute value
    virtual void abs();

    /**
     * retrieve the max component as well as the index of the max component
     */
    //virtual real_type max( size_type& index ) const = 0;

    /**
     * @returns the \f$l_1\f$-norm of the vector, i.e.
     * the sum of the absolute values.
     */
    virtual real_type l1Norm () const = 0;

    /**
     * @returns the \f$l_2\f$-norm of the vector, i.e.
     * the square root of the sum of the
     * squares of the elements.
     */
    virtual real_type l2Norm () const = 0;

    /**
     * @returns the maximum absolute value of the
     * elements of this vector, which is the
     * \f$l_\infty\f$-norm of a vector.
     */
    virtual real_type linftyNorm () const = 0;

    /**
     * @returns dimension of the vector. This
     * function was formerly called \p n(), but
     * was renamed to get the \p Vector<T> class
     * closer to the C++ standard library's
     * \p std::vector container.
     */
    virtual size_type size () const
    {
        return M_map->nDof();
    }

    /**
     * @returns the local size of the vector
     * (index_stop-index_start)
     */
    virtual size_type localSize() const
    {
        return M_map->nLocalDofWithGhost();
    }

    /**
     * @returns the index of the first vector element
     * actually stored on this processor.  Hint: the
     * minimum for this index is \p 0.
     */
    virtual size_type firstLocalIndex() const
    {
        return M_map->minMyGID();
    }

    /**
     * @returns the index+1 of the last vector element
     * actually stored on this processor.  Hint: the
     * maximum for this index is \p size().
     */
    virtual size_type lastLocalIndex() const
    {
        return M_map->maxMyGID()+1;
    }

    virtual bool localIndexIsGhost(size_type localDof) const
    {
        return M_map->dofGlobalProcessIsGhost(localDof);
    }

    /**
     * \return the communicator
     */
    WorldComm const& comm() const
    {
        return M_map->comm();
    }
    /**
     * \return the communicator
     */
    worldcomm_ptr_t const& commPtr() const
    {
        return M_map->worldCommPtr();
    }
    /**
     * \return the communicator
     */
    worldcomm_ptr_t const& worldCommPtr() const
    {
        return M_map->worldCommPtr();
    }
    /**
     * \return the communicator
     */
    worldcomm_ptr_t & worldCommPtr() 
    {
        return M_map->worldCommPtr();
    }

    /**
     * Access components, returns \p U(i).
     */
    virtual T operator() ( const size_type i ) const = 0;

    virtual T& operator() ( const size_type i ) = 0;

    /**
     * Addition operator.
     * Fast equivalent to \p U.add(1, V).
     */
    virtual Vector<T> & operator += ( const Vector<value_type> &V ) = 0;

    /**
     * Subtraction operator.
     * Fast equivalent to \p U.add(-1, V).
     */
    virtual Vector<T> & operator -= ( const Vector<value_type> &V ) = 0;

    /**
     * v(i) = value
     */
    virtual void set ( const size_type i, const value_type& value ) = 0;

    /**
     * v([i1,i2,...,in]) = [value1,...,valuen]
     */
    virtual void setVector ( int* i, int n, value_type* v ) = 0;

    /**
     * v(i) += value
     */
    virtual void add ( const size_type i, const value_type& value ) = 0;

    /**
     * v([i1,i2,...,in]) += [value1,...,valuen]
     */
    virtual void addVector ( int* i, int n, value_type* v, size_type K, size_type K2 ) = 0;

    /**
     * \f$U(0-DIM)+=s\f$.
     * Addition of \p s to all components. Note
     * that \p s is a scalar and not a vector.
     */
    virtual void add ( const value_type& s ) = 0;

    /**
     * \f$U+=V\f$:
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    virtual void add ( const Vector<value_type>& V ) = 0;

    /**
     * \f$U+=a*V\f$.
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    virtual void add ( const value_type& a, const Vector<value_type>& v ) = 0;

    /**
     * \f$U+=a*V\f$.
     * Simple vector addition, equal to the
     * \p operator +=.
     */
    void add ( const value_type& a, const std::shared_ptr<Vector<value_type> >& v )
    {
        add( a, *v );
    }

    /**
     * @brief add a vector to a vector
     * 
     * @param a vector of coefficients to multiply the vectors
     * @param v vector of vectors to add
     */
    virtual void add ( const eigen_vector_type<Eigen::Dynamic,value_type>& a, const std::vector<vector_ptrtype>& v ) {}

    /**
     * @brief add to a vector  a given number of vectors
     * 
     * @param a 
     * @param v 
     * @param N 
     */
    virtual void add ( const eigen_vector_type<Eigen::Dynamic,value_type>& a, const std::vector<vector_ptrtype>& v, int N ) 
    {
        if ( N >= v.size() )
            throw std::invalid_argument( fmt::format("invalid number of vectors {} should be less than {}", N, v.size() ) ); 
        if ( a.size() < N )
            throw std::invalid_argument( fmt::format("invalid number of coefficients {} should be greater of equal to {}", a.size(), N ) );
        auto sv = v| ranges::views::take(N);
        this->add( a, {sv.begin(), sv.end()} );
    }

    /**
     * \f$ U+=v \f$ where v is a ublas::vector<T>
     * and you
     * want to specify WHERE to add it
     */
    virtual void addVector ( const std::vector<T>& v,
                             const std::vector<size_type>& dof_indices ) = 0;

    /**
     * \f$U+=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to add
     * the Vector<T> V
     */
    virtual void addVector ( const Vector<T>& V,
                             const std::vector<size_type>& dof_indices ) = 0;

    /**
     * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    virtual void addVector ( const Vector<T>& V_in,
                             const MatrixSparse<T>& A_in ) = 0;

    /**
     * \f$U+=A*V\f$, add the product of a \p SparseMatrix \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const std::shared_ptr<Vector<T> >& V_in,
                     const std::shared_ptr<MatrixSparse<T> >& A_in )
    {
        addVector( *V_in, *A_in );
    }

    /**
     * \f$U+=A*V\f$, add the product of a \p MatrixShell \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const Vector<T>& V_in,
                     const MatrixShell<T>& A_in );

    /**
     * \f$U+=A*V\f$, add the product of a \p MatrixShell \p A
     * and a \p Vector \p V to \p this, where \p this=U.
     */
    void addVector ( const std::shared_ptr<Vector<T> >& V_in,
                     const std::shared_ptr<MatrixShell<T> >& A_in );

#if 0
    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to add
     * the DenseVector<T> V
     */
    virtual void addVector ( const DenseVector<T>& V,
                             const std::vector<size_type>& dof_indices ) = 0;
#endif

    virtual value_type dot( Vector<T> const& v ) const = 0;
    virtual value_type dot( std::shared_ptr<Vector<T> > const& v ) const { return dot( *v ); }

    /**
     * @brief multiple dot products
     * 
     * @return eigen_vector_type<Eigen::Dynamic,value_type> 
     */
    virtual eigen_vector_type<Eigen::Dynamic,value_type> mDot( std::vector<vector_ptrtype> const& vs ) const  { return eigen_vector_type<Eigen::Dynamic,value_type>(); }

    virtual eigen_vector_type<Eigen::Dynamic,value_type> mDot( std::vector<vector_ptrtype> const& vs, int N ) const  
    {
        if ( N >= vs.size() )
            throw std::invalid_argument( fmt::format( "invalid number of vectors {} should be less than {}", N, vs.size() ) );
        auto sv = vs | ranges::views::take( N );
        return this->mDot( { sv.begin(), sv.end() } );
    }
   

    /**
     * \f$ U=v \f$ where v is a DenseVector<T>
     * and you want to specify WHERE to insert it
     */
    virtual void insert ( const std::vector<T>& v,
                          const std::vector<size_type>& dof_indices ) = 0;

    /**
     * \f$U=V\f$, where U and V are type
     * Vector<T> and you
     * want to specify WHERE to insert
     * the Vector<T> V
     */
    virtual void insert ( const Vector<T>& V,
                          const std::vector<size_type>& dof_indices ) = 0;

    /**
     * \f$ U+=V \f$ where U and V are type
     * DenseVector<T> and you
     * want to specify WHERE to insert
     * the DenseVector<T> V
     */
    virtual void insert ( const ublas::vector<T>& V,
                          const std::vector<size_type>& dof_indices ) = 0;

    /**
     * Scale each element of the
     * vector by the given factor.
     */
    virtual void scale ( const T factor ) = 0;

    /**
     * @returns \p -1 when \p this is equivalent to \p other_vector,
     * up to the given \p threshold.  When differences occur,
     * the return value contains the first index where
     * the difference exceeded the threshold.  When
     * no threshold is given, the \p Application \p TOLERANCE
     * is used.
     */
    virtual int compare ( const Vector<T> &other_vector,
                          const real_type threshold = 1e-10 ) const;



    /**
     * Prints the contents of the vector to the screen.
     */
    virtual void print( std::ostream& os=std::cout ) const;

    /**
     * Same as above but allows you to use stream syntax.
     */
    friend std::ostream& operator << ( std::ostream& os, const Vector<T,size_type>& v )
    {
        v.print( os );
        return os;
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    virtual void printMatlab( const std::string name="NULL", bool renumber = false ) const = 0;

    /**
     * Creates the subvector "subvector" from the indices in the
     * "rows" array.  Similar to the create_submatrix routine for
     * the SparseMatrix class, it is currently only implemented for
     * PetscVectors.
     */
    virtual void createSubvector( Vector<T>& ,
                                  const std::vector<size_type>& ) const
    {
        std::cerr << "ERROR: Not Implemented in base class yet!" << std::endl;
        FEELPP_ASSERT( 0 ).error( "invalid call" );
    }

    /**
     * copy vector entries in subvector ( subvector is already built from a createSubVector)
     * row indices given in the "rows" entries.
     */
    virtual
    void
    updateSubVector( std::shared_ptr<Vector<T> > & subvector,
                     std::vector<size_type> const& rows,
                     bool init=true )
        {
            CHECK( false ) << "invalid call : Not Implemented in base class";
        }

    /**
     * Creates the subvector "subvector" from the indices in the
     * "rows" array.  Similar to the create_submatrix routine for
     * the SparseMatrix class, it is currently only implemented for
     * PetscVectors.
     */
    virtual
    std::shared_ptr<Vector<T> >
    createSubVector( std::vector<size_type> const& rows,
                     bool checkAndFixRange=true ) const
        {
            CHECK( false ) << "invalid call : Not Implemented in base class";
            std::shared_ptr<Vector<T> > res;
            return res;
        }

    virtual void save( std::string const& filename="default_archive_name", std::string const& format="binary" )
    {}

    virtual void load( std::string const& filename="default_archive_name", std::string const& format="binary" )
    {}

protected:

    /**
     * Flag to see if the Numeric
     * assemble routines have been called yet
     */
    bool M_is_closed;

    /**
     * Flag to tell if init
     * has been called yet
     */
    bool M_is_initialized;

    /**
     * data distribution map of the vector over the processors
     */
    datamap_ptrtype M_map;
};

typedef Vector<double> vector_type;
typedef std::shared_ptr<vector_type> vector_ptrtype;
typedef std::unique_ptr<vector_type> vector_uptrtype;

/*----------------------- Inline functions ----------------------------------*/

/**
 * Computes the inner product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the inner product of \p v1 and \p v2
 */
template <typename T>
typename type_traits<T>::real_type
inner_product( Vector<T> const& v1, Vector<T> const& v2 )
{
    FEELPP_ASSERT( v1.localSize() == v2.localSize() &&
                   v1.size() == v2.size() )
    ( v1.localSize() )( v2.localSize() )
    ( v1.size() )( v2.size() ).error( "incompatible vector sizes" );

    return v1.dot( v2 );
}
/**
 * Computes the inner product of two vectors and eventually in parallel
 * \param v1 vector (eventually distributed)
 * \param v2 vector (eventually distributed)
 *
 * \return the inner product of \p v1 and \p v2
 */
template <typename T>
FEELPP_EXPORT typename type_traits<T>::real_type
inner_product( std::shared_ptr<Vector<T> > const& v1,
               std::shared_ptr<Vector<T> > const& v2 )
{
    return inner_product( *v1, *v2 );
}

template <typename T>
FEELPP_EXPORT typename type_traits<T>::real_type
dot( std::shared_ptr<Vector<T> > const& v1,
     std::shared_ptr<Vector<T> > const& v2 )
{
    return inner_product( *v1, *v2 );
}
template <typename T>
FEELPP_EXPORT typename type_traits<T>::real_type
dot( Vector<T> const& v1,
     Vector<T> const& v2 )
{
    return inner_product( v1, v2 );
}


namespace detail
{
template <class VectorType>
struct is_vector_ptr : mpl::false_ {};

template <class VectorType>
struct is_vector_ptr<std::shared_ptr<VectorType> >
        :
        boost::is_base_of<Vector<typename VectorType::value_type>,
        VectorType>
{};

template <typename T,typename SizeT=uint32_type>
struct FEELPP_EXPORT syncOperator
{
    using size_type = SizeT;
    typedef std::set<std::pair< rank_type, T > > storage_ghostdof_type;
    syncOperator() {}
    syncOperator( std::map<size_type, std::set<rank_type> > const& m )
        :
        M_activeDofClusterUsedByProc( m )
        {}
    syncOperator( syncOperator const& obj ) = default;

    virtual T operator()( size_type gcdof, rank_type activeProcId, T activeDofValue, storage_ghostdof_type const& ghostDofs ) const = 0;
    virtual bool hasOperator() const = 0;
    std::map<size_type, std::set<rank_type> > const& activeDofClusterUsedByProc() const { return M_activeDofClusterUsedByProc; }
    void setActiveDofClusterUsedByProc( std::map<size_type, std::set<rank_type> > const& m ) { M_activeDofClusterUsedByProc = m; }
private :
    std::map<size_type, std::set<rank_type> > M_activeDofClusterUsedByProc;

};

} // namespace detail

template <typename T,typename SizeT>
FEELPP_EXPORT void
sync( Vector<T,SizeT> & v, std::string const& opSyncStr = "=" );

template <typename T, typename SizeT>
FEELPP_EXPORT void
sync( Vector<T,SizeT> & v, std::string const& opSyncStr, std::set<SizeT> const& dofGlobalProcessPresent );

template <typename T,typename SizeT>
FEELPP_EXPORT void
sync( Vector<T,SizeT> & v, Feel::detail::syncOperator<T,SizeT> const& opSync );


} // Feel

#if 0
BOOST_SERIALIZATION_ASSUME_ABSTRACT(Feel::Vector)

namespace boost {
    namespace serialization {

    template<typename T, class Archive>
    FEELPP_EXPORT void save(Archive & ar, const Feel::Vector<T> & v, const unsigned int version)
    {
        Feel::DataMap const& map = v.map();
        ar & BOOST_SERIALIZATION_NVP(map);
    }
    template<typename T, class Archive>
    FEELPP_EXPORT void load(Archive & ar, Feel::Vector<T> & v, const unsigned int version)
    {
        Feel::DataMap map;
        ar & BOOST_SERIALIZATION_NVP(map);
        auto mapPtr = std::make_shared<Feel::DataMap>(map);
        v.init(mapPtr);
    }
    template<typename T, class Archive>
    FEELPP_EXPORT void serialize(Archive & ar, Feel::Vector<T> & v, const unsigned int version)
    {
        split_free( ar, v, version );
    }
    } // namespace serialization
} // namespace boost
#endif

#endif  // #ifdef __numeric_vector_h__
