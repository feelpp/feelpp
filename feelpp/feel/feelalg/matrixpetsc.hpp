/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file matrixpetsc.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#ifndef __MatrixPetsc_H
#define __MatrixPetsc_H 1

#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <feel/feelconfig.h>

#if defined(FEELPP_HAS_PETSC_H)

#include <Eigen/Core>

#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/graphcsr.hpp>


extern "C"
{
#include <petscmat.h>
}

#ifndef USE_COMPLEX_NUMBERS
extern "C" {
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
}
#else
# include <petscversion.h>
# if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 1)
#   include <petscsles.h>
# else
#   include <petscksp.h>
# endif
#endif



namespace Feel
{
template<typename T> class VectorPetsc;
template<typename T> class VectorPetscMPI;
/**
 * \class MatrixPetsc
 * \brief Wrapper for petsc matrices
 *
 * Petsc matrix. Provides a nice interface to the
 * Petsc C-based data structures for parallel,
 * sparse matrices.
 *
 * @author Benjamin S. Kirk, 2002
 * @author Christophe Prud'homme
 * @see
 */
template<typename T>
class FEELPP_EXPORT MatrixPetsc : public MatrixSparse<T>
{
    typedef MatrixSparse<T> super;
public:
    /** @name Typedefs
     */
    //@{

    static const bool is_row_major = true;

    //@{

    /** @name Typedefs
     */
    //@{

    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;
    using size_type = typename super::size_type;
    typedef std::vector<std::set<size_type> > pattern_type;

    typedef typename super::graph_type graph_type;
    typedef typename super::graph_ptrtype graph_ptrtype;

    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;

    typedef typename super::indexsplit_type indexsplit_type;
    typedef typename super::indexsplit_ptrtype indexsplit_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Constructor; initializes the matrix to
     * be empty, without any structure, i.e.
     * the matrix is not usable at all. This
     * constructor is therefore only useful
     * for matrices which are members of a
     * class. All other matrices should be
     * created at a point in the data flow
     * where all necessary information is
     * available.
     *
     * You have to initialize
     * the matrix before usage with
     * \p init(...).
     */
    explicit MatrixPetsc( worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );
    MatrixPetsc( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol );
    MatrixPetsc( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, worldcomm_ptr_t const& worldComm );


    /**
     * Constructor.  Creates a PetscMatrix assuming you already
     * have a valid Mat object.  In this case, m is NOT destroyed
     * by the PetscMatrix destructor when this object goes out of scope.
     * This allows ownership of m to remain with the original creator,
     * and to simply provide additional functionality with the PetscMatrix.
     */
    MatrixPetsc ( Mat m );
    MatrixPetsc ( Mat m, datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, bool destroyMatOnExit=false );
    MatrixPetsc ( MatrixSparse<value_type> const& M, IS& isrow, IS& iscol );
    MatrixPetsc ( MatrixSparse<value_type> const& M, std::vector<int> const& rowIndex, std::vector<int> const& colIndex );
    /**
     * Destructor. Free all memory, but do not
     * release the memory of the sparsity
     * structure.
     */
    ~MatrixPetsc() override;

    //@}

    /** @name Operator overloads
     */
    //@{

    /**
     * Return the value of the entry \p (i,j).  This may be an
     * expensive operation and you should always take care where to
     * call this function.  In order to avoid abuse, this function
     * throws an exception if the required element does not exist in
     * the matrix.
     *
     * In case you want a function that returns zero instead (for
     * entries that are not in the sparsity pattern of the matrix),
     * use the \p el function.
     */
    value_type operator () ( const size_type i,
                             const size_type j ) const override;

    /**
     * create a PETSc matrix which is a copy of \p M
     * \param M the matrix to copy
     */
    MatrixPetsc& operator=( MatrixSparse<value_type> const& M ) override;

    using clone_ptrtype = typename super::clone_ptrtype;
    /**
     * @brief clone the petsc matrix
     * 
     * @return clone_ptrtype 
     */
    clone_ptrtype clone() const override;

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns \p m, the row-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size1 () const override;

    /**
     * @returns \p n, the column-dimension of
     * the matrix where the marix is \f$ M \times N \f$.
     */
    size_type size2 () const override;

    /**
     * @returns \p (m,n) the shape of the matrix
     */
    std::pair<size_t,size_t> shape () const { return std::pair<size_t,size_t> (size1() , size2()); }

    /**
     * return row_start, the index of the first
     * matrix row stored on this processor
     */
    size_type rowStart () const override;

    /**
     * return row_stop, the index of the last
     * matrix row (+1) stored on this processor
     */
    size_type rowStop () const override;

    //!
    //! @return the number of non-zero entries
    //!
    size_type nnz() const override;

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Initialize a Petsc matrix that is of global
     * dimension \f$ m \times  n \f$ with local dimensions
     * \f$ m_l \times n_l \f$.  \p nnz is the number of on-processor
     * nonzeros per row (defaults to 30).
     * \p noz is the number of on-processor
     * nonzeros per row (defaults to 30).
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                const size_type nnz=30,
                const size_type noz=10 ) override;

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                graph_ptrtype const& graph ) override;

    /**
     *
     */
    void setIndexSplit( indexsplit_ptrtype const& indexSplit ) override;

    /**
     * reinitialize the matrix
     */
    void resize( size_type m, size_type n, bool /*preserve*/ )
    {
        this->setInitialized( false );
        init( m, n, m, n );
    }

    void fill( pattern_type const& patt )
    {}
    /**
     * Release all memory and return
     * to a state just like after
     * having called the default
     * constructor.
     */
    void clear () override;

    /**
     * Set all entries to 0. This method retains
     * sparsity structure.
     */
    void zero () override;

    /**
     * Set all entries to 0 in the range
     * [start1-stop1,start2-stop2]. This method retains sparsity
     * structure.
     */
    void zero ( size_type start1, size_type stop1, size_type start2, size_type stop2 ) override;

    /**
     * Call the Petsc assemble routines.
     * sends necessary messages to other
     * processors
     */
    void close () const override;

    /**
     * see if Petsc matrix has been closed
     * and fully assembled yet
     */
    bool closed() const override;

    /**
     * Return the l1-norm of the matrix, that is
     * \f$|M|_1=max_{all columns j}\sum_{all rows i} |M_ij|\f$, (max. sum of columns).
     *
     * This is the natural matrix norm that is compatible to the
     * l1-norm for vectors, i.e.  \f$|Mv|_1\leq |M|_1 |v|_1\f$.
     */
    real_type l1Norm () const override;

    /**
     * Return the linfty-norm of the matrix, that is
     *
     * \f$|M|_\infty=max_{all rows i}\sum_{all columns j} |M_ij|\f$,
     *
     * (max. sum of rows).
     * This is the natural matrix norm that is
     * compatible to the linfty-norm of vectors, i.e.
     * \f$|Mv|_\infty \leq |M|_\infty |v|_\infty\f$.
     */
    real_type linftyNorm () const override;

    /**
     * Set the element \p (i,j) to \p value.
     * Throws an error if the entry does
     * not exist. Still, it is allowed to store
     * zero values in non-existent fields.
     */
    void set ( const size_type i,
               const size_type j,
               const value_type& value ) override;

    /**
     * Add \p value to the element
     * \p (i,j).  Throws an error if
     * the entry does not
     * exist. Still, it is allowed to
     * store zero values in
     * non-existent fields.
     */
    void add ( const size_type i,
               const size_type j,
               const value_type& value ) override;

    /**
     * Add the full matrix to the
     * Petsc matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix ( const ublas::matrix<value_type> &dm,
                     const std::vector<size_type> &rows,
                     const std::vector<size_type> &cols ) override;

    /**
     * Add the full matrix to the
     * Petsc matrix.  This is useful
     * for adding an element matrix
     * at assembly time
     */
    void addMatrix ( int* rows, int nrows,
                     int* cols, int ncols,
                     value_type* data,
                     size_type K = 0,
                     size_type K2 = invalid_v<size_type>) override;

    /**
     * Same, but assumes the row and column maps are the same.
     * Thus the matrix \p dm must be square.
     */
    void addMatrix ( const ublas::matrix<value_type> &dm,
                     const std::vector<size_type> &dof_indices ) override
    {
        this->addMatrix ( dm, dof_indices, dof_indices );
    }

    /**
     * Add a Sparse matrix \p X, scaled with \p a, to \p this,
     * stores the result in \p this:
     * \f$\texttt{this} = a*X + \texttt{this} \f$.
     * Use this with caution, the sparse matrices need to have the
     * same nonzero pattern, otherwise \p PETSc will crash!
     * It is advisable to not only allocate appropriate memory with
     * \p init() , but also explicitly zero the terms of \p this
     * whenever you add a non-zero value to \p X.  Note: \p X will
     * be closed, if not already done, before performing any work.
     */
    void addMatrix ( const T a, MatrixSparse<T> const&X, Feel::MatrixStructure matStruc = Feel::SAME_NONZERO_PATTERN ) override;

    /**
     * set diagonal entries from vector
     */
    void setDiagonal( const Vector<T>& vecDiag ) override;

    /**
     * add diagonal entries from vector
     */
    void addDiagonal( const Vector<T>& vecDiag ) override;

    /**
     * Multiplies the matrix with \p arg and stores the result in \p
     * dest.
     */
    void multVector( const Vector<T>& arg, Vector<T>& dest, bool transpose ) const override;

    /**
     * Multiply this by a Sparse matrix \p In,
     * stores the result in \p Res:
     * \f$ Res = \texttt{this}*In \f$.
     */
    void matMatMult ( MatrixSparse<T> const& In, MatrixSparse<T> &Res ) const override;

    /**
     * Creates the matrix product C = P^T * A * P with A the current matrix
     */
    void PtAP( MatrixSparse<value_type> const& P, MatrixSparse<value_type> & C ) const override;

    /**
     * Creates the matrix product C = P * A * P^T with A the current matrix
     */
    void PAPt( MatrixSparse<value_type> const& P, MatrixSparse<value_type> & C ) const override;

    /**
     * scale the matrix by the factor \p a
     * \param a scaling factor
     */
    void scale( T const a ) override;

    /**
     * Copies the diagonal part of the matrix into \p dest.
     */
    void diagonal ( Vector<value_type>& dest ) const override;

    /**
     * Return copy vector of the diagonal part of the matrix.
     */
    std::shared_ptr<Vector<T> > diagonal() const override;

    /**
     * Returns the transpose of a matrix
     *
     * \param Mt the matrix transposed
     * \param options options for tranpose
     */
    void transpose( MatrixSparse<value_type>& Mt, size_type options ) const override;

    /**
     * Returns the transpose of a matrix
     *
     * \param options options for tranpose
     */
    std::shared_ptr<MatrixSparse<T> > transpose( size_type options ) const override;

    /**
    * Returns the symmetric part of the matrix
    */
    void symmetricPart( MatrixSparse<value_type>& Ms ) const override;

    /**
     * Returns the raw PETSc matrix context pointer.  Note this is generally
     * not required in user-level code. Just don't do anything crazy like
     * calling MatDestroy()!
     */
    Mat mat () const
    {
        FEELPP_ASSERT ( M_mat != NULL ).error( "null petsc matrix" );
        return M_mat;
    }
    Mat& mat ()
    {
        FEELPP_ASSERT ( M_mat != NULL ).warn( "null petsc matrix" );
        return M_mat;
    }

    /**
     * Print the contents of the matrix in Matlab's
     * sparse matrix format. Optionally prints the
     * matrix to the file named \p name.  If \p name
     * is not specified it is dumped to the screen.
     */
    void printMatlab( const std::string name="NULL" ) const override;

    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     */
    std::shared_ptr<MatrixSparse<T> >
    createSubMatrix( std::vector<size_type> const& rows,
                     std::vector<size_type> const& cols,
                     bool useSameDataMap=false,
                     bool checkAndFixRange=true ) const override;

    /**
     * copy matrix entries in submatrix ( submatrix is already built from a createSubMatrix)
     * row and column indices given in the "rows" and "cols" entries.
     */
    void
    updateSubMatrix( std::shared_ptr<MatrixSparse<T> > & submatrix,
                     std::vector<size_type> const& rows,
                     std::vector<size_type> const& cols, bool doClose = true ) override;


    /**
     * This function creates a matrix called "submatrix" which is defined
     * by the row and column indices given in the "rows" and "cols" entries.
     * Currently this operation is only defined for the PetscMatrix type.
     */
    void createSubmatrix( MatrixSparse<T>& submatrix,
                          const std::vector<size_type>& rows,
                          const std::vector<size_type>& cols ) const override;

    /**
     * \return \f$ v^T M u \f$
     */
    real_type
    energy( Vector<value_type> const& __v,
            Vector<value_type> const& __u,
            bool transpose = false ) const override;

    /**
     * eliminate row without change pattern, and put 1 on the diagonal
     * entry
     *
     *\warning if the matrix was symmetric before this operation, it
     * won't be afterwards. So use the proper solver (nonsymmetric)
     */
    void zeroRows( std::vector<int> const& rows, Vector<value_type> const& values, Vector<value_type>& rhs, Context const& on_context, value_type value_on_diagonal ) override;

    /**
     * update a block matrix
     */
    void updateBlockMat( std::shared_ptr<MatrixSparse<T> > const& m, std::vector<size_type> const& start_i, std::vector<size_type> const& start_j ) override;

    void updatePCFieldSplit( PC & pc, indexsplit_ptrtype const& is );
    void updatePCFieldSplit( PC & pc );

    std::vector<IS> const& petscSplitIS() const { return M_petscIS; }
    std::map<PC*,bool > & mapSplitPC() { return M_mapPC; }

    std::vector<PetscInt> ia() { return M_ia; }
    std::vector<PetscInt> ja() { return M_ja; }


    bool isSymmetric ( bool check = false ) const override;

    bool isTransposeOf ( MatrixSparse<T> &Trans ) const override;

    //@}

    /*
     * Set zero entries diagonal if missing : only for PETSC!
     */
    void zeroEntriesDiagonal();

    //!
    //! get some matrix information use MatInfo data structure from Petsc
    //!
    void getMatInfo(std::vector<double> &) override;

    //!
    //!
    //!
    void threshold( void ) override;

    void save( std::string const& filename="default_archive_name", std::string const& format="binary" ) override;

    void load( std::string const& filename="default_archive_name", std::string const& format="binary" ) override;


private:

    // disable
    MatrixPetsc( MatrixPetsc const & );

    void getSubMatrixPetsc( std::vector<size_type> const& rows,
                            std::vector<size_type> const& cols,
                            Mat &submat, bool doClose = true ) const;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version );

protected:

    /**
     * Petsc matrix datatype to store values
     */
    Mat M_mat;



private:


    std::vector<IS> M_petscIS;

    std::map<PC*,bool > M_mapPC;

    /**
     * This boolean value should only be set to false
     * for the constructor which takes a PETSc Mat object.
     */
    const bool M_destroy_mat_on_exit;
    std::vector<PetscInt> M_ia,M_ja;

    MatInfo M_info;
};



template<typename T>
class MatrixPetscMPI : public MatrixPetsc<T>
{
    typedef MatrixPetsc<T> super;

public :
    using size_type = typename super::size_type;
    typedef typename super::graph_type graph_type;
    typedef typename super::graph_ptrtype graph_ptrtype;
    typedef typename super::value_type value_type;
    typedef typename super::real_type real_type;

    typedef typename super::datamap_type datamap_type;
    typedef typename super::datamap_ptrtype datamap_ptrtype;

    explicit MatrixPetscMPI( worldcomm_ptr_t const& worldComm=Environment::worldCommPtr() );
    MatrixPetscMPI( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol );
    MatrixPetscMPI( datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, worldcomm_ptr_t const& worldComm );

    MatrixPetscMPI( Mat m, datamap_ptrtype const& dmRow, datamap_ptrtype const& dmCol, bool initLocalToGlobalMapping=false, bool destroyMatOnExit=false );

    ~MatrixPetscMPI() override
    {
        this->clear();
    }

    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                const size_type nnz=30,
                const size_type noz=10 ) override;

    /**
     * Initialize using sparsity structure computed by \p dof_map.
     */
    void init ( const size_type m,
                const size_type n,
                const size_type m_l,
                const size_type n_l,
                graph_ptrtype const& graph ) override;

    /**
     * define in petsc matrix the dof mapping between
     * numbering between local process to global world
     */
    void initLocalToGlobalMapping();


    void set( const size_type i,
              const size_type j,
              const value_type& value ) override;

    void add ( const size_type i,
               const size_type j,
               const value_type& value ) override;

    void addMatrix( const ublas::matrix<value_type>& dm,
                    const std::vector<size_type>& rows,
                    const std::vector<size_type>& cols ) override;

    void addMatrix( int* rows, int nrows,
                    int* cols, int ncols,
                    value_type* data, size_type K = 0, size_type K2 = invalid_v<size_type> ) override;

    //void addMatrix( const T a, MatrixSparse<T> const&X );


    void zero() override;
    void zero( size_type start1, size_type stop1, size_type start2, size_type stop2 ) override;
    //void zeroEntriesDiagonal();
    void zeroRows( std::vector<int> const& rows,
                   Vector<value_type> const& values,
                   Vector<value_type>& rhs,
                   Context const& on_context,
                   value_type value_on_diagonal ) override;

private :

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version )
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(super);
    }

    void addMatrixSameNonZeroPattern( const T a, MatrixSparse<T> &X );
};

} // Feel


#endif /* FEELPP_HAS_PETSC */
#endif /* __MatrixPetsc_H */
