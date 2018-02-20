/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-01-16

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
/**
   \file preconditioner.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-01-16
 */
#ifndef FEELPP_PRECONDITIONER_HPP
#define FEELPP_PRECONDITIONER_HPP 1

#include <boost/parameter.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/traits.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vector.hpp>

#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/nullspace.hpp>

namespace Feel
{
class BackendBase;
template<typename T> class Backend;
typedef Backend<double> backend_type;
typedef boost::shared_ptr<Backend<double> > backend_ptrtype;

template<typename T> class OperatorPCDBase;

/**
 * \class Preconditioner
 * \brief base class for preconditioner
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename T>
class Preconditioner
{
public:


    /** @name Enums
     */
    //@{

    /**
     * preconditioner side
     */
    enum Side
    {
        LEFT=0, // default
        RIGHT=1,
        SYMMETRIC=2
    };

    //@}

    /** @name Typedefs
     */
    //@{
    typedef T value_type;
    typedef Preconditioner<T> preconditioner_type;
    typedef boost::shared_ptr<preconditioner_type > preconditioner_ptrtype;

    typedef boost::shared_ptr<MatrixSparse<T> > sparse_matrix_ptrtype;
    typedef boost::shared_ptr<Vector<T> > vector_ptrtype;

    typedef boost::shared_ptr<OperatorPCDBase<T> > operator_pcdbase_ptrtype;
    
    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    Preconditioner( std::string const& name = "", WorldComm const& worldComm=Environment::worldComm() );

    //! copy constructor
    Preconditioner( Preconditioner const & o )
    :
    M_name(),
    M_worldComm( o.M_worldComm ),
    M_matrix( o.M_matrix ),
    M_side( o.M_side ),
    M_preconditioner_type( o.M_preconditioner_type ),
    M_matSolverPackage_type( o.M_matSolverPackage_type ),
    M_prec_matrix_structure ( o.M_prec_matrix_structure ),
    M_is_initialized( o.M_is_initialized ),
    M_mat_has_changed( o.M_mat_has_changed ),
    M_nearNullSpace( o.M_nearNullSpace )
        {}

    //! destructor
    virtual ~Preconditioner();

    static preconditioner_ptrtype build( 
            std::string const& name = "", 
#if FEELPP_HAS_PETSC
            BackendType = BACKEND_PETSC, 
#else
            BackendType = BACKEND_NONE,
#endif
            WorldComm const& worldComm=Environment::worldComm() );

    /**
     * Initialize data structures if not done so already.
     */
    virtual void init () {};

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Preconditioner& operator=( Preconditioner const & o ) 
        {
            if ( this != &o )
            {
                M_name = o.M_name;
                M_worldComm = o.M_worldComm;
                M_matrix = o.M_matrix;
                M_side = o.M_side;
                M_is_initialized = o.M_is_initialized;
                M_matSolverPackage_type = o.M_matSolverPackage_type;
                M_prec_matrix_structure = o.M_prec_matrix_structure;
                M_preconditioner_type = o.M_preconditioner_type;
                M_mat_has_changed = o.M_mat_has_changed;
            }

            return *this;
        }

    void operator()()
        {
            this->clear();
        }
    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @returns true if the data structures are
     * initialized, false otherwise.
     */
    bool initialized () const
        {
            return M_is_initialized;
        }

    WorldComm const& worldComm() const { return M_worldComm; }

    /**
     * View preconditioner context
     */
    virtual void view() const {};

    /**
     * Computes the preconditioned vector "y" based on input "x".
     * Usually by solving Py=x to get the action of P^-1 x.
     */
    virtual void apply( const Vector<T> & x, Vector<T> & y ) const = 0;

    /**
     * Computes the preconditioned vector "y" based on input "x".
     * Usually by solving Py=x to get the action of P^-1 x.
     */
    void apply( vector_ptrtype const& x, vector_ptrtype& y ) const
        {
            this->apply( *x, *y );
        }

    /**
     * Release all memory and clear data structures.
     */
    virtual void clear () {}



    /**
     * Returns the type of preconditioner to use.
     */
    PreconditionerType type () const
        {
            return M_preconditioner_type;
        }

    virtual std::string name() const { return M_name; }

    sparse_matrix_ptrtype const& matrix() const { return M_matrix; }

    /**
     * Return true if the preconditioner will be reuse
     */
    bool reusePrec() const { return M_prec_matrix_structure == MatrixStructure::SAME_PRECONDITIONER; }

    /**
     * @return the side of the system to which the preconditioner applies
     */
    Side side() const { return M_side; }

    bool hasNearNullSpace( std::set<int> const& splitIds ) const { return M_nearNullSpace.find(splitIds) != M_nearNullSpace.end(); }
    boost::shared_ptr<NullSpace<value_type> > const& nearNullSpace( std::set<int> const& splitIds ) const
    {
        CHECK( this->hasNearNullSpace( splitIds ) ) << " near null space not given for index split ";
        return M_nearNullSpace.find(splitIds)->second;
    }

    bool hasAuxiliaryVector( std::string const& key ) const { return M_auxiliaryVector.find( key ) != M_auxiliaryVector.end(); }
    vector_ptrtype const& auxiliaryVector( std::string const& key ) const
    {
        CHECK( this->hasAuxiliaryVector( key ) ) << " auxiliary vector not given for this key : " << key ;
        return M_auxiliaryVector.find( key )->second;
    }

    bool hasAuxiliarySparseMatrix( std::string const& key ) const { return M_auxiliarySparseMatrix.find( key ) != M_auxiliarySparseMatrix.end(); }
    sparse_matrix_ptrtype const& auxiliarySparseMatrix( std::string const& key ) const
    {
        CHECK( this->hasAuxiliarySparseMatrix( key ) ) << " auxiliary sparse matrix not given for this key : " << key ;
        return M_auxiliarySparseMatrix.find( key )->second;
    }

    bool hasInHousePreconditioners( std::string const& key ) const { return M_inHousePreconditioners.find( key ) != M_inHousePreconditioners.end(); }
    preconditioner_ptrtype const& inHousePreconditioners( std::string const& key ) const
    {
        CHECK( this->hasInHousePreconditioners( key ) ) << " in house preconditioner not given for this key : " << key ;
        return M_inHousePreconditioners.find(key)->second;
    }
    preconditioner_ptrtype & inHousePreconditioners( std::string const& key )
    {
        CHECK( this->hasInHousePreconditioners( key ) ) << " in house preconditioner not given for this key : " << key ;
        return M_inHousePreconditioners[key];
    }

    bool hasOperatorPCD( std::string const& key ) const { return M_operatorPCD.find( key ) != M_operatorPCD.end(); }
    operator_pcdbase_ptrtype const& operatorPCD( std::string const& key ) const
        {
            CHECK( this->hasOperatorPCD( key ) ) << " operator PCD not given for this key : " << key ;
            return M_operatorPCD.find(key)->second;
        }

    //@}

    /** @name  Mutators
     */
    //@{

    virtual void setName( std::string const& n ) { M_name = n; }

    /**
     * Sets the matrix P to be preconditioned.
     */
    void setMatrix( sparse_matrix_ptrtype  mat );

    /**
     * Sets the type of preconditioner to use.
     */
    void setType ( const PreconditionerType pct );

    /**
     * the software that is used to perform the factorization
     */
    void setMatSolverPackageType( const MatSolverPackageType mspt );

    /**
     * information about the preconditioner matrix structure during successive linear solves
     */
    virtual void setPrecMatrixStructure( MatrixStructure mstruct  );

    /**
     * set the side \p s of the linear system to which the preconditioner applies
     */
    void setSide( Side s ) { M_side = s; }

    void attachNearNullSpace( int k, boost::shared_ptr<NullSpace<value_type> > const& nearNullSpace )
    {
        std::set<int> splitIds; splitIds.insert( k );
        this->attachNearNullSpace( splitIds, nearNullSpace );
    }
    void attachNearNullSpace( std::set<int> const& splitIds, boost::shared_ptr<NullSpace<value_type> > const& nearNullSpace )
    {
        M_nearNullSpace[splitIds] = nearNullSpace;
    }

    void attachAuxiliaryVector( std::string const& key,vector_ptrtype const& vec )
    {
        M_auxiliaryVector[key] = vec;
    }

    void attachAuxiliarySparseMatrix( std::string const& key,sparse_matrix_ptrtype const& mat )
    {
        M_auxiliarySparseMatrix[key] = mat;
    }

    void attachInHousePreconditioners( std::string const& key, preconditioner_ptrtype const& pc )
    {
        M_inHousePreconditioners[key] = pc;
    }

    void attachOperatorPCD( std::string const& key, operator_pcdbase_ptrtype const& opPCD )
    {
        M_operatorPCD[key] = opPCD;
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}


protected:

    /**
     * name of the preconditioner
     */
    std::string M_name;

    /**
     * Communicator
     */
    WorldComm M_worldComm;

    /**
     * The matrix P... ie the matrix to be preconditioned.
     * This is often the actual system matrix of a linear sytem.
     */
    sparse_matrix_ptrtype  M_matrix;

    /**
     * side of the preconditioner
     */
    Side M_side;
    
    /**
     * Enum stating with type of preconditioner to use.
     */
    PreconditionerType M_preconditioner_type;

    /**
     * Enum the software that is used to perform the factorization
     */
    MatSolverPackageType M_matSolverPackage_type;

    /**
     * Enum that indicating information about the preconditioner matrix structure during successive linear solves
     */
    MatrixStructure M_prec_matrix_structure;

    /**
     * Flag indicating if the data structures have been initialized.
     */
    bool M_is_initialized, M_mat_has_changed;

    /**
     *  Near Null Space for Field Split
     */
    std::map<std::set<int>,boost::shared_ptr<NullSpace<value_type> > >  M_nearNullSpace;

    std::map<std::string,sparse_matrix_ptrtype> M_auxiliarySparseMatrix;
    std::map<std::string,vector_ptrtype> M_auxiliaryVector;

    std::map<std::string,preconditioner_ptrtype> M_inHousePreconditioners;

    std::map<std::string,operator_pcdbase_ptrtype> M_operatorPCD;
};

typedef Preconditioner<double> preconditioner_type;
typedef boost::shared_ptr<Preconditioner<double> > preconditioner_ptrtype;


template <typename T>
FEELPP_STRONG_INLINE
Preconditioner<T>::Preconditioner ( std::string const& name, WorldComm const& worldComm )
:
M_name(name),
M_worldComm(worldComm),
M_matrix(),
M_side( LEFT ),
M_preconditioner_type   ( ILU_PRECOND ),
#if FEELPP_HAS_PETSC
M_matSolverPackage_type ( MATSOLVER_PETSC ),
#else
M_matSolverPackage_type ( MATSOLVER_NONE ),
#endif
M_prec_matrix_structure ( MatrixStructure::SAME_NONZERO_PATTERN ),
M_is_initialized        ( false ),
M_mat_has_changed       ( false )
{
}



template <typename T>
FEELPP_STRONG_INLINE
Preconditioner<T>::~Preconditioner ()
{
    this->clear ();
}

typedef Preconditioner<double> preconditioner_type;
typedef boost::shared_ptr<preconditioner_type> preconditioner_ptrtype;


template<typename Args>
struct compute_prec_return
{
    typedef typename parameter::value_type<Args, tag::backend>::type::element_type::value_type value_type;
    typedef boost::shared_ptr<Preconditioner<value_type>> type;
};

BOOST_PARAMETER_FUNCTION( ( boost::shared_ptr<Preconditioner<double> > ),
                          preconditioner,
                          tag,
                          ( required
                            ( pc,( PreconditionerType ) )
                            ( backend, *(boost::is_convertible<mpl::_, boost::shared_ptr<BackendBase>>) ) )
                          ( optional
                            ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                            ( matrix,( d_sparse_matrix_ptrtype ),d_sparse_matrix_ptrtype() )
                            ( pcfactormatsolverpackage,( MatSolverPackageType ), MATSOLVER_DEFAULT )
                            ( rebuild,      (bool), false )
                            )
                          )
{
    using value_type = typename compute_prec_return<Args>::value_type;
    preconditioner_ptrtype p = Preconditioner<value_type>::build( prefix, backend->type(), backend->comm() );
    p->setType( pc );
    p->setMatSolverPackageType( pcfactormatsolverpackage );

    if ( matrix )
    {
        p->setMatrix( matrix );
    }
    return p;
}

/**
 * FEELPP_INSTANTIATE_PRECONDITIONER is never defined except in preconditioner.cpp
 * where we do the instantiate. This allows to reduce the Preconditioner
 * instantiation to the strict minimum
 */
#if !defined( FEELPP_INSTANTIATE_PRECONDITIONER )
extern template class Preconditioner<double>;
extern template class Preconditioner<std::complex<double>>;
#endif


} // Feel
#endif /* __Preconditioner_H */
