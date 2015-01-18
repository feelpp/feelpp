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
#ifndef __Preconditioner_H
#define __Preconditioner_H 1

#include <boost/parameter.hpp>
#include <feel/feelcore/singleton.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vector.hpp>

#include <feel/feelalg/enums.hpp>

namespace Feel
{
template<typename T> class Backend;
typedef Backend<double> backend_type;
typedef boost::shared_ptr<Backend<double> > backend_ptrtype;

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

    typedef Preconditioner<T> preconditioner_type;
    typedef boost::shared_ptr<Preconditioner<T> > preconditioner_ptrtype;

    typedef boost::shared_ptr<MatrixSparse<T> > sparse_matrix_ptrtype;
    typedef boost::shared_ptr<Vector<T> > vector_ptrtype;

    
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
    M_mat_has_changed( o.M_mat_has_changed )
        {}

    //! destructor
    ~Preconditioner();

    static preconditioner_ptrtype build( std::string const& name = "", BackendType = BACKEND_PETSC, WorldComm const& worldComm=Environment::worldComm() );

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
     * @return the side of the system to which the preconditioner applies
     */
    Side side() const { return M_side; }
    
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
     * Enum statitng with type of preconditioner to use.
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
M_matSolverPackage_type ( MATSOLVER_PETSC ),
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

namespace detail
{
class PreconditionerManagerImpl:
        public std::map<std::pair<backend_ptrtype,std::string>, preconditioner_ptrtype >,
        public boost::noncopyable
{
public:
    typedef preconditioner_ptrtype value_type;
    typedef std::pair<backend_ptrtype,std::string> key_type;
    typedef std::map<key_type, value_type> preconditioner_manager_type;

};
typedef Feel::Singleton<PreconditionerManagerImpl> PreconditionerManager;

struct PreconditionerManagerDeleterImpl
{
    void operator()() const
        {
            VLOG(2) << "[PreconditionerManagerDeleter] clear PreconditionerManager Singleton: " << Feel::detail::PreconditionerManager::instance().size() << "\n";
            Feel::detail::PreconditionerManager::instance().clear();
            VLOG(2) << "[PreconditionerManagerDeleter] clear PreconditionerManager done\n";
        }
};
typedef Feel::Singleton<PreconditionerManagerDeleterImpl> PreconditionerManagerDeleter;
} // detail


BOOST_PARAMETER_MEMBER_FUNCTION( ( boost::shared_ptr<Preconditioner<double> > ),
                                 preconditioner,
                                 tag,
                                 ( required
                                   ( pc,( PreconditionerType ) )
                                   ( backend, (backend_ptrtype) ) )
                                 ( optional
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( matrix,( d_sparse_matrix_ptrtype ),d_sparse_matrix_ptrtype() )

                                   ( pcfactormatsolverpackage,( MatSolverPackageType ), MATSOLVER_DEFAULT )
                                   ( rebuild,      (bool), false )
                                     )
    )
{
    // register the PreconditionerManager into Feel::Environment so that it gets the
    // PreconditionerManager is cleared up when the Environment is deleted
    static bool observed=false;
    if ( !observed )
    {
        Environment::addDeleteObserver( Feel::detail::PreconditionerManagerDeleter::instance() );
        observed = true;
    }


    Feel::detail::ignore_unused_variable_warning( args );

    auto git = Feel::detail::PreconditionerManager::instance().find( std::make_pair( backend, prefix ) );

    if (  git != Feel::detail::PreconditionerManager::instance().end() && ( rebuild == false ) )
    {
        VLOG(2) << "[preconditioner] found preconditioner name=" << prefix << " rebuild=" << rebuild << "\n";
        return git->second;
    }

    else
    {

        preconditioner_ptrtype p = Preconditioner<double>::build( prefix, backend->type(), backend->comm() );
        p->setType( pc );
        p->setMatSolverPackageType( pcfactormatsolverpackage );

        if ( matrix )
        {
            p->setMatrix( matrix );
        }
        VLOG(2) << "storing preconditionerin singleton" << "\n";
        Feel::detail::PreconditionerManager::instance().operator[]( std::make_pair( backend, prefix ) ) = p;
        backend->addDeleteObserver( p );
        return p;
    }
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
