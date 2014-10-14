/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 30 Sep 2014
 
 Copyright (C) 2014 Feel++ Consortium
 
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

#ifndef FEELPP_FEELALG_OPERATOR_HPP
#define FEELPP_FEELALG_OPERATOR_HPP 1

#include <feel/feelcore/worldcomm.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/vector.hpp>

namespace Feel
{
/**
 * Base class for Operators
 */
class OperatorBase
{
  public:
 
    OperatorBase( WorldComm const& comm, std::string label, bool use_transpose  ) 
        : 
        M_label( label ),
        M_comm( comm ), 
        M_use_transpose( use_transpose ) 
        {}
    OperatorBase( std::string label, bool use_transpose = false  ) 
        : 
        M_label( label ),
        M_comm( Environment::worldComm() ), 
        M_use_transpose( use_transpose ) 
        {}
    
    virtual ~OperatorBase() {};
    
    virtual int setUseTranspose(bool UseTranspose)  { M_use_transpose = UseTranspose; }
    
    virtual int apply(const vector_ptrtype& X, vector_ptrtype& Y) const { return apply( *X, *Y ); }
    virtual int apply(const vector_type& X, vector_type& Y) const;
    
    virtual int applyInverse(const vector_ptrtype& X, vector_ptrtype& Y) const { return apply( *X, *Y ); }
    virtual int applyInverse(const vector_type& X, vector_type& Y) const = 0;
 
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
     */ 
    virtual double normInf() const = 0;
    
    /**
     * \return the label of the operator
     */
    virtual std::string const& label() const { return M_label; }

    /**
     * st the label
     */
    virtual void setLabel( std::string const& l ) { M_label = l; }

    /**
     * \return true is transposed should be used, false otherwise
     */
    virtual bool useTranspose() const { return M_use_transpose; }
 
    /**
     * \return  true if Operator can compute the inf norm, false otherwise
     */
    virtual bool hasNormInf() const = 0;
    
    /**
     * \return the WorldComm of the Operator
     */
    virtual const WorldComm& comm() const { return M_comm; }
 
protected:
    std::string M_label;
    WorldComm M_comm;
    bool M_use_transpose;
};

template<typename T>
class OperatorMatrix : public OperatorBase
{
public:
    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;
    
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    typedef OperatorBase prec_type;
    typedef boost::shared_ptr<prec_type> prec_ptrtype;

    OperatorMatrix( sparse_matrix_ptrtype const& F, std::string _label, bool transpose = 0 )
        :
        OperatorBase( F->comm(), _label, transpose ),
        M_F( F ),
        M_Prec( ),
        M_hasInverse( 1 ),
        M_hasApply( 1 )
    {
        DVLOG(2) << "Create operator " << this->label() << " ...\n";
    }

    OperatorMatrix( sparse_matrix_ptrtype const& F,
                    std::string _label,
                    prec_ptrtype  Prec, bool transpose = 0 )
        :
        OperatorBase( F->comm(), _label, transpose ),
        M_F( F ),
        M_Prec( Prec ),
        M_hasInverse( 1 ),
        M_hasApply( 1 )
    {
        DVLOG(2) << "Create operator " << this->label() << " ...\n";

    }

    OperatorMatrix( const OperatorMatrix& tc )
        :
        OperatorBase( tc ),
        M_F( tc.M_F ),
        M_Prec( tc.M_Prec ),
        M_hasInverse( tc.M_hasInverse ),
        M_hasApply( tc.M_hasApply )
    {
        DVLOG(2) << "Copy operator " << this->label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_hasInverse;
    }

    bool hasApply() const
    {
        return M_hasApply;
    }


    int apply( const vector_type& X, vector_type& Y ) const
    {
        M_F->multVector( X, Y );
        return !hasApply();
    }
    
    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        CHECK( hasInverse() ) << "Operator " << this->label() << "cannot be inverted.";

        auto r = backend(_prefix=this->label())->solve( _matrix=M_F, _rhs=X, _solution=Y );
        return r->isConverged();
    }

    value_type normInf() const
    {
        return M_F->linftyNorm();
    }

    ~OperatorMatrix()
    {
        DVLOG(2) << "Destroyed matrix operator: " << this->label() << " ...\n";
    };

private:

    sparse_matrix_ptrtype M_F;
    prec_ptrtype M_Prec;

    bool M_hasInverse, M_hasApply;
};

/**
 * Operator class to model an inverse operator
 */
template< typename operator_type >
class OperatorInverse : public OperatorBase
{
    typedef OperatorBase super;
public:

    typedef boost::shared_ptr<operator_type> operator_ptrtype;
    typedef typename operator_type::value_type value_type;

    // This constructor implements the F^-1 operator
    OperatorInverse( operator_ptrtype& F, std::string label )
        :
        OperatorBase( F->comm(), label, F->useTranspose() ),
        M_F( F )
    {
        this->setName();
        DVLOG(2) << "Create inverse operator " << this->label() << "...\n";
    }

    OperatorInverse( const OperatorInverse& tc )
        :
        super(tc),
        M_F( tc.M_F )
    {
        DVLOG(2) << "Copy inverse operator " << this->label() << "...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasApply();
    }

    bool hasApply() const
    {
        return M_F->hasInverse();
    }

    int apply( const vector_type & X, vector_type & Y ) const
    {
        DVLOG(2) << "apply matrix " << label() << "\n";

        CHECK( hasApply() ) << "This operator" << this->label() << " cannot be applied.";
        M_F->applyInverse( X,Y );

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        DVLOG(2) << "applyInverse matrix " << label() << "\n";

        CHECK( hasInverse() ) << "This operator" << this->label() << " cannot be inverted.";

        M_F->apply( X,Y );

        return !hasInverse();
    }

    value_type normInf() const
    {
        return( false );
    }

    ~OperatorInverse()
    {
        DVLOG(2) << "Destroyed inverse operator: " << this->label() << " ...\n";
    };

private:

    operator_ptrtype M_F;

    void setName()
    {
        std::string L = M_F->label();
        L.append( ")" );
        std::string temp( "inv(" );
        temp.append( L );

        this->setLabel(temp);
    }
};






template< typename op1_type, typename op2_type >
class OperatorCompose : public OperatorBase
{
    typedef OperatorBase super;
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef boost::shared_ptr<op2_type> op2_ptrtype;
    typedef typename op2_type::value_type value_type;
    // This constructor implements the (F o G) operator

    OperatorCompose()
        :
        M_F(),
        M_G()
    {
        std::string t( M_F->label() );
        std::string u( M_G->label() );

        t.append( "*" );
        t.append( u );
        this->setLabel(t);
    }

    OperatorCompose( op1_ptrtype& F, op2_ptrtype& G )
        :
        M_F( F ),
        M_G( G )
    {
        std::string t( F->label() );
        std::string u( G->label() );

        t.append( "*" );
        t.append( u );
        this->setLabel(t);
        DVLOG(2) << "Create operator " << label() << " ...\n";
    }

    OperatorCompose( const OperatorCompose& tc )
        :
        super(tc),
        M_F( tc.M_F ),        
        M_G( tc.M_G )
    {
        DVLOG(2) << "Copy operator " << label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasInverse() * M_G->hasInverse();
    }

    bool hasApply() const
    {
        return M_F->hasApply() * M_G->hasApply();
    }

    int apply( const vector_type & X, vector_type & Y ) const
    {
        CHECK( hasApply() ) << "This operator " << this->label() << " cannot be applied.";

        DVLOG(2) << "apply operator " << label() << " ...\n";

        vector_ptrtype Z;

        M_G->apply( X,Z );
        M_F->apply( Z,Y );

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        CHECK( hasInverse() ) << "This operator " << this->label() << " cannot be inverted.";

        DVLOG(2) << "apply Inverse operator " << label() << " ...\n";

        vector_ptrtype Z = X.clone();

        M_F->applyInverse( X,Z );
        M_G->applyInverse( Z,Y );

        return hasInverse();
    }


    value_type normInf() const
    {
        return 0;
    }

    ~OperatorCompose()
    {
        DVLOG(2) << "Destroyed compose operator: " << this->label() << " ...\n";
    };

private:

    op1_ptrtype M_F;
    op2_ptrtype M_G;

};

/**
 * Scaling Operator class
 */
template< typename op1_type>
class OperatorScale : public OperatorBase
{
    typedef OperatorBase super;
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef typename op1_type::value_type value_type;

    // This constructor implements the (\alpha F) operator
    OperatorScale()
        :
        M_F(),
        M_alpha( 0 )
    {
    }

    OperatorScale( op1_ptrtype& F )
        :
        M_F( F ),
        M_alpha( 1 )
    {
        std::string temp( "alpha" );
        std::string t = M_F->label();

        temp.append( "." );
        temp.append( t );
        
        this->setLabel(temp);

        DVLOG(2) << "Create scale operator " << label() << " ...\n";
    }

    OperatorScale( op1_ptrtype& F, value_type alpha )
        :
        M_F( F ),
        M_alpha( alpha )
    {
        std::string temp( "alpha" );
        std::string t = M_F->label();

        temp.append( "." );
        temp.append( t );

        this->setLabel(temp);

        DVLOG(2) << "Create scale operator " << label() << " ...\n";
    }

    OperatorScale( const OperatorScale& tc )
        :
        super( tc ),
        M_F( tc.M_F ),
        M_alpha( tc.M_alpha )
    {
        DVLOG(2) << "Copy scale operator " << label() << " ...\n";
    }

    bool hasInverse() const
    {
        return M_F->hasApply()*( M_alpha != 0 );
    }

    bool hasApply() const
    {
        return M_F->hasApply();
    }



    int apply( const vector_type& X, vector_type & Y ) const
    {
        DVLOG(2) << "apply scale operator " << label() << "\n";

        M_F->apply( X,Y );

        Y.scale( M_alpha );

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        DVLOG(2) << "applyInverse scale operator " << label() << "\n";

        FEELPP_ASSERT( hasInverse() && ( M_alpha != 0 ) ).error( "This operator cannot be inverted." );

        vector_ptrtype Z =  X.clone();
        Z->scale( 1./M_alpha );

        M_F->applyInverse( Z,Y );

        return !hasInverse();
    }

    value_type normInf() const
        {
            return 0;
        }

    ~OperatorScale()
    {
        //DVLOG(2) << "Destroyed scale operator: " << label() << " ...\n";
    };

private:

    op1_ptrtype M_F;

    value_type M_alpha;
};






template< typename operator_type >
class OperatorFree : public OperatorBase
{
    typedef OperatorBase super;
public:

    typedef boost::shared_ptr<operator_type> operator_ptrtype;
    typedef typename operator_type::value_type value_type;
    typedef OperatorBase prec_type;
    typedef boost::shared_ptr<prec_type> prec_ptrtype;

    OperatorFree( operator_ptrtype F )
        :
        super( F->label(), F->useTranspose() ),
        M_op( F ),
        M_hasInverse( 0 ),
        M_hasApply( F->hasApply() )
    {
        DVLOG(2) << "Create operator " << label() << " ...\n";
    }

    OperatorFree( operator_ptrtype F, prec_ptrtype Prec )
        :
        super(F->label(), F->useTranspose()),
        M_op( F ),
        M_hasInverse( 1 ),
        M_hasApply( F->hasApply() ),
        M_Prec( Prec )
    {
        DVLOG(2) << "Create operator " << label() << " ...\n";
    }


    OperatorFree( const OperatorFree& tc )
        :
        super( tc ),
        M_op( tc.M_op ),
        M_hasInverse( tc.M_hasInverse ),
        M_hasApply( tc.M_hasApply ),
        M_Prec( tc.M_Prec )
    {
        DVLOG(2) << "Copy operator " << label() << " ...\n";
    }


    bool hasInverse() const
    {
        return M_hasInverse;
    }

    bool hasApply() const
    {
        return M_hasApply;
    }


    int apply( const vector_ptrtype & X, vector_ptrtype & Y ) const
    {
        DVLOG(2) << "apply operator " << label() << "\n";
        M_op->apply( X,Y );
        DVLOG(2) << "Finished apply operator " << label() << "\n";

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        DVLOG(2) << "applyInverse operator " << label() << "\n";

        FEELPP_ASSERT( hasInverse() ).error( "This operator cannot be inverted." );

#if 0
        std::pair<unsigned int, real_type> result = M_Solver->solve( M_op, M_Prec, Y, X, M_tol, M_maxiter );
#endif
#warning TODO : solve system

        DVLOG(2) << "Finished applyInverse operator " << label() << "\n";
        return !hasInverse();
    }

    value_type NormInf() const
    {
        return 0;
    }

    bool HasNormInf () const
    {
        return( true );
    }

    ~OperatorFree()
    {
        DVLOG(2) << "Destroyed operator: " << label() << " ...\n";
    };

private:

    operator_ptrtype M_op;


    bool M_hasInverse, M_hasApply;

    prec_ptrtype M_Prec;
};








} // Feel

#endif /* FEELPP_FEELAG_OPERATOR_HPP */
