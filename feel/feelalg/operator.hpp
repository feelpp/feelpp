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
template<typename T>
class OperatorBase
{
  public:
    typedef T value_type;
    typedef DataMap datamap_type;
    typedef boost::shared_ptr<DataMap> datamap_ptrtype;

    OperatorBase() 
        : 
        M_domain_map(),
        M_image_map(),
        M_label(""),
        M_use_transpose( false ) ,
        M_has_norminf( false )
        {}

    OperatorBase( datamap_ptrtype const& map, std::string label, bool use_transpose, bool has_norminf  ) 
        : 
        M_domain_map( map ),
        M_image_map( map ),
        M_label( label ),
        M_use_transpose( use_transpose ) ,
        M_has_norminf( has_norminf )
        {}
    OperatorBase( datamap_ptrtype  const& dmap, datamap_ptrtype const& imap, std::string label, bool use_transpose, bool has_norminf  ) 
        : 
        M_domain_map( dmap ),
        M_image_map( imap ),
        M_label( label ),
        M_use_transpose( use_transpose ) ,
        M_has_norminf( has_norminf )
        {}

#if 0
    OperatorBase( std::string label, bool use_transpose = false, bool has_norminf = false  ) 
        : 
        M_label( label ),
        M_use_transpose( use_transpose ) ,
        M_has_norminf( has_norminf )
        {}
#endif
    virtual ~OperatorBase() {};
    
    virtual void setUseTranspose(bool UseTranspose)  { M_use_transpose = UseTranspose; }
    
    virtual int apply(const vector_ptrtype& X, vector_ptrtype& Y) const { return apply( *X, *Y ); }
    virtual int apply(const vector_type& X, vector_type& Y) const = 0;
    
    virtual int applyInverse(const vector_ptrtype& X, vector_ptrtype& Y) const { return apply( *X, *Y ); }
    virtual int applyInverse(const vector_type& X, vector_type& Y) const = 0;
 
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
     */ 
    virtual double normInf() const { return 0; };
    
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
    virtual bool hasNormInf() const { return M_has_norminf; };
    
    /**
     * set if the operator can compute its inf norm
     */
    void setHasNormInf( bool b ) { M_has_norminf = b; }

    /**
     * domain map
     */
    DataMap const& domainMap() const { return *M_domain_map; }

    /**
     * domain map
     */
    datamap_ptrtype const& domainMapPtr() const { return M_domain_map; }

    /**
     * domain map
     */
    datamap_ptrtype domainMapPtr() { return M_domain_map; }

    /**
     * image map
     */
    DataMap const& imageMap() const { return *M_image_map; }

    /**
     * image map
     */
    datamap_ptrtype const& imageMapPtr() const { return M_image_map; }

    /**
     * image map
     */
    datamap_ptrtype imageMapPtr() { return M_image_map; }


    /**
     * \return the WorldComm of the Operator
     */
    virtual const WorldComm& comm() const { return M_domain_map->worldComm(); }
 
protected:
    // domain and image map 
    datamap_ptrtype M_domain_map, M_image_map;
    std::string M_label;
    bool M_use_transpose;
    bool M_has_norminf;
};

template<typename T>
class OperatorMatrix : public OperatorBase<T>
{
public:
    typedef T value_type;
    typedef typename type_traits<T>::real_type real_type;
    
    typedef MatrixSparse<value_type> sparse_matrix_type;
    typedef boost::shared_ptr<sparse_matrix_type> sparse_matrix_ptrtype;

    typedef Vector<value_type> vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;


    OperatorMatrix()
        :
        OperatorBase<T>()
        {}
    OperatorMatrix( sparse_matrix_ptrtype const& F, std::string _label, bool transpose = 0 )
        :
        OperatorBase<T>( F->mapColPtr(), F->mapRowPtr(), _label, transpose, true ),
        M_F( F ),
        M_hasInverse( 1 ),
        M_hasApply( 1 )
    {
        auto b = backend(_name=this->label(),_rebuild=true);
        LOG(INFO) << "Create operator " << this->label() << " ...\n";
    }

    OperatorMatrix( const OperatorMatrix& tc )
        :
        OperatorBase<T>( tc ),
        M_F( tc.M_F ),
        M_hasInverse( tc.M_hasInverse ),
        M_hasApply( tc.M_hasApply )
    {
        LOG(INFO) << "Copy operator " << this->label() << " ...\n";
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
        LOG(INFO) << "OperatorMatrix: apply(X,Y)";
        M_F->multVector( X, Y );
        Y.close();
        return !hasApply();
    }
    
    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        CHECK( hasInverse() ) << "Operator " << this->label() << "cannot be inverted.";
        LOG(INFO) << "OperatorMatrix: applyInverse(X,Y)";
        auto xx = backend(_name=this->label())->newVector( X.mapPtr() );
        *xx = X;
        xx->close();
        auto yy = backend(_name=this->label())->newVector( Y.mapPtr() );
        //auto r = backend(_name=this->label())->solve( _matrix=M_F, _rhs=X.shared_from_this(), _solution=Y.shared_from_this() );
        auto r = backend(_name=this->label())->solve( _matrix=M_F, _rhs=xx, _solution=yy );
        Y=*yy;
        Y.close();
        return r.isConverged();
    }


    value_type normInf() const
    {
        return M_F->linftyNorm();
    }

    virtual ~OperatorMatrix()
    {
        LOG(INFO) << "Destroyed matrix operator: " << this->label() << " ...\n";
    };

private:

    sparse_matrix_ptrtype M_F;

    bool M_hasInverse, M_hasApply;
};
/**
 * \param M matrix
 * \oaram l label of the operator
 * \param transpose boolean to say wether we want the matrix or its transpose
 * \return the Operator associated to the matrix \p M
 */
template<typename MatrixType>
boost::shared_ptr<OperatorMatrix<typename MatrixType::value_type>>
op( boost::shared_ptr<MatrixType> M, std::string label, bool transpose = false )
{
    return boost::make_shared<OperatorMatrix<typename MatrixType::value_type>>(M,label,transpose) ;
}
/**
 * Operator class to model an inverse operator
 */
template< typename operator_type >
class OperatorInverse : public OperatorBase<typename operator_type::value_type>
{
    typedef OperatorBase<typename operator_type::value_type> super;
public:

    typedef boost::shared_ptr<operator_type> operator_ptrtype;
    typedef typename operator_type::value_type value_type;

    // This constructor implements the F^-1 operator
    OperatorInverse( operator_ptrtype& F )
        :
        super( F->imageMapPtr(), F->domainMapPtr(), F->label(), F->useTranspose(), false ),
        M_F( F )
    {
        this->setName();
        LOG(INFO) << "Create inverse operator " << this->label() << "...\n";
    }

    OperatorInverse( const OperatorInverse& tc )
        :
        super(tc),
        M_F( tc.M_F )
    {
        LOG(INFO) << "Copy inverse operator " << this->label() << "...\n";
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
        LOG(INFO) << "OperatorInverse: apply matrix " << this->label() << "\n";

        CHECK( hasApply() ) << "This operator " << this->label() << " cannot be applied.";
        M_F->applyInverse( X,Y );
        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        LOG(INFO) << "OperatorInverse: applyInverse matrix " << this->label() << "\n";

        CHECK( hasInverse() ) << "This operator" << this->label() << " cannot be inverted.";
        CHECK(M_F) << "Invalid operator " << this->label() << " to inverse";
        M_F->apply( X,Y );

        return !hasInverse();
    }

    value_type normInf() const
    {
        return( false );
    }

    ~OperatorInverse()
    {
        LOG(INFO) << "Destroyed inverse operator: " << this->label() << " ...\n";
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


/**
 * \param op an operator
 * \oaram l label of the operator
 * \param transpose boolean to say wether we want the matrix or its transpose
 * \return the Operator associated to the matrix \p M
 */
template<typename OpType>
boost::shared_ptr<OperatorInverse<OpType>>
inv( boost::shared_ptr<OpType>  M )
{
    return boost::make_shared<OperatorInverse<OpType>>(M) ;
}


/**
 * implement  \f$f \circ g \f$ operator
 */
template< typename op1_type, typename op2_type >
class OperatorCompose : public OperatorBase<typename op1_type::value_type>
{
    typedef OperatorBase<typename op1_type::value_type> super;
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef boost::shared_ptr<op2_type> op2_ptrtype;
    typedef typename op2_type::value_type value_type;
    // This constructor implements the (F o G) operator


    OperatorCompose( op1_ptrtype F, op2_ptrtype G )
        :
        super(G->domainMapPtr(),F->imageMapPtr(), F->label(),F->useTranspose(),false),
        M_F( F ),
        M_G( G ),
        M_ZG ( backend()->newVector( M_G->imageMapPtr() ) ),
        M_ZF ( backend()->newVector( M_F->domainMapPtr() ) )
    {
        // TODO:We should ensure here that the domain map of F corresponds to the
        // image map of G.
        
        std::string t( F->label() );
        std::string u( G->label() );

        t.append( "*" );
        t.append( u );
        this->setLabel(t);
        LOG(INFO) << "Create operator " << this->label() << " ...\n";
    }

    OperatorCompose( const OperatorCompose& tc )
        :
        super(tc),
        M_F( tc.M_F ),        
        M_G( tc.M_G ),
        M_ZG( tc.M_ZG ),
        M_ZF( tc.M_ZF )
    {
        LOG(INFO) << "Copy operator " << this->label() << " ...\n";
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

        LOG(INFO) << "OperatorCompose: apply operator " << this->label() << " ...\n";

        
        
        LOG(INFO) << "  - apply operator " << M_G->label() << " ...\n";
        M_G->apply( X,*M_ZG );
        LOG(INFO) << "  - apply operator " << M_F->label() << " ...\n";
        M_F->apply( *M_ZG,Y );

        LOG(INFO) << "OperatorCompose apply operator " << this->label() << " done.\n";

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        CHECK( hasInverse() ) << "This operator " << this->label() << " cannot be inverted.";

        LOG(INFO) << "OperatorCompose apply operator " << this->label() << " ...\n";

        LOG(INFO) << "  - apply operator " << M_F->label() << " ...\n";
        M_F->applyInverse( X,*M_ZF );
        LOG(INFO) << "  - apply operator " << M_G->label() << " ...\n";
        M_G->applyInverse( *M_ZF,Y );
        LOG(INFO) << "OperatorCompose applyInverse operator " << this->label() << " done.\n";
        return hasInverse();
    }


    value_type normInf() const
    {
        return 0;
    }

    ~OperatorCompose()
    {
        LOG(INFO) << "Destroyed compose operator: " << this->label() << " ...\n";
    };

private:

    op1_ptrtype M_F;
    op2_ptrtype M_G;
    vector_ptrtype M_ZG;
    vector_ptrtype M_ZF;
};

/**
 * \param op1 an operator
 * \param op2 an operator
 * \return the operator which is the composition of op1 with op2
 */
template<typename Op1Type, typename Op2Type>
boost::shared_ptr<OperatorCompose<Op1Type,Op2Type>>
compose( boost::shared_ptr<Op1Type>  op1,  boost::shared_ptr<Op2Type>  op2  )
{
    return boost::make_shared<OperatorCompose<Op1Type,Op2Type>>(op1,op2) ;
}

/**
 * Scaling Operator class
 */
template< typename op1_type>
class OperatorScale : public OperatorBase<typename op1_type::value_type>
{
    typedef OperatorBase<typename op1_type::value_type> super;
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef typename op1_type::value_type value_type;


    OperatorScale()
        :
        OperatorBase<value_type>()
        {}
    OperatorScale( op1_ptrtype& F )
        :
        super( F->domainMapPtr(), F->imageMapPtr(), F->label(), F->useTranspose(), F->hasNormInf() ),
        M_F( F ),
        M_alpha( 1 )
    {
        std::string temp( "alpha" );
        std::string t = M_F->label();

        temp.append( "." );
        temp.append( t );
        
        this->setLabel(temp);

        LOG(INFO) << "Create scale operator " << this->label() << " ...\n";
    }

    OperatorScale( op1_ptrtype& F, value_type alpha )
        :
        super( F->domainMapPtr(), F->imageMapPtr(), "", false, false ),
        M_F( F ),
        M_alpha( alpha )
    {
        std::string temp( "alpha" );
        std::string t = M_F->label();

        temp.append( "." );
        temp.append( t );

        this->setLabel(temp);

        LOG(INFO) << "Create scale operator " << this->label() << " ...\n";
    }

    OperatorScale( const OperatorScale& tc )
        :
        super( tc ),
        M_F( tc.M_F ),
        M_alpha( tc.M_alpha )
    {
        LOG(INFO) << "Copy scale operator " << this->label() << " ...\n";
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
        LOG(INFO) << "apply scale operator " << this->label() << "\n";

        M_F->apply( X,Y );

        Y.scale( M_alpha );
        Y.close();

        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        LOG(INFO) << "applyInverse scale operator " << this->label() << "\n";

        FEELPP_ASSERT( hasInverse() && ( M_alpha != 0 ) ).error( "This operator cannot be inverted." );

        
        M_F->applyInverse( X,Y );
        Y.scale( 1./M_alpha );
        Y.close();

        return !hasInverse();
    }

    value_type normInf() const
        {
            return 0;
        }

    ~OperatorScale()
    {
        //LOG(INFO) << "Destroyed scale operator: " << label() << " ...\n";
    };

private:

    op1_ptrtype M_F;

    value_type M_alpha;
};

/**
 * \param op an operator
 * \param s a scalar
 * \return the operator which is the scaling of op by s
 */
template<typename OpType>
boost::shared_ptr<OperatorScale<OpType>>
scale( boost::shared_ptr<OpType> M, typename OpType::value_type s )
{
    return boost::make_shared<OperatorScale<OpType>>(M,s) ;
}


/**
 * Diagonal of an operator
 */
template< typename op1_type>
class OperatorDiag : public OperatorBase<typename op1_type::value_type>
{
    typedef OperatorBase<typename op1_type::value_type> super;
public:

    typedef boost::shared_ptr<op1_type> op1_ptrtype;
    typedef typename op1_type::value_type value_type;


    OperatorDiag()
        :
        OperatorBase<value_type>()
        {}
    OperatorDiag( op1_ptrtype& F )
        :
        super( F->domainMapPtr(), F->imageMapPtr(), F->label(), F->useTranspose(), F->hasNormInf() ),
        M_D( backend()->newVector(F->domainMapPtr()) )
    {
        std::string temp( "diag(" );
        std::string t = F->label();


        temp.append( t );
        temp.append( ")" );        
        this->setLabel(temp);

        LOG(INFO) << "Create diag operator " << this->label() << " ...\n";
    }

    OperatorDiag( const OperatorDiag& tc )
        :
        super( tc ),
        M_D( tc.M_D )
    {
        LOG(INFO) << "Copy diag operator " << this->label() << " ...\n";
    }

    bool hasInverse() const
    {
        return true;
    }

    bool hasApply() const
    {
        return true;
    }



    int apply( const vector_type& X, vector_type & Y ) const
    {
        LOG(INFO) << "apply diag operator " << this->label() << "\n";
        Y.pointwiseMult( X, *M_D );
        
        return !hasApply();
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const
    {
        LOG(INFO) << "applyInverse diag operator " << this->label() << "\n";

        Y.pointwiseDivide( X, *M_D);
        
        return !hasInverse();
    }

    value_type normInf() const
        {
            return 0;
        }

    ~OperatorDiag()
    {
        //LOG(INFO) << "Destroyed scale operator: " << label() << " ...\n";
    };

private:

    vector_ptrtype M_D;
};

/**
 * \param op an operator
 * \param s a scalar
 * \return the operator which is the scaling of op by s
 */
template<typename OpType>
boost::shared_ptr<OperatorDiag<OpType>>
diag( boost::shared_ptr<OpType> M )
{
    return boost::make_shared<OperatorDiag<OpType>>(M) ;
}

} // Feel

#endif /* FEELPP_FEELAG_OPERATOR_HPP */
