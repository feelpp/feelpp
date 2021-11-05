/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-04-26

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file operatorlinearfree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-26
 */
#ifndef FEELPP_DISCR_OPERATORLINEARCOMPOSITE_H
#define FEELPP_DISCR_OPERATORLINEARCOMPOSITE_H 1

#include <feel/feelcore/pslogger.hpp>
#include <feel/feeldiscr/operatorlinear.hpp>

namespace Feel
{

template<class DomainSpace, class DualImageSpace>
class OperatorLinearComposite : public OperatorLinear<DomainSpace, DualImageSpace>
{

    typedef OperatorLinear<DomainSpace,DualImageSpace> super;
    typedef std::shared_ptr<super> super_ptrtype;

public :

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::dual_image_space_type  dual_image_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename super::dual_image_space_ptrtype  dual_image_space_ptrtype;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;

    typedef typename super::matrix_type matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;
    typedef typename super::domain_element_type domain_element_type;
    typedef typename super::dual_image_element_type dual_image_element_type;

    typedef typename super::domain_element_range_type domain_element_range_type;
    typedef typename super::domain_element_slice_type domain_element_slice_type;
    typedef typename super::dual_image_element_range_type dual_image_element_range_type;
    typedef typename super::dual_image_element_slice_type dual_image_element_slice_type;

    typedef OperatorLinearComposite<DomainSpace, DualImageSpace> this_type;
    typedef std::shared_ptr<this_type> this_ptrtype;

    OperatorLinearComposite (domain_space_ptrtype     domainSpace,
                             dual_image_space_ptrtype dualImageSpace,
                             backend_ptrtype          backend,
                             size_type                pattern=Pattern::COUPLED)
        :
        super( domainSpace,dualImageSpace,backend , false ),
        M_backend( backend ),
        M_pattern( pattern )
    {}

    ~OperatorLinearComposite() override {}


    //if we have a list of operators
    //i.e. \sum_{q=0}^Q Aq(.,.)
    void addElement( super_ptrtype const& op )
    {
        int size = M_operators1.size();
        M_operators1.insert( std::make_pair( size , op ) ) ;
    }


    void addElement( int q, super_ptrtype const& op )
    {
        M_operators1.insert( std::make_pair( q , op ) ) ;
    }

    void addList( std::vector< super_ptrtype > const& vec )
    {
        int q_max=vec.size();
        for(int q=0; q<q_max; q++)
            this->addElement( q , vec[q] );
    }

    void displayOperatorsNames()
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();

        LOG( INFO ) << " the composite operator linear "<<this->name()<<" has following operators : ";

        if( size1  > 0 )
        {
            auto end = M_operators1.end();
            for(auto it=M_operators1.begin(); it!=end; it++)
                LOG(INFO)<<it->second->name();
        }
        else
        {
            auto end = M_operators2.end();
            for(auto it=M_operators2.begin(); it!=end; it++)
                LOG(INFO)<<it->second->name();
        }
    }

    int size()
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();
        return size1+size2;
    }




    //return number of terms
    //Q-terms (i.e. without using EIM in a CRB context)
    //and M-terms (i.e. using EIM in a CRB context)
    //Q-terms is an int
    //and then for each term (from 0 to Q-terms) we count number of sub-terms
    //finally we fill a vector to store them
    //vector V of size Q-terms
    //V[1] = number of sub-terms associated to q=1
    //if V[3]=2 that means there is 2 sub-terms associated to q=3
    std::vector<int> countAllContributions()
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();

        if( size1 == 0 )
        {
            //first : count Q
            int Q=1;
            if( size2 == 0 )
                Q=0;
            //initialization
            auto it_=M_operators2.begin();
            auto tuple_=it_->first;
            int old_q = tuple_.template get<0>();
            //loop over all operators
            auto end = M_operators2.end();
            for(auto it=M_operators2.begin(); it!=end; it++)
            {
                auto tuple = it->first;
                int q = tuple.template get<0>();
                if (q!=old_q)
                {
                    Q++;
                    old_q=q;
                }
            }
            std::vector<int> V(Q);
            //now count sub-terms
            int count=0;
            old_q = tuple_.template get<0>();
            for(auto it=M_operators2.begin(); it!=end; it++)
            {
                auto tuple = it->first;
                int q = tuple.template get<0>();
                if (q!=old_q)
                {
                    V[old_q]=count;
                    count=1;
                    old_q=q;
                }
                else
                {
                    count++;
                }
            }
            V[old_q]=count;
            return V;
        }//end of case operator2
        else
        {
            int Q = this->size();
            std::vector<int> V(Q);
            auto it_=M_operators1.begin();
            auto tuple_=it_->first;
            auto end = M_operators1.end();
            for(auto it=M_operators1.begin(); it!=end; it++)
            {
                int q = it->first;
                V[q]=1;
            }
            return V;
        }//end of case operator1
    }

    //for a given index q, return the number
    //of operators that have q in the tuple
    //Warning : it works only with M_operator2 !
    int subSize(int q)
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();
        CHECK( size1 == 0 )<<"the function subSize can only be called when using operators with 2 index\n";

        int count=0;
        auto end = M_operators2.end();
        for(auto it=M_operators2.begin(); it!=end; it++)
        {
            auto tuple = it->first;
            int q_ = tuple.template get<0>();
            int m_ = tuple.template get<1>();
            if (q_==q)
                count++;
        }
        return count;
    }

    //if we have a list of list of operators
    //i.e. \sum_{q=0}^Q \sum_{m=0}^M A_{qm}(.,.)
    void addElement(  boost::tuple<int,int> const& qm , super_ptrtype const& op )
    {
        M_operators2.insert( std::make_pair( qm , op ) );
    }

    void addList( std::vector< std::vector< super_ptrtype > > const& vec )
    {
        int q_max = vec.size();
        for(int q=0; q<q_max; q++)
        {
            int m_max = vec[q].size();
            for(int m=0; m<m_max; m++)
            {
                auto tuple = boost::make_tuple( q , m );
                this->addElement( tuple , vec[q][m] );
            }
        }
    }

    void setScalars( std::vector<double> scalars )
    {
        M_scalars1 = scalars;
    }

    void setScalars( std::vector< std::vector< double > > scalars )
    {
        int size1 = M_operators1.size();
        if( size1 == 0 )
        {
            M_scalars2 = scalars;
        }
        else
        {
            int qsize=scalars.size();
            std::vector< double > new_scalars( qsize );
            for(int q=0; q<qsize; q++)
            {
                int msize=scalars[q].size();
                CHECK(msize==1)<<"Error ! You should use a vector of double to call setScalars(), or use a vector of vector of matrices.\n";
                new_scalars[q]=scalars[q][0];
            }
            M_scalars1=new_scalars;
        }
    }

    void sumAllMatrices(matrix_ptrtype & matrix, bool use_scalar_one=false ) const
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();
        bool size_error=false;
        if( size1 > 0 && size2 > 0 )
            size_error=true;
        if( (size1 + size2) == 0 )
            size_error=true;

        FEELPP_ASSERT( !size_error )( size1 )( size2 ).error( "OperatorLinearComposite has no elements, or both maps have elements" );

        if( size1 > 0 )
            sumAllMatrices1( matrix, use_scalar_one );
        else
            sumAllMatrices2( matrix, use_scalar_one );
    }

    //return the sum of matrices given
    //by all opertors in M_vectors
    //arguments : a vector of scalars and a bool ( use scalar=1 if true )
    void sumAllMatrices1( matrix_ptrtype & matrix, bool use_scalar_one=false ) const
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();

        FEELPP_ASSERT( size1 > 0 )( size1 )( size2 ).error( "OperatorLinearComposite has no elements" );

        matrix->zero();
        auto temp_matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );

        auto end = M_operators1.end();
        for(auto it=M_operators1.begin(); it!=end; ++it)
        {
            int position = it->first;
            double scalar=1;
            if( ! use_scalar_one )
                scalar = M_scalars1[position];
            it->second->matPtr(temp_matrix);
            matrix->addMatrix( scalar , temp_matrix );
        }


    }//sumAllMatrices

    //return the sum of matrices given
    //by all opertors in M_vectors
    //arguments : a vector of vector of scalars and a bool ( use scalar=1 if true )
    void sumAllMatrices2( matrix_ptrtype & matrix, bool use_scalar_one=false ) const
    {
        int size1 = M_operators1.size();
        int size2 = M_operators2.size();

        FEELPP_ASSERT( size2 > 0 )( size1 )( size2 ).error( "OperatorLinearComposite has no elements" );

        matrix->zero();
        auto temp_matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );

        auto end = M_operators2.end();
        for(auto it=M_operators2.begin(); it!=end; ++it)
        {
            auto position = it->first;
            int q = position.template get<0>();
            int m = position.template get<1>();
            double scalar=1;
            if( ! use_scalar_one )
                scalar = M_scalars2[q][m];
            it->second->matPtr(temp_matrix);
            matrix->addMatrix( scalar , temp_matrix );
        }

        temp_matrix.reset();
    }//sumAllMatrices


    //Access to a specific matrix
    void matrixPtr(int q , matrix_ptrtype& matrix)
    {
        int q_max = M_operators1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "OperatorLinearComposite has not enough elements" );
        //fill matrix
        M_operators1.at(q)->matPtr(matrix);
    }

    void matrixPtr(int q , matrix_ptrtype& matrix) const
    {
        int q_max = M_operators1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "OperatorLinearComposite has not enough elements" );
        //fill matrix
        M_operators1.at(q)->matPtr(matrix);
    }

    //Access to a specific matrix
    void matrixPtr(int q, int m , matrix_ptrtype& matrix)
    {
        //fill matrix
        auto tuple = boost::make_tuple(q,m);
        M_operators2.at(tuple)->matPtr(matrix);
    }

    void matrixPtr(int q, int m , matrix_ptrtype& matrix) const
    {
        auto tuple = boost::make_tuple(q,m);
        M_operators2.at(tuple)->matPtr(matrix);
    }


    //Acces to operator
    super_ptrtype & operatorlinear(int q)
    {
        int q_max = M_operators1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "OperatorLinearComposite has not enough elements" );
        return M_operators1.at(q);
    }

    super_ptrtype & operatorlinear(int q, int m)
    {
        auto tuple = boost::make_tuple(q,m);
        return M_operators2.at(tuple);
    }


    void
    apply( const domain_element_type& de,
           image_element_type&        ie ) const override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.space() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    double
    energy( const typename domain_space_type::element_type & de,
            const typename dual_image_space_type::element_type & ie ) const override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        *_v2 = ie;
        vector_ptrtype _v3( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v3 );
        return inner_product( _v2, _v3 );
    }


    void
    apply( const typename domain_space_type::element_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices(  matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    void
    apply( const domain_element_range_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices(matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const typename domain_space_type::element_type & de,
           dual_image_element_range_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_range_type & de,
           dual_image_element_range_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_slice_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    void
    apply( const typename domain_space_type::element_type & de,
           dual_image_element_slice_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( /*const*/ domain_element_slice_type /*&*/ de,
                     dual_image_element_slice_type /*&*/ ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }



    void
    apply( const domain_element_range_type & de,
           dual_image_element_slice_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_slice_type & de,
           dual_image_element_range_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    //! apply the inverse of the operator: \f$de = O^{-1} ie\f$
    void
    applyInverse( domain_element_type&      de,
                  const image_element_type& ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace(), _trial=this->domainSpace(), _pattern=M_pattern );
        sumAllMatrices( matrix, true );
        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        vector_ptrtype _v2( M_backend->newVector( _test=ie.space() ) );
        *_v2 = ie.container();
        M_backend->solve( matrix, matrix, _v1, _v2 );
        de = *_v1;
    }

    this_type& operator=( this_type const& m )
    {
        M_backend = m.M_backend;
        M_pattern = m.M_pattern;
        M_operators1 = m.M_operators1;
        M_operators2 = m.M_operators2;
        M_scalars1 = m.M_scalars1;
        M_scalars2 = m.M_scalars2;
        return *this;
    }



private :

    std::map< int, super_ptrtype > M_operators1;
    std::map< boost::tuple<int,int> , super_ptrtype > M_operators2;
    backend_ptrtype M_backend;
    size_type M_pattern;
    std::vector< double > M_scalars1;
    std::vector< std::vector<double> > M_scalars2;
};//class OperatorLinearComposite


template <typename ... Ts>
auto opLinearComposite( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.get(_domainSpace);
    auto && imageSpace = args.get(_imageSpace);
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    size_type pattern = args.get_else( _pattern, Pattern::COUPLED );

    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    return std::make_shared<OperatorLinearComposite<domain_space_type, image_space_type>>( domainSpace,imageSpace,backend,pattern );
}


}//namespace Feel
#endif

