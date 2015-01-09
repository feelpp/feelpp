/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-04-26

  Copyright (C) 2013 Feel++ Consortium

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
   \file fsfunctionallinearfree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-26
 */
#ifndef __FSFUNCTIONALLINEARCOMPOSITE_H
#define __FSFUNCTIONALLINEARCOMPOSITE_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/fsfunctional.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{

template<class Space>
class FsFunctionalLinearComposite : public FsFunctionalLinear<Space>
{
public:

    // -- TYPEDEFS --
    typedef FsFunctionalLinearComposite<Space> this_type;
    typedef FsFunctionalLinear<Space> super_type;
    typedef boost::shared_ptr<super_type> super_ptrtype;

    typedef Space space_type;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    typedef typename space_type::value_type value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    FsFunctionalLinearComposite( space_ptrtype space ) :
        super_type( space ),
        M_backend( backend_type::build( BACKEND_PETSC ) )
    {}

    FsFunctionalLinearComposite( space_ptrtype space, backend_ptrtype backend ) :
        super_type( space ),
        M_backend( backend )
    {}

    virtual ~FsFunctionalLinearComposite() {}

    int size()
    {
        int size1 = M_functionals1.size();
        int size2 = M_functionals2.size();
        return size1+size2;
    }


    std::vector<int> countAllContributions()
    {
        int size1 = M_functionals1.size();
        int size2 = M_functionals2.size();

        if( size1==0 )
        {
            //first : count Q
            int Q=1;
            if( size2 == 0 )
                Q=0;
            //initialization
            auto it_=M_functionals2.begin();
            auto tuple_=it_->first;
            int old_q = tuple_.template get<0>();
            //loop over all operators
            auto end = M_functionals2.end();
            for(auto it=M_functionals2.begin(); it!=end; it++)
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
            for(auto it=M_functionals2.begin(); it!=end; it++)
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
        }//end of case functional2
        else
        {
            int Q = this->size();
            std::vector<int> V(Q);
            auto it_=M_functionals1.begin();
            auto tuple_=it_->first;
            auto end = M_functionals1.end();
            for(auto it=M_functionals1.begin(); it!=end; it++)
            {
                int q = it->first;
                V[q]=1;
            }
            return V;
        }//end of case functional1
    }


    //if we have a list of functionals
    //i.e. \sum_{q=0}^Q Fq(.)
    void addElement( super_ptrtype const& functional )
    {
        int size = M_functionals1.size();
        M_functionals1.insert( std::make_pair( size , functional ) ) ;
    }

    void addElement( int q, super_ptrtype const& functional )
    {
        M_functionals1.insert( std::make_pair( q , functional ) ) ;
    }

    void addList( std::vector< super_ptrtype > const& vec )
    {
        int q_max=vec.size();
        for(int q=0; q<q_max; q++)
            this->addElement( q , vec[q] );
    }


    //if we have a list of list of functionals
    //i.e. \sum_{q=0}^Q \sum_{m=0}^M F_{qm}(.)
    void addElement(  boost::tuple<int,int> const& qm , super_ptrtype const& functional )
    {
        M_functionals2.insert( std::make_pair( qm , functional ) );
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
        int size1 = M_functionals1.size();
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
                CHECK(msize==1)<<"Error ! You should use a vector of double to call setScalars() in your case.\n";
                new_scalars[q]=scalars[q][0];
            }
            M_scalars1=new_scalars;
        }

    }

    this_type& operator=( this_type const& m )
    {
        M_backend = m.M_backend;
        M_functionals1 = m.M_functionals1;
        M_functionals1 = m.M_functionals2;
        M_scalars1 = m.M_scalars1;
        M_scalars2 = m.M_scalars2;
        return *this;
    }


    void sumAllVectors( vector_ptrtype & vector, bool use_scalar_one=false ) const
    {
        int size1 = M_functionals1.size();
        int size2 = M_functionals2.size();
        bool size_error=false;
        if( size1 > 0 && size2 > 0 )
            size_error=true;
        if( (size1 + size2) == 0 )
            size_error=true;

        FEELPP_ASSERT( !size_error )( size1 )( size2 ).error( "FsFunctionalLinearComposite has no elements, or both maps have elements" );

        if( size1 > 0 )
            sumAllVectors1( vector, use_scalar_one );
        else
            sumAllVectors2( vector, use_scalar_one );
    }


    void sumAllVectors1( vector_ptrtype & vector, bool use_scalar_one=false ) const
    {
        int size1 = M_functionals1.size();
        int size2 = M_functionals2.size();

        FEELPP_ASSERT( size1 > 0 )( size1 )( size2 ).error( "FsFunctionalLinearComposite has no elements" );

        vector->zero();
        auto temp_vector = M_backend->newVector( this->space() );
        auto end = M_functionals1.end();
        for(auto it=M_functionals1.begin(); it!=end; ++it)
        {
            int position = it->first;
            double scalar=1;
            if( ! use_scalar_one )
                scalar = M_scalars1[position];
            it->second->containerPtr(temp_vector);
            vector->add( scalar , temp_vector );
        }

    }//sumAllVectors

    void sumAllVectors2( vector_ptrtype & vector, bool use_scalar_one=false ) const
    {
        int size1 = M_functionals1.size();
        int size2 = M_functionals2.size();

        FEELPP_ASSERT( size2 > 0 )( size1 )( size2 ).error( "FsFunctionalLinearComposite has no elements" );

        vector->zero();
        auto temp_vector = M_backend->newVector( this->space() );
        auto end = M_functionals2.end();
        for(auto it=M_functionals2.begin(); it!=end; ++it)
        {
            auto position = it->first;
            int q = position.template get<0>();
            int m = position.template get<1>();
            double scalar=1;
            if( ! use_scalar_one )
                scalar = M_scalars2[q][m];
            it->second->containerPtr(temp_vector);
            vector->add( scalar , temp_vector );
        }

    }//sumAllVectors


    super_ptrtype& functionallinear( int q )
    {
        int q_max = M_functionals1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "FsFunctionalLinearComposite has not enough elements" );
        return M_functionals1.template at(q);
    }

    super_ptrtype& functionallinear( int q , int m)
    {
        auto tuple = boost::make_tuple( q , m );
        return M_functionals2.template at(tuple);
    }


    //Access to a specific vector
    void vecPtr( int q , vector_ptrtype& vector )
    {
        int q_max = M_functionals1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "FsFunctionalLinearComposite has not enough elements" );
        //auto vector = M_backend->newVector( this->space() );
        M_functionals1.template at(q)->containerPtr(vector);
        //return vector;
    }

    void vecPtr( int q , vector_ptrtype& vector ) const
    {
        int q_max = M_functionals1.size();
        FEELPP_ASSERT( q < q_max )( q_max )( q ).error( "FsFunctionalLinearComposite has not enough elements" );
        M_functionals1.template at(q)->containerPtr(vector);
    }

    //Access to a specific vector
    void vecPtr( int q , int m , vector_ptrtype& vector )
    {
        auto tuple = boost::make_tuple( q , m );
        M_functionals2.template at(tuple)->containerPtr(vector);
    }

    void vecPtr( int q , int m , vector_ptrtype& vector ) const
    {
        auto tuple = boost::make_tuple( q , m );
        M_functionals2.template at(tuple)->containerPtr(vector);
    }


    //return the sum of all vectors
    virtual void containerPtr( vector_ptrtype & vector_to_fill )
    {
        //vector_to_fill = sumAllVectors(true);
        sumAllVectors( vector_to_fill , true );
        vector_to_fill->close();
    }

    // apply the functional
    virtual value_type
    operator()( const element_type& x ) const
    {
        //auto vector = sumAllVectors( true );
        auto vector = M_backend->newVector( this->space() );
        sumAllVectors( vector, true );
        vector->close();
        return M_backend->dot( *vector, x.container() );
    }


private:

    std::map< int, super_ptrtype > M_functionals1;
    std::map< boost::tuple<int,int> , super_ptrtype > M_functionals2;
    backend_ptrtype M_backend;
    std::vector< double > M_scalars1;
    std::vector< std::vector<double> > M_scalars2;
};

namespace detail
{

template<typename Args>
struct compute_functionalLinearComposite_return
{
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::space>::type>::type::element_type space_type;

    typedef FsFunctionalLinearComposite<space_type> type;
    typedef boost::shared_ptr<FsFunctionalLinearComposite<space_type> > ptrtype;
};
}

BOOST_PARAMETER_FUNCTION(
    ( typename Feel::detail::compute_functionalLinearComposite_return<Args>::ptrtype ), // 1. return type
    functionalLinearComposite,                        // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required
      ( space,    *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
    ) // required
    ( optional
      ( backend,        *, backend() )
    ) // optionnal
)
{

    Feel::detail::ignore_unused_variable_warning( args );
    typedef typename Feel::detail::compute_functionalLinearComposite_return<Args>::type functionalcomposite_type;
    typedef typename Feel::detail::compute_functionalLinearComposite_return<Args>::ptrtype functionalcomposite_ptrtype;
    return functionalcomposite_ptrtype ( new functionalcomposite_type( space , backend ) );

} // functionalLinearComposite

}//Feel


#endif /* _FSFUNCTIONALLINEARCOMPOSITE_HPP_ */
