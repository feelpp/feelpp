/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): JB Wahl <wahl.jb@gmail.com>
 Date: 24 nov. 2015

 Copyright (C) 2015 Feel++ Consortium

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
#include <feel/feelts/cnab2.hpp>

#ifndef FEELPP_CNAB2X2_HPP
#define FEELPP_CNAB2X2_HPP 1

#define CN(X) fusion::at_c<X>( M_cnab2_vec )

namespace Feel{

template<typename Element1, typename Element2>
class CNAB2x2 : public TimeSet
{
public :
    typedef typename fusion::vector<Element1,Element2> elements_vector_t;
    typedef typename fusion::vector< CNAB2<Element1>,CNAB2<Element2> > cnab2_vector_t;

    template <int N>
    using field_t = typename mpl::at_c<elements_vector_t,N>::type;


    CNAB2x2() = default;

    template< typename ElementType1, typename ElementType2 >
    CNAB2x2( ElementType1 element1, ElementType2 element2 ) :
        steady( boption( _name="cnab2.steady" ) ),
        keps( doption( _name="cnab2.keps" ) ),
        steady_tol( doption( _name="cnab2.steady-tol" ) ),
        T( doption(_name="cnab2.time-final") ),
        M_cnab2_vec( {element1,element2} ),
        nstar( ioption(_name="cnab2.nstar") ),
        M_rejectedError(0),
        M_rejectedSolver(0)
        {}

    template<int N>
    field_t<N> const& rateOfChange() const
        { return CN(N).rateOfChange(); }

    template<int N>
    field_t<N> const& field() const
        { return CN(N).field(); }

    template<int N>
    field_t<N> const& extrapolateVelocity()
        { return CN(N).extrapolateVelocity(); }


    template<typename FT1, typename FT2>
    void addStep( double k1, std::pair<FT1,FT1> elements1, std::pair<FT2,FT2> elements2 )
        {
            this->push_back( k1, this->size()?t()+k1:k1 );
            CN(0).addStep( k1, elements1.first, elements1.second );
            CN(1).addStep( k1, elements2.first, elements2.second );
        }

    template<typename ElementType1, typename ElementType2, typename BCType1, typename BCType2>
    void start( ElementType1& element1, ElementType2& element2, BCType1& bc1, BCType2& bc2 )
        {
            next( element1, element2, bc1, bc2 );
        }

    template<typename ElementType1, typename ElementType2, typename BCType1, typename BCType2>
    bool next( ElementType1& element1, ElementType2& element2, BCType1& bc1, BCType2& bc2 , bool is_converged=true );

    bool isFinished() const
        {
            bool b1 = CN(0).isFinished();
            bool b2 =  CN(1).isFinished();
            bool is_finished = b1 && b2;
            if (is_finished)
                displayStats();
            return is_finished;
        }

    void displayStats() const
        {
            if (Environment::isMasterRank())
            {
                std::cout << "CNAB2x2 stats :"
                          << " -> Number of rejected TS (error)= " << M_rejectedError
                          << " -> Number of rejected TS (solver)= " << M_rejectedSolver
                          << std::endl;
            }
        }

    std::pair<bool,double> computeStep( double step, double err1, double err2 )
        {
            double e = std::sqrt( err1*err1 + err2*err2 );
            return std::make_pair( !(e > std::pow(1./0.7,3.)*keps), step*std::pow(keps/e,1./3.) );
        }

private :
    bool steady;
    double keps,steady_tol, T;

    cnab2_vector_t M_cnab2_vec;

    int nstar;
    int M_rejectedError;
    int M_rejectedSolver;

}; // class CNAB2x2


template<typename Element1, typename Element2>
template<typename ElementType1, typename ElementType2, typename BCType1, typename BCType2>
bool
CNAB2x2<Element1,Element2>::next( ElementType1& element1, ElementType2& element2, BCType1& bc1, BCType2& bc2, bool is_converged )

{
    if ( is_converged )
    {
        if ( Environment::isMasterRank() )
            std::cout << "trying next step (index:" << this->index() << "), with kn1=" << this->k() << " kn=" << this->kprev(1) << " at t=" << this->t() << std::endl;

        double err1 = CN(0).computeError( element1 );
        double err2 = CN(1).computeError( element2 );

        auto ktry = this->computeStep( k(), err1, err2 );

        if ( Environment::isMasterRank() )
            std::cout << " --> k_n="<< this->kprev(1) << ", k_{n+1}=" << this->k()
                      << ", err1="<<err1 <<", err2="<<err2 << ", keps=" <<keps
                      << ", accepted="<<ktry.first << ", proposed kn2="<<ktry.second
                      << std::endl;

        if ( ktry.first )
        {
            double k1 = ktry.second;
            if ( (index() > 0) && (index() % nstar == 0) )
            {
                double tstar = tprev(1);
                double kstar = kprev(1);

                kprev(1) = 0.5*kstar;
                tprev(1) = tstar-kprev(1);

                kstar = k();
                t() = tstar+.5*k();
                k() = kprev(1) + 0.5*kstar;

                CN(0).averaging(bc1, k1 );
                CN(1).averaging(bc2, k1 );
            }
            else
            {
                CN(0).addStep( k1 );
                CN(1).addStep( k1 );
            }
            this->push_back( k1, this->t()+k1 );
        }
        else
        {
            CN(0).eraseStep(ktry.second);
            CN(1).eraseStep(ktry.second);
            this->k() = ktry.second;
            this->t() = tprev(1)+this->k();
            M_rejectedError++;
        }
        is_converged =  ktry.first;

    }
    else
    {
        double new_k = this->kprev(1);
        CN(0).eraseStep(new_k);
        CN(1).eraseStep(new_k);
        this->k() = new_k;
        this->t() = tprev(1)+this->k();

        M_rejectedSolver++;
    }

    return is_converged;
}





} // namespace Feel

#endif
