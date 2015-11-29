/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 28 sept. 2015

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
#ifndef FEELPP_CNAB2_HPP
#define FEELPP_CNAB2_HPP 1

#include <vector>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelts/timeset.hpp>


namespace Feel {

template<typename FieldType>
class CNAB2 : public TimeSet
{
public:
    using super = TimeSet;
    using field_t = FieldType;
    using mesh_t = typename field_t::functionspace_type::mesh_type;
    using mesh_ptr_t = typename field_t::functionspace_type::mesh_ptrtype;

    CNAB2() = default;
    template<typename FT>
    CNAB2( FT const& f );


    field_t const& rateOfChange() const { return acc; }
    field_t const& field() const { return u; }

    template<typename FT>
    void addStep( double k1, FT const& f_u, FT const& f_acc )
    {
        this->push_back( k1, this->size()?t()+k1:k1 );
        up = u;
        accp=acc;
        u = f_u;
        acc = f_acc;
    }

    /**
     * @return the l2 error of the increment relative to the latest solution
     */
    double l2Error() const
        {
            return normL2(_range=elements(mesh), _expr=idv(u)-idv(up) )
                /normL2(_range=elements(mesh), _expr=idv(u)) /k();
        }

    /**
     * @return true if the adaptive algorithm is finished, false otherwise
     */
    bool isFinished() const
    {
        double e = 1;
        if ( steady )
            e = l2Error();
        if ( Environment::isMasterRank() )
        {
            if ( steady )
                std::cout << "checking for steady state ||u-up||_2 = "
                          << e << " tolerance: " << steady_tol << std::endl;
        }
        return steady ? (e<steady_tol && index()>10) : !(t()<T);
    }
    /**
     * compute new time step using time \p step and error \p err
     * @return a pair of bool and double, bool states whether the step is accepted (true)   * or not (false) and the double stores the new timestep in case of acceptance or
     * rejection
     */
    std::pair<bool,double> computeStep( double step, double err )
        {
            return std::make_pair( !(err > std::pow(1./0.7,3.)*keps), step*std::pow(keps/err,1./3.) );
        };

    void acceptedStep( double k )
        {
            this->addStep( k, u_try, acc_try );
        }

    void eraseStep( double k )
        {
            this->k() = k;
            this->t() = tprev(1)+this->k();
        }

    template<typename FT>
    void start( FT const& d );

    template<typename FT>
    bool next( FT const& d, bool is_converged=true );

    template<typename FT>
    std::pair<double,bool> computeError( FT const& d );


    field_t const& extrapolateVelocity();

    bool steady;
    double keps,k0,T, steady_tol;
    int nstar;

    mesh_ptr_t mesh;
    field_t u, up, acc, accp, d, w, uab2;
    field_t u_try, acc_try, accstar, ustar;
};

template<typename FieldType>
template<typename FT>
CNAB2<FieldType>::CNAB2( FT const& f )
    :
    steady( boption( _name="cnab2.steady" ) ),
    keps( doption( _name="cnab2.keps" ) ),
    k0( doption(_name="cnab2.time-initial") ),
    T( doption(_name="cnab2.time-final") ),
    steady_tol( doption( _name="cnab2.steady-tol" ) ),
    nstar( ioption(_name="cnab2.nstar") ),
    mesh( f.functionSpace()->mesh() ),
    u ( f.functionSpace() ),
    up ( f.functionSpace() ),
    acc ( f.functionSpace() ),
    accp ( f.functionSpace() ),
    d ( f.functionSpace() ),
    w ( f.functionSpace() ),
    uab2 ( f.functionSpace() ),
    u_try ( f.functionSpace() ),
    acc_try ( f.functionSpace() ),
    accstar( f.functionSpace() ),
    ustar( f.functionSpace() )
{
    if ( Environment::isMasterRank() )
        std::cout << "--> k0=" << k0 << " T=" << T << " n*=" << nstar << " keps=" << keps  << std::endl;

}
template<typename FieldType>
template<typename FT>
void
CNAB2<FieldType>::start( FT const& f_d )
{
    next( f_d );
}

template<typename FieldType>
template<typename FT>
std::pair<double,bool>
CNAB2<FieldType>::computeError( FT const& fd )
{
    bool averaging=false;
    if ( (index() > 0) && (index() % nstar == 0) )
    {
        double tstar = tprev(1);
        ustar.on(_range=elements(mesh), _expr=idv(u));
        accstar.on(_range=elements(mesh), _expr=idv(acc));
        tprev(1) = tprev(2) + kprev(1)/2;
        u.on(_range=elements(mesh), _expr=0.5*(idv(up)+idv(ustar)));
        acc.on(_range=elements(mesh), _expr=0.5*(idv(accp)+idv(accstar)));

        // t_{n+1}
        t() = tstar + k()/2;
        u_try.on(_range=elements(mesh), _expr=idv(ustar)+k()*idv(fd)/2);
        acc_try.on(_range=elements(mesh), _expr=idv(fd));
        k() = t()-tprev(1);
        kprev(1) = tprev(1)-tprev(2);
        if ( Environment::isMasterRank() )
            std::cout << " --> averaging t_{n}=" << tprev(1) << " t={n+1}=" << t() << std::endl;
        averaging = true;
    }
    else
    {
        // d^n is obtained, now updated u^n+1 and acc^n+1
        u_try.on(_range=elements(mesh), _expr=idv(u)+k()*idv(fd));
        acc_try.on(_range=elements(mesh), _expr=2*idv(fd)-idv(acc));
    }
    // compute AB2 velocity
    uab2.on( _range=elements(mesh), _expr=idv(u)+(k()/2.)*( (2+k()/kprev(1))*idv(acc)-(k()/kprev(1))*idv(accp) ) );

    // compute error
    return std::make_pair(normL2( _range=elements(mesh), _expr=(idv(u_try)-idv(uab2)))/(3.*(1+kprev(1)/k())), averaging);

}


template<typename FieldType>
template<typename FT>
bool
CNAB2<FieldType>::next( FT const& fd, bool is_converged )
{
    if ( is_converged )
    {
        if ( Environment::isMasterRank() )
            std::cout << "trying next step (index:" << index() << ")with kn1" << k() << " kn=" << kprev(1) << " at t=" << t() << std::endl;
        // every nstar iteration do averaging to avoid ringing (and time step stagnation)

        double err = computeError( fd ).first;

        double l2_u = normL2( _range=elements(mesh), _expr=idv(u_try));
        double l2_uab2 = normL2( _range=elements(mesh), _expr=idv(uab2));

        auto ktry = this->computeStep( k(), err );

        if ( Environment::isMasterRank() )
            std::cout << " --> k_n="<< kprev(1) << " k_{n+1}=" << k()
                      << " err=" << err << "(" << l2_u << "," << l2_uab2 << "," << (3.*(1+kprev(1)/k()))
                      << ")"  << " keps=" << keps
                      << " accepted= " << ktry.first << " proposed kn2=" << ktry.second
                      << std::endl;

        if ( ktry.first )
        {
            // accept the tried step and move forward
            this->addStep( ktry.second, u_try, acc_try );
        }
        else
        {
            // do not accept current kn1, overwrite it with a smaller step
            // adjust time to new smaller time step
            this->k() = ktry.second;
            this->t() = tprev(1)+this->k();
        }

        is_converged = is_converged && ktry.first;
    }
    else
    {
        // in this case solver did not converged and so we reduce the time step

        double new_k = this->k()/2. ;
        this->k() = new_k;
        this->t() = tprev(1)+this->k();
    }



    return is_converged;
}

template<typename FieldType>
typename CNAB2<FieldType>::field_t const&
CNAB2<FieldType>::extrapolateVelocity()
{
    // compute convection velocity (2nd order accurate)
    w.on(_range=elements(mesh), _expr=(1+k()/kprev(1))*idv(u) - (k()/kprev(1))*idv(up) );
    return w;
}


} // namespace Feel
#endif
