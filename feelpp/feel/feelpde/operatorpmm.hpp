/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 Date: 02 Oct 2023

 Copyright (C) 2014-2016 Feel++ Consortium

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
#ifndef FEELPP_OPERATORPMM_HPP
#define FEELPP_OPERATORPMM_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pdh.hpp>
#include <feel/feelpde/boundaryconditions.hpp>

namespace Feel
{
template<typename T>
class OperatorPMM : public OperatorBase<T>
{
public:
    using super = OperatorBase<T>;
    using value_type = T;
    using value_t = T;
    using type = OperatorPMM<T>;
    typedef std::shared_ptr<type> ptrtype;

    using sparse_matrix_type = backend_type::sparse_matrix_type;
    using sparse_matrix_ptrtype = backend_type::sparse_matrix_ptrtype;

    using vector_type = backend_type::vector_type;
    using vector_ptrtype = backend_type::vector_ptrtype;

    using op_mat_type = OperatorMatrix<value_type>;
    using op_mat_ptrtype = std::shared_ptr<op_mat_type>;

    using op_inv_type = OperatorInverse<op_mat_type>;
    using op_inv_ptrtype = std::shared_ptr<op_inv_type>;

    using op_type = OperatorBase<value_type>;
    using op_ptrtype = std::shared_ptr<op_type>;

    static inline const uint16_type Dim = space_type::nDim;


    OperatorPMM( space_ptrtype Qh,
                 backend_ptrtype b,
                 BoundaryConditions const& bcFlags,
                 std::string const& p,
                 bool acc = false, bool applyInPETSc = false );

    OperatorPMM( const OperatorPMM& tc ) = default;
    OperatorPMM( OperatorPMM&& tc ) = default;
    OperatorPMM& operator=( const OperatorPMM& tc ) = default;
    OperatorPMM& operator=( OperatorPMM&& tc ) = default;

    void initialize();

    void setProblemType( std::string prob_type )
        {
            M_prob_type = prob_type;
        }

    std::string problemType() const
        {
            return M_prob_type;
        }

    ~OperatorPMM() {};

    sparse_matrix_ptrtype pressureMassMatrix() const override { return M_mass; }

    void setParameterValues( std::map<std::string,double> const& pv );

    int apply(const vector_type& X, vector_type& Y) const override;
    int applyInverse(const vector_type& X, vector_type& Y) const override;


private:
    backend_ptrtype M_b;
    space_ptrtype M_Xh;
    pressure_space_ptrtype M_Qh;

    pressure_element_type p;

    sparse_matrix_ptrtype M_mass;

    op_mat_ptrtype massOp;

    std::string M_prefix;

    op_ptrtype precOp;

    std::string M_prob_type;

    bool M_applyInPETSc;

    void assembleMass();
};




template < typename space_type>
OperatorPMM<space_type>::OperatorPMM( space_ptrtype Qh,
                                      backend_ptrtype b,
                                      std::string const& prefix,
                                      bool applyInPETSc )
    :
    super( Qh->template functionSpace<1>()->mapPtr(), "PMM", false, false ),
    M_b( b),
    M_Xh( Qh ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    p( M_Qh, "p" ),
    M_mass( backend()->newMatrix(_test=M_Qh,_trial=M_Qh) ),
    M_prefix( prefix ),
    M_applyInPETSc( applyInPETSc )
{
    initialize();

    this->assembleMass();
}
template < typename space_type>
void
OperatorPMM<space_type>::initialize()
{
}

template < typename space_type>
void
OperatorPMM<space_type>::setParameterValues( std::map<std::string,double> const& pv )
{
}

template < typename space_type>
void
OperatorPMM<space_type>::assembleMass()
{
    tic();
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(p) );
    M_mass->close();
    massOp = op( M_mass, "Mp" );
    toc("OperatorPMM::mass assembly",FLAGS_v>0);
}

template < typename space_type>
int
OperatorPMM<space_type>::apply(const vector_type& X, vector_type& Y) const
{
    return precOp->apply( X, Y );
}
template < typename space_type>
int
OperatorPMM<space_type>::applyInverse(const vector_type& X, vector_type& Y) const
{
    return precOp->applyInverse( X, Y );
}

} // Feel

#endif
