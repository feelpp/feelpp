/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
      Date: 30 Jul 2024

 Copyright (C) 2014-2024 Feel++ Consortium

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
#include <feel/feelpde/operatorpmmbase.hpp>

namespace Feel
{

template <typename space_type>
class OperatorPMM : public OperatorPMMBase<typename space_type::value_type>
{
    typedef OperatorPMMBase<typename space_type::value_type> super;

  public:
    typedef OperatorPMM<space_type> type;
    typedef std::shared_ptr<type> ptrtype;

    typedef typename space_type::value_type value_type;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef std::shared_ptr<space_type> space_ptrtype;

    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    using pressure_element_type = element_type;
    typedef OperatorMatrix<value_type> op_mat_type;
    typedef std::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef std::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorBase<value_type> op_type;
    typedef std::shared_ptr<op_type> op_ptrtype;

    static inline const uint16_type Dim = space_type::nDim;
    static inline const uint16_type pOrder = space_type::basis_type::nOrder;

    OperatorPMM( space_ptrtype Qh,
                 backend_ptrtype b,
                 std::string const& p,
                 bool applyInPETSc = false );

    OperatorPMM( const OperatorPMM& tc ) = default;
    OperatorPMM( OperatorPMM&& tc ) = default;
    OperatorPMM& operator=( const OperatorPMM& tc ) = default;
    OperatorPMM& operator=( OperatorPMM&& tc ) = default;

    void initialize();

    template <typename ExprMu>
    void update( ExprMu const& expr_mu );

    ~OperatorPMM() {};

    sparse_matrix_ptrtype pressureMassMatrix() const override { return M_mass; }
    void setParameterValues( std::map<std::string, double> const& pv );

    int apply( const vector_type& X, vector_type& Y ) const override;
    int applyInverse( const vector_type& X, vector_type& Y ) const override;

  private:
    backend_ptrtype M_b;
    space_ptrtype M_Qh;
    pressure_element_type p;

    sparse_matrix_ptrtype M_mass;
    vector_ptrtype M_rhs;

    op_mat_ptrtype massOp;

    std::string M_prefix;

    op_ptrtype precOp;

    bool M_applyInPETSc;

    void assembleMass();
};

template <typename space_type>
OperatorPMM<space_type>::OperatorPMM( space_ptrtype Qh,
                                      backend_ptrtype b,
                                      std::string const& p,
                                      bool applyInPETSc )
    : super( Qh->mapPtr(), "PMM", false, false ),
      M_b( b ),
      M_Qh( Qh ),
      p( M_Qh, "p" ),
      M_mass( backend()->newMatrix( _test = M_Qh, _trial = M_Qh ) ),
      M_rhs( backend()->newVector( M_Qh ) ),
      M_prefix( p ),
      M_applyInPETSc( applyInPETSc )
{
    initialize();

    this->assembleMass();
}
template <typename space_type>
void OperatorPMM<space_type>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template <typename space_type>
void OperatorPMM<space_type>::setParameterValues( std::map<std::string, double> const& pv )
{
}

template <typename space_type>
template <typename ExprMu>
void OperatorPMM<space_type>::update( ExprMu const& expr_mu )
{
    tic();
    if ( !M_applyInPETSc && !precOp )
    {
        // S = F G^-1 M
        LOG( INFO ) << "[OperatorPMM] setting pmm operator...\n";
        precOp = massOp;
        LOG( INFO ) << "[OperatorPMM] setting pmm operator done.\n";
        // init_G = true;
    }
    toc( "Operator::PMM update", FLAGS_v > 0 );
}

template <typename space_type>
void OperatorPMM<space_type>::assembleMass()
{
    tic();
    auto m = form2( _test = M_Qh, _trial = M_Qh, _matrix = M_mass );
    m = integrate( elements( M_Qh->mesh() ), idt( p ) * id( p ) );
    M_mass->close();
    if ( !M_applyInPETSc )
        massOp = op( M_mass, "Mp" );
    toc( "OperatorPMM::mass assembly", FLAGS_v > 0 );
}


template <typename space_type>
int OperatorPMM<space_type>::apply( const vector_type& X, vector_type& Y ) const
{
    return precOp->apply( X, Y );
}
template <typename space_type>
int OperatorPMM<space_type>::applyInverse( const vector_type& X, vector_type& Y ) const
{
    return precOp->applyInverse( X, Y );
}

namespace Alternatives
{

template <typename SpacePressureType>
class OperatorPMM : public OperatorPMMBase<typename SpacePressureType::value_type>
{
public:

    static inline const int Dim = SpacePressureType::nDim;

    using value_type = typename SpacePressureType::value_type;
    using super = OperatorPMMBase<value_type>;
    using type = OperatorPMM<SpacePressureType>;
    using space_pressure_type = SpacePressureType;
    typedef std::shared_ptr<type> ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef std::shared_ptr<space_pressure_type> space_pressure_ptrtype;

    using mesh_type = typename space_pressure_type::mesh_type;
    using mesh_ptrtype = typename space_pressure_type::mesh_ptrtype;
    using element_pressure_type = typename space_pressure_type::element_type;
    using element_pressure_ptrtype = typename space_pressure_type::element_ptrtype;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef std::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef std::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorBase<value_type> op_type;
    typedef std::shared_ptr<op_type> op_ptrtype;

    using range_type = Range<mesh_type, MESH_ELEMENTS>;
    using range_faces_type = Range<mesh_type, MESH_FACES>;


    OperatorPMM( space_pressure_ptrtype Ph,
                 backend_ptrtype b,
                 std::string const& p,
                 bool applyInPETSc = false );

    OperatorPMM( const OperatorPMM& tc ) = default;
    OperatorPMM( OperatorPMM&& tc ) = default;
    OperatorPMM& operator=( const OperatorPMM& tc ) = default;
    OperatorPMM& operator=( OperatorPMM&& tc ) = default;
    ~OperatorPMM() {};

    void initialize();
    void assemble();

    void updateStart();

    void updateFinish();

    sparse_matrix_ptrtype pressureMassMatrix() const override { return M_mass; }

    int apply( const vector_type& X, vector_type& Y ) const override;
    int applyInverse( const vector_type& X, vector_type& Y ) const override;

private:
    void assembleMass();
private:
    backend_ptrtype M_b;
    space_pressure_ptrtype M_Ph;
    element_pressure_ptrtype M_p;

    sparse_matrix_ptrtype M_mass;

    op_mat_ptrtype massOp;

    op_ptrtype precOp;

    bool M_applyInPETSc;


};

template <typename SpacePressureType>
OperatorPMM<SpacePressureType>::OperatorPMM( space_pressure_ptrtype Ph,
                                                                backend_ptrtype b,
                                                                std::string const& prefix,
                                                                bool applyInPETSc )
    : super( Ph->mapPtr(), "PMM", false, false ),
      M_b( b ),
      M_Ph( Ph ),
      M_applyInPETSc( applyInPETSc )
{
}

template <typename SpacePressureType>
void OperatorPMM<SpacePressureType>::initialize()
{
    M_p = M_Ph->elementPtr();
    M_mass = M_b->newMatrix( _test = M_Ph, _trial = M_Ph );
    this->assemble();
}

template <typename SpacePressureType>
void OperatorPMM<SpacePressureType>::assemble()
{
    this->assembleMass();
}

template <typename SpacePressureType>
void OperatorPMM<SpacePressureType>::updateStart()
{
}

template <typename SpacePressureType>
void OperatorPMM<SpacePressureType>::updateFinish()
{
    tic();
    if ( !M_applyInPETSc && !precOp )
    {
        // S = F G^-1 M
        LOG( INFO ) << "[OperatorPMM] setting pmm operator...\n";
        precOp = massOp;
        LOG( INFO ) << "[OperatorPMM] setting pmm operator done.\n";
    }
    toc("Operator::PMM updateFinish", FLAGS_v > 0);
}

template <typename SpacePressureType>
void OperatorPMM<SpacePressureType>::assembleMass()
{
    tic();
    auto rangeElt = M_Ph->template rangeElements<0>();
    auto m = form2( _test = M_Ph, _trial = M_Ph, _matrix = M_mass );
    m = integrate( _range = rangeElt, _expr = idt( M_p ) * id( M_p ) );
    M_mass->close();
    if ( !M_applyInPETSc )
        massOp = op( M_mass, "Mp" );
    toc( "OperatorPMM::mass assembly", FLAGS_v > 0 );
}

template <typename SpacePressureType>
int OperatorPMM<SpacePressureType>::apply( const vector_type& X, vector_type& Y ) const
{
    return precOp->apply( X, Y );
}
template <typename SpacePressureType>
int OperatorPMM<SpacePressureType>::applyInverse( const vector_type& X, vector_type& Y ) const
{
    return precOp->applyInverse( X, Y );
}

} // namespace Alternatives

} // namespace Feel

#endif
