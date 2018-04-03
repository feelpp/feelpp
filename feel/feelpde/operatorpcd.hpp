/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
            Goncalo Pena  <gpena@mat.uc.pt>
 Date: 02 Oct 2014

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
#ifndef FEELPP_OPERATORPCD_HPP
#define FEELPP_OPERATORPCD_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feeldiscr/pdh.hpp>

namespace Feel
{

template<typename space_type, typename PropertiesSpaceType = Pdh_type<typename space_type::mesh_type,0>>
class OperatorPCD : public OperatorBase<typename space_type::value_type>
{
    typedef OperatorBase<typename space_type::value_type> super;
public:

    typedef OperatorPCD<space_type, PropertiesSpaceType> type;
    typedef boost::shared_ptr<type> ptrtype;

    typedef typename space_type::value_type value_type;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;

    using properties_space_type = PropertiesSpaceType;
    using properties_space_ptrtype = boost::shared_ptr<properties_space_type>;
    using property_type = typename properties_space_type::element_type;
    typedef boost::shared_ptr<property_type> property_ptrtype;

    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type velocity_space_type;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype pressure_space_ptrtype;
    typedef typename velocity_space_type::element_type velocity_element_type;
    typedef typename pressure_space_type::element_type pressure_element_type;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef boost::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef boost::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorCompose<op_inv_type, op_mat_type> comp1_type;
    typedef boost::shared_ptr<comp1_type> comp1_ptrtype;

    typedef OperatorCompose<op_mat_type, comp1_type> comp2_type;
    typedef boost::shared_ptr<comp2_type> comp2_ptrtype;

    typedef super op_type;
    typedef boost::shared_ptr<op_type> op_ptrtype;


    static const uint16_type Dim = space_type::nDim;
    static const uint16_type pOrder = pressure_space_type::basis_type::nOrder;

    OperatorPCD( space_ptrtype Qh,
                 sparse_matrix_ptrtype A,
                 backend_ptrtype b,
                 BoundaryConditions const& bcFlags,
                 std::string const& p,
                 property_ptrtype mu,
                 property_ptrtype rho,
                 property_ptrtype alpha,
                 bool acc);

    OperatorPCD( const OperatorPCD& tc ) = default;
    OperatorPCD( OperatorPCD&& tc ) = default;
    OperatorPCD& operator=( const OperatorPCD& tc ) = default;
    OperatorPCD& operator=( OperatorPCD&& tc ) = default;

    void initialize();

    template < typename ExprConvection, typename ExprBC >
    void update( ExprConvection const& expr_b, ExprBC const& ebc, bool hasConvection=true , double tn = 0., double tn1 = 0 );

    void setProblemType( std::string prob_type )
        {
            M_prob_type = prob_type;
        }

    std::string problemType() const
        {
            return M_prob_type;
        }

    ~OperatorPCD() {};

    int apply(const vector_type& X, vector_type& Y) const;
    int applyInverse(const vector_type& X, vector_type& Y) const;


    // void setAlpha( property_ptrtype a ) { M_alpha = a; }
    // void setMu( property_ptrtype m ) { M_mu = m; }
    // void setRho( property_ptrtype r ) { M_rho = r; }
private:
    backend_ptrtype M_b;
    space_ptrtype M_Xh;
    velocity_space_ptrtype M_Vh;
    pressure_space_ptrtype M_Qh;

    velocity_element_type u, v;
    pressure_element_type p, q;

    sparse_matrix_ptrtype G;
    form2_type<pressure_space_type,pressure_space_type> form2_conv;
    sparse_matrix_ptrtype M_mass, M_diff, M_conv, M_B, M_massv_inv, M_Bt;
    vector_ptrtype rhs;

    op_mat_ptrtype massOp, diffOp, convOp;

    BoundaryConditions M_bcFlags;
    std::string M_prefix;

    op_ptrtype precOp;

    property_ptrtype M_mu, M_rho, M_alpha;

    std::string M_prob_type;

    bool M_accel;

    void assembleMass();

    void assembleDiffusion();

    void applyBC( sparse_matrix_ptrtype& A );


};




template < typename space_type, typename PropertiesSpaceType>
OperatorPCD<space_type, PropertiesSpaceType>::OperatorPCD( space_ptrtype Qh,
                                      sparse_matrix_ptrtype Bt,
                                      backend_ptrtype b,
                                      BoundaryConditions const& bcFlags,
                                      std::string const& p,
                                      property_ptrtype mu,
                                      property_ptrtype rho,
                                      property_ptrtype alpha,
                                      bool acc)
    :
    super( Qh->template functionSpace<1>()->mapPtr(), "PCD", false, false ),
    M_b( b),
    M_Xh( Qh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    G( backend()->newMatrix(M_Qh, M_Qh) ),
    form2_conv( "pcdconv", M_Qh, M_Qh, G ),
    M_mass( backend()->newMatrix(M_Qh, M_Qh) ),
    M_diff( backend()->newMatrix(M_Qh, M_Qh) ),
    M_conv( backend()->newMatrix(M_Qh, M_Qh) ),
    M_B( backend()->newMatrix(_trial=M_Vh, _test=M_Qh) ),
    M_massv_inv( backend()->newMatrix(_trial=M_Vh, _test=M_Vh) ),
    M_Bt( Bt ),
    rhs( backend()->newVector( M_Qh ) ),
    M_bcFlags( bcFlags ),
    M_prefix( p ),
    M_mu( mu ),
    M_rho( rho ),
    M_alpha( alpha ),
    M_accel( acc )
{
    initialize();

    LOG(INFO) << "[Pressure Correction Diffusion Operator] Constructor: using mu=" << M_mu << "\n";
    LOG(INFO) << "[Pressure Correction Diffusion Operator] Constructor: using rho=" << M_rho << "\n";
    LOG(INFO) << "[Pressure Correction Diffusion Operator] Constructor: using alpha=" << M_alpha << "\n";

    this->assembleMass();
    this->assembleDiffusion();
}
template < typename space_type, typename PropertiesSpaceType>
void
OperatorPCD<space_type, PropertiesSpaceType>::initialize()
{
    rhs->zero();
    rhs->close();
}

template < typename space_type, typename PropertiesSpaceType>
template < typename ExprConvection, typename ExprBC >
void
OperatorPCD<space_type, PropertiesSpaceType>::update( ExprConvection const& expr_b,
                                 ExprBC const& ebc,
                                 bool hasConvection,
                                 double tn,
                                 double tn1 )
{
    tic();

    double time_step = M_accel?tn1-tn:1;
    tic();
    form2_conv = integrate( _range=elements(M_Qh->mesh()), _expr=idv(*M_mu)*gradt(p)*trans(grad(q)));
    toc("OperatorPCD::update apply diffusion",FLAGS_v>0);

    if ( hasConvection )
    {
        tic();
        form2_conv += integrate( _range=elements(M_Qh->mesh()), _expr=(trans(expr_b)*trans(gradt(p)))*id(q));
        toc("OperatorPCD::update apply convection",FLAGS_v>0);
    }

    if ( soption("blockns.pcd.inflow") == "Robin" )
    {
        tic();
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
            LOG(INFO) << "Setting Robin condition on " << dir.marker();
            if ( ebc.find( dir.marker() ) != ebc.end() )
            {
                if ( M_accel )
                {
                    auto en = ebc.find(dir.marker())->second;
                    en.setParameterValues( { { "t", tn } } );
                    auto en1 = ebc.find(dir.marker())->second;
                    en1.setParameterValues( { { "t", tn1 } } );

                    form2_conv += integrate( _range=markedfaces(M_Qh->mesh(), dir.meshMarkers()),
                                             _expr=-idv(*M_rho)*trans((en1-en)/time_step)*N()*idt(p)*id(q));
                }
                else
                {
                    form2_conv += integrate( _range=markedfaces(M_Qh->mesh(), dir.meshMarkers()),
                                             _expr=-idv(*M_rho)*trans(ebc.find(dir.marker())->second)*N()*idt(p)*id(q));
                }
            }
        }
        toc("OperatorPCD::update apply Robin",FLAGS_v>0);
    }

    //auto alphaMax = normLinf( _range=elements(M_Qh->mesh()), _expr=idv(M_alpha), _pset=_Q<5>() );
    //if ( alphaMax.value() > 1e-15 )
    if ( M_alpha->linftyNorm() > 1e-15 )
    {
        LOG(INFO) << "[OperatorPCD] Add mass matrix...\n";
        tic();
        form2_conv += integrate( _range=elements(M_Qh->mesh()), _expr=idv(*M_alpha)/time_step*idt(p)*id(q) );
        toc("OperatorPCD::update apply mass",FLAGS_v>0);
    }

    G->close();
    tic();
    if ( soption("blockns.pcd.inflow") != "Robin" )
        for( auto dir : M_bcFlags[M_prefix]["Dirichlet"])
        {
            form2_conv += on( markedfaces(M_Qh->mesh(),dir.meshMarkers()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }

    // on neumann boundary on velocity, apply Dirichlet condition on pressure
    if ( soption("blockns.pcd.outflow") == "Dirichlet" )
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
            form2_conv += on( markedfaces(M_Qh->mesh(),cond.meshMarkers()), _element=p, _rhs=rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
    toc("OperatorPCD::update apply on()",FLAGS_v>0);

    this->applyBC(G);

    static bool init_G = false;

    //if ( !init_G )
    {
        // S = F G^-1 M
        LOG(INFO) << "[OperatorPCD] setting pcd operator...\n";
        if ( ioption("blockns.pcd.order") == 1 )
            precOp = compose( massOp, compose(inv(op(G,"Fp")),diffOp) );
        else
            precOp = compose( diffOp, compose(inv(op(G,"Fp")),massOp) );
        LOG(INFO) << "[OperatorPCD] setting pcd operator done.\n";
        init_G = true;
    }
    toc("Operator::PCD update",FLAGS_v>0);
}





template < typename space_type, typename PropertiesSpaceType>
void
OperatorPCD<space_type, PropertiesSpaceType>::assembleMass()
{
    tic();
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(q) );
    M_mass->close();
    massOp = op( M_mass, "Mp" );
    toc("OperatorPCD::mass assembly",FLAGS_v>0);
}

template < typename space_type, typename PropertiesSpaceType>
void
OperatorPCD<space_type, PropertiesSpaceType>::assembleDiffusion()
{
    tic();
    if ( soption("blockns.pcd.diffusion") == "Laplacian" )
    {
        auto d = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_diff );
        d = integrate( _range=elements(M_Qh->mesh()), _expr=gradt(p)*trans(grad(q)));
        LOG(INFO) << "blockns.pcd.diffusion is Laplacian";
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
            LOG(INFO) << "Diffusion Setting Dirichlet condition on pressure on " << cond.marker();
            if ( boption("blockns.weakdir" ) )
                d+= integrate( markedfaces(M_Qh->mesh(),cond.meshMarkers()),
                               _expr=-gradt(p)*N()*id(p)-grad(p)*N()*idt(p)+doption("penaldir")*idt(p)*id(p)/hFace() );
            else
                d += on( markedfaces(M_Qh->mesh(),cond.meshMarkers()), _element=p, _rhs=rhs,
                         _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
        M_diff->close();
     }
    if ( soption("blockns.pcd.diffusion") == "BTBt" )
    {
        tic();
        std::vector<size_type> M_Vh_indices( M_Vh->nLocalDofWithGhost() );
        std::vector<size_type> M_Qh_indices( M_Qh->nLocalDofWithGhost() );
        std::iota( M_Vh_indices.begin(), M_Vh_indices.end(), 0 );
        //std::iota( M_Qh_indices.begin(), M_Qh_indices.end(),
        //M_Vh->nLocalDofWithGhost() );
        //M_B = M_A->createSubMatrix( M_Qh_indices, M_Vh_indices );
        toc(" - OperatorPCD Extracted B" );
        tic();
        auto m = form2( _test=M_Vh, _trial=M_Vh );
        m = integrate( elements(M_Vh->mesh()), trans(idt(u))*id(v) );
        m.matrixPtr()->close();
        toc(" - OperatorPCD Velocity Mass Matrix" );
        tic();
        auto d = M_b->newVector( M_Vh );
        M_b->diag( m.matrixPtr(), d );
        d->reciprocal();
        M_b->diag( d, M_massv_inv );
        M_massv_inv->close();
        toc(" - OperatorPCD inverse diagonal mass matrix extracted" );
        tic();
        M_b->PtAP( M_massv_inv, M_Bt, M_diff );
        M_diff->close();
        toc(" - OperatorPCD B T^-1 B^T built");
        if ( Environment::numberOfProcessors() == 1 )
            M_diff->printMatlab( "BTBt.m" );
    }

    diffOp = op( M_diff, "Ap" );
    toc("OperatorPCD::diffusion assembly",FLAGS_v>0);
}


template < typename space_type, typename PropertiesSpaceType>
void
OperatorPCD<space_type, PropertiesSpaceType>::applyBC( sparse_matrix_ptrtype& A )
{
}

template < typename space_type, typename PropertiesSpaceType>
int
OperatorPCD<space_type, PropertiesSpaceType>::apply(const vector_type& X, vector_type& Y) const
{
    return precOp->apply( X, Y );
}
template < typename space_type, typename PropertiesSpaceType>
int
OperatorPCD<space_type, PropertiesSpaceType>::applyInverse(const vector_type& X, vector_type& Y) const
{
    return precOp->applyInverse( X, Y );
}


} // Feel

#endif
