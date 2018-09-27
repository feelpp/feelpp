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


#include <feel/feelpde/operatorpcdbase.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/pdh.hpp>

namespace Feel
{

template<typename space_type>
class OperatorPCD : public OperatorPCDBase<typename space_type::value_type>
{
    typedef OperatorPCDBase<typename space_type::value_type> super;
public:

    typedef OperatorPCD<space_type> type;
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
    typedef typename space_type::template sub_functionspace<0>::type velocity_space_type;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype pressure_space_ptrtype;
    typedef typename velocity_space_type::element_type velocity_element_type;
    typedef typename pressure_space_type::element_type pressure_element_type;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef std::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef OperatorInverse<op_mat_type> op_inv_type;
    typedef std::shared_ptr<op_inv_type> op_inv_ptrtype;

    typedef OperatorCompose<op_inv_type, op_mat_type> comp1_type;
    typedef std::shared_ptr<comp1_type> comp1_ptrtype;

    typedef OperatorCompose<op_mat_type, comp1_type> comp2_type;
    typedef std::shared_ptr<comp2_type> comp2_ptrtype;

    typedef OperatorBase<value_type> op_type;
    typedef std::shared_ptr<op_type> op_ptrtype;


    static const uint16_type Dim = space_type::nDim;
    static const uint16_type pOrder = pressure_space_type::basis_type::nOrder;

    OperatorPCD( space_ptrtype Qh,
                 backend_ptrtype b,
                 BoundaryConditions const& bcFlags,
                 std::string const& p,
                 bool acc = false, bool applyInPETSc = false );

    OperatorPCD( const OperatorPCD& tc ) = default;
    OperatorPCD( OperatorPCD&& tc ) = default;
    OperatorPCD& operator=( const OperatorPCD& tc ) = default;
    OperatorPCD& operator=( OperatorPCD&& tc ) = default;

    void initialize();

    template < typename ExprRho, typename ExprMu, typename ExprConvection, typename ExprAlpha >
    void update( ExprRho const& expr_rho, ExprMu const& expr_mu, ExprConvection const& expr_b, ExprAlpha const& expr_alpha,
                 bool hasConvection=true , bool hasAlpha=true, double tn = 0., double tn1 = 0 );
    template < typename ExprRho, typename ExprMu, typename ExprConvection, typename ExprBC, typename ExprAlpha >
    void update( ExprRho const& expr_rho, ExprMu const& expr_mu, ExprConvection const& expr_b, ExprBC const& ebc, ExprAlpha const& expr_alpha,
                 bool hasConvection=true , bool hasAlpha=true, double tn = 0., double tn1 = 0 );

    void setProblemType( std::string prob_type )
        {
            M_prob_type = prob_type;
        }

    std::string problemType() const
        {
            return M_prob_type;
        }

    ~OperatorPCD() {};

    sparse_matrix_ptrtype pressureMassMatrix() const override { return M_mass; }
    sparse_matrix_ptrtype pressureLaplacianMatrix() const override { return M_diff; }
    sparse_matrix_ptrtype pressureDiffusionConvectionMatrix() const override { return M_conv; }
    sparse_matrix_ptrtype velocityMassMatrix() const override { return M_massv; }

    BoundaryConditions const& bcFlags() const { return M_bcFlags; }

    int pcdOrder() const override { return M_pcdOrder; }
    std::string const& pcdDiffusionType() const override { return M_pcdDiffusionType; }

    void setBt( sparse_matrix_ptrtype Bt ) { M_Bt = Bt; }
    void setParameterValues( std::map<std::string,double> const& pv );

    int apply(const vector_type& X, vector_type& Y) const override;
    int applyInverse(const vector_type& X, vector_type& Y) const override;


private:
    backend_ptrtype M_b;
    space_ptrtype M_Xh;
    velocity_space_ptrtype M_Vh;
    pressure_space_ptrtype M_Qh;

    velocity_element_type u;
    pressure_element_type p;

    sparse_matrix_ptrtype M_mass, M_diff, M_conv, M_massv, M_massv_inv, M_Bt;
    vector_ptrtype M_rhs;

    op_mat_ptrtype massOp, diffOp, convOp;

    int M_pcdOrder;
    std::string M_pcdDiffusionType;
    bool M_pcdDiffusionLaplacianWeakDir;
    double M_pcdDiffusionLaplacianWeakDirParameter;
    std::string M_bcInflowType, M_bcOutflowType;
    BoundaryConditions M_bcFlags;
    map_vector_field<Dim,1,2> M_bcDirichlet;
    std::string M_prefix;

    op_ptrtype precOp;

    std::string M_prob_type;

    bool M_accel;
    bool M_applyInPETSc;

    void assembleMass();

    void assembleDiffusion();

};




template < typename space_type>
OperatorPCD<space_type>::OperatorPCD( space_ptrtype Qh,
                                      backend_ptrtype b,
                                      BoundaryConditions const& bcFlags,
                                      std::string const& p,
                                      bool acc, bool applyInPETSc )
    :
    super( Qh->template functionSpace<1>()->mapPtr(), "PCD", false, false ),
    M_b( b),
    M_Xh( Qh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    u( M_Vh, "u" ),
    p( M_Qh, "p" ),
    M_mass( backend()->newMatrix(_test=M_Qh,_trial=M_Qh) ),
    M_diff( backend()->newMatrix(_test=M_Qh,_trial=M_Qh) ),
    M_conv( backend()->newMatrix(_test=M_Qh,_trial=M_Qh) ),
    M_rhs( backend()->newVector( M_Qh ) ),
    M_pcdOrder( ioption("blockns.pcd.order") ),
    M_pcdDiffusionType( soption("blockns.pcd.diffusion") ),
    M_pcdDiffusionLaplacianWeakDir( boption("blockns.weakdir" ) ),
    M_pcdDiffusionLaplacianWeakDirParameter( doption("blockns.weakdir.penaldir") ),
    M_bcInflowType( soption("blockns.pcd.inflow") ),
    M_bcOutflowType( soption("blockns.pcd.outflow") ),
    M_bcFlags( bcFlags ),
    M_prefix( p ),
    M_accel( acc ),
    M_applyInPETSc( applyInPETSc )
{
    M_bcDirichlet = M_bcFlags.template getVectorFields<Dim>( std::string(M_prefix), "Dirichlet" );

    initialize();

    this->assembleMass();
    this->assembleDiffusion();
}
template < typename space_type>
void
OperatorPCD<space_type>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template < typename space_type>
void
OperatorPCD<space_type>::setParameterValues( std::map<std::string,double> const& pv )
{
    M_bcDirichlet.setParameterValues( pv );
}

template < typename space_type>
template < typename ExprRho, typename ExprMu, typename ExprConvection, typename ExprAlpha >
void
OperatorPCD<space_type>::update( ExprRho const& expr_rho, ExprMu const& expr_mu,
                                 ExprConvection const& expr_b,
                                 ExprAlpha const& expr_alpha,
                                 bool hasConvection, bool hasAlpha,
                                 double tn, double tn1 )
{
    this->update( expr_rho,expr_mu,expr_b,M_bcDirichlet,expr_alpha,hasConvection, hasAlpha,tn,tn1 );
}

template < typename space_type>
template < typename ExprRho, typename ExprMu, typename ExprConvection, typename ExprBC, typename ExprAlpha >
void
OperatorPCD<space_type>::update( ExprRho const& expr_rho, ExprMu const& expr_mu,
                                 ExprConvection const& expr_b,
                                 ExprBC const& ebc,
                                 ExprAlpha const& expr_alpha,
                                 bool hasConvection, bool hasAlpha,
                                 double tn, double tn1 )
{
    tic();

    auto form2_conv = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_conv );

    double time_step = M_accel?tn1-tn:1;
    tic();
    form2_conv = integrate( _range=elements(M_Qh->mesh()), _expr=expr_mu*gradt(p)*trans(grad(p)));
    toc("OperatorPCD::update apply diffusion",FLAGS_v>0);

    if ( hasConvection )
    {
        tic();
        form2_conv += integrate( _range=elements(M_Qh->mesh()), _expr=(trans(expr_b)*trans(gradt(p)))*id(p));
        toc("OperatorPCD::update apply convection",FLAGS_v>0);
    }

    if ( hasAlpha )
    {
        LOG(INFO) << "[OperatorPCD] Add mass matrix...\n";
        tic();
        form2_conv += integrate( _range=elements(M_Qh->mesh()), _expr=expr_alpha/time_step*idt(p)*id(p) );
        toc("OperatorPCD::update apply mass",FLAGS_v>0);
    }

    if ( M_bcInflowType == "Robin" )
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
                                             _expr=-expr_rho*trans((en1-en)/time_step)*N()*idt(p)*id(p));
                }
                else
                {
                    form2_conv += integrate( _range=markedfaces(M_Qh->mesh(), dir.meshMarkers()),
                                             _expr=-expr_rho*trans(ebc.find(dir.marker())->second)*N()*idt(p)*id(p));
                }
            }
        }
        toc("OperatorPCD::update apply Robin",FLAGS_v>0);
    }


    tic();
    std::set<std::string> markersApplyStrongDirichlet;
    if ( M_bcInflowType != "Robin" )
        for( auto const& dir : M_bcFlags[M_prefix]["Dirichlet"])
            markersApplyStrongDirichlet.insert( dir.meshMarkers().begin(),dir.meshMarkers().end() );
    // on neumann boundary on velocity, apply Dirichlet condition on pressure
    if ( M_bcOutflowType == "Dirichlet" )
        for( auto const& cond : M_bcFlags[M_prefix]["Neumann"])
            markersApplyStrongDirichlet.insert( cond.meshMarkers().begin(),cond.meshMarkers().end() );

    if ( !markersApplyStrongDirichlet.empty() )
        form2_conv += on( _range=markedfaces(M_Qh->mesh(),markersApplyStrongDirichlet),
                          _element=p, _rhs=M_rhs, _expr=cst(0.), _type="elimination_keep_diagonal" );
    toc("OperatorPCD::update apply on()",FLAGS_v>0);

    //this->applyBC(G);
    M_conv->close();
    //static bool init_G = false;

    //if ( !init_G )
    if ( !M_applyInPETSc && !precOp )
    {
        // S = F G^-1 M
        LOG(INFO) << "[OperatorPCD] setting pcd operator...\n";
        if ( M_pcdOrder == 1 )
            precOp = compose( massOp, compose(inv(op(M_conv,"Fp")),diffOp) );
        else
            precOp = compose( diffOp, compose(inv(op(M_conv,"Fp")),massOp) );
        LOG(INFO) << "[OperatorPCD] setting pcd operator done.\n";
        //init_G = true;
    }
    toc("Operator::PCD update",FLAGS_v>0);
}





template < typename space_type>
void
OperatorPCD<space_type>::assembleMass()
{
    tic();
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(p) );
    M_mass->close();
    if ( !M_applyInPETSc )
        massOp = op( M_mass, "Mp" );
    toc("OperatorPCD::mass assembly",FLAGS_v>0);
}

template < typename space_type>
void
OperatorPCD<space_type>::assembleDiffusion()
{
    tic();
    if ( M_pcdDiffusionType == "Laplacian" )
    {
        auto d = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_diff );
        d = integrate( _range=elements(M_Qh->mesh()), _expr=gradt(p)*trans(grad(p)));
        LOG(INFO) << "blockns.pcd.diffusion is Laplacian";
        for( auto cond : M_bcFlags[M_prefix]["Neumann"])
        {
            LOG(INFO) << "Diffusion Setting Dirichlet condition on pressure on " << cond.marker();
            if ( M_pcdDiffusionLaplacianWeakDir )
                d+= integrate( markedfaces(M_Qh->mesh(),cond.meshMarkers()),
                               _expr=-gradt(p)*N()*id(p)-grad(p)*N()*idt(p)+M_pcdDiffusionLaplacianWeakDirParameter*idt(p)*id(p)/hFace() );
            else
                d += on( markedfaces(M_Qh->mesh(),cond.meshMarkers()), _element=p, _rhs=M_rhs,
                         _expr=cst(0.), _type="elimination_keep_diagonal" );
        }
        M_diff->close();
     }
    if ( M_pcdDiffusionType == "BTBt" )
    {
        if ( M_applyInPETSc )
        {
            if ( !M_massv )
                M_massv = backend()->newMatrix(_trial=M_Vh, _test=M_Vh);
            auto m = form2( _test=M_Vh, _trial=M_Vh,_matrix=M_massv );
            m = integrate( _range=elements(M_Vh->mesh()), _expr=inner(idt(u),id(u)) );
            M_massv->close();
        }
        else
        {
            if ( !M_massv_inv )
                M_massv_inv = backend()->newMatrix(_trial=M_Vh, _test=M_Vh);
            tic();
            auto m = form2( _test=M_Vh, _trial=M_Vh );
            m = integrate( _range=elements(M_Vh->mesh()), _expr=inner(idt(u),id(u)) );
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
            M_diff->clear(); // stencil will change
            M_b->PtAP( M_massv_inv, M_Bt, M_diff );
            M_diff->close();
            toc(" - OperatorPCD B T^-1 B^T built");
            //if ( Environment::numberOfProcessors() == 1 )
            //    M_diff->printMatlab( "BTBt.m" );
        }
    }
    if ( !M_applyInPETSc )
        diffOp = op( M_diff, "Ap" );
    toc("OperatorPCD::diffusion assembly",FLAGS_v>0);
}

template < typename space_type>
int
OperatorPCD<space_type>::apply(const vector_type& X, vector_type& Y) const
{
    return precOp->apply( X, Y );
}
template < typename space_type>
int
OperatorPCD<space_type>::applyInverse(const vector_type& X, vector_type& Y) const
{
    return precOp->applyInverse( X, Y );
}


} // Feel

#endif
