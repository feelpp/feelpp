/* -*- mode: c++ -*-

 This file is part of the Feel library

 Author(s): JB Wahl <wahl.jb@gmail.com>
 Date: 2018-04-04

 Copyright (C) 2018 Feel++ Consortium

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
 * @file   crbmodelaero.hpp
 * @author JB Wahl
 * @date   2018
 */
#ifndef CRBMODELAERO_H
#define CRBMODELAERO_H

#include <feel/feelmor/crbmodelsaddlepoint.hpp>

namespace Feel
{
template<typename ModelType>
class CRBModelAero :
        public CRBModelSaddlePoint<ModelType>
{
    typedef CRBModelSaddlePoint<ModelType> super;
public :
    typedef ModelType model_type;
    typedef std::shared_ptr<ModelType> model_ptrtype;

    typedef typename model_type::value_type value_type;
    typedef typename model_type::parameter_type parameter_type;

    typedef typename ModelType::mesh_type mesh_type;
    typedef typename ModelType::mesh_ptrtype mesh_ptrtype;

    typedef typename model_type::space_type space_type;
    typedef typename model_type::element_type element_type;
    typedef typename model_type::beta_vector_type beta_vector_type;

    CRBModelAero( crb::stage stage, int level = 0 ) :
        CRBModelAero( std::make_shared<model_type>(), stage, level )
        {}

    CRBModelAero( model_ptrtype const& model , crb::stage stage, int level = 0 ) :
        super ( model, stage, level ),
        M_use_psit( boption(_prefix=this->M_prefix,_name="crb.aero.use-psit")),
        M_delta( doption(_prefix=this->M_prefix,_name="crb.aero.psit.delta0" ) ),
        M_rez(-1),
        M_QTri( 0 ),
        M_set_0_mean_pressure( boption(_prefix=this->M_prefix,_name="crb.aero.fix-mean-pressure") )
        {
            this->M_addSupremizer = boption(_prefix=this->M_prefix,_name="crb.aero.add-supremizer");
            if ( stage == crb::stage::offline )
                this->initTrilinear();
        }
    element_type solve( parameter_type const& mu ) override
        { return this->M_model->solve(mu); }

    void initTrilinear();
    element_type offlineSolveAD( parameter_type const& mu );
    int QTri() { return M_QTri; }
    sparse_matrix_ptrtype assembleTrilinearOperator( element_type const& u, int const& q )
        { return this->M_model->assembleTrilinearOperator( u,q ); }
    std::vector<double> computeBetaTri( parameter_type const& mu )
        { return this->M_model->computeBetaTri(mu); }

    bool isTrilinear() const override
        { return true; }

    int sizeOfBilinearJ() const
        { return this->M_model->sizeOfBilinearJ()==-1 ? this->Qa():this->M_model->sizeOfBilinearJ(); }
    int sizeOfLinearR() const
        { return this->M_model->sizeOfLinearR()==-1 ? this->Ql(0):this->M_model->sizeOfLinearR(); }



private:
    void updateJacobianAD( const vector_ptrtype& X, sparse_matrix_ptrtype & J , const parameter_type & mu);
    void updateResidualAD( const vector_ptrtype& X, vector_ptrtype& R , const parameter_type & mu);
    void postSolve( vector_ptrtype rhs, vector_ptrtype sol );
    double updateR( vector_ptrtype const& X, parameter_type const& mu )
        {
            auto R = backend()->newVector(this->functionSpace());
            updateResidualAD(X, R, mu );
            return R->l2Norm();
        }
    void updatePsiT( const vector_ptrtype& X, parameter_type const& mu );
private:
    bool M_use_psit;
    double M_delta, M_rez;
    int M_QTri;
    sparse_matrix_ptrtype M_J;
    vector_ptrtype M_R;

    bool M_set_0_mean_pressure;
    vector_ptrtype M_mean_p_vec;

    using super::M_Jqm;
    using super::M_Rqm;
}; // class CRBModelAero

template<typename ModelType>
void
CRBModelAero<ModelType>::initTrilinear()
{
    LOG(INFO) <<"CRB Model Aero : init Trilinear\n";
    auto mu = this->M_model->parameterSpace()->element();
    auto beta_tri = this->M_model->computeBetaTri( mu );
    M_QTri = beta_tri.size();
}

template<typename ModelType>
typename CRBModelAero<ModelType>::element_type
CRBModelAero<ModelType>::offlineSolveAD( parameter_type const& mu )
{
    if ( !M_J )
        M_J = this->newMatrix();
    if ( !M_R )
        M_R = this->newVector();
    auto initialguess = this->assembleInitialGuess( mu );

    if ( M_set_0_mean_pressure )
    {
        auto mesh = this->functionSpace()->mesh();
        auto q = this->functionSpace()->template functionSpace<1>()->element();
        M_mean_p_vec = form1_mean( _test=this->functionSpace()->template functionSpace<1>(),
                                   _range=elements(mesh), _expr=id(q) ).vectorPtr();
    }

    boost::tie( boost::tuples::ignore , M_Jqm, M_Rqm ) = this->computeAffineDecomposition();

    this->M_backend_primal = backend( _name="backend-primal", _rebuild=true );

    this->M_backend_primal->nlSolver()->jacobian = boost::bind( &CRBModelAero<ModelType>::updateJacobianAD, boost::ref( *this ), _1, _2, mu );
    this->M_backend_primal->nlSolver()->residual = boost::bind( &CRBModelAero<ModelType>::updateResidualAD, boost::ref( *this ), _1, _2, mu );
    auto post_solve = boost::bind( &CRBModelAero<ModelType>::postSolve, boost::ref( *this ), _1, _2 );

    if ( M_use_psit )
    {
        M_rez = -1;
        M_delta = doption(_prefix=this->M_prefix,_name="crb.aero.psit.delta0");
    }

    auto U = this->functionSpace()->element();
    this->M_backend_primal->nlSolve( _jacobian=M_J, _residual=M_R, _solution=U, _post=post_solve );

    return U;
}

template<typename ModelType>
void
CRBModelAero<ModelType>::updateJacobianAD( const vector_ptrtype& X, sparse_matrix_ptrtype & J , const parameter_type & mu)
{
    J->zero();
    auto Xh = this->functionSpace();
    auto U = Xh->element();
    U = *X;

    beta_vector_type betaJqm;
    boost::tie( boost::tuples::ignore, betaJqm, boost::tuples::ignore ) = this->computeBetaQm( U , mu , 0 );
    auto betaTri = this->M_model->computeBetaTri( mu );

    // trilinear part
    for ( int q=0; q<betaTri.size(); q ++ )
    {
        auto JTri = this->M_model->assembleTrilinearJ( U, q );
        J->addMatrix( betaTri[q], JTri );
    }

    // rest : bilinear, non-linear+eim
    for ( size_type q = 0; q < this->Qa(); ++q )
    {
        for(int m=0; m<this->mMaxA(q); m++)
        {
            J->addMatrix( betaJqm[q][m], M_Jqm[q][m] );
        }
    }

    if ( M_use_psit )
    {
        auto u = U.template element<0>();
        auto T = U.template element<2>();
        updatePsiT( X, mu );
        auto norm_u=max(1e-12,sqrt(inner(idv(u),idv(u))));
        auto range = u.functionSpace()->template rangeElements<0>();
        form2( _test=Xh, _trial=Xh, _matrix=J )
            += integrate( range, norm_u/hMin()/M_delta*(inner(idt(u),id(u)) + inner(id(T),idt(T))) );
    }
    J->close();
}


template<typename ModelType>
void
CRBModelAero<ModelType>::updateResidualAD( const vector_ptrtype& X, vector_ptrtype& R , const parameter_type & mu)
{
    R->zero();
    auto temp = backend()->newMatrix( this->functionSpace(),this->functionSpace() );
    auto U = this->functionSpace()->element();
    U = *X;
    beta_vector_type betaJqm;
    std::vector<beta_vector_type> betaRqm;
    boost::tie( boost::tuples::ignore, betaJqm, betaRqm ) = this->computeBetaQm( U, mu );

    auto betaTri = this->M_model->computeBetaTri( mu );

    // trilinear part
    for ( int q=0; q<M_QTri; q ++ )
    {
        auto RTri = this->M_model->assembleTrilinearR( U, q );
        R->add( betaTri[q], RTri );
    }

    // bilinear part
    temp->zero();
    for ( size_type q = 0; q < this->sizeOfBilinearJ(); ++q )
    {
        temp->addMatrix( betaJqm[q][0], M_Jqm[q][0] );
    }
    R->addVector( X, temp );

    // linear part and non-linear+eim
    for ( size_type q = 0; q < this->Ql( 0 ); ++q )
    {
        for( int m=0; m<this->mMaxF(0,q); m++ )
        {
            R->add( betaRqm[0][q][m] , M_Rqm[0][q][m] );
        }

    }
    R->close();
}

template<typename ModelType>
void
CRBModelAero<ModelType>::updatePsiT( const vector_ptrtype& X, parameter_type const& mu )
{
    double new_rez = updateR( X, mu );
    if ( M_rez==-1 )
        M_rez=new_rez;
    M_delta = M_delta*M_rez/new_rez;
    M_rez = new_rez;
}

template<typename ModelType>
void
CRBModelAero<ModelType>::postSolve( vector_ptrtype rhs, vector_ptrtype sol )
{
    if( M_set_0_mean_pressure )
    {
        auto U = this->functionSpace()->element( sol, 0 );
        auto p = U.template element<1>();
        double meanPressureCurrent = inner_product( *M_mean_p_vec, p );
        p.add( -meanPressureCurrent );
    }
}


} // namespace Feel

#endif
