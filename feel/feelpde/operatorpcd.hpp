/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
            Goncalo Pena  <gpena@mat.uc.pt>
 Date: 02 Oct 2014

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
#ifndef FEELPP_OPERATORPCD_HPP
#define FEELPP_OPERATORPCD_HPP 1


#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>

namespace Feel
{

template<typename space_type>
class OperatorPCD : public OperatorBase<typename space_type::value_type>
{
    typedef OperatorBase<typename space_type::value_type> super;
public:

    typedef OperatorPCD<space_type> type;
    typedef boost::shared_ptr<type> ptrtype;

    typedef typename space_type::value_type value_type;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type velocity_space_type;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype pressure_space_ptrtype;
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

    OperatorPCD(){}

    OperatorPCD( pressure_space_ptrtype Qh,
                 std::map< std::string, std::set<flag_type> > bcFlags,
                 double nu,
                 double alpha );

    OperatorPCD( const OperatorPCD& tc );

    void initialize();

    template < typename ExprConvection, typename ExprBC >
    void update( ExprConvection const& expr_b, ExprBC const& ebc );

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

private:

    pressure_space_ptrtype M_Qh;

    pressure_element_type p, q;

    sparse_matrix_ptrtype M_mass, M_diff, M_conv, G;
    vector_ptrtype rhs;

    op_mat_ptrtype massOp, diffOp, convOp;

    std::map< std::string, std::set<flag_type> > M_bcFlags;

    op_ptrtype precOp;

    double M_nu, M_alpha;

    std::set<flag_type>::const_iterator inflowIter;

    std::string M_prob_type;

    void assembleMass();

    void assembleDiffusion();

    void assembleConvection();

    void applyBC( sparse_matrix_ptrtype& A );


};




template < typename space_type>
OperatorPCD<space_type>::OperatorPCD( pressure_space_ptrtype Qh,
                                                      std::map< std::string, std::set<flag_type> > bcFlags,
                                                      double nu,
                                                      double alpha )
    :
    super( Qh->map(), "PCD", false, false ),
    M_Qh( Qh ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_mass( backend()->newMatrix(Qh, Qh) ),
    M_diff( backend()->newMatrix(Qh, Qh) ),
    M_conv( backend()->newMatrix(Qh, Qh) ),
    G( backend()->newMatrix(Qh, Qh) ),
    rhs( backend()->newVector( Qh ) ),
    M_bcFlags( bcFlags ),
    M_nu( nu ),
    M_alpha( alpha )
{
    initialize();

    LOG(INFO) << "[Pressure Correction Diffusion Operator] Constructor: using nu=" << M_nu << "\n";
    LOG(INFO) << "[Pressure Correction Diffusion Operator] Constructor: using alpha=" << M_alpha << "\n";

    if ( alpha == 0 )
        this->setProblemType( "steady" );
    else
        this->setProblemType( "unsteady" );

    this->assembleMass();
    this->assembleDiffusion();
    this->assembleConvection();
}



template < typename space_type>
OperatorPCD<space_type>::OperatorPCD( const OperatorPCD& tc )
    :
    super(tc),
    M_Qh( tc.M_Qh ),
    p( tc.p ),
    q( tc.q ),
    M_mass( tc.M_mass ),
    M_diff( tc.M_diff ),
    M_conv( tc.M_conv ),
    G( tc.G ),
    rhs( tc.rhs ),
    massOp( tc.massOp ),
    diffOp( tc.diffOp ),
    convOp( tc.convOp ),
    precOp( tc.precOp ),
    M_bcFlags( tc.M_bcFlags ),
    M_nu( tc.M_nu ),
    M_alpha( tc.M_alpha ),
    inflowIter( tc.inflowIter ),
    M_prob_type( tc.M_prob_type )
{
    //LOG(INFO) << "Call for OperatorPCD copy constructor...\n";
}


template < typename space_type>
void
OperatorPCD<space_type>::initialize()
{
    rhs->zero();
    rhs->close();
}

template < typename space_type>
template < typename ExprConvection, typename ExprBC >
void
OperatorPCD<space_type>::update( ExprConvection const& expr_b,
                                 ExprBC const& ebc )
{

    auto conv  = form2( _test=M_Qh, _trial=M_Qh, _matrix=G );
    G->zero();

    conv = integrate( elements(M_Qh->mesh()), (trans(expr_b)*trans(gradt(p)))*id(q));

    for( auto dir : M_bcFlags["Dirichlet"])
    {
        std::string m = M_Qh->mesh()->markerName(dir);
        if (( soption("btpcd.pcd.inflow") == "Robin" ) &&
            ( ( m == "inlet" )  || ( m == "inflow") ) )
        {
            conv += integrate( _range=markedfaces(M_Qh->mesh(), dir), _expr=trans(ebc.find(M_Qh->mesh()->markerName(dir))->second)*N()*idt(p)*id(q));
        }
    }

    G->close();
    LOG(INFO) << "[OperatorPCD] Add diffusion matrix...\n";
    G->addMatrix( M_nu, M_diff );



    if ( this->problemType() == "unsteady" )
    {
        LOG(INFO) << "[OperatorPCD] Add mass matrix...\n";
        G->addMatrix( M_alpha, M_mass );
    }

    this->applyBC(G);

    static bool init_G = false;

    if ( !init_G )
    {
        // S = F G^-1 M
        LOG(INFO) << "[OperatorPCD] setting pcd operator...\n";
        if ( ioption("btpcd.pcd.order") == 1 )
            precOp = compose( massOp, compose(inv(op(G,"Ap")),diffOp) );
        else
            precOp = compose( diffOp, compose(inv(op(G,"Ap")),massOp) );
        LOG(INFO) << "[OperatorPCD] setting pcd operator done.\n";
        init_G = true;
    }
}





template < typename space_type>
void
OperatorPCD<space_type>::assembleMass()
{
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(q) );
    M_mass->close();
    massOp = op( M_mass, "Mp" );
}

template < typename space_type>
void
OperatorPCD<space_type>::assembleDiffusion()
{
    auto d = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_diff );
    d = integrate( elements(M_Qh->mesh()), trace(trans(gradt(p))*grad(q)));

    this->applyBC(M_diff);

    diffOp = op( M_diff, "Fp" );
}

template < typename space_type>
void
OperatorPCD<space_type>::assembleConvection()
{
    /*
     form2( M_Qh, M_Qh, M_conv, _init=true ) =
     integrate( elements(M_Qh->mesh()), typename MyIm<2*pOrder>::type(),
     trace(trans(gradt(p))*grad(q))
     );

     M_conv->close();

     this->applyBC(M_conv);
     */

}

template < typename space_type>
void
OperatorPCD<space_type>::applyBC( sparse_matrix_ptrtype& A )
{
    auto a = form2( _test=M_Qh, _trial=M_Qh, _matrix=A );
    for( auto dir : M_bcFlags["Dirichlet"])
    {
        std::string m = M_Qh->mesh()->markerName(dir);
        if (( soption("btpcd.pcd.inflow") == "Robin" ) &&
            ( ( m == "inlet" )  || ( m == "inflow") ) )
            continue;

        a += on( markedfaces(M_Qh->mesh(),dir), _element=p, _rhs=rhs, _expr=cst(0.) );
    }
    // on neumann boundary on velocity, apply Dirichlet condition on pressure
    if ( soption("btpcd.pcd.outflow") == "Dirichlet" )
        for( auto dir : M_bcFlags["Neumann"])
        {
            std::string m = M_Qh->mesh()->markerName(dir);
            if ( (m=="outlet") || (m == "outflow") )
                a += on( markedfaces(M_Qh->mesh(),dir), _element=p, _rhs=rhs, _expr=cst(0.) );
        }
    rhs->close();
    A->close();
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
