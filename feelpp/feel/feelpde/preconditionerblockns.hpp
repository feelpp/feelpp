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
#ifndef FEELPP_PRECONDITIONERBLOCKNS_HPP
#define FEELPP_PRECONDITIONERBLOCKNS_HPP 1

#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelpde/operatorpcd.hpp>
#include <feel/feelpde/boundaryconditions.hpp>
#include <feel/feeldiscr/pdh.hpp>

namespace Feel
{
template< typename SpaceType, typename PropertiesSpaceType = Pdh_type<typename SpaceType::mesh_type,0> >
class PreconditionerBlockNS : public Preconditioner<typename SpaceType::value_type>
{
    typedef Preconditioner<typename SpaceType::value_type> super;
public:

    enum Type
    {
        PCD = 0, // pressure convection diffusion
        PCD_ACCELERATION = 1, // pressure convection diffusion
        PMM = 2, // pressure mass matrix
        SIMPLE=3 //
    };
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    using space_type = SpaceType;
    typedef std::shared_ptr<space_type> space_ptrtype;

    using properties_space_type = PropertiesSpaceType;
    using properties_space_ptrtype = std::shared_ptr<properties_space_type>;
    using property_type = typename properties_space_type::element_type;
    typedef std::shared_ptr<property_type> property_ptrtype;

    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type velocity_space_type;
    typedef typename space_type::template sub_functionspace<1>::type pressure_space_type;
    typedef typename space_type::template sub_functionspace<0>::ptrtype velocity_space_ptrtype;
    typedef typename space_type::template sub_functionspace<1>::ptrtype pressure_space_ptrtype;

    typedef typename velocity_space_type::element_type velocity_element_type;
    typedef typename pressure_space_type::element_type pressure_element_type;

    typedef typename space_type::value_type value_type;

    static inline const uint16_type Dim = space_type::nDim;
    static inline const uint16_type uOrder = velocity_space_type::basis_type::nOrder;
    static inline const uint16_type pOrder = pressure_space_type::basis_type::nOrder;

    typedef OperatorMatrix<value_type> op_mat_type;
    typedef std::shared_ptr<op_mat_type> op_mat_ptrtype;

    typedef typename OperatorPCD<space_type>::type op_pcd_type;
    typedef typename OperatorPCD<space_type>::ptrtype op_pcd_ptrtype;
    typedef OperatorBase<value_type> op_type;
    typedef std::shared_ptr<op_type> op_ptrtype;

    /**
     * \param Xh velocity/pressure space type
     * \param bcFlags the boundary conditions flags
     * \param A the full matrix \f$(F B^T;B 0)\f$
     * \param mu viscosity
     * \param rho density
     * \param alpha mass term
     */
    PreconditionerBlockNS( std::string t,
                           space_ptrtype Xh,
                           BoundaryConditions bcFlags,
                           std::string const& s,
                           sparse_matrix_ptrtype A,
                           double mu,
                           double rho,
                           double alpha = 0 );

    template<typename MuExprT, typename RhoExprT, typename AlphaExprT>
    PreconditionerBlockNS( std::string t,
                           space_ptrtype Xh,
                           properties_space_ptrtype Ph,
                           BoundaryConditions bcFlags,
                           std::string const& s,
                           sparse_matrix_ptrtype A,
                           MuExprT mu,
                           RhoExprT rho,
                           AlphaExprT alpha );

    Type type() const { return M_type; }
    std::string typeStr() const
        {
            switch ( M_type )
            {
            default:
            case  PCD:
                return "PCD" ;
            case PCD_ACCELERATION:
                return "PCD_ACCELERATION";
            case PMM:
                return "PMM";
            case SIMPLE:
                return "SIMPLE";
            }
        }
    void setType( std::string t );

    BoundaryConditions const& bcFlags() const { return M_bcFlags; }
    void setParameterValues( std::map<std::string,double> const& pv ) { M_bcExprParameterValues=pv; }

    void initialize();

    void assembleHelmholtz( double mu, double rho, double alpha = 0 );

    void assembleDivergence();

    void assembleGradient();

    void assembleSchurApp( double mu, double rho, double alpha = 0 );

#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( sparse_matrix_ptrtype ),
                                     update,
                                     tag,
                                     ( required
                                       ( matrix,*( boost::is_convertible<mpl::_,std::shared_ptr<MatrixSparse> > ) )
                                       ( bc, *) )
                                     ( optional
                                       ( diffusion,*( boost::is_convertible<mpl::_,ExprBase> ), cst(M_mu) )
                                       ( density,*( boost::is_convertible<mpl::_,ExprBase> ), cst(M_rho) )
                                       ( convection,*( boost::is_convertible<mpl::_,ExprBase> ), zero<Dim>() )
                                       ( tn,  *( std::is_arithmetic<mpl::_> ), 0. )
                                       ( tn1,  *( std::is_arithmetic<mpl::_> ), 0. )
                                     ) )
    {
        this->setDiffusion( diffusion );
        this->setDensity( diffusion );
        this->setConvection( convection);
        this->setBC( bc );
        this->setTimes( std::pair<double,double>( tn1, tn) )
        this->update( matrix );
    }
#endif

#if 0
    template<typename E>
    void setDiffusion( E && e )
        {
            M_muh.on( _range=elements(mesh), _expr=e );
        }
    template<typename E>
    void setDensity( E && e )
        {
            M_rhoh.on( _range=elements(mesh), _expr=e );
        }
    template<typename E>
    void setConvection( E && e )
        {
            M_convh.on( _range=elements(mesh), _expr=e );
        }
#endif
    template< typename Expr_convection, typename Expr_diffusion, typename Expr_bc >
    void updateImpl( sparse_matrix_ptrtype A,
                     Expr_diffusion const& expr_d,
                     Expr_convection const& expr_b,
                     Expr_bc const& g,
                     double tn=0, double tn1 = 0 )
     {
         CHECK(0) << "Invalid call";
     }

    template< typename Expr_convection, typename Expr_bc, typename MuExprT, typename RhoExprT, typename AlphaExprT >
    void update( sparse_matrix_ptrtype A, Expr_convection const& expr_b,
                 Expr_bc const& g, MuExprT mu, RhoExprT rho, AlphaExprT alpha,
                 bool hasConvection=true,
                 double tn=0, double tn1 = 0 );

    template< typename Expr_convection, typename Expr_bc >
    void update( sparse_matrix_ptrtype A, Expr_convection const& expr_b,
                 Expr_bc const& g, bool hasConvection=true,
                 double tn=0, double tn1 = 0 );

    template< typename Expr_convection >
    void update( sparse_matrix_ptrtype A, Expr_convection const& expr_b,
                 bool hasConvection=true, double tn=0, double tn1 = 0  );

    void update( sparse_matrix_ptrtype A );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        applyInverse( X, Y );
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;
    int guess( vector_type& U ) const;

    // other function
    int SetUseTranspose( bool UseTranspose)
    {
        return(false);
    }

    double NormInf() const
    {
        return 0;
    }

    const char * Label () const
    {
        return("Triangular Blockns Preconditioner");
    }

    bool UseTranspose() const
    {
        return(false);
    }

    bool HasNormInf () const
    {
        return(false);
    }

    virtual ~PreconditionerBlockNS(){};

    template<typename ExprType>
    void setMu( ExprType m ) { M_mu->on( _range=elements( M_mu->mesh() ), _expr=m ); M_updatePMM = true; }
    template<typename ExprType>
    void setRho( ExprType r ) { M_rho->on( _range=elements( M_rho->mesh() ), _expr=r ); M_updatePMM = true; }
    template<typename ExprType>
    void setAlpha( ExprType a ) { M_alpha->on( _range=elements( M_alpha->mesh() ), _expr=a ); }
private:
    void createSubMatrices();
    void updatePMM();
private:

    Type M_type;

    backend_ptrtype M_b;

    space_ptrtype M_Xh;

    velocity_space_ptrtype M_Vh;
    pressure_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;

    properties_space_ptrtype M_Ph;

    property_ptrtype M_mu, M_alpha, M_rho;

    sparse_matrix_ptrtype /*M_helm, G, M_div*/ M_F, /*M_B,*/ M_Bt, M_mass/*, M_massv_inv*/;
    bool M_updatePMM;
    op_mat_ptrtype divOp, helmOp;

    mutable vector_ptrtype M_rhs, M_aux, M_vin,M_pin, M_vout, M_pout;

    mutable element_type U;
    velocity_element_type u, v;
    pressure_element_type p, q;

    op_pcd_ptrtype pcdOp;
    op_ptrtype pm;

    BoundaryConditions M_bcFlags;
    std::map<std::string,double> M_bcExprParameterValues;
    std::string M_prefix;
    op_mat_ptrtype precHelm;
};


template < typename SpaceType, typename PropertiesSpaceType >
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::
PreconditionerBlockNS( std::string t,
                       space_ptrtype Xh,
                       BoundaryConditions bcFlags,
                       std::string const& s,
                       sparse_matrix_ptrtype A,
                       double mu, double rho, double alpha )
    :
    M_type( PCD ),
    M_b(backend()),
    M_Xh( Xh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_Vh_indices( A->mapRow().dofIdToContainerId( 0 ) ),
    M_Qh_indices( A->mapRow().dofIdToContainerId( 1 ) ),
    M_Ph( properties_space_type::New(Xh->mesh()) ),
    M_mass( M_b->newMatrix( _trial=M_Qh, _test=M_Qh) ),
    M_updatePMM( false ),
    M_rhs( M_b->newVector( M_Vh )  ),
    M_aux( M_b->newVector( M_Vh )  ),
    M_vin( M_b->newVector( M_Vh )  ),
    M_pin( M_b->newVector( M_Qh )  ),
    M_vout( M_b->newVector( M_Vh )  ),
    M_pout( M_b->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_bcFlags( bcFlags ),
    M_prefix( s )
{
    M_mu = std::make_shared<property_type>( M_Ph->element(cst(mu), "mu") );
    M_alpha = std::make_shared<property_type>( M_Ph->element(cst(alpha), "alpha") );
    M_rho = std::make_shared<property_type>( M_Ph->element(cst(rho), "rho") );

    tic();
    LOG(INFO) << "[PreconditionerBlockNS] setup starts";
    this->setMatrix( A );

    this->createSubMatrices();


    initialize();

    this->setType ( t );
    toc( "[PreconditionerBlockNS] setup done ", FLAGS_v > 0 );
}

template < typename SpaceType, typename PropertiesSpaceType >
template<typename MuExprT, typename RhoExprT, typename AlphaExprT>
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::
PreconditionerBlockNS( std::string t,
                       space_ptrtype Xh,
                       properties_space_ptrtype Ph,
                       BoundaryConditions bcFlags,
                       std::string const& s,
                       sparse_matrix_ptrtype A,
                       MuExprT mu,
                       RhoExprT rho,
                       AlphaExprT alpha )
    :
    M_type( PCD ),
    M_b(backend()),
    M_Xh( Xh ),
    M_Vh( M_Xh->template functionSpace<0>() ),
    M_Qh( M_Xh->template functionSpace<1>() ),
    M_Vh_indices( A->mapRow().dofIdToContainerId( 0 ) ),
    M_Qh_indices( A->mapRow().dofIdToContainerId( 1 ) ),
    M_Ph( Ph ),
    M_mass( M_b->newMatrix( _trial=M_Qh, _test=M_Qh) ),
    M_updatePMM( false ),
    M_rhs( M_b->newVector( M_Vh )  ),
    M_aux( M_b->newVector( M_Vh )  ),
    M_vin( M_b->newVector( M_Vh )  ),
    M_pin( M_b->newVector( M_Qh )  ),
    M_vout( M_b->newVector( M_Vh )  ),
    M_pout( M_b->newVector( M_Qh )  ),
    U( M_Xh, "U" ),
    u( M_Vh, "u" ),
    v( M_Vh, "v" ),
    p( M_Qh, "p" ),
    q( M_Qh, "q" ),
    M_bcFlags( bcFlags ),
    M_prefix( s )
{
    M_mu = std::make_shared<property_type>( M_Ph->element(mu, "mu") );
    M_alpha = std::make_shared<property_type>( M_Ph->element(alpha, "alpha") );
    M_rho = std::make_shared<property_type>( M_Ph->element(rho, "rho") );

    tic();
    LOG(INFO) << "[PreconditionerBlockNS] setup starts";
    this->setMatrix( A );

    this->createSubMatrices();


    initialize();

    this->setType ( t );
    toc( "[PreconditionerBlockNS] setup done ", FLAGS_v > 0 );
}

template < typename SpaceType, typename PropertiesSpaceType >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::initialize()
{
    M_rhs->zero();
    M_rhs->close();
}

template < typename SpaceType, typename PropertiesSpaceType >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::createSubMatrices()
{
    tic();
    if ( !M_F )
    {
        tic();
        M_F = this->matrix()->createSubMatrix( M_Vh_indices, M_Vh_indices, true );
        M_F->mapRowPtr()->setIndexSplit( M_Vh->dof()->indexSplit() );
        if ( M_Vh->dof()->hasIndexSplitWithComponents() )
            M_F->mapRowPtr()->setIndexSplitWithComponents( M_Vh->dof()->indexSplitWithComponents() );
        //M_B = this->matrix()->createSubMatrix( M_Qh_indices, M_Vh_indices );
        tic();
        M_Bt = this->matrix()->createSubMatrix( M_Vh_indices, M_Qh_indices );
        toc("submatrix B^T",FLAGS_v>0);
        helmOp = op( M_F, "Fu" );
        helmOp->setCloseMatrixRhs( false );
        divOp = op( M_Bt, "Bt");
        toc("create submatrix", FLAGS_v>0);
    }
    else
    {
        tic();
        this->matrix()->updateSubMatrix( M_F, M_Vh_indices, M_Vh_indices );
        //this->matrix()->updateSubMatrix( M_B, M_Qh_indices, M_Vh_indices );
        tic();
        this->matrix()->updateSubMatrix( M_Bt, M_Vh_indices, M_Qh_indices );
        toc("update submatrix B^T",FLAGS_v>0);
        toc("update submatrix",FLAGS_v>0);
    }
    toc( "PreconditionerBlockNS::createSubMatrix(Fu,B^T)", FLAGS_v > 0 );
}

template < typename SpaceType, typename PropertiesSpaceType >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::updatePMM()
{
    LOG(INFO) << "Updating PMM preconditioner...";
    auto m = form2( _test=M_Qh, _trial=M_Qh, _matrix=M_mass );
    m = integrate( elements(M_Qh->mesh()), idt(p)*id(q)/idv(*M_mu) );
    M_mass->close();
    LOG(INFO) << "Updating PMM preconditioner done.";
}

template < typename SpaceType, typename PropertiesSpaceType >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::setType( std::string t )
{
    if ( t == "PCD") M_type = PCD;
    if ( t == "PCD_ACCELERATION") M_type = PCD_ACCELERATION;
    if ( t == "PMM") M_type = PMM;
    if ( t == "SIMPLE") M_type = SIMPLE;

    LOG(INFO) << "setting preconditioner " << t << " type: " << M_type;
    switch( M_type )
    {
    case PCD:
    case PCD_ACCELERATION:
        tic();
        pcdOp = std::make_shared<op_pcd_type>( M_Xh, M_b, M_bcFlags, M_prefix, M_type==PCD_ACCELERATION );
        pcdOp->setBt( M_Bt );
        this->setSide( super::RIGHT );

        toc( "Preconditioner::setType " + typeStr(), FLAGS_v > 0 );

        break;
    case PMM:
    {
        tic();
        M_updatePMM = true;

        if ( boption( "blockns.pmm.diag" ) )
        {
            pm = diag( op( M_mass, "Mp" ) );
        }
        else
        {
            pm = op( M_mass, "Mp" );
        }
        this->setSide( super::RIGHT );
        toc( "Preconditioner::setType PMM", FLAGS_v > 0 );
    }
    break;
    case SIMPLE:
        break;
    }
}

template < typename SpaceType, typename PropertiesSpaceType >
template< typename Expr_convection, typename Expr_bc, typename MuExprT, typename RhoExprT, typename AlphaExprT >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::
update( sparse_matrix_ptrtype A, Expr_convection const& expr_b,
        Expr_bc const& g, 
        MuExprT mu, 
        RhoExprT rho, 
        AlphaExprT alpha,
        bool hasConvection,
        double tn, double tn1 )
{
    tic();
    this->setMatrix( A );
    this->createSubMatrices();

    this->setMu(mu);
    this->setRho(rho);
    this->setAlpha(alpha);

    if ( type() == PCD || type() == PCD_ACCELERATION )
    {
        tic();
        bool hasAlpha = ( M_alpha->linftyNorm() > 1e-15 );
        pcdOp->update( idv(M_rho), idv(M_mu), expr_b, g, idv(M_alpha), hasConvection, hasAlpha, tn, tn1 );
        toc( "Preconditioner::update "+ typeStr(), FLAGS_v > 0 );
    }
    if ( type() == PMM && M_updatePMM )
    {
        tic();
        updatePMM();
        toc( "Preconditioner::update "+ typeStr(), FLAGS_v > 0 );
    }
    toc( "Preconditioner::update", FLAGS_v > 0 );
}

template < typename SpaceType, typename PropertiesSpaceType >
template< typename Expr_convection, typename Expr_bc >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::
update( sparse_matrix_ptrtype A,
        Expr_convection const& expr_b,
        Expr_bc const& g,
        bool hasConvection,
        double tn, double tn1 )
{
    tic();
    this->setMatrix( A );
    this->createSubMatrices();
    if ( type() == PCD || type() == PCD_ACCELERATION )
    {
        tic();
        bool hasAlpha = ( M_alpha->linftyNorm() > 1e-15 );
        pcdOp->update( idv(M_rho), idv(M_mu), expr_b, g, idv(M_alpha), hasConvection, hasAlpha, tn, tn1 );
        toc( "Preconditioner::update "+ typeStr(), FLAGS_v > 0 );
    }
    if ( type() == PMM && M_updatePMM )
    {
        tic();
        updatePMM();
        toc( "Preconditioner::update "+ typeStr(), FLAGS_v > 0 );
    }

    toc( "Preconditioner::update", FLAGS_v > 0 );
}
template < typename SpaceType, typename PropertiesSpaceType >
template< typename Expr_convection >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::
update( sparse_matrix_ptrtype A,
        Expr_convection const& expr_b,
        bool hasConvection,
        double tn, double tn1 )
{
    map_vector_field<Dim,1,2> m_dirichlet { M_bcFlags.template getVectorFields<Dim> ( std::string(M_prefix), "Dirichlet" ) };
    if ( !M_bcExprParameterValues.empty() )
        m_dirichlet.setParameterValues( M_bcExprParameterValues );
    this->update( A, expr_b, m_dirichlet, hasConvection, tn, tn1 );
}
template < typename SpaceType, typename PropertiesSpaceType >
void
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::update( sparse_matrix_ptrtype A )
{
    this->update( A, zero<Dim,1>(), false, 0., 0. );
}


template < typename SpaceType, typename PropertiesSpaceType >
int
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
#if 0
    tic();
    U = X;
    U.close();
    LOG(INFO) << "Create velocity/pressure component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();
    *M_aux = *M_vin;
    M_aux->close();

    if ( this->type() == PMM )
    {
        LOG(INFO) << "Applying PMM:  pressure mass matrix";
        tic();
        pm->applyInverse( *M_pin, *M_pout );
        M_pout->scale(-1);
        M_pout->close();
        toc("PreconditionerBlockNS::applyInverse PMM::Q^-1",FLAGS_v>0);
        LOG(INFO) << "Applying PMM done";
    }
    if ( this->type() == PCD || this->type() == PCD_ACCELERATION )
    {
        if ( boption("blockns.pcd") )
        {
            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion...\n";
            CHECK(pcdOp) << "Invalid PCD operator\n";
            CHECK(M_aux) << "Invalid aux vector\n";
            CHECK(M_pout) << "Invalid aux vector\n";
            tic();
            pcdOp->applyInverse( *M_pin, *M_pout );
            M_pout->scale(-1);
            M_pout->close();
            toc("PreconditionerBlockNS::applyInverse PCD::S^-1",FLAGS_v>0);

            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion done\n";
        }
        else
        {
            *M_pout = *M_pin;
            M_pout->close();
        }
    }
    LOG(INFO) << "pressure/velocity blockns : apply divergence...\n";
    tic();
    divOp->apply( *M_pout, *M_vout );


    M_aux->add( -1.0, *M_vout );
    M_aux->close();
    toc("PreconditionerBlockNS::applyInverse apply B^T",FLAGS_v>0);

    if ( boption("blockns.cd") )
    {

        LOG(INFO) << "velocity blockns : apply inverse convection diffusion...\n";
        tic();
        helmOp->applyInverse(*M_aux, *M_vout);
        toc("PreconditionerBlockNS::applyInverse Fu^-1",FLAGS_v>0);
    }
    else
    {
        *M_vout = *M_vin;
        M_vout->close();
    }


    LOG(INFO) << "Update output velocity/pressure...\n";
    tic();
    U.template element<0>() = *M_vout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();
    toc("PreconditionerBlockNS::applyInverse update solution",FLAGS_v>0);
    toc("PreconditionerBlockNS::applyInverse", FLAGS_v > 0 );
    return 0;
#else
    tic();
    LOG(INFO) << "Create velocity/pressure component...\n";
    auto XUblas = M_Xh->element( X );
    auto vinUblas = XUblas.template element<0>();
    auto pinUblas = XUblas.template element<1>();
    auto YUblas = M_Xh->element( Y );
    auto voutUblas = YUblas.template element<0>();
    auto poutUblas = YUblas.template element<1>();

    auto vin = M_b->toBackendVectorPtr( vinUblas );
    auto pin = M_b->toBackendVectorPtr( pinUblas );
    auto vout = M_b->toBackendVectorPtr( voutUblas );
    auto pout = M_b->toBackendVectorPtr( poutUblas );

    *M_aux = *vin;

    if ( this->type() == PMM )
    {
        LOG(INFO) << "Applying PMM:  pressure mass matrix";
        tic();
        pm->applyInverse( *pin, *pout );
        pout->scale(-1);
        toc("PreconditionerBlockNS::applyInverse PMM::Q^-1",FLAGS_v>0);
        LOG(INFO) << "Applying PMM done";
    }
    if ( this->type() == PCD || this->type() == PCD_ACCELERATION )
    {
        if ( boption("blockns.pcd") )
        {
            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion...\n";
            CHECK(pcdOp) << "Invalid PCD operator\n";
            CHECK(M_aux) << "Invalid aux vector\n";
            CHECK(pout) << "Invalid aux vector\n";
            tic();
            pcdOp->applyInverse( *pin, *pout );
            pout->scale(-1);
            toc("PreconditionerBlockNS::applyInverse PCD::S^-1",FLAGS_v>0);

            LOG(INFO) << "pressure blockns: Solve for the pressure convection diffusion done\n";
        }
        else
        {
            *pout = *pin;
        }
    }
    LOG(INFO) << "pressure/velocity blockns : apply divergence...\n";
    tic();
    divOp->apply( *pout, *vout );


    M_aux->add( -1.0, *vout );
    toc("PreconditionerBlockNS::applyInverse apply B^T",FLAGS_v>0);

    if ( boption("blockns.cd") )
    {

        LOG(INFO) << "velocity blockns : apply inverse convection diffusion...\n";
        tic();
        helmOp->applyInverse(*M_aux, *vout);
        toc("PreconditionerBlockNS::applyInverse Fu^-1",FLAGS_v>0);
    }
    else
    {
        *vout = *vin;
    }

    toc("PreconditionerBlockNS::applyInverse", FLAGS_v > 0 );
    return 0;
#endif
}

template < typename SpaceType, typename PropertiesSpaceType >
int
PreconditionerBlockNS<SpaceType,PropertiesSpaceType>::guess ( vector_type& Y ) const
{
    U = Y;
    U.close();

    LOG(INFO) << "Create velocity/pressure component...\n";
    *M_vin = U.template element<0>();
    M_vin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    LOG(INFO) << "pressure/velocity blockns : apply divergence...\n";
    divOp->apply( *M_pout, *M_vin );

    M_aux->zero();
    M_aux->add( -1.0, *M_vin );
    M_aux->close();

    LOG(INFO) << "velocity blockns : apply inverse convection diffusion...\n";
    helmOp->applyInverse(*M_aux, *M_vin);
    LOG(INFO) << "Update output velocity/pressure...\n";

    U.template element<0>() = *M_vin;
    U.template element<1>() = *M_pin;
    U.close();
    Y=U;
    Y.close();

    return 0;
}
namespace meta
{
template< typename SpaceType, typename PropertiesSpaceType >
struct blockns
{
    typedef PreconditionerBlockNS<SpaceType,PropertiesSpaceType> type;
    typedef std::shared_ptr<type> ptrtype;
};
}
BOOST_PARAMETER_MEMBER_FUNCTION( ( typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type,
                                                          typename parameter::value_type<Args, tag::properties_space>::type::element_type>::ptrtype ),
                                 blockns,
                                 tag,
                                 ( required
                                   ( space, *)
                                   ( matrix, *)
                                   ( type, *)
                                   )
                                 ( optional
                                   ( properties_space, *, (Pdh<0>( space->mesh() )))
                                   ( prefix, *( boost::is_convertible<mpl::_,std::string> ), "" )
                                   ( bc, *, (BoundaryConditions ()) )
                                   ( mu,  *, cst(doption("mu")) )
                                   ( rho,  *, cst(doption("rho")) )
                                   ( alpha, *, cst(0.) )
                                   )
                                 )
{
    typedef typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type,
                                   typename parameter::value_type<Args, tag::properties_space>::type::element_type>::ptrtype pblockns_t;
    typedef typename meta::blockns<typename parameter::value_type<Args, tag::space>::type::element_type,
                                   typename parameter::value_type<Args, tag::properties_space>::type::element_type>::type blockns_t;
    pblockns_t p( new blockns_t( type, space, properties_space, bc, prefix, matrix, mu, rho, alpha ) );
    return p;
} // btcpd
} // Feel
#endif
