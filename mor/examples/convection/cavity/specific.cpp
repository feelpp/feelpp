/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4 */

#include "../convection.hpp"


std::string ConvectionCrb::name()
{
    std::stringstream app_name;
    app_name<< "cavity"
            << CONVECTION_DIM
            << "d";
    return app_name.str();
}

/**
 * \return the number of terms of the affine decompostion of the left hand side
 * bilinear form
 */
int ConvectionCrb::Qa() const { return 3;}

/**
 * \return the number of terms of the affine decomposition of the left hand side
 * trilinear form
 */
int ConvectionCrb::QaTri() const { return 1; }

/**
 * \param l the index of output
 * \return number of terms  in affine decomposition of the \p q th output term
 */
int ConvectionCrb::Ql( int l ) const
{
    return 1;
}

int ConvectionCrb::mMaxA( int q )
{
    if ( q < 3 )
        return 1;
    else
        throw std::logic_error( "[Model] ERROR : try to acces to mMaxA(q) with a bad value of q");
}

int ConvectionCrb::mMaxF( int output_index, int q)
{
    int dumy=0;
    if( output_index < 2 )
    {
        if ( q < Ql(output_index) )
            return 1;
        else
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxF(output_index,q) with a bad value of q");
    }
    if( output_index >= 2 )
            throw std::logic_error( "[Model] ERROR : try to acces to mMaxF(output_index,q) with a bad value of output_index");
    return dumy;
}

/**
 * \return the number of outputs associated to the model,
 * output number 0 is the RHS of the primal problem
 */
int ConvectionCrb::Nl() const
{
    return 2;
}


ConvectionCrb::betaqm_type ConvectionCrb::computeBetaQm( parameter_type const& mu )
{
    M_betaAqm.resize( Qa() );
    M_betaAqm[0].resize(1);
    M_betaAqm[1].resize(1);
    M_betaAqm[2].resize(1);
    M_betaAqm[0][0] = 1/math::sqrt( mu( 0 ) );
    M_betaAqm[1][0] = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) );
    M_betaAqm[2][0] = 1;

    M_betaFqm.resize( Nl() );

    M_betaFqm[0].resize( Ql( 0 ) );
    M_betaFqm[0][0].resize(1);
    M_betaFqm[0][0][0] = 1/( mu( 1 )*math::sqrt( mu ( 0 ) ) );

    M_betaFqm[1].resize( Ql( 1 ) );
    M_betaFqm[1][0].resize(1);
    M_betaFqm[1][0][0] = 1;

    return boost::make_tuple( M_betaAqm, M_betaFqm );
}

/**
 * \brief solve the model for parameter \p mu
 * \param mu the model parameter
 * \param T solution
 */
void
ConvectionCrb::solve( parameter_type const& mu, element_ptrtype& T )
{
    T->zero();
    using namespace vf;

    M_backend->nlSolver()->jacobian =
        std::bind( &self_type::updateJ, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );
    M_backend->nlSolver()->residual =
        std::bind( &self_type::updateR, std::ref( *this ), std::placeholders::_1, std::placeholders::_2 );

    vector_ptrtype R( M_backend->newVector( Xh ) );
    sparse_matrix_ptrtype J( M_backend->newMatrix( Xh,Xh ) );

    double current_gr=mu(0);
    double current_pr=mu(1);
    bool use_continuity = boption( "use_continuity");
    int N=1;

    if( use_continuity )
        N = 0.5*std::ceil( std::max( 2.0, std::max( std::log(mu(0)), std::log(mu(1)) ) ) );

    A_OUT << " Resolution for gr=" << mu(0) << ", pr="<< mu(1) << std::endl;
    A_OUT << " N = " << N << std::endl;
    for ( int i = 1; i <= N; ++i )
    {
        if( use_continuity )
        {
            current_gr = math::exp( i*math::log(mu(0))/N );
            current_pr = math::exp( i*math::log(mu(1))/N );
            A_OUT << " Continuity iteration : "<< i <<"/"<< N << std::endl;
        }

        auto current_mu = Dmu->element();
        current_mu << current_gr, current_pr;
        this->computeBetaQm( current_mu );
        this->update( current_mu );

        M_backend->nlSolve (_solution=T, _jacobian=J, _residual=R );
    }
}


/**
 * \param output_index number of output to evaluate
 * \param mu parameter
 */
double
ConvectionCrb::output( int output_index,
                       parameter_type const& mu,
                       element_type& unknown, bool need_to_solve)
{
    using namespace vf;

    auto mesh = Xh->mesh();
    element_type U = Xh->element( "U" );
    U = *pT;

    auto u = U.template element<0>(); // velocity
    auto p = U.template element<1>(); // pressure
    auto xi = U.template element<2>(); // lagrange multiplier
    auto t = U.template element<3>(); // temperature

    A_OUT << "====================================================\n"
          << "                     MEASURES\n";
    // domain measure
    double meas = integrate( elements(mesh), cst(1.0) ).evaluate()(0,0);
    A_OUT << "measure of the domain = " << meas << std::endl;

    // mean value of div(u)
    double mean_div_u = integrate( elements(mesh), divv(u) ).evaluate()(0,0);
    A_OUT << "mean div(u) = " << mean_div_u << std::endl;

    // L2 norme of div(u) (should be 0)
    double div_u_error_L2 = integrate( elements(mesh), divv(u)*divv(u) ).evaluate()(0,0);
    A_OUT << "||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    double mean_T = integrate( elements(mesh), idv(t) ).evaluate()(0,0) ;
    A_OUT << "mean temperature = " << mean_T << std::endl;


    // Evaluation of the specific output
    double output = 0.0;
    switch ( output_index )
    {
    // RHS (compliant)
    case 0 :
    {
        break;
    }
    // output 1 : Nusselt number
    case 1 :
    {
        output = integrate( markedfaces( mesh,"wallRight" ) ,
                            idv( t ) ).evaluate()( 0,0 ) ;
        A_OUT << "output 1 : Nu = " << output << std::endl;
        break;
    }
    default :
    {
        A_OUT << "ERROR : you tried to evaluate an output out of range\n";
    }
    }
    A_OUT << "====================================================\n";
    return output;
}

namespace Feel
{
    FEELPP_CRBTRILINEAR_PLUGIN( ConvectionCrb, cavity2d )
}
