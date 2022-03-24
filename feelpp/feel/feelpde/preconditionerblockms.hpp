/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
   -*- vim: set ft=cpp fenc=utf-8 sw=4 ts=4 sts=4 tw=80 et cin cino=N-s,c0,(0,W4,g0:

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
#ifndef FEELPP_PRECONDITIONERBlockMS_HPP
#define FEELPP_PRECONDITIONERBlockMS_HPP 1

#if FEELPP_DIM == 2
#define curl_op curlx
#define curlt_op curlxt
#define curlv_op curlxv
#else
#define curl_op curl
#define curlt_op curlt
#define curlv_op curlv
#endif

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/operator.hpp>
#include <feel/feelalg/preconditioner.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feelalg/backendpetsc.hpp>

#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelpde/boundaryconditions.hpp>

namespace Feel
{
template< typename space_type >
class PreconditionerBlockMS : public Preconditioner<typename space_type::value_type, typename space_type::size_type>
{
    typedef Preconditioner<typename space_type::value_type,typename space_type::size_type> super;
public:
    static const uint16_type Dim = space_type::nDim;
    typedef typename space_type::value_type value_type;
    typedef typename space_type::size_type size_type;

    typedef typename backend_type::solve_return_type solve_return_type;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;

    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::indexsplit_ptrtype  indexsplit_ptrtype;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::mesh_ptrtype mesh_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef typename space_type::template sub_functionspace<0>::type potential_space_type;
    typedef typename space_type::template sub_functionspace<1>::type lagrange_space_type;
    typedef typename std::shared_ptr<potential_space_type> potential_space_ptrtype;
    typedef typename std::shared_ptr<lagrange_space_type> lagrange_space_ptrtype;

    // Qh3
    typedef Lagrange<1,Vectorial> lag_v_type;
    typedef FunctionSpace<mesh_type, bases< lag_v_type >> lag_v_space_type;
    typedef std::shared_ptr<lag_v_space_type> lag_v_space_ptrtype;

    // Potential
    typedef typename potential_space_type::element_type potential_element_type;
    typedef typename lagrange_space_type::element_type lagrange_element_type;

    // Grad operator
    typedef Grad_t<lagrange_space_type, potential_space_type> grad_type;

    /// Forms type
    typedef typename Feel::meta::LinearForm<potential_space_type>::type form1_potential_type;
    typedef typename Feel::meta::BilinearForm<potential_space_type, potential_space_type>::type form2_potential_type;

    typedef typename Feel::meta::LinearForm<lagrange_space_type>::type form1_lagrange_type;
    typedef typename Feel::meta::BilinearForm<lagrange_space_type, lagrange_space_type>::type form2_lagrange_type;

    typedef typename Feel::meta::LinearForm<lag_v_space_type>::type form1_lag_v_type;
    typedef typename Feel::meta::BilinearForm<lag_v_space_type, lag_v_space_type>::type form2_lag_v_type;

    /**
     * \param t Kind of prec (Simple or AFP)
     * \param Xh potential/lagrange space type
     * \param Mh Permeability space type
     * \param bcFlags the boundary conditions flags
     * \param s name of backend
     * \param A the full matrix
     */
    PreconditionerBlockMS( space_ptrtype Xh,
                           ModelProperties model,
                           std::string const& s,
                           sparse_matrix_ptrtype A,
                           value_type relax=1);

    void init( void );
    void initAMS( void );

    void apply( const vector_type & X, vector_type & Y ) const
    {
        this->applyInverse(X,Y);
    }

    int applyInverse ( const vector_type& X, vector_type& Y ) const;

    void attacheAAlpha(sparse_matrix_ptrtype _a) {M_a_alpha = _a;}
    void attacheABeta (sparse_matrix_ptrtype _b) {M_a_beta  = _b;}
    virtual ~PreconditionerBlockMS(){};
#if 0
    void printMatSize(int i, std::ostream &os)
    {
        if(i == 1)
            os << M_11->graph()->size();
        else if(i==2)
            os << M_L->graph()->size();
    }
    void printIter(int i, std::ostream & os){
        os << "<ul>";
        if(i == 1)
            for(auto i : NbIter1)
                os << "<li>" << i.nIterations() << " -- " << i.residual() << " -- " << i.isConverged() << "</li>";
        else if(i == 2)
            for(auto i : NbIter2)
                os << "<li>" << i.nIterations() << " -- " << i.residual() << " -- " << i.isConverged() << "</li>";
        else
            os << "<li>**error**</li>" << std::endl;
        os << "</ul>";
    }
    void iterMinMaxMean(void)
    {
        std::vector<double> tmp;
        auto minmax_1 = std::minmax_element(NbIter1.begin(), NbIter1.end(),[] (solve_return_type const& lhs, solve_return_type const& rhs) {return lhs.nIterations() < rhs.nIterations();});
        auto minmax_2 = std::minmax_element(NbIter2.begin(), NbIter2.end(),[] (solve_return_type const& lhs, solve_return_type const& rhs) {return lhs.nIterations() < rhs.nIterations();});
        double sum1 = std::accumulate(NbIter1.begin(), NbIter1.end(), 0.0, [] (double d, solve_return_type rhs) { return d+rhs.nIterations();});
        double mean1 = sum1 / NbIter1.size();
        double sum2 = std::accumulate(NbIter2.begin(), NbIter2.end(), 0.0, [] (double d, solve_return_type rhs) { return d+rhs.nIterations();});
        double mean2 = sum2 / NbIter2.size();

        tmp.push_back(minmax_1.first->nIterations());
        tmp.push_back(minmax_1.second->nIterations());
        tmp.push_back(mean1);
        M_minMaxMean.push_back(tmp);
        tmp.clear();
        tmp.push_back(minmax_2.first->nIterations());
        tmp.push_back(minmax_2.second->nIterations());
        tmp.push_back(mean2);
        M_minMaxMean.push_back(tmp);
    }
    double printMinMaxMean(int i, int j)
    {
        return M_minMaxMean.at(i).at(j);
    }
#endif
private:
    backend_ptrtype M_backend;

#if 0
    // for behavior analysis purpose
    mutable std::vector<solve_return_type> NbIter1, NbIter2;
    std::vector<std::vector<double>> M_minMaxMean;
#endif

    space_ptrtype M_Xh;

    potential_space_ptrtype M_Vh;
    lagrange_space_ptrtype M_Qh;
    std::vector<size_type> M_Vh_indices;
    std::vector<size_type> M_Qh_indices;

    // The two blocks: rhs and unknows
    mutable vector_ptrtype M_uin;
    mutable vector_ptrtype M_uout;
    mutable vector_ptrtype M_pin;
    mutable vector_ptrtype M_pout;

    mutable element_type U;

    sparse_matrix_ptrtype M_11;
    sparse_matrix_ptrtype M_mass;
    sparse_matrix_ptrtype M_L;

    /// Warning: at this point the permittivity is set to one for the domain
    value_type M_er; // permittivity

    ModelProperties M_model;

    std::string M_prefix;
    std::string M_prefix_11;
    std::string M_prefix_22;

    potential_element_type u;
    potential_element_type ozz;
    potential_element_type zoz;
    potential_element_type zzo;
    mutable vector_ptrtype M_ozz;
    mutable vector_ptrtype M_zoz;
    mutable vector_ptrtype M_zzo;

    lagrange_element_type X;
    lagrange_element_type Y;
    lagrange_element_type Z;
    mutable vector_ptrtype M_X;
    mutable vector_ptrtype M_Y;
    mutable vector_ptrtype M_Z;
    lagrange_element_type phi;

    grad_type M_grad;

    sparse_matrix_ptrtype M_a_alpha;
    sparse_matrix_ptrtype M_a_beta;

    value_type M_relax;
};

    template < typename space_type >
PreconditionerBlockMS<space_type>::PreconditionerBlockMS(space_ptrtype Xh,             // (u)x(p)
                                                         ModelProperties model,        // model
                                                         std::string const& p,         // prefix
                                                         sparse_matrix_ptrtype AA, value_type relax )    // The matrix
    :
        M_backend(backend()),           // the backend associated to the PC
        M_Xh( Xh ),
        M_Vh( Xh->template functionSpace<0>() ), // Potential
        M_Qh( Xh->template functionSpace<1>() ), // Lagrange
        M_Vh_indices( AA->mapRow().dofIdToContainerId( 0 ) ),
        M_Qh_indices( AA->mapRow().dofIdToContainerId( 1 ) ),
        M_uin( M_backend->newVector( M_Vh )  ),
        M_uout( M_backend->newVector( M_Vh )  ),
        M_pin( M_backend->newVector( M_Qh )  ),
        M_pout( M_backend->newVector( M_Qh )  ),
        U( M_Xh, "U" ),
        M_mass(M_backend->newMatrix(_test=M_Vh,_trial=M_Vh)),
        M_L(M_backend->newMatrix(_test=M_Qh,_trial=M_Qh)),
        M_er( 1. ),
        M_model( model ),
        M_prefix( p ),
        M_prefix_11( p+".11" ),
        M_prefix_22( p+".22" ),
        u(M_Vh, "u"),
        ozz(M_Vh, "ozz"),
        zoz(M_Vh, "zoz"),
        zzo(M_Vh, "zzo"),
        M_ozz(M_backend->newVector( M_Vh )),
        M_zoz(M_backend->newVector( M_Vh )),
        M_zzo(M_backend->newVector( M_Vh )),
        X(M_Qh, "X"),
        Y(M_Qh, "Y"),
        Z(M_Qh, "Z"),
        M_X(M_backend->newVector( M_Qh )),
        M_Y(M_backend->newVector( M_Qh )),
        M_Z(M_backend->newVector( M_Qh )),
        phi(M_Qh, "phi"),
        M_relax(relax)
{
    tic();
    LOG(INFO) << "[PreconditionerBlockMS] setup starts";
    //Feel::cout << "[PreconditionerBlockMS] relax (="<<doption(_prefix=M_prefix_11, _name="relax")<<") parameter is now provided thanks to the options" << std::endl;
    M_relax = doption(_prefix=M_prefix_11, _name="relax");
    this->setMatrix( AA );
    this->setName(M_prefix);

    M_11 = AA->createSubMatrix( M_Vh_indices, M_Vh_indices, true, true);

    /* Boundary conditions */
    CHECK( false ) << "TODO fix bc";
    BoundaryConditions M_bc;// = M_model.boundaryConditions();
    //map_vector_field<FEELPP_DIM,1,2> m_dirichlet_u { M_bc.getVectorFields<FEELPP_DIM> ( "u", "Dirichlet" ) };
    map_scalar_field<2> m_dirichlet_p { M_bc.getScalarFields<2> ( "phi", "Dirichlet" ) };

    //map_vector_field<FEELPP_DIM,1,2> m_weak_u { M_bc.getVectorFields<FEELPP_DIM> ( "u", "Weakdir" ) };
    map_scalar_field<2> m_weak_p { M_bc.getScalarFields<2> ( "phi", "Weakdir" ) };

    /* Compute the mass matrix (needed in first block, constant) */
    auto f2A = form2(_test=M_Vh, _trial=M_Vh, _matrix=M_mass);
    auto f1A = form1(_test=M_Vh);
    f2A = integrate(_range=elements(M_Vh->mesh()), _expr=inner(idt(u),id(u))); // M
    // the BC are applied during the init() function
    //for(auto const & it : m_dirichlet_u )
    //{
    //    LOG(INFO) << "Applying " << it.second << " on " << it.first << " for "<<M_prefix_11<<"\n";
    //    f2A += on(_range=markedfaces(M_Vh->mesh(),it.first), _expr=it.second,_rhs=f1A, _element=u, _type="elimination_symmetric");
    //}

    /* Compute the L (= er * grad grad) matrix (the second block) */
    auto f2L = form2(_test=M_Qh,_trial=M_Qh, _matrix=M_L);
#if 0
    //If you want to manage the relative permittivity materials per material,
    //here is the entry to deal with.
    for(auto it : M_model.materials() )
    {
        f2L += integrate(_range=markedelements(M_Qh->mesh(),marker(it)), _expr=M_er*inner(gradt(phi), grad(phi)));
    }
#else
    f2L += integrate(_range=elements(M_Qh->mesh()), _expr=M_er*inner(gradt(phi), grad(phi)));
#endif
    auto f1LQ = form1(_test=M_Qh);

    for(auto const & it : m_weak_p)
    {
        LOG(INFO) << "Applying (weak) " << it.second.first << " on " << it.first << " for "<<M_prefix_22<<"\n";
        f2L += integrate(_range=markedfaces(M_Qh->mesh(),it.second.second),
                         _expr=M_er*inner(trans(gradt(phi)),N())*id(phi)
                         + (doption(_prefix=M_prefix,_name="penaldir")/hFace())*idt(phi)*id(phi)
                         );
    }

    for(auto const & it : m_dirichlet_p)
    {
        LOG(INFO) << "Applying (on)" << it.second << " on " << it.first << " for "<<M_prefix_22<<"\n";
        f2L += on(_range=markedfaces(M_Qh->mesh(),it.second.second),_element=phi, _expr=it.second.first, _rhs=f1LQ, _type="elimination_symmetric");
    }

    init();
    toc( "[PreconditionerBlockMS] setup done ", FLAGS_v > 0 );
}

template < typename space_type >
void
PreconditionerBlockMS<space_type>::init( void )
{
    //Feel::cout << "Init preconditioner blockms\n";
    LOG(INFO) << "Init preconditioner blockms...\n";
    tic();
    CHECK( false ) << "TODO fix bc";
    BoundaryConditions M_bc;// = M_model.boundaryConditions();

    LOG(INFO) << "Create sub Matrix\n";
    map_vector_field<FEELPP_DIM,1,2> m_dirichlet_u { M_bc.getVectorFields<FEELPP_DIM> ( "u", "Dirichlet" )};
    //map_scalar_field<2> m_dirichlet_p { M_bc.getScalarFields<2> ( "phi", "Dirichlet" ) };
    map_vector_field<FEELPP_DIM,1,2> m_weak_u { M_bc.getVectorFields<FEELPP_DIM> ( "u", "Weakdir" ) };
    //map_scalar_field<2> m_weak_p { M_bc.getScalarFields<2> ( "phi", "Weakdir" ) };

    /*
     * AA = [[ A - k^2 M, B^t],
     *      [ B        , 0  ]]
     * We need to extract A-k^2 M and add it M to form A+(1-k^2) M = A+g M
     */
    // Is the zero() necessary ?
    M_11->zero();
    this->matrix()->updateSubMatrix(M_11, M_Vh_indices, M_Vh_indices, false); // M_11 = A-k^2 M
    LOG(INFO) << "Use relax = " << M_relax << std::endl;
    M_11->addMatrix(M_relax,M_mass);                            // A-k^2 M + M_relax*M = A+(M_relax-k^2) M
    auto f2A = form2(_test=M_Vh, _trial=M_Vh,_matrix=M_11);
    auto f1A = form1(_test=M_Vh);
    for(auto const & it : m_weak_u )
    {
        LOG(INFO) << "Applying (weak) " << it.second.first << " on " << it.first << " for "<<M_prefix_11<<"\n";
        f2A += integrate(_range=markedfaces(M_Vh->mesh(),it.second.second),
            _expr=-(1./M_relax)*trans(curlt_op(u))*(cross(N(),id(u)) )
            - (1./M_relax)*trans(curl_op(u))*(cross(N(),idt(u)) )
            + doption(_prefix=M_prefix,_name="penaldir")/(hFace()*M_relax)*inner(cross(idt(u),N()),cross(id(u),N())) );
    }
    for(auto const & it : m_dirichlet_u )
    {
        LOG(INFO) << "Applying (on) " << it.second.first << " on " << it.first << " for "<<M_prefix_11<<"\n";
        f2A += on(_range=markedfaces(M_Vh->mesh(),it.second.second), _expr=it.second.first,_rhs=f1A, _element=u, _type="elimination_symmetric");
    }

    /* 
     * Rebuilding sub-backend
     */
    backend(_name=M_prefix_11, _rebuild=true);
    backend(_name=M_prefix_22, _rebuild=true);
    // We have to set the G, Px,Py,Pz or X,Y,Z matrices to AMS
    if(soption(_name="pc-type", _prefix=M_prefix_11) == "ams")
    {
#if FEELPP_DIM == 3
    initAMS();
    {
        if(boption(_name="setAlphaBeta",_prefix=M_prefix_11))
        {
            auto prec = preconditioner(_pc=pcTypeConvertStrToEnum(soption(M_prefix_11+".pc-type")),
                                       _backend=backend(_name=M_prefix_11),
                                       _prefix=M_prefix_11,
                                       _matrix=M_11
                                      );
            prec->setMatrix(M_11);
            prec->attachAuxiliarySparseMatrix("a_alpha",M_a_alpha);
            prec->attachAuxiliarySparseMatrix("a_beta",M_a_beta);
        }
    }
#else
    std::cerr << "ams preconditioner is not interfaced in two dimensions\n";
#endif
    }
    toc("[PreconditionerBlockMS] Init",FLAGS_v>0);
    LOG(INFO) << "Init done\n";
}

template < typename space_type >
void
PreconditionerBlockMS<space_type>::initAMS( void )
    {
        M_grad  = Grad( _domainSpace=M_Qh, _imageSpace=M_Vh);

        // This preconditioner is linked to that backend : the backend will
        // automatically use the preconditioner.
        auto prec = preconditioner(_pc=pcTypeConvertStrToEnum(soption(M_prefix_11+".pc-type")),
                                   _backend=backend(_name=M_prefix_11),
                                   _prefix=M_prefix_11,
                                   _matrix=M_11
                                  );
        backend(_name=M_prefix_11)->attachPreconditioner( prec );
        prec->setMatrix(M_11);
        prec->attachAuxiliarySparseMatrix("G",M_grad.matPtr());
        if(boption(M_prefix_11+".useEdge"))
        {
            LOG(INFO) << "[ AMS ] : using SetConstantEdgeVector \n";
            ozz.on(_range=elements(M_Vh->mesh()),_expr=vec(cst(1),cst(0),cst(0)));
            zoz.on(_range=elements(M_Vh->mesh()),_expr=vec(cst(0),cst(1),cst(0)));
            zzo.on(_range=elements(M_Vh->mesh()),_expr=vec(cst(0),cst(0),cst(1)));
            *M_ozz = ozz; M_ozz->close();
            *M_zoz = zoz; M_zoz->close();
            *M_zzo = zzo; M_zzo->close();

            prec->attachAuxiliaryVector("Px",M_ozz);
            prec->attachAuxiliaryVector("Py",M_zoz);
            prec->attachAuxiliaryVector("Pz",M_zzo);
        }
        else
        {
            LOG(INFO) << "[ AMS ] : using SetCoordinates \n";
            X.on(_range=elements(M_Vh->mesh()),_expr=Px());
            Y.on(_range=elements(M_Vh->mesh()),_expr=Py());
            Z.on(_range=elements(M_Vh->mesh()),_expr=Pz());
            *M_X = X; M_X->close();
            *M_Y = Y; M_Y->close();
            *M_Z = Z; M_Z->close();
            prec->attachAuxiliaryVector("X",M_X);
            prec->attachAuxiliaryVector("Y",M_Y);
            prec->attachAuxiliaryVector("Z",M_Z);
        }
    }


template < typename space_type >
int
PreconditionerBlockMS<space_type>::applyInverse ( const vector_type& X, vector_type& Y ) const
{
    tic();
    U = X;
    U.close();
    *M_uin = U.template element<0>();
    M_uin->close();
    *M_pin = U.template element<1>();
    M_pin->close();

    // Solve eq (12)
    // solve here eq 15 : Pm v = c
    backend(_name=M_prefix_11)->solve(_matrix=M_11,
                                      _rhs=M_uin,
                                      _solution=M_uout
                                     ) ;
    M_uout->close();

    // solve here eq 16
    backend(_name=M_prefix_22)->solve(_matrix=M_L,
                                      _rhs=M_pin,
                                      _solution=M_pout
                                     );
    M_pout->close();

    U.template element<0>() = *M_uout;
    U.template element<1>() = *M_pout;
    U.close();
    Y=U;
    Y.close();
    toc("[PreconditionerBlockMS] applyInverse update solution",FLAGS_v>0);
    return 0;
}

template <typename ... Ts>
auto blockms( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && space = args.get(_space);
    auto && matrix = args.get(_matrix);
    auto && model = args.get(_model);
    std::string const& prefix = args.get_else(_prefix, "blockms" );

    using space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(space)>>>;
    return std::make_shared<PreconditionerBlockMS<space_type>>( space, model, prefix, matrix );
}


} // Feel
#endif
