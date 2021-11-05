/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 This file is part of the Feel++ library
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 24 Jul 2016
 Copyright (C) 2016 Feel++ Consortium
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
#ifndef FEELPP_BLOCKFORMS_HPP
#define FEELPP_BLOCKFORMS_HPP 1

#include <feel/feelvf/form.hpp>
#include <feel/feelalg/vectorblock.hpp>
#include <feel/feelalg/matrixcondensed.hpp>
#include <feel/feelalg/vectorcondensed.hpp>
#include <feel/feeldiscr/csrgraphblocks.hpp>
#include <feel/feeldiscr/product.hpp>


namespace Feel {

//!
//! forward declarations of @c BlockBilinearForm and @c blocform2()
//!
template<typename PS>
class BlockBilinearForm;


template<typename PS>
BlockBilinearForm<PS>
blockform2( PS&& ps );

template<typename PS, typename BackendT>
BlockBilinearForm<PS>
blockform2( PS&& ps, BackendT&& b );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
BlockBilinearForm<PS>
blockform2( PS&& ps, solve::strategy s, BackendT&& b,
            size_type pattern = Pattern::COUPLED,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
BlockBilinearForm<PS>
blockform2( PS&& ps, solve::strategy s, BackendT&& b,
            std::vector<size_type> const& patterns,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b,
            size_type pattern = Pattern::COUPLED,
            RangeMapT r = stencilRangeMap() );

template<typename PS, typename BackendT, typename RangeMapT = StencilRangeMap0Type>
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b,
            std::vector<size_type> const& patterns,
            RangeMapT r = stencilRangeMap() );

template<typename PS,typename T>
BlockBilinearForm<PS>
blockform2( PS&& ps, condensed_matrix_ptr_t<T> & m );

//!
//! forward declarations of @c BlockLinearForm and @c blockform1()
//!
template<typename PS>
class BlockLinearForm;

template<typename PS>
BlockLinearForm<PS>
blockform1( PS&& ps );

template<typename PS,typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, BackendT&& b );

template<typename PS,typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, solve::strategy s, BackendT&& b );


template<typename PS,typename T>
BlockLinearForm<PS>
blockform1( PS&& ps, condensed_vector_ptr_t<T> v );


/**
 * Handles bilinear form over a product of spaces
 */
template<typename PS>
class BlockBilinearForm
{
public :
    using value_type = typename Feel::decay_type<PS>::value_type;
    using condensed_matrix_type = MatrixCondensed<value_type>;
    using size_type = typename condensed_matrix_type::size_type;
    using condensed_matrix_ptrtype = std::shared_ptr<condensed_matrix_type>;
    using product_space_t = decay_type<PS>;

    BlockBilinearForm() = default;
    BlockBilinearForm( BlockBilinearForm const& ) = default;
    BlockBilinearForm( BlockBilinearForm && ) = default;
    
    template<typename T>
    BlockBilinearForm( T&& ps, std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value>* = nullptr)
        :
        M_ps(std::forward<T>(ps)),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_ps, Pattern::COUPLED), backend(), false ) )
        {}

    template<typename T>
    BlockBilinearForm( T&& ps, std::enable_if_t<std::is_base_of<ProductSpaceBase,decay_type<T>>::value>* = nullptr)
        :
        M_ps(std::forward<T>(ps)),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_ps, Pattern::COUPLED), backend(), false ) )
        {}    
    
    template<typename T,typename BackendT, typename RangeMapT>
    BlockBilinearForm( T&& ps, BackendT&& b, RangeMapT r = stencilRangeMap(),
                       std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value && std::is_base_of<BackendBase,decay_type<BackendT>>::value>* = nullptr )
        :
        M_ps(std::forward<T>(ps)),
        M_matrix( std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_ps, Pattern::COUPLED, r), std::forward<BackendT>(b), false ) )
        {}

    template<typename T, typename BackendT, typename RangeMapT>
    BlockBilinearForm( T&& ps, solve::strategy s, BackendT&& b, size_type pattern = Pattern::COUPLED, RangeMapT r = stencilRangeMap(),
                       std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value && std::is_base_of<BackendBase,decay_type<BackendT>>::value>* = nullptr )
        :
        M_ps(std::forward<T>(ps)),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:pattern,r),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {}
    template<typename T, typename BackendT>
    BlockBilinearForm( T&& ps, solve::strategy s, BackendT&& b, std::vector<size_type> const& patterns,
                       std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value && std::is_base_of<BackendBase,decay_type<BackendT>>::value>* = nullptr )
        :
        M_ps(std::forward<T>(ps)),
        M_matrix( std::make_shared<condensed_matrix_type>( s,
                                                             csrGraphBlocks(M_ps, (s>=solve::strategy::static_condensation)?pattern::toZero(patterns):patterns),
                                                             std::forward<BackendT>(b),
                                                             (s>=solve::strategy::static_condensation)?false:true )  )
        {}

    BlockBilinearForm(product_space_t&& ps, condensed_matrix_ptrtype & m)
        :
        M_ps(ps),
        M_matrix(m)
        {}

    BlockBilinearForm(product_space_t const& ps, condensed_matrix_ptrtype & m)
        :
        M_ps(ps),
        M_matrix(m)
        {}
    BlockBilinearForm& operator=( BlockBilinearForm const& a )
        {
            if ( this == &a )
                return *this;

            bool same_spaces = (M_ps == a.M_ps);
            M_ps = a.M_ps;
            if ( !this->isMatrixAllocated() || !same_spaces )
            {
                this->allocateMatrix( a.M_matrix->solveStrategy(), a.M_matrix->backend() );
                M_matrix->setBackend( a.M_matrix->backend()->clone() );
            }
            M_matrix->zero();
            M_matrix->addMatrix( 1.,(MatrixSparse<value_type> const&)*a.M_matrix->getSparseMatrix() );
            
            return *this;
        }
    BlockBilinearForm& operator=( BlockBilinearForm && a ) = default;
    BlockBilinearForm& operator+=( BlockBilinearForm& a )
        {
            if ( this == &a )
            {
                M_matrix->scale( 2.0 );
                return *this;
            }

            M_matrix->addMatrix( 1.0, a.M_matrix );

            return *this;
        }

    //!
    //! allocate algebraic representation of the bilinear form
    //! @param s the solve strategy (monolithic, static condensation or local)
    //! @param b the algebraic backend to use (petsc or eigen)
    //!
    void allocateMatrix( solve::strategy s = solve::strategy::monolithic, backend_ptrtype const& b = backend() )
        {
            M_matrix = std::make_shared<condensed_matrix_type>( s, csrGraphBlocks(M_ps, (s==solve::strategy::static_condensation)?Pattern::ZERO:Pattern::COUPLED), b, (s==solve::strategy::static_condensation)?false:true );
        }
    //!
    //! @return true if allocated, false otherwise
    //!
    bool isMatrixAllocated() const
        {
            return (bool)M_matrix;
        }
#if 0
    template<typename N1,typename N2>
    decltype(auto) operator()( N1 n1, N2 n2 )
        {
            cout << "filling out matrix block (" << n1 << "," << n2 << ")\n";

        }
#endif
    template<typename N1, typename N2>
    decltype(auto) operator()( N1 n1, N2 n2, int s1 = 0, int s2 = 0 )
        {
            int n = 0;
            auto&& spaces=remove_shared_ptr_f( M_ps );
            #if 0
            hana::if_( hana::bool_<Feel::is_shared_ptr_v<PS>>{},
                                     []( auto&& x ) { return *x; },
                                     []( auto&& x ) { return x; } )(M_ps);
            #endif
            auto test_space = hana::at( spaces.tupleSpaces(), n1 );
            auto trial_space = hana::at( spaces.tupleSpaces(), n2 );

            return hana::eval_if(std::is_base_of<ProductSpaceBase,decay_type<decltype(test_space)>>{},
                                 [&]( auto _ ) { return hana::eval_if( std::is_base_of<ProductSpaceBase,decay_type<decltype(trial_space)>>{},
                                                                 [&] (auto _) {
                                                                     LOG(INFO) << "filling out dyn matrix block (" << int(n1) + s1 << "," << int(n2)+s2  << ")\n";
                                                                     return form2(_test=(*_(test_space))[s1],_trial=(*_(trial_space))[s2],
                                                                                  _name="bilinearform.a"s+"("+std::to_string(int(n1)+s1)+","s+std::to_string(int(n2)+s2)+")"s,
                                                                                  _matrix=M_matrix->block(int(n1)+s1,int(n2)+s2), _rowstart=int(n1)+s1, _colstart=int(n2)+s2 );
                                                                 },
                                                                 [&] (auto _){
                                                                     LOG(INFO) << "filling out dyn matrix block (" << int(n1) + s1 << "," << int(n2)+s2  << ")\n";
                                                                     return form2(_test=(*_(test_space))[s1],_trial=_(trial_space),
                                                                                  _name="bilinearform.a"s+"("+std::to_string(int(n1)+s1)+","s+std::to_string(int(n2))+")"s,
                                                                                  _matrix=M_matrix->block(int(n1)+s1,int(n2)), _rowstart=int(n1)+s1, _colstart=int(n2) );
                                                                 }); },
                                 [&]( auto _ ) { return hana::eval_if( std::is_base_of<ProductSpaceBase,decay_type<decltype(trial_space)>>{},
                                                                 [&] (auto _) {
                                                                     LOG(INFO) << "filling out dyn matrix block (" << int(n1) + s1 << "," << int(n2)+s2  << ")\n";
                                                                     return form2(_test=_(test_space),_trial=(*_(trial_space))[s2],
                                                                                  _name="bilinearform.a"s+"("+std::to_string(int(n1))+","s+std::to_string(int(n2)+s2)+")"s,
                                                                                  _matrix=M_matrix->block(int(n1),int(n2)+s2), _rowstart=int(n1), _colstart=int(n2)+s2 );
                                                                 },
                                                                 [&] (auto _){
                                                                     LOG(INFO) << "filling out dyn matrix block (" << int(n1) + s1 << "," << int(n2)+s2  << ")\n";
                                                                     return form2(_test=_(test_space),_trial=_(trial_space),
                                                                                  _name="bilinearform.a"s+"("+std::to_string(int(n1))+","s+std::to_string(int(n2))+")"s,
                                                                                  _matrix=M_matrix->block(int(n1),int(n2)), _rowstart=int(n1), _colstart=int(n2) );
                                                                 }); });


        }

    decltype(auto) operator()( int n1, int n2 )
        {
            cout << "filling out matrix block (" << n1 << "," << n2 << ")\n";
            return form2(_test=M_ps[n1],_trial=M_ps[n2], _matrix=M_matrix, _rowstart=int(n1), _colstart=int(n2) );
        }

    template<typename T>
    void setFunctionSpace( T&& ps )
        {
            M_ps = std::forward<T>(ps);
        }
    template<typename BackendT>
    void setStrategy( BackendT&& b )
        {
            
            M_matrix = std::make_shared<condensed_matrix_type>( csrGraphBlocks(M_ps, Pattern::COUPLED), std::forward<BackendT>(b), false );
        }
    template<typename BackendT>
    void setStrategy( solve::strategy s, BackendT&& b, size_type pattern = Pattern::COUPLED )
        {
            M_matrix =  std::make_shared<condensed_matrix_type>( s,
                                                                   csrGraphBlocks(M_ps, (s>=solve::strategy::static_condensation)?Pattern::ZERO:pattern),
                                                                   std::forward<BackendT>(b),
                                                                   (s>=solve::strategy::static_condensation)?false:true );
        }
    template<typename BackendT>
    void setStrategy( solve::strategy s, BackendT&& b, std::vector<size_type> const& patterns )
        {
            M_matrix =  std::make_shared<condensed_matrix_type>( s,
                                                                   csrGraphBlocks(M_ps, (s>=solve::strategy::static_condensation)?pattern::toZero(patterns):patterns),
                                                                   std::forward<BackendT>(b),
                                                                   (s>=solve::strategy::static_condensation)?false:true );
        }
    void close()
        {
            M_matrix->close();
        }
    void zero()
        {
            M_matrix->zero();
        }
    void zero(int n1, int n2 )
        {
            M_matrix->zero( n1, n2 );
        }
    void transpose(int n1, int n2 )
        {
            M_matrix->transposeBlock( n1, n2 );
        }
    //!
    //! @return the number of non-zero entries in matrix representation
    //!
    std::size_t nnz() const { return M_matrix->nnz(); }
    
    void syncLocalMatrix()
        {
            int s = M_ps.numberOfSpaces();
            int n = 0;
            auto pst = M_ps.tupleSpaces();
            auto cp = hana::cartesian_product( hana::make_tuple( pst, pst ) );
            int nstatic = hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(hana::back(pst))>>{},
                                    [s] (auto&& x ) { return s-hana::back(std::forward<decltype(x)>(x))->numberOfSpaces()+1; },
                                    [s] (auto&& x ) { return s; } )( pst );
            hana::for_each( cp, [&]( auto const& e )
                            {
                                int r = n/nstatic;
                                int c = n%nstatic;
                                auto test_space = e[0_c];
                                auto trial_space = e[1_c];
                                DVLOG( 1 ) << "syncLocalMatrix("<< r << ","<< c <<")dim "<< test_space->mesh()->dimension()<<" , " << trial_space->mesh()->dimension() <<"\n";
                                M_matrix->sc(r,c)->syncLocalMatrix( test_space,trial_space );
                                ++n;
                            });
        }
    sparse_matrix_ptrtype matrixPtr() const { return M_matrix; }
    sparse_matrix_ptrtype matrixPtr() { return M_matrix; }
    using pre_solve_type = typename Backend<value_type>::pre_solve_type;
    using post_solve_type = typename Backend<value_type>::post_solve_type;
#if 0
    BOOST_PARAMETER_MEMBER_FUNCTION( ( typename Backend<double>::solve_return_type ),
                                     solve,
                                     tag,
                                     ( required
                                       ( in_out( solution ),* )
                                       ( rhs, * ) )
                                     ( optional
                                       ( condense,           ( bool ), false )
                                       ( condenser,     (*), condenser_poisson() )
                                       ( local,          ( bool ), false )
                                       ( name,           ( std::string ), "" )
                                       ( kind,           ( std::string ), soption(_prefix=name,_name="backend") )
                                       ( rebuild,        ( bool ), boption(_prefix=name,_name="backend.rebuild") )
                                       ( pre, (pre_solve_type), pre_solve_type() )
                                       ( post, (post_solve_type), post_solve_type() )
                                       ) )
#endif
    template <typename ... Ts>
    typename Backend<double>::solve_return_type solve( Ts && ... v )
        {
            auto args = NA::make_arguments( std::forward<Ts>(v)... );
            auto && solution = args.get(_solution);
            auto && rhs = args.get(_rhs);
            bool condense = args.get_else(_condense,false);
            auto && condenser = args.get_else(_condenser,condenser_poisson() );
            bool local = args.get_else(_local,false);
            std::string const& name = args.get_else(_name,"" );
            std::string const& kind = args.get_else_invocable(_kind,[&name](){ return soption(_prefix=name,_name="backend"); } );
            bool rebuild = args.get_else_invocable(_rebuild,[&name](){ return boption(_prefix=name,_name="backend.rebuild"); } );
            pre_solve_type pre = args.get_else(_pre, pre_solve_type() );
            post_solve_type post = args.get_else(_post,post_solve_type() );

            if constexpr ( std::is_same_v<decay_type<decltype(condenser)>,condenser_stokes> )
                return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
            else
            {
                if ( condense )
                    return solveImplCondense( M_ps, solution, rhs, name, kind, rebuild, pre, post, condenser );
                if ( local )
                {
                    return solveImplLocal( M_ps, solution, rhs, name, kind, rebuild, pre, post );
                }
                return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
            }
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplLocal( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                    bool rebuild, pre_solve_type pre, post_solve_type post,
                    std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            auto sc = M_matrix->sc();
            tic();
            cout << " . starting local Solve" << std::endl;
            auto vsc = rhs.vectorPtr()->sc();
            sc->localSolve ( vsc, solution);
            cout << " . local Solve done" << std::endl;
            toc("blockform.local.localsolve",FLAGS_v>0);
            typename Backend<double>::solve_return_type r;
            return r;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplLocal( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                    bool rebuild, pre_solve_type pre, post_solve_type post,
                    std::enable_if_t<std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            typename Backend<double>::solve_return_type r;
            return r;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t, typename CT>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post, CT ct,
                       std::enable_if_t<std::is_base_of_v<ProductSpaceBase,decay_type<PS_t>> && is_condenser_v<CT>>* = nullptr )
        {
            return solveImpl( solution, rhs, name, kind, rebuild, pre, post );
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t, typename CT>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post, CT ct,
                       std::enable_if_t< hana::Foldable<typename decay_type<PS_t>::tuple_spaces_type>::value &&
                                         !std::is_base_of_v<ProductSpaceBase,decay_type<PS_t>> &&
                                         is_condenser_v<CT>>* = nullptr )
        {
            return solveImplCondense( ps, solution, rhs, name, kind, rebuild, pre, post,hana::integral_constant<int,decltype(hana::size( M_ps.tupleSpaces() ))::value>() );
        }
    template <typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImpl( Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind = "petsc",
               bool rebuild = false, pre_solve_type pre = pre_solve_type(), post_solve_type post = post_solve_type() )
        {
            auto U = backend()->newBlockVector(_block=solution, _copy_values=false);
            tic();
            auto r1 = backend( _name=name, _kind=kind, _rebuild=rebuild,
                               _worldcomm=Environment::worldCommPtr() )->solve( _matrix=M_matrix->getSparseMatrix(),
                                                                             _rhs=rhs.vectorPtr()->getVector(),
                                                                             _solution=U,
                                                                             _pre=pre,
                                                                             _post=post
                                                                             );
            toc("blockform.monolithic",FLAGS_v>0);
            if ( Environment::isSequential() && boption("exporter.matlab") )
            {
                M_matrix->getSparseMatrix()->printMatlab("A.m");
                rhs.vectorPtr()->getVector()->printMatlab("b.m");
            } 
            solution.localize(U);
            return r1;
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,1>,
                       std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            return typename Backend<double>::solve_return_type{};
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t> 
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,2>,
                       std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            return typename Backend<double>::solve_return_type{};
        }
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,3>,
                       std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
#if 1
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();
            tic();
            auto psS = product( e3.functionSpace() );
            toc("blockform.sc.space",FLAGS_v>0);
            tic();
            auto S = blockform2( psS, solve::strategy::monolithic, backend(), Pattern::HDG  );
            toc("blockform.sc.bilinearform",FLAGS_v>0);
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto V = blockform1( psS, solve::strategy::monolithic, backend() );
            
            tic();
            this->syncLocalMatrix();
            toc("blockform.sc.sync", FLAGS_v>0);
            tic();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", FLAGS_v>0);
            S.close();V.close();
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;
            auto U = psS.element();

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );
            //auto r = backend(_name=prefixvm(name,"sc"),_rebuild=rebuild)->solve( _matrix=S.matrixPtr(), _rhs=V.vectorPtr(), _solution=e3);
            solution(2_c)=U(0_c);
            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", FLAGS_v>0);

#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->localSolve ( rhs.vectorPtr()->sc(), solution);
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",FLAGS_v>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m");
#endif


            return r;
#else
            return {};
#endif
        }

    //!
    //! solve using static condensation in the case of 2 trace spaces
    //!
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,4>,
                       std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();

            auto Th = product2( M_ps[3_c], M_ps[2_c] );
            auto S = blockform2(Th, solve::strategy::monolithic, backend(), Pattern::HDG);
            auto V = blockform1(Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", FLAGS_v>0);
            S.close();V.close();
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );

            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", FLAGS_v>0);

            solution(2_c)=U(0_c);
            for( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution(3_c,i)=U(1_c,i);
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_ps[3_c]->numberOfSpaces());
            sc->localSolve ( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",FLAGS_v>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m"); 
#endif

            return r;
        }
 
    //!
    //! solve using static condensation in the case of 2 trace spaces
    //!
    template <typename PS_t, typename Solution_t, typename Rhs_t>
    typename Backend<double>::solve_return_type
    solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                       bool rebuild, pre_solve_type pre, post_solve_type post , hana::integral_constant<int,5>,
                       std::enable_if_t<!std::is_base_of<ProductSpaceBase,decay_type<PS_t>>::value>* = nullptr )
        {
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution(2_c);
            auto& e2 = solution(1_c);
            auto& e1 = solution(0_c);

            auto sc = M_matrix->sc();

            auto Th = product2( M_ps[3_c], M_ps[4_c], ps[2_c] );
            std::vector<size_type> patterns = {Pattern::HDG,Pattern::HDG,Pattern::ZERO,
                                               Pattern::HDG,Pattern::HDG,Pattern::COUPLED,
                                               Pattern::ZERO,Pattern::COUPLED,Pattern::COUPLED};
            auto S = blockform2(Th, solve::strategy::monolithic, backend(), patterns);
            auto V = blockform1(Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense ( rhs.vectorPtr()->sc(), solution, S, V );
            toc("blockform.sc.condense", FLAGS_v>0);
            S.close();V.close();
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution=U, _rhs=V, _name=prefixvm(name,"sc"),_rebuild=rebuild );//, _condense=true );

            cout << " . Solve done" << std::endl;
            toc("blockform.sc.solve", FLAGS_v>0);

            solution(2_c)=U(0_c);
            for( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution(3_c,i)=U(1_c,i);
            for( int i = 0; i < Th[2_c]->numberOfSpaces(); ++i )
                solution(4_c,i) = U(2_c,i);
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_ps[3_c]->numberOfSpaces());
            sc->localSolve ( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc("blockform.sc.localsolve",FLAGS_v>0);
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m"); 
#endif

            return r;
        }

        //!
        //! solve using static condensation in the case of 2 trace spaces
        //!
        template <typename PS_t, typename Solution_t, typename Rhs_t>
        typename Backend<double>::solve_return_type
        solveImplCondense( PS_t& ps, Solution_t& solution, Rhs_t const& rhs, std::string const& name, std::string const& kind,
                           bool rebuild, pre_solve_type pre, post_solve_type post, hana::integral_constant<int, 3*5>,
                           std::enable_if_t<!std::is_base_of<ProductSpaceBase, decay_type<PS_t>>::value>* = nullptr )
        {
#if 0            
            //auto& e4 = solution(3_c,0);
            auto& e3 = solution( 2_c );
            auto& e2 = solution( 1_c );
            auto& e1 = solution( 0_c );

            auto sc = M_matrix->sc();

            auto Th = product2( M_ps[3_c], M_ps[4_c], ps[2_c] );
            std::vector<size_type> patterns = { Pattern::HDG, Pattern::HDG, Pattern::ZERO,
                                                Pattern::HDG, Pattern::HDG, Pattern::COUPLED,
                                                Pattern::ZERO, Pattern::COUPLED, Pattern::COUPLED };
            auto S = blockform2( Th, solve::strategy::monolithic, backend(), patterns );
            auto V = blockform1( Th, solve::strategy::monolithic, backend() );
            //MatSetOption ( dynamic_cast<MatrixPetsc<double>*>(S.matrixPtr().get())->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE );
            auto U = Th.element();
            //e3.printMatlab("phat.m");
            tic();
            this->syncLocalMatrix();
            sc->condense( rhs.vectorPtr()->sc(), solution, S, V );
            toc( "blockform.sc.condense", FLAGS_v > 0 );
            S.close();
            V.close();
            cout << " . Condensation done" << std::endl;
            tic();
            cout << " . starting Solve" << std::endl;

            auto r = S.solve( _solution = U, _rhs = V, _name = prefixvm( name, "sc" ), _rebuild = rebuild ); //, _condense=true );

            cout << " . Solve done" << std::endl;
            toc( "blockform.sc.solve", FLAGS_v > 0 );

            solution( 2_c ) = U( 0_c );
            for ( int i = 0; i < Th[1_c]->numberOfSpaces(); ++i )
                solution( 3_c, i ) = U( 1_c, i );
            for ( int i = 0; i < Th[2_c]->numberOfSpaces(); ++i )
                solution( 4_c, i ) = U( 2_c, i );
#if 0
            S.matrixPtr()->printMatlab("S.m");
            V.vectorPtr()->printMatlab("g.m");
            e3.printMatlab("phat1.m");
            e1.printMatlab("u.m");
            e2.printMatlab("p.m");
#endif
            tic();
            cout << " . starting local Solve" << std::endl;
            sc->setDim4( M_ps[3_c]->numberOfSpaces() );
            sc->localSolve( rhs.vectorPtr()->sc(), solution );
            cout << " . local Solve done" << std::endl;
            toc( "blockform.sc.localsolve", FLAGS_v > 0 );
#if 0
            e1.printMatlab("u1.m");
            e2.printMatlab("p1.m");
#endif

            return r;
#else
            return typename Backend<double>::solve_return_type{};
#endif
        }
        product_space_t M_ps;
        condensed_matrix_ptrtype M_matrix;
};

template<typename PS>
BlockBilinearForm<PS>
blockform2( PS&& ps )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps) );
}

template<typename PS, typename BackendT>
BlockBilinearForm<PS>
blockform2( PS&& ps, BackendT&& b )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps), std::forward<BackendT>(b) );
}

template<typename PS, typename BackendT, typename RangeMapT>
BlockBilinearForm<PS>
blockform2( PS && ps, solve::strategy s, BackendT&& b, size_type pattern, RangeMapT r )
{
    return BlockBilinearForm<PS>( std::forward<PS>( ps ), s, std::forward<BackendT>(b), pattern, r );
}

template<typename PS, typename BackendT, typename RangeMapT>
BlockBilinearForm<PS>
blockform2( PS && ps, solve::strategy s, BackendT&& b, std::vector<size_type>  const& patterns, RangeMapT r )
{
    return BlockBilinearForm<PS>( std::forward<PS>( ps ), s, std::forward<BackendT>(b), patterns, r );
}

template<typename PS, typename BackendT, typename RangeMapT>
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b, size_type pattern, RangeMapT r )
{
    return BlockBilinearForm<PS>( ps, s, std::forward<BackendT>(b), pattern, r );
}

template<typename PS, typename BackendT, typename RangeMapT>
BlockBilinearForm<PS>
blockform2( PS const& ps, solve::strategy s, BackendT&& b, std::vector<size_type> const& patterns, RangeMapT r )
{
    return BlockBilinearForm<PS>( ps, s, std::forward<BackendT>(b), patterns, r );
}

template<typename PS,typename T>
BlockBilinearForm<PS>
blockform2( PS&& ps, condensed_matrix_ptr_t<T> & m )
{
    return BlockBilinearForm<PS>( std::forward<PS>(ps), m );
}
/**
 * Handles linear form over a product of spaces
 */
template<typename PS>
class BlockLinearForm
{
public :
    using value_type = typename decay_type<PS>::value_type;
    using product_space_t = decay_type<PS>;
    using condensed_vector_type = VectorCondensed<value_type>;
    using condensed_vector_ptrtype = std::shared_ptr<condensed_vector_type>;
    using vector_ptrtype = condensed_vector_ptrtype;
    
    BlockLinearForm() = default;
    BlockLinearForm( BlockLinearForm const& ) = default;

    template<typename T, typename BackendT>
    BlockLinearForm( T&& ps, solve::strategy s, BackendT&& b, std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value>* = nullptr )
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(s, blockVector(M_ps), std::forward<BackendT>(b), false))
        {}
    template<typename T>
    BlockLinearForm(T&& ps, std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value>* = nullptr)
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), backend(), false))
        {}
    template<typename T>
    BlockLinearForm(T&& ps, std::enable_if_t<std::is_base_of<ProductSpaceBase,decay_type<T>>::value>* = nullptr)
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), backend(), false))
        {}    
    template<typename T, typename BackendT>
    BlockLinearForm(T&& ps, BackendT&& b, std::enable_if_t<std::is_base_of<ProductSpacesBase,decay_type<T>>::value>* = nullptr )
        :
        M_ps(std::forward<T>(ps)),
        M_vector(std::make_shared<condensed_vector_type>(blockVector(M_ps), std::forward<BackendT>(b), false))
        {}
    template<typename T>
    BlockLinearForm(T&& ps, condensed_vector_ptrtype v )
        :
        M_ps(std::forward<T>(ps)),
        M_vector(v)
        {}
    BlockLinearForm( BlockLinearForm&& ) = default;
    BlockLinearForm& operator=( BlockLinearForm && lf ) = default;
    BlockLinearForm& operator=( BlockLinearForm const& lf )
        {
            if ( this == &lf )
                return *this;

            bool same_spaces = ( M_ps == lf.M_ps );
            M_ps = lf.M_ps;
            if ( !M_vector || !same_spaces )
            {
                M_vector = std::make_shared<condensed_vector_type>( lf.M_vector->solveStrategy(), blockVector(M_ps), lf.M_vector->backend(), false );
                M_vector->setBackend( lf.M_vector->backend()->clone() );
            }
            M_vector->zero();
            M_vector->add( 1., *lf.M_vector->getVector() );
            
            return *this;
        }

#if 0
    template<typename N1>
    decltype(auto) operator()( N1 n1 )
        {
            cout << "filling out vector block (" << n1 << ")\n";
            return form1(_test=M_ps[n1],_vector=M_vector, _rowstart=int(n1) );
        }
#endif
    template<typename N1>
    decltype(auto) operator()( N1 n1, int s = 0 )
        {
            int n = 0;
            auto&& spaces=hana::if_( hana::bool_<Feel::is_shared_ptr_v<PS>>{},
                                     []( auto&& x ) { return *x; },
                                     []( auto&& x ) { return x; } )(M_ps);
            auto space = hana::at( spaces.tupleSpaces(), n1 );

            return hana::eval_if(std::is_base_of<ProductSpaceBase,decay_type<decltype(space)>>{},
                                 [&] (auto _) {
                                     VLOG(2) << "filling out dyn vector block (" << int(n1) + s  << ") condense=" << M_vector->staticCondensation() << "\n";
                                     return form1(_test=(*_(space))[s],_vector=M_vector->block(int(n1)+s), _rowstart=int(n1)+s );
                                 },
                                 [&] (auto _){
                                     VLOG(2) << "filling out vector block (" << n1  << ") condense=" << M_vector->staticCondensation() << "\n";
                                     return form1(_test=_(space),_vector=M_vector->block(int(n1)), _rowstart=int(n1) );
                                 });
        }

    decltype(auto) operator()( int n1 )
        {
            VLOG(2) << "filling out vector block (" << n1 << ") condense=" << M_vector->staticCondensation() << "\n";
            return form1(_test=M_ps[n1],_vector=M_vector->block(int(n1)), _rowstart=int(n1) );
        }
    template<typename T>
    void setFunctionSpace( T&& ps )
        {
            M_ps = std::forward<T>(ps);
        }
    template<typename BackendT>
    void setStrategy( BackendT&& b )
        {
            M_vector = std::make_shared<condensed_vector_type>(blockVector(M_ps), std::forward<BackendT>(b), false);
            
            
        }
    template<typename BackendT>
    void setStrategy( solve::strategy s, BackendT&& b )
        {
            M_vector = std::make_shared<condensed_vector_type>(s, blockVector(M_ps), std::forward<BackendT>(b), false);
        }
    void close()
        {
            M_vector->close();
        }
    product_space_t functionSpace() const { return M_ps; }

    condensed_vector_ptrtype const& vectorPtr() const { return M_vector; }
    condensed_vector_ptrtype vectorPtr() { return M_vector; }

    /**
     * set linear form to 0
     */
    void zero() { M_vector->zero(); }

    //!
    //! zero block @param n1
    //!
    void zeroBlock( int n1 ) { M_vector->zeroBlock( n1 ); }

    BlockLinearForm& operator+=( BlockLinearForm const& l )
        {
            if ( this == &l )
            {
                M_vector->scale( 2. );
                return *this;
            }

            *M_vector += *l.M_vector;

            return *this;
        }
    product_space_t M_ps;
    condensed_vector_ptrtype M_vector;
};

template<typename PS>
using blockform1_t = BlockLinearForm<PS>;
template<typename PS>
using blockform2_t = BlockBilinearForm<PS>;


template<typename PS>
BlockLinearForm<PS>
blockform1( PS&& ps )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps) );
}

template<typename PS, typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, BackendT&& b )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), std::forward<BackendT>(b) );
}

template<typename PS, typename BackendT>
BlockLinearForm<PS>
blockform1( PS&& ps, solve::strategy s, BackendT&& b )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), s, std::forward<BackendT>(b) );
}

template<typename PS,typename T>
BlockLinearForm<PS>
blockform1( PS&& ps, condensed_vector_ptr_t<T> v )
{
    return BlockLinearForm<PS>( std::forward<PS>(ps), v );
}


}
#endif
