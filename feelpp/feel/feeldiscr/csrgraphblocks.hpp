
#pragma once

#include <feel/feeldiscr/stencil.hpp>

namespace Feel
{
/**
 * Build blocks of CSR graphs with function spaces \p
 * (args1,args2,argn). Variadic template is used to handle an arbitrary number of
 * function spaces.
 * The blocks are organized then matrix wise with the stencil associated of pairs of function spaces in \p (arg1,...,argn)
 *
 */
template<typename PS, typename RangeMapT = StencilRangeMap0Type, 
         typename = std::enable_if_t<std::is_base_of_v<ProductSpacesBase,std::remove_reference_t<PS>>
                                    && std::is_base_of_v<StencilRangeMapTypeBase,RangeMapT>> >
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                uint32_type pattern = Pattern::COUPLED,
                RangeMapT range = stencilRangeMap() )
{
    int s = ps.numberOfSpaces();
    BlocksBaseGraphCSR g( s, s );

    int n = 0;
    auto pst = ps.tupleSpaces();
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

                        hana::if_(std::is_base_of<ProductSpaceBase,decay_type<decltype(test_space)>>{},
                                  [&]( auto&& xx, auto&& yy ) { return hana::if_( std::is_base_of<ProductSpaceBase,decay_type<decltype(yy)>>{},
                                                                                [&] (auto&&x,auto&& y) {
                                                                                    for( int i = 0; i < x->numberOfSpaces(); ++i)
                                                                                        for( int j = 0; j < y->numberOfSpaces(); ++j)
                                                                                        {
                                                                                            LOG(INFO) << "filling out stencil (" << r+i << "," << c+j << ")\n";
                                                                                            g( r+i, c+j ) =
                                                                                                stencil( _test=(*x)[i],
                                                                                                         _trial=(*y)[j],
                                                                                                         _pattern=pattern,
                                                                                                         _range=range,
                                                                                                     _diag_is_nonzero=false, _close=false)->graph();
                                                                                        }
                                                                                },
                                                                                [&] (auto&&x,auto && y){
                                                                                    for( int i = 0; i < x->numberOfSpaces(); ++i)
                                                                                    {
                                                                                        LOG(INFO) << "filling out stencil (" << r+i << "," << c << ")\n";
                                                                                        g( r+i, c ) =
                                                                                            stencil( _test=(*x)[i],
                                                                                                     _trial=y,
                                                                                                     _pattern=pattern,
                                                                                                     _range=range,
                                                                                                     _diag_is_nonzero=false, _close=false)->graph();
                                                                                    }
                                                                                })(xx,yy); },
                                  [&]( auto &&xx, auto &&yy ) { return hana::if_( std::is_base_of<ProductSpaceBase,decay_type<decltype(yy)>>{},
                                                                                      [&] (auto &&x, auto&& y) {
                                                                                          for( int i = 0; i < y->numberOfSpaces(); ++i)
                                                                                          {
                                                                                              LOG(INFO) << "filling out stencil (" << r << "," << c+i << ")\n";
                                                                                              g( r, c+i ) =
                                                                                                  stencil( _test=x,
                                                                                                           _trial=(*y)[i],
                                                                                                           _pattern=pattern,
                                                                                                           _range=range,
                                                                                                           _diag_is_nonzero=false, _close=false)->graph();
                                                                                          }
                                                                                      },
                                                                                      [&] (auto && x, auto &&y){
                                                                                          LOG(INFO) << "filling out stencil (" << r << "," << c << ")\n";
                                                                                          if ( r != c )
                                                                                            g( r, c ) =
                                                                                                  stencil( _test=x,
                                                                                                           _trial=y,
                                                                                                           _pattern=pattern,
                                                                                                           _range=range,
                                                                                                            _diag_is_nonzero=false, _close=false)->graph();
                                                                                          else
                                                                                            g( r, c ) =
                                                                                                  stencil( _test=x,
                                                                                                           _trial=y,
                                                                                                           _pattern=pattern,
                                                                                                           _diag_is_nonzero=false, _close=false)->graph();
                                                                                    
                                                                                      })(xx,yy); })(test_space,trial_space);



                        ++n;
                    });
    return g;
}

template<typename PS>
BlocksBaseGraphCSR
csrGraphBlocks( PS&& ps,
                uint32_type pattern = Pattern::COUPLED,
                std::enable_if_t<std::is_base_of<ProductSpaceBase,std::remove_reference_t<PS>>::value>* = nullptr )
{
    int s = ps.numberOfSpaces();
    BlocksBaseGraphCSR g( s, s );


    for( int i = 0; i < ps.numberOfSpaces(); ++i )
        for( int j = 0; j < ps.numberOfSpaces(); ++j )
        {
            g( i, j ) = stencil( _test=ps[i],_trial=ps[j], _pattern=pattern, _diag_is_nonzero=false, _close=false)->graph();
            LOG(INFO) << "filling out stencil (" << i << "," << j << ")\n";
        }
    return g;
}
}