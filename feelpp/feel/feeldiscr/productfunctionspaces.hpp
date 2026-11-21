/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- */

#if !defined(FEELPP_FEELDISCR_PRODUCTFUNCTIONSPACES_HPP)
#define FEELPP_FEELDISCR_PRODUCTFUNCTIONSPACES_HPP 1

#include <tuple>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/vector.hpp>

#include <feel/feeldiscr/dofcomposite.hpp>
#include <feel/feeldiscr/product.hpp>

namespace Feel
{

template<typename ProductSpacesType>
class ProductFunctionSpaces;

template<typename ProductSpacesType>
using LegacyCompositeSpaceAdapter = ProductFunctionSpaces<ProductSpacesType>;

template<typename... SpaceList>
class ProductFunctionSpaces<ProductSpaces<SpaceList...>>
    : public FunctionSpaceBase,
      public ProductSpacesBase
{
public:
    static_assert( sizeof...(SpaceList) > 1,
                   "ProductFunctionSpaces is a compatibility facade for composite spaces" );

    using product_space_type = ProductSpaces<SpaceList...>;
    using functionspace_type = ProductFunctionSpaces<product_space_type>;
    using space_type = functionspace_type;
    using pointer_type = std::shared_ptr<functionspace_type>;
    using ptrtype = pointer_type;
    using functionspace_ptrtype = pointer_type;
    using tuple_spaces_type = typename product_space_type::tuple_spaces_type;
    using functionspace_vector_type = boost::fusion::vector<SpaceList...>;
    using first_space_ptrtype = std::tuple_element_t<0, std::tuple<SpaceList...>>;
    using first_space_type = typename first_space_ptrtype::element_type;
    using value_type = typename first_space_type::value_type;
    using mesh_type = typename first_space_type::mesh_type;
    using mesh_ptrtype = typename first_space_type::mesh_ptrtype;
    using index_type = typename first_space_type::index_type;
    using size_type = typename first_space_type::size_type;
    using datamap_ptrtype = std::shared_ptr<DataMap<>>;
    using dof_type = DofComposite;
    using dof_ptrtype = std::shared_ptr<dof_type>;

    static inline const bool is_composite = true;
    static constexpr bool legacy_composite_enabled = false;
    static constexpr bool is_legacy_composite = false;
    static constexpr bool uses_internal_composite = false;
    static constexpr bool is_product_backed_composite = true;
    static constexpr uint16_type nSpaces = sizeof...(SpaceList);

    template<int I>
    struct sub_functionspace
    {
        static_assert( I >= 0 && I < static_cast<int>( nSpaces ), "invalid product function-space index" );
        using ptrtype = std::tuple_element_t<static_cast<std::size_t>( I ), std::tuple<SpaceList...>>;
        using type = typename ptrtype::element_type;
    };
    template<int I> using sub_functionspace_type = typename sub_functionspace<I>::type;
    template<int I> using sub_functionspace_ptrtype = typename sub_functionspace<I>::ptrtype;

    class Element;
    using element_type = Element;
    using element_ptrtype = std::shared_ptr<element_type>;

    explicit ProductFunctionSpaces( SpaceList... spaces )
        :
        FunctionSpaceBase( "ProductFunctionSpaces", firstWorldComm( spaces... ) ),
        M_productSpaces( spaces... ),
        M_functionspaces( spaces... ),
        M_worldsComm( makeWorldsComm( nSpaces, firstWorldComm( spaces... ) ) ),
        M_dof( this->makeDofComposite() )
    {
    }

    static ptrtype New( SpaceList... spaces )
    {
        return std::make_shared<functionspace_type>( spaces... );
    }

    int numberOfSpaces() const { return M_productSpaces.numberOfSpaces(); }
    uint16_type nSubFunctionSpace() const { return nSpaces; }

    product_space_type& productSpace() { return M_productSpaces; }
    product_space_type const& productSpace() const { return M_productSpaces; }

    tuple_spaces_type const& tupleSpaces() const { return M_productSpaces.tupleSpaces(); }
    tuple_spaces_type& tupleSpaces() { return M_productSpaces.tupleSpaces(); }

    functionspace_vector_type const& functionSpaces() const { return M_functionspaces; }
    functionspace_vector_type& functionSpaces() { return M_functionspaces; }

    template<int I>
    decltype(auto) space() const
    {
        return M_productSpaces.template space<I>();
    }

    template<int I>
    decltype(auto) space()
    {
        return M_productSpaces.template space<I>();
    }

    template<typename N>
    decltype(auto) space( N const& n ) const
    {
        return M_productSpaces.space( n );
    }

    template<typename N>
    decltype(auto) space( N const& n )
    {
        return M_productSpaces.space( n );
    }

    template<int I>
    typename boost::fusion::result_of::at_c<functionspace_vector_type const, I>::type
    functionSpace() const
    {
        return boost::fusion::at_c<I>( M_functionspaces );
    }

    template<int I>
    typename boost::fusion::result_of::at_c<functionspace_vector_type, I>::type
    functionSpace()
    {
        return boost::fusion::at_c<I>( M_functionspaces );
    }

    mesh_ptrtype mesh() const { return this->template functionSpace<0>()->mesh(); }
    template<int I>
    decltype(auto) mesh() const { return this->template functionSpace<I>()->mesh(); }

    worldscomm_ptr_t& worldsComm() { return M_worldsComm; }
    worldscomm_ptr_t const& worldsComm() const { return M_worldsComm; }

    WorldComm& worldComm() { return this->template functionSpace<0>()->worldComm(); }
    WorldComm const& worldComm() const { return this->template functionSpace<0>()->worldComm(); }

    size_type nDof() const { return M_productSpaces.nDof(); }
    size_type nLocalDof() const { return M_productSpaces.nLocalDof(); }
    size_type nLocalDofWithGhost() const { return this->nLocalDof(); }
    size_type nLocalDofWithoutGhost() const
    {
        return this->sumSpaces( []( auto const& space ) { return space->nLocalDofWithoutGhost(); } );
    }
    size_type nLocalDofWithGhost( rank_type proc ) const
    {
        return this->sumSpaces( [proc]( auto const& space ) { return space->nLocalDofWithGhost( proc ); } );
    }
    size_type nLocalDofWithoutGhost( rank_type proc ) const
    {
        return this->sumSpaces( [proc]( auto const& space ) { return space->nLocalDofWithoutGhost( proc ); } );
    }

    size_type nDofStart( size_type i = 0 ) const { return M_productSpaces.nDofStart( i ); }
    size_type nLocalDofStart( size_type i = 0 ) const { return M_productSpaces.nLocalDofStart( i ); }
    size_type nLocalDofWithGhostStart( size_type i = 0 ) const { return this->nLocalDofStart( i ); }
    size_type nLocalDofWithoutGhostStart( size_type i = 0 ) const
    {
        return this->sumSpacesUntil( i, []( auto const& space ) { return space->nLocalDofWithoutGhost(); } );
    }
    size_type nLocalDofWithGhostOnProcStart( rank_type proc, size_type i = 0 ) const
    {
        return this->sumSpacesUntil( i, [proc]( auto const& space ) { return space->nLocalDofWithGhost( proc ); } );
    }
    size_type nLocalDofWithoutGhostOnProcStart( rank_type proc, size_type i = 0 ) const
    {
        return this->sumSpacesUntil( i, [proc]( auto const& space ) { return space->nLocalDofWithoutGhost( proc ); } );
    }

    size_type blockDofStart( size_type i = 0 ) const { return M_productSpaces.blockDofStart( i ); }
    size_type blockLocalDofStart( size_type i = 0 ) const { return M_productSpaces.blockLocalDofStart( i ); }
    datamap_ptrtype blockMapPtr( size_type i ) const { return M_productSpaces.blockMapPtr( i ); }

    DataMap<> const& map() const { return *M_dof; }
    datamap_ptrtype mapPtr() const override { return M_dof; }
    dof_ptrtype const& dof() const { return M_dof; }
    dof_ptrtype const& dofOn() const { return M_dof; }
    dof_ptrtype const& dofOnOff() const { return M_dof; }

    element_type element( std::string const& name = "u" )
    {
        return element_type( *this, name );
    }

    element_ptrtype elementPtr( std::string const& name = "u" )
    {
        return std::make_shared<element_type>( *this, name );
    }

    class Element : public FunctionSpaceBase::ElementBase
    {
    public:
        using functionspace_type = ProductFunctionSpaces<ProductSpaces<SpaceList...>>;
        using product_element_type = typename functionspace_type::product_space_type::element_type;
        using value_type = typename functionspace_type::value_type;
        using size_type = typename functionspace_type::size_type;

        template<int I>
        struct sub_element
        {
            using type = std::remove_reference_t<decltype( std::declval<product_element_type&>()[hana::int_c<I>] )>;
            using ptrtype = std::shared_ptr<type>;
        };
        template<int I> using sub_element_type = typename sub_element<I>::type;
        template<int I> using sub_element_ptrtype = typename sub_element<I>::ptrtype;

        Element() = default;
        Element( functionspace_type& space, std::string name = "u" )
            :
            M_space( &space ),
            M_productElement( space.productSpace().element() ),
            M_name( std::move( name ) )
        {
        }

        size_type nDof() const { return this->functionSpace().nDof(); }
        size_type size() const { return this->nDof(); }
        std::string const& name() const { return M_name; }

        functionspace_type& functionSpace()
        {
            CHECK( M_space ) << "invalid product function-space element";
            return *M_space;
        }
        functionspace_type const& functionSpace() const
        {
            CHECK( M_space ) << "invalid product function-space element";
            return *M_space;
        }

        functionspace_vector_type const& functionSpaces() const
        {
            return this->functionSpace().functionSpaces();
        }

        template<int I>
        decltype(auto) functionSpace() const
        {
            return this->functionSpace().template functionSpace<I>();
        }

        template<int I>
        decltype(auto) element()
        {
            return M_productElement[hana::int_c<I>];
        }

        template<int I>
        decltype(auto) element() const
        {
            return M_productElement[hana::int_c<I>];
        }

        template<typename N>
        decltype(auto) operator[]( N const& n )
        {
            return M_productElement[n];
        }

        template<typename N>
        decltype(auto) operator[]( N const& n ) const
        {
            return M_productElement[n];
        }

        product_element_type& productElement() { return M_productElement; }
        product_element_type const& productElement() const { return M_productElement; }

    private:
        functionspace_type* M_space = nullptr;
        product_element_type M_productElement;
        std::string M_name;
    };

private:
    template<typename First, typename... Rest>
    static worldcomm_ptr_t firstWorldComm( First const& first, Rest const&... )
    {
        return first->worldCommPtr();
    }

    dof_ptrtype makeDofComposite() const
    {
        std::vector<datamap_ptrtype> subdm;
        subdm.reserve( nSpaces );
        boost::fusion::for_each( M_functionspaces,
                                 [&subdm]( auto const& space )
                                 {
                                     subdm.push_back( space->mapPtr() );
                                 } );
        return std::make_shared<dof_type>( subdm, this->worldCommPtr() );
    }

    template<typename F>
    size_type sumSpaces( F&& f ) const
    {
        size_type n = 0;
        boost::fusion::for_each( M_functionspaces,
                                 [&n,&f]( auto const& space )
                                 {
                                     n += f( space );
                                 } );
        return n;
    }

    template<typename F>
    size_type sumSpacesUntil( size_type i, F&& f ) const
    {
        CHECK( i <= nSpaces ) << "invalid block index " << i << " for product function space with " << nSpaces << " spaces";
        size_type n = 0;
        size_type c = 0;
        boost::fusion::for_each( M_functionspaces,
                                 [&n,&c,i,&f]( auto const& space )
                                 {
                                     if ( c < i )
                                         n += f( space );
                                     ++c;
                                 } );
        return n;
    }

private:
    product_space_type M_productSpaces;
    functionspace_vector_type M_functionspaces;
    worldscomm_ptr_t M_worldsComm;
    dof_ptrtype M_dof;
};

template<typename... SpaceList>
using product_function_spaces_t = ProductFunctionSpaces<ProductSpaces<SpaceList...>>;

template<typename... SpaceList>
using product_function_spaces_ptr_t = std::shared_ptr<product_function_spaces_t<SpaceList...>>;

template<typename... SpaceList>
product_function_spaces_t<SpaceList...>
productFunctionSpaces( SpaceList... spaces )
{
    return product_function_spaces_t<SpaceList...>( spaces... );
}

template<typename... SpaceList>
product_function_spaces_t<SpaceList...>
productFunctionSpaces( ProductSpaces<SpaceList...> const& spaces )
{
    return hana::unpack( spaces.tupleSpaces(),
                         []( auto const&... space )
                         {
                             return productFunctionSpaces( space... );
                         } );
}

template<typename... SpaceList>
product_function_spaces_ptr_t<SpaceList...>
productFunctionSpacesPtr( SpaceList... spaces )
{
    return std::make_shared<product_function_spaces_t<SpaceList...>>( spaces... );
}

template<typename... SpaceList>
product_function_spaces_ptr_t<SpaceList...>
productFunctionSpacesPtr( ProductSpaces<SpaceList...> const& spaces )
{
    return hana::unpack( spaces.tupleSpaces(),
                         []( auto const&... space )
                         {
                             return productFunctionSpacesPtr( space... );
                         } );
}

} // namespace Feel

#endif /* FEELPP_FEELDISCR_PRODUCTFUNCTIONSPACES_HPP */
