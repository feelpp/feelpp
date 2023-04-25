/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4 
 */

#ifndef FEELPP_TOOLBOXES_MODELCORE_CACHEDMODELFIELDS_H
#define FEELPP_TOOLBOXES_MODELCORE_CACHEDMODELFIELDS_H 1

#include <memory>
#include <functional>
#include <type_traits>

namespace Feel
{
namespace FeelModels
{

template <typename FieldType>
struct ReturnValueUpdatePolicy
{
    typedef FieldType field_type;
    typedef std::shared_ptr<field_type> field_ptrtype;
    using update_function_type = std::function<field_type()>;

    template<typename UpdateFunctionType>
    ReturnValueUpdatePolicy( UpdateFunctionType && f ) :
        M_updateFunction( std::forward<UpdateFunctionType>(f) )
    {}

    void update( field_ptrtype & fieldPtr )
    {
        if( fieldPtr ) // if field ptr is not empty, we update the pointed value
            *fieldPtr = M_updateFunction();
        else // otherwise we reset the pointer
            fieldPtr.reset( new field_type( M_updateFunction() ) );
    }

protected:
    update_function_type M_updateFunction;
};

template <typename FieldType>
struct InplaceUpdatePolicy
{
    typedef FieldType field_type;
    typedef std::shared_ptr<field_type> field_ptrtype;
    using update_function_type = std::function<void(field_ptrtype &)>;

    template<typename UpdateFunctionType>
    InplaceUpdatePolicy( UpdateFunctionType && f ) :
        M_updateFunction( std::forward<UpdateFunctionType>(f) )
    {}

    void update( field_ptrtype & fieldPtr )
    {
        M_updateFunction( fieldPtr );
    }

protected:
    update_function_type M_updateFunction;
};

template <typename FieldType, template<typename> class UpdatePolicy = ReturnValueUpdatePolicy>
struct CachedModelField : UpdatePolicy<FieldType>
{
    typedef CachedModelField<FieldType, UpdatePolicy> self_type;
    typedef UpdatePolicy<FieldType> update_policy;

    typedef FieldType field_type;
    typedef std::shared_ptr<field_type> field_ptrtype;
    typedef typename field_type::functionspace_type functionspace_type;
    typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;

    template<typename UpdateFunctionType>
    CachedModelField( UpdateFunctionType && f ) :
        update_policy( std::forward<UpdateFunctionType>(f) ),
        M_fieldPtr( new field_type ),
        M_doUpdateField( true )
    {}
    template<typename UpdateFunctionType>
    CachedModelField( field_ptrtype const& fieldPtr, UpdateFunctionType && f ) :
        update_policy( std::forward<UpdateFunctionType>(f) ),
        M_fieldPtr( fieldPtr ),
        M_doUpdateField( true )
    {}

    CachedModelField( CachedModelField const& ) = default;
    CachedModelField( CachedModelField && ) = default;

    field_ptrtype fieldPtr( bool update = true ) const 
    {
        if( update )
            const_cast<self_type*>(this)->update();
        return M_fieldPtr;
    }
    field_type const& field( bool update = true ) const 
    {
        if( update )
            const_cast<self_type*>(this)->update();
        return *M_fieldPtr;
    }

    void setDoUpdate( bool doUpdate ) { M_doUpdateField = doUpdate; }

    void update()
    {
        if( M_doUpdateField )
        {
            update_policy::update( M_fieldPtr );
        }

        M_doUpdateField = false;
    }

private:
    field_ptrtype M_fieldPtr;
    bool M_doUpdateField {true};
};

template<typename T>
struct is_cached_model_field : std::false_type {};

template<typename T, template<typename> class U>
struct is_cached_model_field<CachedModelField<T,U>> : std::true_type {};

template<typename T>
inline constexpr bool is_cached_model_field_v = is_cached_model_field<T>::value;

template<typename T>
struct raw_field_type { typedef T type; };
template<typename T, template<typename> class U>
struct raw_field_type<CachedModelField<T,U>> : raw_field_type<T> {};
template<typename T>
struct raw_field_type<std::shared_ptr<T>>: raw_field_type<T> {};

template<typename T>
using raw_field_t = typename raw_field_type<T>::type;


} // namespace FeelModels
} // namespace Feel

#endif
