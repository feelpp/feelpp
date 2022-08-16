/**
 * @file viewfactorbase.hpp
 * @author Christophe Prud'homme
 * @brief
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#pragma once

#include <feel/feelalg/glas.hpp>

namespace Feel
{
/**
 * @brief base class for View Factor computation
 * 
 * 
 */
template <typename MeshType>
class ViewFactorBase
{
public:
    using value_type = double;
    using mesh_t = MeshType;
    using mesh_ptrtype = std::shared_ptr<mesh_t>;
    using mesh_ptr_t = mesh_ptrtype;
    using mesh_const_ptrtype = std::shared_ptr<mesh_t const>;

    ViewFactorBase() = default;
    ViewFactorBase( mesh_ptr_t const& mesh, nl::json const& specs )
        : mesh_( mesh ),
          j_( specs ),
          list_of_bdys_( specs["markers"].get<std::vector<std::string>>() ),
          vf_( list_of_bdys_.size(), list_of_bdys_.size() ),
          areas_( list_of_bdys_.size() )
    {}
    ViewFactorBase( const ViewFactorBase& ) = default;
    ViewFactorBase( ViewFactorBase&& ) = default;
    ViewFactorBase& operator=( const ViewFactorBase& ) = default;
    ViewFactorBase& operator=( ViewFactorBase&& ) = default;
    virtual ~ViewFactorBase() = default;
//    virtual void init( std::vector<std::string> const& list_of_bdys ) { list_of_bdys_ = list_of_bdys; }

    /// this function computes the deviation from reciprocity defined as Fij - Aj/Ai * Fji
    value_type devReciprocity( unsigned int i, unsigned int j ) const;

    /// this function computes the maximum absolute value of the deviation from reciprocity
    Real maxDevReciprocity() const;

    /**
     * @brief compute the view factor matrix
     * 
     */
    virtual void compute() = 0;

    /**
     * @brief get the view factor matrix
     * 
     */
    eigen_matrix_xx_type<value_type> const& viewFactors() const { return vf_; }

    /**
     * @brief get the area matrix
     * 
     */
    eigen_vector_x_col_type<value_type> const& areas() const { return areas_; }
    
  protected:

    /// the mesh
    mesh_ptr_t mesh_;

    /// the json object
    nl::json j_;

    //! list of boundary markers where the view factor is computed
    std::vector<std::string> list_of_bdys_;

    ///! view factors
    eigen_matrix_xx_type<value_type> vf_;

    //! areas of the markedfaces
    eigen_vector_x_col_type<value_type> areas_;

    //! tolerance
    double vf_tol_ = 1e-6;

    //! normalize view factors so that vf_.array().rowwise().sum() == 1
    bool vf_normalize_ = true;

    inline static constexpr value_type exponent_ = (mesh_t::nDim==2)?1:2;
    inline static constexpr value_type divisor_ = (mesh_t::nDim==2)?2:M_PI;
};

template <typename MeshType>
typename ViewFactorBase<MeshType>::value_type
ViewFactorBase<MeshType>::devReciprocity( unsigned int i, unsigned int j ) const
{
    return vf_( i, j ) - areas_( j ) / areas_( i ) * vf_( j, i );
}
template <typename MeshType>
typename ViewFactorBase<MeshType>::value_type
ViewFactorBase<MeshType>::maxDevReciprocity() const
{
    return 1.;//vf_.array().rowwise()
}
} // namespace Feel
