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
class ViewFactorBase
{
public:
    using value_type = double;
    ViewFactorBase() = default;
    ViewFactorBase( std::vector<std::string> const& list_of_bdys ) 
        : 
        list_of_bdys_( list_of_bdys ),
        vf_( list_of_bdys.size(), list_of_bdys.size() ),
        areas_( list_of_bdys.size() )
        {}
    ViewFactorBase( const ViewFactorBase& ) = default;
    ViewFactorBase( ViewFactorBase&& ) = default;
    ViewFactorBase& operator=( const ViewFactorBase& ) = default;
    ViewFactorBase& operator=( ViewFactorBase&& ) = default;
    virtual ~ViewFactorBase() = default;
    virtual void init( std::vector<std::string> const& list_of_bdys ) { list_of_bdys_ = list_of_bdys; }

    /// this function computes the deviation from reciprocity defined as Fij - Aj/Ai * Fji
    value_type devReciprocity( unsigned int i, unsigned int j ) const;

    /// this function computes the maximum absolute value of the deviation from reciprocity
    Real maxDevReciprocity() const;

  private:

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

};

} // namespace Feel
