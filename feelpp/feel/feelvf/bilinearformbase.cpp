/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2017 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#define FEELPP_BILINEARFORMBASE_NOEXTERN 1
#include <feel/feelvf/bilinearformbase.hpp>



namespace Feel {




template<typename T>
BilinearFormBase<T>::BilinearFormBase( BilinearFormBase const& __vf )
    :
    super( __vf ),
    M_name( __vf.M_name ),
    M_pattern( __vf.M_pattern ),
    M_matrix( __vf.M_matrix ),
    M_lb( __vf.M_lb ),
    M_row_startInMatrix( __vf.M_row_startInMatrix ),
    M_col_startInMatrix( __vf.M_col_startInMatrix ),
    M_do_build( __vf.M_do_build ),
    M_do_threshold( __vf.M_do_threshold ),
    M_threshold( __vf.M_threshold ),
    M_dirichletState( __vf.M_dirichletState ? __vf.M_dirichletState :
                                            std::make_shared<deferred_dirichlet_state_type>() ),
    b_mutex()
{
    auto dmTest = M_matrix->mapRowPtr();
    auto dmTrial = M_matrix->mapColPtr();
    this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
    this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );
}


template<typename T>
BilinearFormBase<T>&
BilinearFormBase<T>::operator=( BilinearFormBase const& form )
{
    if ( this != &form )
    {
        this->invalidateMaterializedDeferredDirichlet();
        super::operator=( form );
        M_name = form.M_name;
        M_pattern = form.M_pattern;
        M_matrix->zero();
        M_matrix->addMatrix( 1.0, form.M_matrix );
        M_row_startInMatrix = form.M_row_startInMatrix;
        M_col_startInMatrix = form.M_col_startInMatrix;
        M_lb = form.M_lb;
        if ( !M_dirichletState )
            M_dirichletState = std::make_shared<deferred_dirichlet_state_type>();
        *M_dirichletState = form.M_dirichletState ? *form.M_dirichletState : deferred_dirichlet_state_type{};
        auto dmTest = M_matrix->mapRowPtr();
        auto dmTrial = M_matrix->mapColPtr();
        this->setDofIdToContainerIdTest( dmTest->dofIdToContainerId( M_row_startInMatrix ) );
        this->setDofIdToContainerIdTrial( dmTrial->dofIdToContainerId( M_col_startInMatrix ) );
    }

    return *this;
}
template<typename T>
void
BilinearFormBase<T>::addMatrix( int* rows, int nrows,
                                int* cols, int ncols,
                                value_type* data,
                                size_type K,
                                size_type K2)
{
    std::lock_guard<std::mutex> guard(b_mutex);
    this->invalidateMaterializedDeferredDirichlet();
    M_matrix->addMatrix( rows, nrows, cols, ncols, data, K, K2 );
}

template<typename T>
void
BilinearFormBase<T>::zeroRows( std::vector<int> const& __dofs,
                               Vector<value_type> const&__values,
                               Vector<value_type>& rhs,
                               Feel::Context const& on_context,
                               double value_on_diagonal )
{
    this->invalidateMaterializedDeferredDirichlet();
    M_matrix->zeroRows( __dofs, __values, rhs, on_context, value_on_diagonal );
}


template class BilinearFormBase<double>;
} // Feel
