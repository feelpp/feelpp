/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes
       Date: 2014-11-19

  Copyright (C) 2012 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file exprtensorbase.hpp
   \author Vincent Chabannes
   \date 2013-11-19
 */
#ifndef __FEELPP_MODELS_VF_TENSORBASE_H
#define __FEELPP_MODELS_VF_TENSORBASE_H 1

namespace Feel
{
namespace FeelModels
{
    enum ExprApplyType { EVAL=0,JACOBIAN=1 };

    template<typename Geo_t, typename Basis_i_t, typename Basis_j_t, typename ShapeType, typename ValueType>
    struct tensorBase
    {
    public :
        typedef ValueType value_type;
        typedef ShapeType shape_type;
        typedef typename mpl::if_<fusion::result_of::has_key<Geo_t, Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type key_type;
        typedef typename fusion::result_of::value_at_key<Geo_t,key_type>::type::element_type gmc_type;
        typedef std::shared_ptr<gmc_type> gmc_ptrtype;
        typedef typename gmc_type::gm_type gm_type;

        typedef Basis_i_t map_basis_fec_test_type;
        typedef Basis_j_t map_basis_fec_trial_type;
        //typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_test_type,vf::detail::gmc<0> >,
        //                        mpl::identity<vf::detail::gmc<0> >,
        //                         mpl::identity<vf::detail::gmc<1> > >::type::type key_fec_test_type;
        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_test_type,Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type basis_fec_test_key_type;

        typedef typename mpl::if_<fusion::result_of::has_key<map_basis_fec_trial_type,Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<0> >,
                                  mpl::identity<Feel::vf::detail::gmc<1> > >::type::type basis_fec_trial_key_type;

        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type  basis_fec_test_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_test_type, basis_fec_test_key_type>::type::element_type* basis_fec_test_ptrtype;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type  basis_fec_trial_type;
        typedef typename fusion::result_of::value_at_key<map_basis_fec_trial_type,basis_fec_trial_key_type>::type::element_type* basis_fec_trial_ptrtype;

        typedef Eigen::Matrix<value_type,shape_type::M,shape_type::N> matrix_shape_type;
        typedef boost::multi_array<matrix_shape_type,1> array_shape_type;
        typedef Eigen::Matrix<matrix_shape_type,Eigen::Dynamic,1> new_array_shape_type;
        
        using ret_type = Eigen::Map<matrix_shape_type const>;
        
        // shapes used
        typedef Shape<shape_type::nDim, Scalar, false, false> shape_scalar;
        //typedef Eigen::Matrix<value_type,shape_scalar::M,shape_scalar::N> loc_scalar_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_scalar::M,shape_scalar::N>> loc_scalar_type;
        typedef boost::multi_array<loc_scalar_type,1> array_scalar_type;
        typedef Eigen::Matrix<loc_scalar_type,Eigen::Dynamic,1> new_array_scalar_type;
        
        typedef Shape<shape_type::nDim, Vectorial, false, false> shape_vectorial;
        //typedef Eigen::Matrix<value_type,shape_vectorial::M,shape_vectorial::N> loc_vectorial_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_vectorial::M,shape_vectorial::N>> loc_vectorial_type;
        typedef boost::multi_array<loc_vectorial_type,1> array_vectorial_type;
        typedef Eigen::Matrix<loc_vectorial_type,Eigen::Dynamic,1> new_array_vectorial_type;
        
        typedef Shape<shape_type::nDim, Vectorial, true, false> shape_vectorial_transpose;
        //typedef Eigen::Matrix<value_type,shape_vectorial_transpose::M,shape_vectorial_transpose::N> loc_vectorial_transpose_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_vectorial_transpose::M,shape_vectorial_transpose::N>> loc_vectorial_transpose_type;
        typedef boost::multi_array<loc_vectorial_transpose_type,1> array_vectorial_transpose_type;
        typedef Eigen::Matrix<loc_vectorial_transpose_type,Eigen::Dynamic,1> new_array_vectorial_transpose_type;
        
        typedef Shape<shape_type::nDim, Tensor2, false, false> shape_tensor2;
        //typedef Eigen::Matrix<value_type,shape_tensor2::M,shape_tensor2::N> loc_tensor2_type;
        typedef Eigen::TensorFixedSize<value_type,Eigen::Sizes<shape_tensor2::M,shape_tensor2::N>> loc_tensor2_type;
        typedef boost::multi_array<loc_tensor2_type,1> array_tensor2_type;
        typedef Eigen::Matrix<loc_tensor2_type,Eigen::Dynamic,1> new_array_tensor2_type;

        typedef Eigen::Matrix<value_type,shape_tensor2::M,shape_tensor2::N> loc_matrix_tensor2_type;
        typedef boost::multi_array<loc_matrix_tensor2_type,1> array_matrix_tensor2_type;

        typedef boost::multi_array<value_type,1> array_value_type;


        tensorBase( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu )
            :
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_fecTest( fusion::at_key<basis_fec_test_key_type>( fev ).get() ),
            M_fecTrial( fusion::at_key<basis_fec_trial_key_type>( feu ).get() ),
            M_locMatrixShape( matrix_shape_type::Zero() ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocTensor2()
            {
                M_zeroLocScalar.setZero();
                M_zeroLocVectorial.setZero();
                M_zeroLocTensor2.setZero();
            }
        tensorBase( Geo_t const& geom, Basis_i_t const& fev )
            :
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_fecTest( fusion::at_key<basis_fec_test_key_type>( fev ).get() ),
            M_locMatrixShape( matrix_shape_type::Zero() ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocTensor2()
            {
                M_zeroLocScalar.setZero();
                M_zeroLocVectorial.setZero();
                M_zeroLocTensor2.setZero();
            }
        tensorBase( Geo_t const& geom )
            :
            M_geot( fusion::at_key<key_type>( geom ) ),
            M_locMatrixShape( matrix_shape_type::Zero() ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocTensor2()
            {
                M_zeroLocScalar.setZero();
                M_zeroLocVectorial.setZero();
                M_zeroLocTensor2.setZero();
            }
        tensorBase( tensorBase const& t )
            :
            M_geot( t.M_geot ),
            M_fecTest( t.M_fecTest ),
            M_fecTrial( t.M_fecTrial ),
            M_zeroLocScalar(),
            M_zeroLocVectorial(),
            M_zeroLocTensor2()
        {}
        virtual ~tensorBase() {}

        void setGmc( Geo_t const& geom ) { M_geot = fusion::at_key<key_type>( geom ); }

        gmc_ptrtype const& gmc() const { return M_geot; }
        basis_fec_test_ptrtype const& fecTest() const { return M_fecTest; }
        basis_fec_trial_ptrtype const& fecTrial() const { return M_fecTrial; }
        matrix_shape_type /*const*/& locMatrixShape() const { return M_locMatrixShape; }

        virtual void update( Geo_t const& geom ) { CHECK( false ) << "should be override"; };
        virtual void update( Geo_t const& geom, uint16_type face ) { CHECK( false ) << "should be override"; };
        virtual void update( Geo_t const& geom, Basis_i_t const& fev, Basis_j_t const& feu ) { CHECK( false ) << "should be override"; }
        virtual void update( Geo_t const& geom, Basis_i_t const& fev ) { CHECK( false ) << "should be override"; }



        virtual
        ret_type
        evalijq( uint16_type i, uint16_type j, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return ret_type(M_locMatrixShape.data());
        }

        virtual
        value_type
        evalijq( uint16_type i, uint16_type j, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return value_type(0);
        }

        virtual
        value_type
        evaliq( uint16_type i, uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return value_type(0);
        }
        virtual
        ret_type
        evaliq( uint16_type i, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return ret_type(M_locMatrixShape.data());
        }

        virtual
        value_type
        evalq( uint16_type c1, uint16_type c2, uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return value_type(0);
        }
        virtual
        ret_type
        evalq( uint16_type q ) const
        {
            CHECK( false ) << "not allow\n";
            return ret_type(M_locMatrixShape.data());
        }

    private:
        gmc_ptrtype M_geot;
        basis_fec_test_ptrtype M_fecTest;
        basis_fec_trial_ptrtype M_fecTrial;

    protected :
        mutable matrix_shape_type M_locMatrixShape;
        loc_scalar_type M_zeroLocScalar;
        loc_vectorial_type M_zeroLocVectorial;
        loc_tensor2_type M_zeroLocTensor2;

    };

} // namespace FeelModels
} // namespace Feel

#endif /* __FEELPP_MODELS_VF_TENSORBASE_H */
