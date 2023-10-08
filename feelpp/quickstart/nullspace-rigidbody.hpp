//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 22 Jul 2019
//! @copyright 2019 Feel++ Consortium
//!
#ifndef FEELPP_NULLSPACE_RIGIDBODY_HPP
#define FEELPP_NULLSPACE_RIGIDBODY_HPP 1

namespace Feel
{
template <typename SpaceType>
NullSpace<double> qsNullSpace( SpaceType const& space )
{
    if constexpr( decay_type<SpaceType>::nRealDim == 2 )
    {
        auto mode1 = space->element( oneX() );
        auto mode2 = space->element( oneY() );
        auto mode3 = space->element( vec(Py(),-Px()) );
        NullSpace<double> userNullSpace( { mode1,mode2,mode3 } );
        return userNullSpace;
    }
    else if constexpr ( decay_type<SpaceType>::nRealDim == 3 )
    {
        auto mode1 = space->element( oneX() );
        auto mode2 = space->element( oneY() );
        auto mode3 = space->element( oneZ() );
        auto mode4 = space->element( vec(Py(),-Px(),cst(0.)) );
        auto mode5 = space->element( vec(-Pz(),cst(0.),Px()) );
        auto mode6 = space->element( vec(cst(0.),Pz(),-Py()) );
        NullSpace<double> userNullSpace( { mode1,mode2,mode3,mode4,mode5,mode6 } );
        return userNullSpace;
    }
}

/**
 * @brief compute the rotation matrix for a rigid body
 * 
 * @param rigidRotationAngles vector of angles in Rad
 * @return eigen_matrix_type<RealDim, RealDim> 
 */
template<int RealDim>
eigen_matrix_type<RealDim, RealDim> 
rigidRotationMatrix( eigen_matrix_type<RealDim, 1> const& rigidRotationAngles )
{
    constexpr int nRealDim = RealDim;
    eigen_matrix_type<nRealDim, nRealDim> res;
    if constexpr ( nRealDim == 2 )
    {
        double angle = rigidRotationAngles(0,0);
        res <<   std::cos(angle), -std::sin(angle),
            /**/ std::sin(angle),  std::cos(angle);
    }
    else
    {
        double angleZ = rigidRotationAngles(2);
        double angleY = rigidRotationAngles(1);
        double angleX = rigidRotationAngles(0);
        eigen_matrix_type<3, 3> rotMatZ,rotMatY,rotMatX;
        rotMatZ << std::cos(angleZ), -std::sin(angleZ), 0,
            /**/   std::sin(angleZ),  std::cos(angleZ), 0,
            /**/                  0,                 0, 1;
        rotMatY << std::cos(angleY), 0, std::sin(angleY),
            /**/                  0, 1,                0,
            /**/  -std::sin(angleY), 0, std::cos(angleY);
        rotMatX << 1,                0,                 0,
            /**/   0, std::cos(angleX), -std::sin(angleX),
            /**/   0, std::sin(angleX),  std::cos(angleX);
        res = rotMatZ*rotMatY*rotMatX;
    }
    return res;
}

eigen_matrix_type<2, 2> 
rigidRotationMatrix( eigen_matrix_type<1, 1> const& r )
{
    return rigidRotationMatrix<2>( eigen_matrix_type<2, 1>{r(0,0),r(0,0)} );
}
eigen_matrix_type<2, 2> 
rigidRotationMatrix( double const& r )
{
    return rigidRotationMatrix<2>( eigen_matrix_type<2, 1>{r,r} );
}
}

#endif

