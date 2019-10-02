//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

using namespace Feel;

template<typename SpaceTypeConductor, typename SpaceTypeMag>
class BiotSavartBase
{
public:
    using self_type = BiotSavartBase<SpaceTypeConductor, SpaceTypeMag>;

    using space_conductor_ptrtype = SpaceTypeConductor;
    using space_conductor_type = typename space_conductor_ptrtype::element_type;
    using element_conductor_type = typename space_conductor_type::element_type;
    using space_mag_ptrtype = SpaceTypeMag;
    using space_mag_type = typename space_mag_ptrtype::element_type;
    using element_mag_type = typename space_mag_type::element_type;

    using dof_point_type = boost::tuple<node_type, size_type, uint16_type >;
    using dof_points_type = typename std::vector<dof_point_type>;

    BiotSavartBase(space_conductor_ptrtype const& XhC, space_mag_ptrtype const& XhM)
        : M_XhC(XhC),
          M_XhM(XhM)
        {
            setupCommunicators();
        }

    void setupCommunicators();
    element_mag_type computeMagneticField(element_conductor_type const& j);

private:
    space_conductor_ptrtype M_XhC;
    space_mag_ptrtype M_XhM;
    std::vector< mpi::communicator > M_commsC1M;
    std::map<int, dof_points_type> M_dofMgn;
};

template<typename SpaceTypeConductor, typename SpaceTypeMag>
void BiotSavartBase<SpaceTypeConductor, SpaceTypeMag>::setupCommunicators()
{
    M_dofMgn.clear();
    M_commsC1M.clear();
    std::vector<int> isIn( 2, 0 );
    if( M_XhC->nLocalDof()  > 0 ) // current proc has dofs in Omega_C
        isIn[0] = 1;
    if( M_XhM->nLocalDof() > 0 ) // current proc has dofs in Omega_M
        isIn[1] = 1;

    std::vector<int> isInGlob( 2*Environment::worldComm().size() );
    mpi::all_gather( Environment::worldComm(), isIn.data(), 2, isInGlob );

    std::vector<int> pC,pM,pC1M;
    for(int i=0; i<Environment::worldComm().size(); i++)
    {
        if( isInGlob[2*i] )
            pC.push_back( i );
        if( isInGlob[2*i + 1] )
            pM.push_back( i ); //Check if pM has at least one element for each procs
    }

    // commC = communicator with all procs which in pC
    auto groupC = Environment::worldComm().group().include(pC.begin(), pC.end());
    mpi::communicator commC(Environment::worldComm(), groupC);
    // commM = communicator with all procs which in pM
    auto groupM = Environment::worldComm().group().include(pM.begin(), pM.end());
    mpi::communicator commM(Environment::worldComm(), groupM);

    dof_points_type dofM; //dofs

    for( int i=0; i<pM.size(); i++)
    {
        pC1M.clear();
        pC1M.resize( pC.size() );
        pC1M = pC;
        auto find = std::find( pC1M.begin(), pC1M.end(), pM[i] );
        if( find != pC1M.end() )
            pC1M.erase(find);  //No double definition
        pC1M.insert(pC1M.begin(), pM[i] );

        auto groupC1M = Environment::worldComm().group().include(pC1M.begin(), pC1M.end());
        mpi::communicator commC1M(Environment::worldComm(), groupC1M);
        M_commsC1M.push_back( commC1M );
        M_dofMgn.insert(std::make_pair(i,dofM));
    }

    // dofM = map<int, dof_points_type>
    int dofSize;
    for(int i=0; i<M_commsC1M.size(); i++)
    {
        if( M_commsC1M[i] ) // Current communicator
        {
            if(M_commsC1M[i].rank() == 0 ) // Proc of rank 0 is the Omega_M part
            {
                for ( size_type dof_id = 0; dof_id < M_XhM->nLocalDofWithGhost() ; ++dof_id )
                {
                    auto dofpoint = M_XhM->dof()->dofPoint(dof_id);
                    M_dofMgn.at(i).push_back( dofpoint );
                }
                dofSize = M_dofMgn.at(i).size();
            }

            mpi::broadcast( M_commsC1M[i], dofSize, 0);
            M_dofMgn.at(i).resize( dofSize );
            mpi::broadcast( M_commsC1M[i], M_dofMgn.at(i).data(), dofSize, 0);
        }
    }
}

template<typename SpaceTypeConductor, typename SpaceTypeMag>
typename BiotSavartBase<SpaceTypeConductor, SpaceTypeMag>::element_mag_type
BiotSavartBase<SpaceTypeConductor, SpaceTypeMag>::computeMagneticField(element_conductor_type const& j)
{
    auto mesh = M_XhC->mesh();
    auto B = M_XhM->element();

    auto coeff = 1/(4*M_PI);
    auto mu0 = 4*M_PI*1e-7; //SI unit : H.m-1 = m.kg.s-2.A-2

    std::vector<double> intLocD, intSumD;
    for(int i = 0; i < M_dofMgn.size(); ++i)
    {
        if( M_commsC1M[i] )
        {
            int dofSize = M_dofMgn.at(i).size();
            intLocD.resize( dofSize, 0 );
            intSumD.resize( dofSize );

            if( M_XhC->nLocalDof() > 0 )
            {
                std::vector<Eigen::Matrix<double,3,1>> coords( dofSize/3 );
                for( int d = 0; d < dofSize; d+= 3 )
                {
                    auto dofCoord = M_dofMgn.at(i)[d].template get<0>();
                    Eigen::Matrix<double,3,1> coord;
                    coord << dofCoord[0], dofCoord[1], dofCoord[2];
                    coords[d/3] = coord;
                }

                std::vector<Eigen::Matrix<double,3,1> > mgnFields;
                auto dist = inner( _e1v-P(), _e1v-P(),
                                   mpl::int_<InnerProperties::IS_SAME|InnerProperties::SQRT>() );
                mgnFields = integrate(_range=M_XhC->dof()->meshSupport()->rangeElements(),
                                      _expr=-mu0*coeff*cross(idv(j),_e1v-P())/(dist*dist*dist)
                                      ).template evaluate(coords);

                for( int d = 0; d < dofSize; ++d )
                {
                    auto dofComp = M_dofMgn.at(i)[d].template get<2>();
                    intLocD[d] += mgnFields[d/3](dofComp,0);
                }
            }

            mpi::reduce(M_commsC1M[i], intLocD.data(), dofSize, intSumD.data(), std::plus<double>(), 0);

            if( M_commsC1M[i].rank() == 0 )
            {
                for(int d=0; d<dofSize; d++)
                    B.set(M_dofMgn.at(i)[d].template get<1>(), intSumD[d]);
            }
        }
        Environment::worldComm().barrier();
    }

    B.close();

    return B;
}
