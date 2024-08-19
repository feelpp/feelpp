/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 01 Jun 2016

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
#ifndef FEELPP_MESHSTRUCTURED_HPP
#define FEELPP_MESHSTRUCTURED_HPP 1

#include <feel/feeldiscr/mesh.hpp>

namespace Feel
{

template <typename T>
using holo3_image = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


/**
 * MeshStructuredSetup class define :
 * - the geometry of structured mesh
 * - the discretisation with points and cells (segement,quad or hexa with respect to the dim)
 * - the partitioning strategy
 */
template <int Dim, int RealDim = Dim, typename T = double, typename IndexT = uint32_type>
class MeshStructuredSetup
{
    using self_type = MeshStructuredSetup<Dim,RealDim,T,IndexT>;
public :
    static constexpr uint16_type nDim = Dim;
    static constexpr uint16_type nRealDim = RealDim;
    using value_type = T;
    using index_type = IndexT;
    using node_eigen_type = eigen_vector_type<nRealDim,value_type>;

    enum class DivisionEdgeType { Equidistributed };

    static constexpr int nPointInCell = Dim==1? 2 : ( Dim==2? 4 : 8 );

    static std::shared_ptr<self_type> New( nl::json const& jDesc );

    MeshStructuredSetup( std::string const& shape )
        :
        M_shape( shape )
        {}

    virtual ~MeshStructuredSetup() {}

    std::vector<index_type> const& nPointByAxis() const { return M_nPointByAxis; }

    virtual bool isCartesian() const { return false; }

    virtual node_type pointCoordinates( std::array<index_type,nDim> const& indexes ) const = 0;

    index_type pointId( std::array<index_type,nDim> const& indexes ) const
        {
            if constexpr ( nDim == 1 )
            {
                return indexes[0];
            }
            else if constexpr ( nDim == 2 )
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                return M_nPointByAxis[1]*i + j;
            }
            else
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                index_type k = indexes[2];
                return M_nPointByAxis[1]*i + j + M_nPointByAxis[0]*M_nPointByAxis[1]*k;
            }
        }
    index_type cellId( std::array<index_type,nDim> const& indexes ) const
        {
            if constexpr ( nDim == 1 )
            {
                return indexes[0];
            }
            else if constexpr ( nDim == 2 )
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                return (M_nPointByAxis[1]-1)*i + j;
            }
            else
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                index_type k = indexes[2];
                return (M_nPointByAxis[1]-1)*i + j + (M_nPointByAxis[0]-1)*(M_nPointByAxis[1]-1)*k;
            }
        }

    std::array<index_type,nPointInCell> pointIdsFromCellId( std::array<index_type,nDim> const& indexes ) const
        {
            if constexpr ( nDim == 1 )
            {
                return { this->pointId( {indexes[0]} ), this->pointId( {indexes[0]+1} ) };
            }
            else if constexpr ( nDim == 2 )
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                return {
                    this->pointId( {i+1,j} ),    // 0
                    this->pointId( {i+1,j+1} ),  // 1
                    this->pointId( {i,j+1} ),    // 2
                    this->pointId( {i,j} )       // 3
                        };
            }
            else
            {
                index_type i = indexes[0];
                index_type j = indexes[1];
                index_type k = indexes[2];

                return {
                    this->pointId( {i+1,j,k} ),    // 0
                        this->pointId( {i+1,j+1,k} ),  // 1
                        this->pointId( {i,j+1,k} ),    // 2
                        this->pointId( {i,j,k} ),       // 3

                        this->pointId( {i+1,j,k+1} ),    // 4
                        this->pointId( {i+1,j+1,k+1} ),  // 5
                        this->pointId( {i,j+1,k+1} ),    // 6
                        this->pointId( {i,j,k+1} )       // 7
                        };
            }

        }
protected:
    void updateDiscretisation( nl::json const& jDesc = {} )
        {
            this->M_nPointByAxis.resize(nDim,10);
            if ( jDesc.contains( "Discretisation" ) )
            {
                auto const& jDescDiscr = jDesc.at( "Discretisation" );
                if ( jDescDiscr.contains( "n_points" ) )
                {
                    auto jDescDiscr_npts = jDescDiscr.at( "n_points" );
                    if ( jDescDiscr_npts.is_array() )
                    {
                        int d=0;
                        for ( auto n : jDescDiscr_npts.template get<std::vector<int>>() )
                        {
                            if ( d <nDim )
                                this->M_nPointByAxis[d] = n;
                            ++d;
                        }
                    }
                    else if ( jDescDiscr_npts.is_number_integer() )
                    {
                        int n = jDescDiscr_npts.template get<int>();
                        for (int d=0;d<nDim;++d)
                            this->M_nPointByAxis[d] = n;
                    }
                }

                M_discretisationEdges.resize( M_geometricEdges.size(), std::make_tuple( DivisionEdgeType::Equidistributed, 0. ) );
                if ( jDescDiscr.contains( "division" ) )
                {
                    auto jDescDiscr_div = jDescDiscr.at( "division" );
                    if ( jDescDiscr_div.is_string() && false )
                    {
                        std::string theDivision = jDescDiscr_div.template get<std::string>();
                        // TODO
                    }
                    else
                    {
                        for ( int k=0;k < M_discretisationEdges.size();++k)
                        {
                            std::get<0>( M_discretisationEdges[k] ) = DivisionEdgeType::Equidistributed;
                        }
                    }
                }
            }

            for ( int k=0;k < M_discretisationEdges.size();++k)
            {
                int axis = std::get<0>( M_geometricEdges[k] );
                auto const& pt0 = M_geometricPoints[ std::get<1>( M_geometricEdges[k] ) ];
                auto const& pt1 = M_geometricPoints[ std::get<2>( M_geometricEdges[k] ) ];
                value_type edgeLength = (pt0-pt1).norm();
                if ( std::get<0>( M_discretisationEdges[k] ) == DivisionEdgeType::Equidistributed )
                {
                    std::get<1>( M_discretisationEdges[k] ) = edgeLength/(this->M_nPointByAxis[axis]-1);
                }
            }
        }
protected:
    std::string M_shape;

    std::vector<node_eigen_type> M_geometricPoints;
    std::vector<std::tuple<int,int,int>> M_geometricEdges;// (axis, ptId_0, ptId_1)

    std::vector<index_type> M_nPointByAxis;
    std::vector<std::tuple<DivisionEdgeType,value_type>> M_discretisationEdges; //( type, spacing param)

};

#if 0
template <int RealDim = 2, typename T = double, typename IndexT = uint32_type>
class MeshStructuredSetupQuadrangle : public MeshStructuredSetup<2,RealDim,T,IndexT>
{
};

template <typename T = double, typename IndexT = uint32_type>
class MeshStructuredSetupQuadrilaterallyFacedHexahedron : public MeshStructuredSetup<3,3,T,IndexT>
{
};
#endif

template <int Dim, int RealDim = Dim, typename T = double, typename IndexT = uint32_type>
class MeshStructuredSetupHyperRectangle : public MeshStructuredSetup<Dim,RealDim,T,IndexT>
{
    using super_type = MeshStructuredSetup<Dim,RealDim,T,IndexT>;
    using self_type = MeshStructuredSetupHyperRectangle<Dim,RealDim,T,IndexT>;
public :
    static constexpr uint16_type nDim = super_type::nDim;
    static constexpr uint16_type nRealDim = super_type::nRealDim;
    using value_type = typename super_type::value_type;
    using index_type = typename super_type::index_type;

    using node_eigen_type = eigen_vector_type<nRealDim,value_type>;
    //static constexpr int nSide = Dim==1? 1 : ( Dim==2? 4 : 6 );
    static constexpr int nPoint = Dim==1? 2 : ( Dim==2? 4 : 8 );
    static constexpr int nEdge = Dim==1? 1 : ( Dim==2? 4 : 12 );

    static
    std::shared_ptr<self_type> NewUnitShape( std::array<index_type,nDim> const& nPts )
        {
            nl::json jDesc;
            jDesc.emplace( "Geometry", nl::json( {
                        { "shape", "hyperrectangle" }
                    } ) );
            jDesc.emplace( "Discretisation", nl::json( {
                        { "n_points", nPts },
                        { "division", "equidistributed" }
                    } ) );
            return std::make_shared<self_type>( jDesc );
        }


    MeshStructuredSetupHyperRectangle( nl::json const& jDesc )
        :
        super_type( "hyperrectangle" ),
        M_originPoint( node_eigen_type::Zero() ),
        M_axis( nDim,node_eigen_type::Zero() ),
        M_lengthByAxis( eigen_vector_type<nDim,value_type>::Ones() )
        {
            for (int d=0;d<nDim;++d)
                M_axis[d][d] = 1.0;

            if ( jDesc.contains( "Geometry" ) )
            {
                auto const& jDescGeo =  jDesc.at( "Geometry" );
                if ( jDescGeo.contains( "domain" ) )
                {
                    auto const& jDescGeoDomain = jDescGeo.at( "domain" );
                    if ( jDescGeoDomain.contains( "origin" ) )
                    {
                        auto vorigin = jDescGeoDomain.at( "origin" ).template get<std::vector<value_type>>();
                        CHECK( vorigin.size() == nRealDim ) << "invalid origin entry : bad vector size";
                        for ( int d=0;d<nRealDim;++d )
                            M_originPoint[d] = vorigin[d];
                    }
                    if ( jDescGeoDomain.contains( "length" ) )
                    {
                        auto const& jlength = jDescGeoDomain.at( "length" );
                        if ( jlength.is_array() )
                        {
                            auto vlength = jlength.template get<std::vector<value_type>>();
                            CHECK( vlength.size() == nDim ) << "invalid length entry : bad vector size";
                            for ( int d=0;d<nDim;++d )
                                M_lengthByAxis[d] = vlength[d];
                        }
                        else if ( jlength.is_number_float() )
                        {
                            value_type vlength = jlength.template get<value_type>();
                            for ( int d=0;d<nDim;++d )
                                M_lengthByAxis[d] = vlength;
                        }
                    }
                }
            }

            this->M_geometricPoints.resize( nPoint );
            this->M_geometricEdges.resize( nEdge );

            if constexpr ( nDim == 1 )
            {
                this->M_geometricPoints[0] = M_originPoint;
                this->M_geometricPoints[1] = M_originPoint + M_lengthByAxis[0]*M_axis[0];
                this->M_geometricEdges[0] = std::make_tuple( 0,0,1 );
            }
            else if constexpr ( nDim == 2 )
            {
                this->M_geometricPoints[0] = M_originPoint;
                this->M_geometricPoints[1] = M_originPoint + M_lengthByAxis[0]*M_axis[0];
                this->M_geometricPoints[2] = M_originPoint + M_lengthByAxis[0]*M_axis[0] + M_lengthByAxis[1]*M_axis[1];
                this->M_geometricPoints[3] = M_originPoint + M_lengthByAxis[1]*M_axis[1];

                this->M_geometricEdges[0] = std::make_tuple( 0,0,1 );
                this->M_geometricEdges[1] = std::make_tuple( 1,1,2 );
                this->M_geometricEdges[2] = std::make_tuple( 0,2,3 );
                this->M_geometricEdges[3] = std::make_tuple( 1,3,0 );
            }
            else
            {
                this->M_geometricPoints[0] = M_originPoint;
                this->M_geometricPoints[1] = M_originPoint + M_lengthByAxis[0]*M_axis[0];
                this->M_geometricPoints[2] = M_originPoint + M_lengthByAxis[0]*M_axis[0] + M_lengthByAxis[1]*M_axis[1];
                this->M_geometricPoints[3] = M_originPoint + M_lengthByAxis[1]*M_axis[1];
                // extrude point
                for (int k=0;k<4;++k)
                    this->M_geometricPoints[4+k] = this->M_geometricPoints[k] + M_lengthByAxis[2]*M_axis[2];

                this->M_geometricEdges[0] = std::make_tuple( 0,0,1 );
                this->M_geometricEdges[1] = std::make_tuple( 1,1,2 );
                this->M_geometricEdges[2] = std::make_tuple( 0,2,3 );
                this->M_geometricEdges[3] = std::make_tuple( 1,3,0 );

                this->M_geometricEdges[4] = std::make_tuple( 0,4,5 );
                this->M_geometricEdges[5] = std::make_tuple( 1,5,6 );
                this->M_geometricEdges[6] = std::make_tuple( 0,6,7 );
                this->M_geometricEdges[7] = std::make_tuple( 1,7,4 );

                for (int k=0;k<4;++k)
                    this->M_geometricEdges[8+k] = std::make_tuple( 2,0+k,4+k );
            }

            this->updateDiscretisation( jDesc );
        }

    bool isCartesian() const override
        {
            for ( auto const& [divisionType,param] : this->M_discretisationEdges)
                if ( divisionType != super_type::DivisionEdgeType::Equidistributed )
                    return false;
            return false;
        }

    node_type pointCoordinates( std::array<index_type,nDim> const& indexes ) const override
        {
            node_eigen_type coords = M_originPoint;
            // only good with all equidistributed
            if constexpr ( nDim <= 2 )
            {
                for (int d=0;d<nDim;++d )
                    coords += indexes[d]*std::get<1>(this->M_discretisationEdges[d])*M_axis[d];
            }
            else
            {
                coords += indexes[0]*std::get<1>(this->M_discretisationEdges[0])*M_axis[0];
                coords += indexes[1]*std::get<1>(this->M_discretisationEdges[1])*M_axis[1];
                coords += indexes[2]*std::get<1>(this->M_discretisationEdges[8])*M_axis[2];
            }

            node_type coordsRet( nRealDim );
            for ( int d=0;d<nRealDim;++d )
                coordsRet[d] = coords[d];
            return coordsRet;
        }
private:
    node_eigen_type M_originPoint;
    std::vector<node_eigen_type> M_axis;
    eigen_vector_type<nDim,value_type> M_lengthByAxis;
};


template <int Dim, int RealDim = Dim, typename T = double, typename IndexT = uint32_type>
class MeshStructuredSetupWithCoordinates : public MeshStructuredSetup<Dim,RealDim,T,IndexT>
{
    using super_type = MeshStructuredSetup<Dim,RealDim,T,IndexT>;
public :
    static constexpr uint16_type nDim = super_type::nDim;
    static constexpr uint16_type nRealDim = super_type::nRealDim;
    using value_type = typename super_type::value_type;
    using index_type = typename super_type::index_type;

    MeshStructuredSetupWithCoordinates( holo3_image<float> const& cx, holo3_image<float> const& cy )
        :
        super_type( "user-coordinates" ),
        M_cx( cx ), M_cy( cy )
        {
            this->M_nPointByAxis.resize(nDim,0);
#if 0
            this->M_nPointByAxis[0] = M_cx->cols();
            this->M_nPointByAxis[1] = M_cy->rows();
#else
            this->M_nPointByAxis[0] = M_cx->rows();
            this->M_nPointByAxis[1] = M_cx->cols();
            CHECK(  (M_cx->rows() == M_cy->rows()) &&  (M_cx->cols() == M_cy->cols()) ) << "invalid coords size";
#endif
        }

    node_type pointCoordinates( std::array<index_type,nDim> const& indexes ) const override
        {
            index_type i = indexes[0];
            index_type j = indexes[1];
            node_type coords( nRealDim );
#if 0
            // TO UNDERSTAND : why x and y are inversed?
            coords[0] = (*M_cy)( j, i );
            if constexpr ( nDim > 1 )
                  coords[1] = (*M_cx)( j, i );
#else
            coords[0] = (*M_cx)( i, j );
            if constexpr ( nRealDim > 1 )
                  coords[1] = (*M_cy)( i, j );
#endif
            return coords;
        }

private:
    std::optional<holo3_image<float>> M_cx; // X-coordinates for nodes
    std::optional<holo3_image<float>> M_cy; // Y-coordinates for nodes
};


template <int Dim, int RealDim, typename T, typename IndexT>
std::shared_ptr<MeshStructuredSetup<Dim,RealDim,T,IndexT>>
MeshStructuredSetup<Dim,RealDim,T,IndexT>::New( nl::json const& jDesc )
{
    return std::make_shared<MeshStructuredSetupHyperRectangle<Dim,RealDim,T,IndexT>>( jDesc );
}



//! Structured mesh class
//!
//! A structured mesh is such that points and elements can
//! be located with a simple index tuple (i,j,k)
//!
//! \code
//! TODO add code example here
//! \endcode
//!

struct MeshStructuredBase {};

template <typename GeoShape, typename T, typename IndexT>
class MeshStructured : public Mesh<GeoShape,T,0,IndexT,EnableSharedFromThis>,
                       public MeshStructuredBase//,
                       //public std::enable_shared_from_this<MeshStructured<GeoShape,T,IndexT>>
{

  public:
    using super = Mesh<GeoShape,T,0,IndexT>;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;
    using point_type = typename super::point_type;
    using element_type = typename super::element_type;
    using face_type = typename super::face_type;
    using node_type = typename super::node_type;

    static constexpr uint16_type nDim = super::nDim;
    static constexpr uint16_type nRealDim = super::nRealDim;

    using setup_type = MeshStructuredSetup<nDim,nRealDim,T,index_type>;

    MeshStructured( worldcomm_ptr_t const& wc  = Environment::worldCommPtr() ): super( wc ) { this->setStructureProperty( "00010" ); }
    MeshStructured( MeshStructured const& ) = delete;
    MeshStructured( MeshStructured&& ) = delete;
    MeshStructured& operator=( MeshStructured const& ) = delete;
    MeshStructured& operator=( MeshStructured&& ) = delete;

    //!
    //!
    //!
    MeshStructured( nl::json const& jDesc, worldcomm_ptr_t const& wc = Environment::worldCommPtr() )
        :
        MeshStructured( setup_type::New( jDesc ), wc )
        {}
    MeshStructured( std::array<index_type,nDim> const& indexes, worldcomm_ptr_t const& wc = Environment::worldCommPtr() )
        :
        MeshStructured( MeshStructuredSetupHyperRectangle<nDim,nRealDim,T,index_type>::NewUnitShape( indexes ), wc )
        {}
    MeshStructured( holo3_image<float> const& cx, holo3_image<float> const& cy,
                    worldcomm_ptr_t const& wc = Environment::worldCommPtr() )
        :
        MeshStructured( std::make_shared< MeshStructuredSetupWithCoordinates<nDim,nRealDim,T,index_type> >( cx,cy ), wc )
        {}

    MeshStructured( std::shared_ptr<setup_type> const& setup, worldcomm_ptr_t const& = Environment::worldCommPtr() );

    //! return current shared_ptr of type MeshBase
    std::shared_ptr<MeshBase<IndexT>> shared_from_this_meshbase() override { return std::dynamic_pointer_cast<MeshBase<IndexT>>( this->shared_from_this() );  }

    //! return current shared_ptr of type MeshBase
    std::shared_ptr<const MeshBase<IndexT>> shared_from_this_meshbase() const override { return std::dynamic_pointer_cast<const MeshBase<IndexT>>( this->shared_from_this() );  }


private:

    void generateStructuredMesh();

    // TO CHECK : maybe inline these methods
    void addStructuredPoint( std::array<index_type,nDim> const& indexes, rank_type partId, bool isGhost );
    std::pair<size_type,size_type> addStructuredElement( std::array<index_type,nDim> const& indexes, rank_type processId, rank_type partId,
                                                         std::vector<rank_type> const& neighborPartitionIds );
   void updateGhostCellInfoByUsingNonBlockingComm(
        std::unordered_map<size_type, size_type> const& idStructuredMeshToFeelMesh,
        std::unordered_map<size_type, boost::tuple<size_type, rank_type>> const& mapGhostElt );

private:
    std::shared_ptr<setup_type> M_setup;
};

/**
 * @brief trait type to detect a @p MeshStructured mesh
 * 
 * @tparam MeshT mesh type
 */
template<typename MeshT>
using is_mesh_structured = std::conditional_t<std::is_base_of_v<MeshStructuredBase, MeshT>, std::true_type, std::false_type>;

/**
 * @brief boolean to detect a @p MeshStructured mesh
 * 
 * @tparam MeshT mesh type
 */
template<typename MeshT>
inline constexpr bool is_mesh_structured_v = is_mesh_structured<MeshT>::value;

} // namespace Feel

#endif
