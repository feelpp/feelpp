/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2011-07-17

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
 \file dofrelationshipmap.hpp
 \author Vincent Chabannes <vincent.chabannes@feelpp.org>
 \date 2011-07-17
 */

#ifndef FEELPP_MODELS_DOFRELATIONSHIPMAP_H
#define FEELPP_MODELS_DOFRELATIONSHIPMAP_H 1

#include <feel/feeldiscr/functionspace.hpp>

namespace Feel
{
namespace FeelModels
{

template< class SpaceType1,class SpaceType2 >
class DofRelationshipMap
{

public :

    typedef SpaceType1 functionspace1_type;
    typedef std::shared_ptr<functionspace1_type> functionspace1_ptrtype;
    typedef SpaceType2 functionspace2_type;
    typedef std::shared_ptr<functionspace2_type> functionspace2_ptrtype;

    typedef typename functionspace1_type::mesh_type mesh1_type;
    typedef typename functionspace2_type::mesh_type mesh2_type;

    static const uint16_type nDim = mesh1_type::nDim;
    static const bool is_simplex = mesh1_type::shape_type::is_simplex;

    static const uint16_type nDofPerVertex = functionspace1_type::fe_type::nDofPerVertex;
    static const uint16_type numVertices = mesh1_type::element_type::numVertices;
    static const uint16_type nDofPerEdge = functionspace1_type::fe_type::nDofPerEdge;
    static const uint16_type numEdges = mesh1_type::element_type::numEdges;
    static const uint16_type nDofPerFace = functionspace1_type::fe_type::nDofPerFace;
    static const uint16_type numGeometricFaces = mesh1_type::element_type::numGeometricFaces;
    static const uint16_type nDofPerVolume = functionspace1_type::fe_type::nDofPerVolume;
    static const uint16_type numVolumes = mesh1_type::element_type::numVolumes;

    DofRelationshipMap(functionspace1_ptrtype __Xh1,functionspace2_ptrtype __Xh2 )
        :
        M_Xh1(__Xh1),
        M_Xh2(__Xh2),
        M_dofRelMapRefToHo( M_Xh1->nLocalDof(), invalid_v<size_type> ),
        M_dofRelMapHoToRef( M_Xh2->nLocalDof(), invalid_v<size_type> )
        {
            buidGeoElementMap();
            buildDofRelMap();
        }

    //std::vector<std::pair<size_type,rank_type> > const & geoElementMap() const { return M_geoElementMap; }
    std::map<size_type, std::pair<size_type,rank_type> > const & geoElementMap() const { return M_geoElementMap; }

    std::vector<size_type> const & dofRelMap() const { return this->dofRelMapRefToHo(); }

    std::vector<size_type> const & dofRelMapRefToHo() const { return M_dofRelMapRefToHo; }

    std::vector<size_type> const & dofRelMapHoToRef() const { return M_dofRelMapHoToRef; }

private :

    bool isIdenticalPoints(typename mesh1_type::element_type::point_type::super const & pt1,
                           typename mesh1_type::element_type::point_type::super const & pt2,
                           double tol = 1e-9 ) const;

    void buidGeoElementMap();

    void buildDofRelMap();

    std::vector<uint16_type>
    buildElementaryMapPoints(typename mesh1_type::element_type const & elt,
                             typename mesh2_type::element_type const & eltRef);

    std::vector<boost::tuple<uint16_type,uint16_type> >
    buildElementaryMapEdges(std::vector<uint16_type> const & mapPoint);

    std::vector<uint16_type>//std::vector<boost::tuple<uint> >
    buildElementaryMapFaces(std::vector<uint16_type> const & mapPoint);

    uint16_type
    convertInternalDofInFace(typename mesh1_type::element_type const & elem,uint16_type nface, uint16_type ilocModif,
                             std::vector<boost::tuple<uint16_type,uint16_type> > const & mapEdge,
                             std::vector<uint16_type> const & mapFace);


    std::vector<boost::tuple<uint16_type,uint16_type> >
    mapTrianglePoints2Edge(std::vector<uint16_type> const & mapPoint);
    std::vector<boost::tuple<uint16_type,uint16_type> >
    mapQuadranglePoints2Edge(std::vector<uint16_type> const & mapPoint);
    std::vector<boost::tuple<uint16_type,uint16_type> >
    mapTetraPoints2Edge(std::vector<uint16_type> const & mapPoint);

    std::vector<uint16_type>
    mapTrianglePoints2Face(std::vector<uint16_type> const & mapPoint);
    std::vector<uint16_type>
    mapQuadranglePoints2Face(std::vector<uint16_type> const & mapPoint);
    std::vector<uint16_type>
    mapTetraPoints2Face(std::vector<uint16_type> const & mapPoint);

    template<uint16_type nDofInFace>
    boost::tuple<uint16_type,uint16_type>
    tableInternalDofFace2Edge(uint16_type dofLoc);
    boost::tuple<uint16_type,uint16_type>
    tableInternalDofFace2Edge1(uint16_type dofLoc);
    boost::tuple<uint16_type,uint16_type>
    tableInternalDofFace2Edge3(uint16_type dofLoc);
    boost::tuple<uint16_type,uint16_type>
    tableInternalDofFace2Edge6(uint16_type dofLoc);

    template<uint16_type nDofInFace>
    uint16_type
    tableInternalDofEdge2Face(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc);
    uint16_type
    tableInternalDofEdge2Face1(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc);
    uint16_type
    tableInternalDofEdge2Face3(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc);
    uint16_type
    tableInternalDofEdge2Face6(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc);


private :

    functionspace1_ptrtype M_Xh1;
    functionspace2_ptrtype M_Xh2;

    //std::vector< std::pair<size_type,rank_type> > M_geoElementMap;
    std::map<size_type, std::pair<size_type,rank_type> > M_geoElementMap;

    std::vector<size_type> M_dofRelMapRefToHo;
    std::vector<size_type> M_dofRelMapHoToRef;

};


//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
bool
DofRelationshipMap<SpaceType1,SpaceType2>::isIdenticalPoints( typename mesh1_type::element_type::point_type::super const & pt1,
                                                              typename mesh1_type::element_type::point_type::super const & pt2,
                                                              double tol ) const
{
    bool res;
    if (mesh1_type::nRealDim==1)
        if ( std::abs( pt1(0)-pt2(0) ) < tol) res=true;
        else res=false;
    else if (mesh1_type::nRealDim==2)
        if ( std::abs( pt1(0)-pt2(0) ) < tol && std::abs( pt1(1)-pt2(1) ) < tol ) res=true;
        else res=false;
    else // if (mesh1_type::nRealDim==3)
        if ( std::abs( pt1(0)-pt2(0) ) < tol && std::abs( pt1(1)-pt2(1) ) < tol && std::abs( pt1(2)-pt2(2) ) < tol ) res=true;
        else res=false;

    return res;
}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
void
DofRelationshipMap<SpaceType1,SpaceType2>::buidGeoElementMap()
{
    std::vector<bool> findPtInElem(mesh1_type::element_type::numVertices);
    std::set<size_type> elt1Done;

    CHECK ( M_Xh1->dof()->buildDofTableMPIExtended() == M_Xh2->dof()->buildDofTableMPIExtended() ) << "buildDofTableMPIExtended between space must be equal \n";
#if 0
    bool upExtendedElt = M_Xh1->dof()->buildDofTableMPIExtended();
    EntityProcessType entityProcess = (upExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    auto rangeElt1 = elements( M_Xh1->mesh(), entityProcess );
    auto rangeElt2 = elements( M_Xh2->mesh(), entityProcess );
#else
    auto rangeElt1 = M_Xh1->template meshSupport<0>()->rangeElements();
    auto rangeElt2 = M_Xh2->template meshSupport<0>()->rangeElements();
#endif

    auto dof1 = M_Xh1->dof();
    auto dof2 = M_Xh2->dof();
    bool isPartialSupport1 = dof1->hasMeshSupport() && dof1->meshSupport()->isPartialSupport();

    M_geoElementMap.clear();

    bool useSubMeshRelation = M_Xh1->mesh()->isSubMeshFrom( M_Xh2->mesh() );

    for ( auto const& elt2RefWrapper : rangeElt2 )
    {
        auto const& elt2 = boost::unwrap_ref( elt2RefWrapper );

        if ( useSubMeshRelation )
        {
            size_type elt1Id = M_Xh1->mesh()->meshToSubMesh( elt2.id() );
            if ( elt1Id == invalid_v<size_type> )
                continue;
            if ( isPartialSupport1 && !dof1->isElementDone( elt1Id ) )
                continue;
            M_geoElementMap[elt1Id] = std::make_pair(elt2.id(),elt2.processId());
        }
        else
        {
            bool find=false;
            for ( auto const& elt1RefWrapper : rangeElt1 )
            {
                auto const& elt1 = boost::unwrap_ref( elt1RefWrapper );

                size_type elt1Id = elt1.id();
                if ( M_geoElementMap.find( elt1Id ) != M_geoElementMap.end() )
                    continue;

                bool findSameElt = true;
                for (uint16_type n=0 ; n<numVertices && findSameElt ; ++n)
                {
                    auto const& pt2 = elt2.point(n);
                    bool findSamePt = false;
                    for (uint16_type m=0 ; m<numVertices && !findSamePt ; ++m)
                        if ( this->isIdenticalPoints(elt1.point(m),pt2) )
                            findSamePt = true;
                    if ( !findSamePt )
                        findSameElt = false;
                }

                if ( !findSameElt )
                    continue;


                find=true;
                M_geoElementMap[elt1Id] = std::make_pair(elt2.id(),elt2.processId());
                break;
            }
            CHECK( find ) <<"not found a relation betwwen elt";
        }
    } // end for it


#if 0
    // add also some ghost elt if has extended dof table
    if ( M_Xh1->dof()->buildDofTableMPIExtended() )
    {
        auto face_it = M_Xh2->mesh()->interProcessFaces().first;
        auto const face_en = M_Xh2->mesh()->interProcessFaces().second;
        for ( ; face_it!=face_en ; ++face_it )
        {
            auto const& elt0 = face_it->element0();
            auto const& elt1 = face_it->element1();
            const bool elt0isGhost = elt0.isGhostCell();
            auto const& eltOffProc = (elt0isGhost)?elt0:elt1;


            auto it1 = M_Xh1->mesh()->beginGhostElement();
            auto en1 = M_Xh1->mesh()->endGhostElement();
            bool find=false;
            while (it1!=en1 && !find )
            {
                std::fill ( findPtInElem.begin(), findPtInElem.end(),false);

                for (uint16_type n=0 ; n<numVertices ; ++n)
                    for (uint16_type m=0 ; m<numVertices ; ++m)
                        if (isIdenticalPoints(it1->point(n),eltOffProc.point(m))) findPtInElem[n]=true;

                // All points of faces are the same?
                find=true;
                auto itbool=findPtInElem.begin();
                auto itbool_end=findPtInElem.end();
                for (  ; itbool != itbool_end; ++itbool )
                    find = (find && *itbool);

                if (find) M_geoElementMap[it1->id()]=std::make_pair(eltOffProc.id(),eltOffProc.processId());
                //if (find) M_geoElementMap[it->id()]=itP1->id();
                ++it1;

            }
            CHECK( find ) <<"\nProbleme!!!!!!!!!\n";
        }
    }
#endif


}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
void
DofRelationshipMap<SpaceType1,SpaceType2>::buildDofRelMap()
{
    auto dof1 = M_Xh1->dof();
    auto dof2 = M_Xh2->dof();

    //auto it1 = M_Xh1->mesh()->beginElementWithProcessId( M_Xh1->mesh()->worldComm().localRank() );
    //auto en1 = M_Xh1->mesh()->endElementWithProcessId( M_Xh1->mesh()->worldComm().localRank() );
    auto it1 = M_Xh1->mesh()->beginElement();
    auto const en1 = M_Xh1->mesh()->endElement();
    for ( ; it1 != en1 ; ++it1 )
    {
        auto const& elem1 = it1->second;
        auto itFindGeoElt = M_geoElementMap.find( elem1.id() );
        if ( itFindGeoElt == M_geoElementMap.end() )
            continue;

        size_type eltIdRelated = itFindGeoElt->second.first;
        rank_type procIdRelated = itFindGeoElt->second.second;
        //if ( M_geoElementMap[it1->id()].first == invalid_size_type_value ) continue;

        auto const& elem2 = M_Xh2->mesh()->element( eltIdRelated );

        auto mapPoint=buildElementaryMapPoints(elem1,elem2);
        auto mapEdge=buildElementaryMapEdges(mapPoint);
        auto mapFace=buildElementaryMapFaces(mapPoint);

        for ( uint16_type iloc1 = 0; iloc1 < functionspace1_type::basis_type::nLocalDof; ++iloc1 )
        {
            //we search this num local
            uint16_type iloc2=iloc1;

            if ( iloc1 < nDofPerVertex*numVertices )
                iloc2=mapPoint[iloc1];
            else if (iloc1 < nDofPerVertex*numVertices + nDofPerEdge*numEdges )
                for (uint16_type ndpe=0;ndpe<numEdges;++ndpe) {
                    if ( (iloc1 >= (nDofPerVertex*numVertices + nDofPerEdge*ndpe)) &&
                         (iloc1 < (nDofPerVertex*numVertices + nDofPerEdge*(ndpe+1))) ) {
                        auto numLocEdge = boost::get<0>(mapEdge[ndpe]);
                        auto PermLocEdge = boost::get<1>(mapEdge[ndpe]);
                        if (PermLocEdge==0)
                            iloc2 = ( nDofPerVertex*numVertices + numLocEdge*nDofPerEdge +
                                      (iloc1-nDofPerVertex*numVertices-ndpe*nDofPerEdge) );
                        else
                            iloc2 = ( nDofPerVertex*numVertices + numLocEdge*nDofPerEdge +
                                      nDofPerEdge-1-(iloc1-nDofPerVertex*numVertices-ndpe*nDofPerEdge) );
                    }
                }
            else if (iloc1 < nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*numGeometricFaces )
                for (uint16_type ndpf=0;ndpf<numGeometricFaces;++ndpf) {
                    if ( (iloc1 >= (nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*ndpf )) &&
                         (iloc1 < (nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*(ndpf+1))) ) {
                        auto numLocFace = mapFace[ndpf];
                        auto ilocModif = iloc1 - (nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*ndpf);
                        auto numDofLocInface=convertInternalDofInFace(elem1,ndpf,ilocModif,mapEdge,mapFace);
                        iloc2 = nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*numLocFace + numDofLocInface;
                    }
                }
            else if (iloc1 < nDofPerVertex*numVertices + nDofPerEdge*numEdges + nDofPerFace*numGeometricFaces + nDofPerVolume*numVolumes)
                iloc2=iloc1;// we guess only one dof in volume ( 3d,order=4 only)

            for ( uint16_type comp = 0;comp < functionspace1_type::basis_type::nComponents;++comp )
            {
#if 0
                CHECK( dof1->isElementDone( elem1.id() ) ) << "dofTable1 not build this elt id : " << elem1.id();
                CHECK( dof2->isElementDone( eltIdRelated ) ) << "dofTable2 not build this elt id : " << eltIdRelated;
#endif
                size_type i1 = dof1->localToGlobal( elem1.id(), iloc1 , comp ).index();
                size_type i2 = dof2->localToGlobal( eltIdRelated , iloc2 , comp ).index();

                M_dofRelMapRefToHo[i1]=i2;
                M_dofRelMapHoToRef[i2]=i1;

            }

        }

    } //end it

}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
std::vector<uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::buildElementaryMapPoints(typename mesh1_type::element_type const & elt,
                                                                    typename mesh2_type::element_type const & eltRef)
{
    std::vector<uint16_type> mapPoint(numVertices);

    for (uint16_type n1=0;n1<numVertices;++n1)
        for (uint16_type n2=0;n2<numVertices;++n2)
            if (isIdenticalPoints(elt.point(n1),eltRef.point(n2))) mapPoint[n1]=n2;

    return mapPoint;
}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
std::vector<boost::tuple<uint16_type,uint16_type> >
DofRelationshipMap<SpaceType1,SpaceType2>::buildElementaryMapEdges(std::vector<uint16_type> const & mapPoint)
{

    if (nDim==2)
    {
        if (is_simplex) { return mapTrianglePoints2Edge(mapPoint); }
        else { return mapQuadranglePoints2Edge(mapPoint); }
    }
    else if (nDim==3)
    {
        if (is_simplex) { return mapTetraPoints2Edge(mapPoint); }
    }

}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
std::vector<uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::buildElementaryMapFaces(std::vector<uint16_type> const & mapPoint)
{
    if (nDim==2)
    {
        if (is_simplex) { return mapTrianglePoints2Face(mapPoint); }
        else { return mapQuadranglePoints2Face(mapPoint); }
    }
    else if (nDim==3)
    {
        if (is_simplex) { return mapTetraPoints2Face(mapPoint); }
    }
}

//---------------------------------------------------------------------------------//


template< class SpaceType1,class SpaceType2 >
uint16_type
DofRelationshipMap<SpaceType1,SpaceType2>::convertInternalDofInFace(typename mesh1_type::element_type const & elem,
                                                                    uint16_type nface, uint16_type ilocModif,
                                                                    std::vector<boost::tuple<uint16_type,uint16_type> > const & mapEdge,
                                                                    std::vector<uint16_type> const & mapFace)
{
    auto traitLocal = /*template*/ tableInternalDofFace2Edge<nDofPerFace>(ilocModif);

    auto traitGlobal = boost::make_tuple( elem.f2e(nface, boost::get<0>(traitLocal)),elem.f2e(nface,boost::get<1>(traitLocal)) );

    auto traitGlobal2 = boost::make_tuple( boost::get<0>(mapEdge[boost::get<0>(traitGlobal)]),
                                           boost::get<0>(mapEdge[boost::get<1>(traitGlobal)]) );

    auto traitLocal2 = boost::make_tuple( elem.f2eLoc( mapFace[nface], boost::get<0>(traitGlobal2)),
                                          elem.f2eLoc( mapFace[nface], boost::get<1>(traitGlobal2)) );

    return tableInternalDofEdge2Face<nDofPerFace>(traitLocal2);
}



template< class SpaceType1,class SpaceType2 >
std::vector<boost::tuple<uint16_type,uint16_type> >
DofRelationshipMap<SpaceType1,SpaceType2>::mapTrianglePoints2Edge(std::vector<uint16_type> const & mapPoint)
{
    // 2 pts -> (num glob of edge, permutation)
    std::map<boost::tuple<uint16_type,uint16_type>,boost::tuple<uint16_type,uint16_type> > mapEdgesBasis;
    mapEdgesBasis[ boost::make_tuple(0,1) ] = boost::make_tuple(2,0);
    mapEdgesBasis[ boost::make_tuple(1,0) ] = boost::make_tuple(2,1);

    mapEdgesBasis[ boost::make_tuple(2,0) ] = boost::make_tuple(1,0);
    mapEdgesBasis[ boost::make_tuple(0,2) ] = boost::make_tuple(1,1);

    mapEdgesBasis[ boost::make_tuple(1,2) ] = boost::make_tuple(0,0);
    mapEdgesBasis[ boost::make_tuple(2,1) ] = boost::make_tuple(0,1);


    //uint le num, uint permutation : 0 or 1
    std::vector<boost::tuple<uint16_type,uint16_type> > mapEdges(3);
    mapEdges[0]= mapEdgesBasis[ boost::make_tuple( mapPoint[1],mapPoint[2]) ];
    mapEdges[1]= mapEdgesBasis[ boost::make_tuple( mapPoint[2],mapPoint[0]) ];
    mapEdges[2]= mapEdgesBasis[ boost::make_tuple( mapPoint[0],mapPoint[1]) ];

    return mapEdges;
}

template< class SpaceType1,class SpaceType2 >
std::vector<boost::tuple<uint16_type,uint16_type> >
DofRelationshipMap<SpaceType1,SpaceType2>::mapQuadranglePoints2Edge(std::vector<uint16_type> const & mapPoint)
{
    // 2 pts -> (num glob of edge, permutation)
    std::map<boost::tuple<uint16_type,uint16_type>,boost::tuple<uint16_type,uint16_type> > mapEdgesBasis;
    mapEdgesBasis[ boost::make_tuple(0,1) ] = boost::make_tuple(0,0);
    mapEdgesBasis[ boost::make_tuple(1,0) ] = boost::make_tuple(0,1);

    mapEdgesBasis[ boost::make_tuple(1,2) ] = boost::make_tuple(1,0);
    mapEdgesBasis[ boost::make_tuple(2,1) ] = boost::make_tuple(1,1);

    mapEdgesBasis[ boost::make_tuple(2,3) ] = boost::make_tuple(2,0);
    mapEdgesBasis[ boost::make_tuple(3,2) ] = boost::make_tuple(2,1);

    mapEdgesBasis[ boost::make_tuple(3,0) ] = boost::make_tuple(3,0);
    mapEdgesBasis[ boost::make_tuple(0,3) ] = boost::make_tuple(3,1);

    //uint le num, uint permutation : 0 or 1
    std::vector<boost::tuple<uint16_type,uint16_type> > mapEdges(4);
    mapEdges[0]= mapEdgesBasis[ boost::make_tuple( mapPoint[0],mapPoint[1]) ];
    mapEdges[1]= mapEdgesBasis[ boost::make_tuple( mapPoint[1],mapPoint[2]) ];
    mapEdges[2]= mapEdgesBasis[ boost::make_tuple( mapPoint[2],mapPoint[3]) ];
    mapEdges[3]= mapEdgesBasis[ boost::make_tuple( mapPoint[3],mapPoint[0]) ];

    return mapEdges;
}


template< class SpaceType1,class SpaceType2 >
std::vector<boost::tuple<uint16_type,uint16_type> >
DofRelationshipMap<SpaceType1,SpaceType2>::mapTetraPoints2Edge(std::vector<uint16_type> const & mapPoint)
{
    // 2 pts -> (num glob of edge, permutation)
    std::map<boost::tuple<uint16_type,uint16_type>,boost::tuple<uint16_type,uint16_type> > mapEdgesBasis;
    mapEdgesBasis[ boost::make_tuple(0,1) ] = boost::make_tuple(2,0);
    mapEdgesBasis[ boost::make_tuple(0,2) ] = boost::make_tuple(1,1);
    mapEdgesBasis[ boost::make_tuple(0,3) ] = boost::make_tuple(3,0);

    mapEdgesBasis[ boost::make_tuple(1,0) ] = boost::make_tuple(2,1);
    mapEdgesBasis[ boost::make_tuple(1,2) ] = boost::make_tuple(0,0);
    mapEdgesBasis[ boost::make_tuple(1,3) ] = boost::make_tuple(4,0);

    mapEdgesBasis[ boost::make_tuple(2,0) ] = boost::make_tuple(1,0);
    mapEdgesBasis[ boost::make_tuple(2,1) ] = boost::make_tuple(0,1);
    mapEdgesBasis[ boost::make_tuple(2,3) ] = boost::make_tuple(5,0);

    mapEdgesBasis[ boost::make_tuple(3,0) ] = boost::make_tuple(3,1);
    mapEdgesBasis[ boost::make_tuple(3,1) ] = boost::make_tuple(4,1);
    mapEdgesBasis[ boost::make_tuple(3,2) ] = boost::make_tuple(5,1);


    //uint le num, uint permutation : 0 or 1
    std::vector<boost::tuple<uint16_type,uint16_type> > mapEdges(6);
    mapEdges[0]= mapEdgesBasis[ boost::make_tuple( mapPoint[1],mapPoint[2]) ];
    mapEdges[1]= mapEdgesBasis[ boost::make_tuple( mapPoint[2],mapPoint[0]) ];
    mapEdges[2]= mapEdgesBasis[ boost::make_tuple( mapPoint[0],mapPoint[1]) ];
    mapEdges[3]= mapEdgesBasis[ boost::make_tuple( mapPoint[0],mapPoint[3]) ];
    mapEdges[4]= mapEdgesBasis[ boost::make_tuple( mapPoint[1],mapPoint[3]) ];
    mapEdges[5]= mapEdgesBasis[ boost::make_tuple( mapPoint[2],mapPoint[3]) ];

    return mapEdges;
}

//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
std::vector<uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::mapTrianglePoints2Face(std::vector<uint16_type> const & mapPoint)
{
    std::vector<uint16_type> mapFaces(1);// 1 faces
    mapFaces[0]= 0; // trival : only 1 face

    return mapFaces;
}

template< class SpaceType1,class SpaceType2 >
std::vector<uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::mapQuadranglePoints2Face(std::vector<uint16_type> const & mapPoint)
{
    std::vector<uint16_type> mapFaces(1);// 1 faces
    mapFaces[0]= 0; // trival : only 1 face

    return mapFaces;
}

template< class SpaceType1,class SpaceType2 >
std::vector<uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::mapTetraPoints2Face(std::vector<uint16_type> const & mapPoint)
{
    // 3 pts -> (num glob of face )
    std::map<boost::tuple<uint16_type,uint16_type,uint16_type>,uint16_type> mapFacesBasis;
    mapFacesBasis[ boost::make_tuple(1,2,3) ] = 0;
    mapFacesBasis[ boost::make_tuple(1,3,2) ] = 0;
    mapFacesBasis[ boost::make_tuple(2,3,1) ] = 0;
    mapFacesBasis[ boost::make_tuple(2,1,3) ] = 0;
    mapFacesBasis[ boost::make_tuple(3,1,2) ] = 0;
    mapFacesBasis[ boost::make_tuple(3,2,1) ] = 0;

    mapFacesBasis[ boost::make_tuple(0,3,2) ] = 1;
    mapFacesBasis[ boost::make_tuple(0,2,3) ] = 1;
    mapFacesBasis[ boost::make_tuple(2,0,3) ] = 1;
    mapFacesBasis[ boost::make_tuple(2,3,0) ] = 1;
    mapFacesBasis[ boost::make_tuple(3,2,0) ] = 1;
    mapFacesBasis[ boost::make_tuple(3,0,2) ] = 1;

    mapFacesBasis[ boost::make_tuple(0,1,3) ] = 2;
    mapFacesBasis[ boost::make_tuple(0,3,1) ] = 2;
    mapFacesBasis[ boost::make_tuple(1,3,0) ] = 2;
    mapFacesBasis[ boost::make_tuple(1,0,3) ] = 2;
    mapFacesBasis[ boost::make_tuple(3,0,1) ] = 2;
    mapFacesBasis[ boost::make_tuple(3,1,0) ] = 2;

    mapFacesBasis[ boost::make_tuple(0,1,2) ] = 3;
    mapFacesBasis[ boost::make_tuple(0,2,1) ] = 3;
    mapFacesBasis[ boost::make_tuple(1,2,0) ] = 3;
    mapFacesBasis[ boost::make_tuple(1,0,2) ] = 3;
    mapFacesBasis[ boost::make_tuple(2,0,1) ] = 3;
    mapFacesBasis[ boost::make_tuple(2,1,0) ] = 3;


    std::vector<uint16_type> mapFaces(4);// 4 faces
    mapFaces[0]= mapFacesBasis[ boost::make_tuple( mapPoint[1],mapPoint[2],mapPoint[3]) ];
    mapFaces[1]= mapFacesBasis[ boost::make_tuple( mapPoint[0],mapPoint[2],mapPoint[3]) ];
    mapFaces[2]= mapFacesBasis[ boost::make_tuple( mapPoint[0],mapPoint[1],mapPoint[3]) ];
    mapFaces[3]= mapFacesBasis[ boost::make_tuple( mapPoint[0],mapPoint[1],mapPoint[2]) ];

    return mapFaces;
}

//---------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------//

template< class SpaceType1,class SpaceType2 >
template<uint16_type nDofInFace>
boost::tuple<uint16_type,uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofFace2Edge(uint16_type dofLoc)
{
    if (nDofInFace==1)
        return tableInternalDofFace2Edge1(dofLoc);
    else if (nDofInFace==3)
        return tableInternalDofFace2Edge3(dofLoc);
    else if (nDofInFace==6)
        return tableInternalDofFace2Edge6(dofLoc);
    else // bad case
        return boost::make_tuple(0,0);

}


template< class SpaceType1,class SpaceType2 >
template<uint16_type nDofInFace>
uint16_type
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofEdge2Face(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc)
{
    if (nDofInFace==1)
        return tableInternalDofEdge2Face1(traitsEdgesLoc);
    else if (nDofInFace==3)
        return tableInternalDofEdge2Face3(traitsEdgesLoc);
    else if (nDofInFace==6)
        return tableInternalDofEdge2Face6(traitsEdgesLoc);
    else // bad case
        return 0;
}

//---------------------------------------------------------------------------------//

//simplex or hypercube
template< class SpaceType1,class SpaceType2 >
boost::tuple<uint16_type,uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofFace2Edge1(uint16_type dofLoc)
{
    //valabale pour 3 dof in the face( thetra bien sur)
    std::vector<boost::tuple<uint16_type,uint16_type> > table(1);
    table[0] = boost::make_tuple(0,0);

    return table[dofLoc];
}

//---------------------------------------------------------------------------------//

//simplex or hypercube
template< class SpaceType1,class SpaceType2 >
uint16_type
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofEdge2Face1(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc)
{
    //valabale pour 1 dof in the face( thetra bien sur)
    std::map<boost::tuple<uint16_type,uint16_type>,uint16_type > table;
    table[boost::make_tuple(0,0)]=0;

    return table[traitsEdgesLoc];
}

//---------------------------------------------------------------------------------//

// simplex
template< class SpaceType1,class SpaceType2 >
boost::tuple<uint16_type,uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofFace2Edge3(uint16_type dofLoc)
{
    //valabale pour 3 dof in the face( thetra bien sur)
    std::vector<boost::tuple<uint16_type,uint16_type> > table(3);
    table[0] = boost::make_tuple(1,2);
    table[1] = boost::make_tuple(0,2);
    table[2] = boost::make_tuple(0,1);

    return table[dofLoc];
}

//---------------------------------------------------------------------------------//

// simplex
template< class SpaceType1,class SpaceType2 >
uint16_type
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofEdge2Face3(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc)
{
    //valabale pour 3 dof in the face( thetra bien sur)
    std::map<boost::tuple<uint16_type,uint16_type>,uint16_type > table;
    table[boost::make_tuple(1,2)]=0;
    table[boost::make_tuple(2,1)]=0;
    table[boost::make_tuple(0,2)]=1;
    table[boost::make_tuple(2,0)]=1;
    table[boost::make_tuple(0,1)]=2;
    table[boost::make_tuple(1,0)]=2;

    return table[traitsEdgesLoc];
}


// simplex
template< class SpaceType1,class SpaceType2 >
boost::tuple<uint16_type,uint16_type>
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofFace2Edge6(uint16_type dofLoc)
{
    //valabale pour 3 dof in the face( thetra bien sur)
    std::vector<boost::tuple<uint16_type,uint16_type> > table(6);
    table[0] = boost::make_tuple(1,2);
    table[1] = boost::make_tuple(2,2);
    table[2] = boost::make_tuple(2,0);
    table[3] = boost::make_tuple(1,1);
    table[4] = boost::make_tuple(0,0);
    table[5] = boost::make_tuple(0,1);

    return table[dofLoc];
}

// simplex
template< class SpaceType1,class SpaceType2 >
uint16_type
DofRelationshipMap<SpaceType1,SpaceType2>::tableInternalDofEdge2Face6(boost::tuple<uint16_type,uint16_type> traitsEdgesLoc)
{
    //valabale pour 3 dof in the face( thetra bien sur)
    std::map<boost::tuple<uint16_type,uint16_type>,uint16_type > table;
    table[boost::make_tuple(1,2)]=0;
    table[boost::make_tuple(2,1)]=0;
    table[boost::make_tuple(2,2)]=1;
    table[boost::make_tuple(2,0)]=2;
    table[boost::make_tuple(0,2)]=2;
    table[boost::make_tuple(1,1)]=3;
    table[boost::make_tuple(0,0)]=4;
    table[boost::make_tuple(0,1)]=5;
    table[boost::make_tuple(1,0)]=5;

    return table[traitsEdgesLoc];
}



template< class SpaceType1,class SpaceType2 >
std::shared_ptr< DofRelationshipMap<SpaceType1,SpaceType2> >
dofRelationshipMap( std::shared_ptr<SpaceType1> Xh1, std::shared_ptr<SpaceType2> Xh2 )
{
    std::shared_ptr< DofRelationshipMap<SpaceType1,SpaceType2> > drm( new DofRelationshipMap<SpaceType1,SpaceType2>(Xh1,Xh2) );
    return drm;
}

} // namespace FeelModels
} // namespace Feel


#endif // FEELPP_MODELS_DOFRELATIONSHIPMAP_H
