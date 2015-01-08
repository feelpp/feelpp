/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2008-02-07

 Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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

/**
 \file disttocurve.hpp
 \author Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
 \date 2014-02-02
 */



/*
 This class provides an algorithm to make a distance function from a curve represented by list of points on the elements crossed by the curve. It can be coupled with a fast marching algorithm to make a distance function in the whole domain.
 */


#ifndef DISTTOCURVE_HPP
#define DISTTOCURVE_HPP 1

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <random>
#include <feel/feelvf/vf.hpp>

//#define DISTANCE_FROM_UNORDERED_POINTS 1
/*
 The distance to a set of unordered points is still experimental for now.
 The absolute distance value is always computed correctly.
 The sign assignment is not always working properly:
 - It works for shapes where the lines are not too thick.
 - It fails when there is some line thickness.
 Depending of future applications, the way to track the side of the curve
 might be changed.

 For now, since no application is using this feature and since it is a
 big piece of code, I disable its compilation via the DISTANCE_FROM_UNORDERED_POINTS macro.
 */



namespace Feel
{

template< typename FunctionSpaceP0Type, typename FunctionSpaceP1Type >
class DistToCurve
{

public :

    /* some requirement on the spaces :
     - P0 discontinuous space
     - P1 space  */

    static_assert( FunctionSpaceP0Type::fe_type::nOrder == 0, "FunctionSpaceP0Type needs to be a finite element space of order 0");
    static_assert( FunctionSpaceP0Type::fe_type::continuity_type::is_discontinuous_totally, "P0 space has to be discontinuous");
    static_assert( FunctionSpaceP1Type::fe_type::nOrder == 1, "FunctionSpaceP1Type needs to be a finite element space of order 1");

    typedef DistToCurve<FunctionSpaceP0Type, FunctionSpaceP1Type> self_type;
    typedef boost::shared_ptr< self_type > self_ptrtype;

    typedef typename boost::shared_ptr< FunctionSpaceP0Type > spaceP0_ptrtype;
    typedef typename boost::shared_ptr< FunctionSpaceP1Type > spaceP1_ptrtype;

    typedef typename FunctionSpaceP1Type::element_ptrtype element_ptrtype;
    typedef typename FunctionSpaceP1Type::value_type value_type;
    typedef typename FunctionSpaceP1Type::mesh_type mesh_type;
    typedef typename mesh_type::node_type node_type;
    
    // Used to store a node and its distance to a point
    typedef std::pair< node_type, double > nodeDist_type;

    enum side_type {sideA, sideB, noSide};

    /** @name Constructors, destructor
     */
    //@{
    DistToCurve( spaceP0_ptrtype spaceP0, spaceP1_ptrtype spaceP1)
        :
        firstDof( spaceP1->dof()->firstDofGlobalCluster() ),
        ndofv( FunctionSpaceP1Type::fe_type::nDof ),
        dim(FunctionSpaceP1Type::fe_type::nRealDim),
        periodT2( 0 ),
        M_spaceP0( spaceP0 ),
        M_spaceP1( spaceP1 ),
        M_mesh( spaceP1->mesh() )
        {
            // make the map ghostClusterToProc
            // TODO : should be in datamap.hpp
            if (Environment::worldComm().size() > 1)
                for (int k=0; k<M_spaceP1->nLocalDof(); ++k)
                    if ( M_spaceP1->dof()->dofGlobalProcessIsGhost( k ) )
                        ghostClusterToProc[ processorToCluster( k ) ] = k;


            // create a P0 elt containing the ids of the elements
            ids = M_spaceP0->element();
            for (auto const& it : elements(M_mesh) )
                ids.assign( it.id(), 0, 0, it.id() );

            eltHavingPoints = M_spaceP0->element();

        } //DistToCurve


    static self_ptrtype New( spaceP0_ptrtype spaceP0, spaceP1_ptrtype spaceP1)
        {
            self_ptrtype dtc( new self_type( spaceP0, spaceP1) );
            return dtc;
        }

    //@}

    // TODO : create the same function where the points are read from a file (x,y) with M_dt = 1



    /* 2d and 3d versions*/
    // TExpr = std::function<double(double, double)> or TExpr = std::function<double(double)>
    template <class TExprX, class TExprY, class TExprZ>
    element_ptrtype fromParametrizedCurve( TExprX xexpr, TExprY yexpr, TExprZ zexpr,
                                           double t1Start, double t1End, double dt1,
                                           double t2Start, double t2End, double dt2,
                                           bool broadenCurveForElementDetection = true,
                                           double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2.,
                                           bool exportPoints = false,
                                           std::string exportName=""
                                           )
        {
            clear();

            generatePointsFromParametrization(xexpr, yexpr, zexpr,
                                              t1Start, t1End, dt1,
                                              t2Start, t2End, dt2,
                                              exportPoints, exportName);

            locateElementsCrossedByCurve(broadenCurveForElementDetection, broadenessAmplitude );

            auto shape = makeDistanceFunctionSequential();

            clear(); // no need for the maps used to create the distance function

            if ( Environment::worldComm().size() > 1)
                reduceDistanceFunction( shape );

            return shape;
        }


#if defined( DISTANCE_FROM_UNORDERED_POINTS )

    /* 2d version only */
    element_ptrtype fromParametrizedCurveDisordered(std::tuple< std::function<double(double)>, std::function<double(double)>, double, double > paramFct,
                                                    double dt,
                                                    bool exportPoints = false,
                                                    std::string exportName=""
                                                    )
        {
            clear();

            generatePointsFromParametrization(get<0>(paramFct), get<1>(paramFct),
                                              get<2>(paramFct), get<3>(paramFct), dt,
                                              exportPoints, exportName);

            for (auto const& tnd : tNodeMap )
                allPoints.push_back( tnd.second );

            locateElementsCrossedByUnorderedPoints();

            auto shape = makeDistanceFunctionSequentialFromUnorderedPoints();

            auto shape_unsigned = *shape;

            node_type pt(dim);
            node_type pt2(dim);
            pt(0) = 0.9; pt(1) = 0.9;
            pt2(0) = 0.9; pt2(1) = 0.1;
            std::vector< node_type > lstPoints( {pt, pt2} );

            setInnerRegion( shape, lstPoints );

            clear();

            auto shape_unreduced = *shape;

            if ( Environment::worldComm().size() > 1)
                reduceDistanceFunction( shape );

            return shape;
        }


    // TFilename = boost::filesystem::path or std::string
    template< class TFilename, class TList = std::vector<node_type> >
    element_ptrtype fromCoordinateFile( TFilename filename, TList insidePoints = std::vector<node_type>() )
        {
            clear();

            readPointsFromFile( filename );
            locateElementsCrossedByUnorderedPoints();
            auto shape = makeDistanceFunctionSequentialFromUnorderedPoints();

            auto shape_unsigned = *shape;

            //            setInnerRegion( shape, insidePoints, true );

            clear();

            if ( Environment::worldComm().size() > 1)
                reduceDistanceFunction( shape );


            auto mark2 = vf::project(M_spaceP0, marked2elements(M_mesh, 1), cst(1) );
            auto exp = exporter(_mesh=M_mesh, _name="disttocurvehpp");
            exp->step(0)->add("shape_unsigned", shape_unsigned);
            exp->step(0)->add("shape_signed", *shape);
            exp->step(0)->add("mark2hpp", mark2);
            exp->save();

            return shape;

        } //fromCoordinateFile

#endif


    /* 2d version only */
    element_ptrtype fromParametrizedCurve( std::function<double(double)> xexpr,
                                           std::function<double(double)> yexpr,
                                           double t1Start, double t1End, double dt1,
                                           bool broadenCurveForElementDetection = true,
                                           double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2.,
                                           bool exportPoints = false,
                                           std::string exportName=""
                                           )
        {
            clear();

            generatePointsFromParametrization(xexpr, yexpr,
                                              t1Start, t1End, dt1,
                                              exportPoints, exportName);

            locateElementsCrossedByCurve(broadenCurveForElementDetection, broadenessAmplitude );

            auto shape = makeDistanceFunctionSequential();

            clear(); // no need for the maps used to create the distance function

            auto shapeunsigned = *shape;

            if ( Environment::worldComm().size() > 1)
                reduceDistanceFunction( shape );

            return shape;
        }





    // same function as the previous one but extract the functions x(t), y(t) from a tuple
    element_ptrtype fromParametrizedCurve( std::tuple< std::function<double(double)>, std::function<double(double)> > paramFct,
                                           double tStart, double tEnd, double dt,
                                           bool broadenCurveForElementDetection = true,
                                           double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2.,
                                           bool exportPoints = false,
                                           std::string exportName=""
                                           )
        {
            return this->fromParametrizedCurve( get<0>(paramFct),
                                                get<1>(paramFct),
                                                tStart,  tEnd,  dt,
                                                broadenCurveForElementDetection,
                                                broadenessAmplitude,
                                                exportPoints,
                                                exportName );
        }

    // same function as the previous one but extract the functions x(t), y(t), tStart, tEnd, from a tuple ( can use the functions given in curveparametrizations.hpp )
    element_ptrtype fromParametrizedCurve( std::tuple< std::function<double(double)>, std::function<double(double)>, double, double > paramFct,
                                           double dt,
                                           bool broadenCurveForElementDetection = true,
                                           double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2.,
                                           bool exportPoints = false,
                                           std::string exportName=""
                                           )
        {
            return this->fromParametrizedCurve( get<0>(paramFct),
                                                get<1>(paramFct),
                                                get<2>(paramFct),
                                                get<3>(paramFct),
                                                dt,
                                                broadenCurveForElementDetection,
                                                broadenessAmplitude,
                                                exportPoints,
                                                exportName );
        }



private :

    // ----------- private attributes ---------
    static constexpr double bigdouble = 1e8;
    const size_type firstDof;
    const uint16_type ndofv;

    const int dim;
    int periodT2;

    // function spaces
    spaceP0_ptrtype M_spaceP0;
    spaceP1_ptrtype M_spaceP1;
    typename FunctionSpaceP1Type::mesh_ptrtype M_mesh;

    typename FunctionSpaceP0Type::element_type ids;

    // store a marker on the last elements having been localized
    typename FunctionSpaceP0Type::element_type eltHavingPoints;

    std::map< size_type, size_type > ghostClusterToProc;

    // ------ for ordered list of points
    // contains parameter number of node, node
    std::map< size_type, node_type > tNodeMap;

    // ------ for unordered list of points
    // contains list of points
    std::vector< node_type > allPoints;
    std::vector< node_type > refPoints;

    // contains as key the id of an element and value the numbers of the nodes (t) which are crossing it
    std::map<size_type, std::set<size_type> > pointsAtIndex;

    // for key=element, store on which side of the curve the dof is
    std::map<size_type, std::array<side_type, FunctionSpaceP1Type::fe_type::nDof > > dofIsOnSide;


    // -------- private methods ---------

    inline size_type clusterToProcessor( size_type dof )
        {return M_spaceP1->dof()->mapGlobalClusterToGlobalProcess( dof - firstDof ); }

    inline size_type processorToCluster( size_type dof )
        {return M_spaceP1->dof()->mapGlobalProcessToGlobalCluster( dof ); }


    double squareDistToPoint(node_type a, node_type b)
        {
            node_type diff = a - b;
            double sdist=0;

            for (int i=0; i<a.size(); ++i)
                sdist += diff[i] * diff[i];

            return sdist;
        }


    typename FunctionSpaceP0Type::element_type getCrossedElements()
        {return eltHavingPoints;}


    void clear()
        {
            tNodeMap.clear();
            pointsAtIndex.clear();
            allPoints.clear();
            dofIsOnSide.clear();
            refPoints.clear();
        }


    /* 2d version */
    template < class TExpr >
    void generatePointsFromParametrization( TExpr xexpr, TExpr yexpr,
                                            double t1Start, double t1End, double dt1,
                                            bool exportPoints = false, std::string exportName="" )
        {
            std::ofstream nodeFile;

            if (exportPoints && (Environment::worldComm().rank() == 0) )
            {
                std::string expName = exportName.empty() ? "nodes.particles" : exportName+".particles";
                nodeFile.open(expName, std::ofstream::out);
            }

            size_type count=0;

            for (double t=t1Start; t<t1End; t+=dt1)
            {
                node_type node(2);
                node[0] = xexpr(t);
                node[1] = yexpr(t);

                tNodeMap[ count ] = node;

                count++;

                if (exportPoints && (Environment::worldComm().rank() == 0) )
                    nodeFile << node[0] << "," << node[1] << ","<< "0" <<std::endl;
            }

            if (exportPoints && (Environment::worldComm().rank() == 0) )
                nodeFile.close();

        } // generatePointsFromParametrization


    /* 3d version*/
    template <class TExprX, class TExprY, class TExprZ>
    void generatePointsFromParametrization( TExprX xexpr, TExprY yexpr, TExprZ zexpr,
                                            double t1Start, double t1End, double dt1,
                                            double t2Start, double t2End, double dt2,
                                            bool exportPoints = false, std::string exportName="" )
        {
            std::ofstream nodeFile;

            if (exportPoints && (Environment::worldComm().rank() == 0) )
            {
                std::string expName = exportName.empty() ? "nodes.particles" : exportName+".particles";
                nodeFile.open(expName, std::ofstream::out);
            }

            int count = 0;

            for (double t1=t1Start; t1<t1End; t1+=dt1)
            {
                periodT2=0;

                for (double t2=t2Start; t2<t2End; t2+=dt2)
                {
                    node_type node(3);
                    node[0] = xexpr(t1,t2);
                    node[1] = yexpr(t1,t2);
                    node[2] = zexpr(t1,t2);

                    tNodeMap[count] = node;

                    if (exportPoints && (Environment::worldComm().rank() == 0) )
                        nodeFile << node[0] << "," << node[1] << ","<< node[2] <<std::endl;
                    periodT2++;
                    count++;
                }
            }

            if (exportPoints && (Environment::worldComm().rank() == 0) )
                nodeFile.close();

        } // generatePointsFromParametrization



    void locateElementsCrossedByCurve(bool randomlyBroadenNodesPositions=false, double randomnessAmplitude = option("gmsh.hsize").as<double>() / 2.)
        {
            // locate the elements crossed by the curve
            // store their ids in a map with the "t" of the nodes being in the element
            // update the marker2 with the elements being crossed
            CHECK( ! tNodeMap.empty() )<<"\n No nodes defining the curve have been loaded.\n";

            auto ctx = M_spaceP0->context();

            std::default_random_engine re( (unsigned int)::time(0) );
            std::uniform_real_distribution<double> smallRd( -randomnessAmplitude, randomnessAmplitude );

            for (auto const& tnode : tNodeMap)
            {
                node_type nodeToAdd = tnode.second;

                if (randomlyBroadenNodesPositions)
                    for (int i=0; i<dim; ++i)
                        nodeToAdd[i] += smallRd(re);

                // the following line takes 99% of the total time of the whole disttocurve algorithm !!!
                // (locates all the points in the mesh)
                ctx.add( nodeToAdd );
            }

            auto allIndexes = ids.evaluate( ctx );

            const int nbPtContext = ctx.nPoints();
            eltHavingPoints.zero();

            auto tnodeit = tNodeMap.begin();
            for (int i=0; i < nbPtContext; ++i )
            {

                if (Environment::worldComm().localRank() == ctx.processorHavingPoint( i ) )
                {
                    const size_type tOfNode = (*tnodeit).first;
                    const size_type index = allIndexes(i);

                    if ( pointsAtIndex.count( index ) )
                        pointsAtIndex[ index ].insert( tOfNode );
                    else
                    {
                        std::set< size_type > newset( {tOfNode} );
                        pointsAtIndex[ index ] = newset;
                    }

                    eltHavingPoints.assign( index, 0, 0, 1);
                }

                tnodeit++;

            }

            M_mesh->updateMarker2( eltHavingPoints );

        } // locateElementsCrossedByCurve




    // make distance function sequential when the points are ordered
    element_ptrtype makeDistanceFunctionSequential()
        {
            auto shape = M_spaceP1->elementPtr();
            *shape = vf::project(M_spaceP1, elements(M_mesh), cst(bigdouble) );


            // squared distance between a point where only its "t" is given, and a node nd2
            auto distToPt = [this] (double t, node_type nd2) -> double
                {
                    node_type nd1 = this->tNodeMap[ t ];
                    node_type diff = nd1 - nd2;
                    double sdist=0;
                    for (int i=0; i<diff.size(); ++i)
                        sdist += diff[i] * diff[i];
                    return sdist;
                };


            auto it_elt = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).first;
            auto en_elt = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).second;

            for(; it_elt!=en_elt; it_elt++)
                for (int j=0; j<ndofv; ++j)
                {
                    double closestDist = bigdouble;
                    double closestPoint = bigdouble;

                    const size_type indexGlobDof = M_spaceP1->dof()->localToGlobal(it_elt->id(), j, 0).index();

                    // coords of the dof
                    const node_type dofCoord = M_spaceP1->dof()->dofPoint( indexGlobDof ).template get<0>();

                    // find the point in the element having the closest distance with the dof. This distance will be the value of shape at this dof (if a smaller distance on the same dof is not found in an other element).
                    //This method assumes that the distance between the points of the curve is very small compared to the size of the mesh
                    for (auto const& pt : pointsAtIndex[ it_elt->id() ] )
                    {
                        const double dtp = distToPt( pt, dofCoord );
                        if( dtp < closestDist )
                        {
                            closestDist = dtp;
                            closestPoint = pt;
                        }
                    }


                    if ( closestDist < (*shape)(indexGlobDof) * (*shape)(indexGlobDof) )
                    {
                        const node_type closestPointCoord = tNodeMap[ closestPoint ];
                        // (tx, ty) = vector tangent to the parametrized curve at the closest point on param curve
                        node_type t1, t2;
                        // try to get the point next to the closest point (in the particular case where closest point is the last point, get the previous one)

                        try
                        { t1 = tNodeMap.at(closestPoint + 1.) - closestPointCoord; }
                        catch (const std::out_of_range& oor)
                        { t1 = closestPointCoord - tNodeMap.at(closestPoint - 1.); }

                        if (dim==3)
                        {
                            try
                            { t2 = tNodeMap.at(closestPoint + periodT2) - closestPointCoord; }
                            catch (const std::out_of_range& oor)
                            { t2 = closestPointCoord - tNodeMap.at(closestPoint - periodT2); }
                        }

                        // v = vector pointing from the closest point on curve to the concerned dof
                        const node_type v = dofCoord - closestPointCoord;

                        // the sign of the distance function is ruled by the vectorial product of the tangent vector and the vector v : sign(v x t)
                        // in 3D, it should be somthing like :  sign( (v x t) . n ) where n is the normal of the param surface pointing outward
                        double signProdVec;

                        if (dim==2)
                            signProdVec = v[0] * t1[1] - v[1] * t1[0] > 0 ? 1 : -1;
                        else
                        {
                            const double nx = t1[1]*t2[2]-t1[2]*t2[1];
                            const double ny = t1[2]*t2[0]-t1[0]*t2[2];
                            const double nz = t1[0]*t2[1]-t1[1]*t2[0];
                            signProdVec = v[0] * nx + v[1] * ny + v[2] * nz > 0 ? 1 : -1;
                        }

                        (*shape)( indexGlobDof ) = std::sqrt(closestDist) * signProdVec;

                    }
                }

            return shape;

        } // makeDistanceFunctionSequential


    void reduceDistanceFunction( element_ptrtype shape )
        {
            // given a distance function made by makeDistanceFunctionSequential which has different values on nodes being at the interface between several subdomain, make a nice, homogeneous distance function (requires several communications though all the proc !)

            auto eltHavingPointP1 = vf::project(M_spaceP1, marked2elements(M_mesh, 1), cst(1) );


            // search for all the dof being marked on at least one proc and being ghost on at least one proc (not necessarily the same proc)
            auto checkMarked = backend()->newVector( M_spaceP1 );
            auto checkGhost = backend()->newVector( M_spaceP1 );
            for (size_type k = 0 ; k < M_spaceP1->nLocalDof() ; ++k )
            {
                checkMarked->add(k,  eltHavingPointP1(k) == 1 ? 1 : 0 );
                checkGhost->add( k, M_spaceP1->dof()->dofGlobalProcessIsGhost( k ) ? 1 : 0 );
            }
            checkMarked->close();
            checkGhost->close();


            // store the value of all the dof being marked and ghost
            std::vector< std::pair< size_type, value_type > > idOnClusterAndValue;
            for (size_type k = 0 ; k < M_spaceP1->nLocalDof() ; ++k )
                if ( (*checkGhost)(k) && (*checkMarked)(k) )
                    idOnClusterAndValue.push_back( { processorToCluster(k), (*shape)(k) } );



            // gather all these values to one single proc
            std::vector< std::vector< std::pair< size_type, value_type > > > allIdOnClusterAndValue;
            mpi::gather( Environment::worldComm().globalComm(),
                         idOnClusterAndValue,
                         allIdOnClusterAndValue, 0 );


            // this proc makes all the work :
            // extract all the values in a map having (key = GlobalIdOnCluster, value = min value of shape )
            std::map< size_type, value_type > idOnClusterAndMinValue;
            if (Environment::worldComm().localRank() == 0)
                for(auto const& v1 : allIdOnClusterAndValue)
                    for (auto const& idValue : v1)
                    {
                        if (idOnClusterAndMinValue.count( idValue.first ) )
                            idOnClusterAndMinValue[ idValue.first ] =
                                std::abs( idValue.second ) < std::abs( idOnClusterAndMinValue[ idValue.first ] )
                                                             ? idValue.second : idOnClusterAndMinValue[ idValue.first ];

                        else
                            idOnClusterAndMinValue.insert( idValue );
                    }

            allIdOnClusterAndValue.clear(); // the info is treated, this vector is not needed anymore
            idOnClusterAndValue.clear();


            // all proc get a copy of the id and the good min value
            mpi::broadcast( Environment::worldComm().globalComm(),
                            idOnClusterAndMinValue, 0);



            // all the elements which are on the proc (ghost or not) have to be set to the min value
            for (auto const& idValue : idOnClusterAndMinValue)
            {
                size_type locId=invalid_size_type_value;

                if ( (M_spaceP1->dof()->dofGlobalClusterIsOnProc( idValue.first )) )
                    locId = clusterToProcessor(idValue.first);

                else if ( ghostClusterToProc.count(idValue.first) )
                    locId = ghostClusterToProc[ idValue.first ];

                if ( locId != invalid_size_type_value )
                    (*shape)(locId) = idValue.second;
            }

        } // reduceDistanceFunction



#if defined( DISTANCE_FROM_UNORDERED_POINTS )

    template <class TFilename>
    void readPointsFromFile( TFilename filename )
        {
            // fill allPoints with the content of the file: filename
            // for now, the format has to be x y separated by spaces (could be changed in the future to support more formats)
            std::ifstream ifile ( filename.c_str(), std::ios_base::in );
            CHECK( ifile.good() )<<"file "<<filename<<" not ok to be read\n";

            // load all the file in a string
            std::string file_str;
            char c;
            file_str.clear();
            while( ifile.get(c) )
                file_str.append(1,c);
            ifile.close();

            boost::regex reg( "^\\s*"
                              "(\\d+\\.\\d+|\\d+\\.|\\d+)"
                              "\\s+(\\d+\\.\\d+|\\d+\\.|\\d+)\\s*$" );

            auto start = file_str.begin();
            auto end = file_str.end();
            boost::match_results< decltype( start ) > what;
            //boost::match_results< boost::smatch > what;

            while( boost::regex_search(start, end, what, reg) )
            {
                double newX = boost::lexical_cast< double >( std::string(what[1].first, what[1].second) );
                double newY = boost::lexical_cast< double >( std::string(what[2].first, what[2].second) );
                start = what[0].second;

                node_type newNd(dim);
                newNd(0)=newX; newNd(1) = newY;
                allPoints.push_back( newNd );
            }

            LOG(INFO)<<"file "<<filename
                     <<"read, "<<allPoints.size()
                     <<" points added to describe the curve\n";


            boost::regex regref( "^r"
                                 "(\\d+\\.\\d+|\\d+\\.|\\d+)"
                                 "\\s+(\\d+\\.\\d+|\\d+\\.|\\d+)\\s*$" );

            start = file_str.begin();
            end = file_str.end();
            while( boost::regex_search(start, end, what, regref) )
            {
                double newX = boost::lexical_cast< double >( std::string(what[1].first, what[1].second) );
                double newY = boost::lexical_cast< double >( std::string(what[2].first, what[2].second) );
                start = what[0].second;

                node_type newNd(dim);
                newNd(0)=newX; newNd(1) = newY;
                refPoints.push_back( newNd );
            }

            LOG(INFO)<<"file "<<filename
                     <<"read, "<<refPoints.size()
                     <<" reference points added to describe the sign of the curve\n";


        } // readPointsFromFile

    void locateElementsCrossedByUnorderedPoints()
        {
            CHECK( ! allPoints.empty() )<<"No points present\n";


            // first create the map containing the number of the element containing each dof
            std::map< size_type, std::set<size_type> > eltsContainingDof;
            auto it = M_mesh->beginElementWithProcessId();
            auto en = M_mesh->endElementWithProcessId();
            for (; it != en; it++)
                for (uint16_type j = 0; j < ndofv; ++j )
                {
                    const size_type index = M_spaceP1->dof()->localToGlobal( *it, j, 0).index();
                    if ( eltsContainingDof.count( index ) )
                        eltsContainingDof[ index ].insert( it->id() );
                    else
                    {
                        std::set< size_type > newindex( {it->id()} );
                        eltsContainingDof[index] = newindex;
                    }
                }


            // localize the points in a narrow band (only on the elements crossed by the points)
            auto ctx = M_spaceP0->context();
            for (auto const& nd : allPoints)
                ctx.add( nd );

            auto allIndexes = ids.evaluate( ctx );

            const int nbPtContext = ctx.nPoints();
            for (size_type i=0; i < nbPtContext; ++i )
                if (Environment::worldComm().localRank() == ctx.processorHavingPoint( i ) )
                {
                    const size_type index = allIndexes(i);
                    if ( pointsAtIndex.count( index ) )
                        pointsAtIndex[ index ].insert( i );
                    else
                    {
                        std::set< size_type > newset( {i} );
                        pointsAtIndex[ index ] = newset;
                    }
                }



            // make sure there are no elements containing less than 2 points
            auto it_pt = pointsAtIndex.begin();
            auto en_pt = pointsAtIndex.end();
            auto eltHavingPoints = M_spaceP0->element();
            std::set< size_type > dofInNarrowBand;

            for( ; it_pt != en_pt; )
            {
                if ( it_pt->second.size() < 2 )
                    it_pt = pointsAtIndex.erase(it_pt);
                else
                {
                    eltHavingPoints.assign(it_pt->first, 0, 0, 1);
                    for (uint16_type j = 0; j < ndofv; ++j )
                    {
                        const size_type index = M_spaceP1->dof()->localToGlobal(it_pt->first, j, 0).index();
                        dofInNarrowBand.insert( index );
                    }

                    ++it_pt;
                }

            }


            // add the elements sharing a dof with the narrow band
            auto widenBand = M_spaceP0->element();
            for (size_type k : dofInNarrowBand )
                for (size_type eltContDof : eltsContainingDof[ k ] )
                {
                    const size_type idx_elt = M_spaceP0->dof()->localToGlobal(eltContDof, 0, 0).index();
                    // if the points was not in the narrow band
                    if ( ! eltHavingPoints[ idx_elt ] )
                    {
                        // add it in the wide band
                        widenBand[ idx_elt ] = 1;

                        // add the points of all the elts in the narrow band sharing the dof
                        for ( size_type elt_sharing_dof : eltsContainingDof[ k ] )
                        {
                            const size_type idx_neighbour = M_spaceP0->dof()->localToGlobal(elt_sharing_dof, 0, 0).index();
                            if ( eltHavingPoints[ idx_neighbour ] )
                            {
                                if ( pointsAtIndex.count(eltContDof) )
                                    pointsAtIndex[eltContDof].insert( pointsAtIndex[elt_sharing_dof].begin(),
                                                                      pointsAtIndex[elt_sharing_dof].end() );
                                else
                                    pointsAtIndex[eltContDof] = pointsAtIndex[elt_sharing_dof];
                            }
                        }
                    }
                }

            eltHavingPoints = vf::project(M_spaceP0, elements(M_mesh),
                                          vf::chi( idv(eltHavingPoints) + idv(widenBand) ) );

            M_mesh->updateMarker2( eltHavingPoints );

            // control loop, every element should have at least two points associated to
            for( auto const & ptsAtId : pointsAtIndex)
                CHECK( ptsAtId.second.size() > 1 )<<"The element at index : "
                                                  <<ptsAtId.first
                                                  <<" has only "
                                                  <<ptsAtId.second.size()
                                                  <<" associated points\n";

        } // locateElementsCrossedByUnorderedPoints


    // Make an unsigned distance function from a set of unordered points.
    element_ptrtype makeDistanceFunctionSequentialFromUnorderedPoints()
        {
            CHECK( dim == 2 )<<"works only in 2d for now\n";

            auto shape = M_spaceP1->elementPtr();
            *shape = vf::project(M_spaceP1, elements(M_mesh), cst(bigdouble) );

            auto it_elt = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).first;
            auto en_elt = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).second;

            for(; it_elt!=en_elt; it_elt++)
            {

                node_type refVec(dim);
                int refSign=0;

                for (int j=0; j<ndofv; ++j)
                {
                    const size_type indexGlobDof = M_spaceP1->dof()->localToGlobal(it_elt->id(), j, 0).index();

                    const node_type dofCoord = M_spaceP1->dof()->dofPoint( indexGlobDof ).get<0>();

                    // store a point and its dist to the considered dof
                    std::vector< nodeDist_type > pointsDistToDof;

                    // compute the distance to the considered dof
                    for (size_type const& idPoint : pointsAtIndex[it_elt->id()] )
                        pointsDistToDof.push_back( { allPoints[ idPoint ],
                                    squareDistToPoint(allPoints[ idPoint ], dofCoord) } );

                    CHECK( pointsDistToDof.size() > 1 )<<"need more than one point in the element\n";

                    // sort the points
                    std::sort( pointsDistToDof.begin(), pointsDistToDof.end(),
                               []( nodeDist_type a, nodeDist_type b){return a.second < b.second;} );

                    const node_type v = pointsDistToDof[0].first - dofCoord;

                    // the reference vector is calculated once by element with respect to the 2 closest points of the first dof
                    if (j==0)
                    {
                        refVec = pointsDistToDof[0].first - pointsDistToDof[1].first;
                        refSign = (refVec[0] * v[1] - refVec[1] * v[0]) > 0 ? 1 : -1;
                        dofIsOnSide[it_elt->id()][0] = sideA;
                    }
                    else
                    {
                        const int signProdVec = (refVec[0] * v[1] - refVec[1] * v[0]) > 0 ? 1 : -1;
                        const side_type side = ( signProdVec == refSign ) ? sideA : sideB;
                        dofIsOnSide[ it_elt->id() ][j] = side;
                    }

                    const double mindist = std::sqrt( pointsDistToDof[0].second );
                    if ( mindist < (*shape)(indexGlobDof) )
                        (*shape)(indexGlobDof) = mindist;
                }
            }


            return shape;

        } //makeDistanceFunctionSequentialFromUnorderedPoints




    // listPoints could be a any container that can be iterated and its elements have operator [0] to [dim]
    // ex vector< node_type > or vector< vector<double> > ...
    template< class TListPoints >
    void
    setInnerRegion(element_ptrtype phi, TListPoints listPoints, bool exportInnerPoints)
        {
            // need the information : for each global index of the marked dofs, the list of the marked elements in which it appears
            std::map< size_type, std::set< size_type > > eltsAtGlobalIndex;
            std::set< size_type > markedDof;

            auto it_marked = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).first;
            const auto en_marked = M_mesh->elementsWithMarker2(1, M_mesh->worldComm().localRank()).second;
            auto checkmarked = backend()->newVector( M_spaceP1 );

            for (; it_marked != en_marked; it_marked++)
                for (uint16_type j = 0; j < ndofv; ++j )
                {
                    const size_type index = M_spaceP1->dof()->localToGlobal( *it_marked, j, 0).index();
                    markedDof.insert( index );
                    checkmarked->add(index, 1);
                    if ( eltsAtGlobalIndex.count( index ) )
                        eltsAtGlobalIndex[ index ].insert( it_marked->id() );
                    else
                    {
                        std::set< size_type > newindex( {it_marked->id()} );
                        eltsAtGlobalIndex[index] = newindex;
                    }
                }

            checkmarked->close();
            for ( size_type k=0; k<M_spaceP1->nLocalDof(); ++k )
                if( (*checkmarked)(k) > 0 )
                    markedDof.insert( k );


            // for all the marked dofs, store the values of the dofs being on a inter-process boundary -> store it in: isOnInterProcessBoundary
            std::set< size_type > isOnInterProcessBoundary;
            if (Environment::worldComm().size()>1)
            {
                auto checkGhost = backend()->newVector( M_spaceP1 );
                for ( size_type k=0; k<M_spaceP1->nLocalDof(); ++k )
                    if ( M_spaceP1->dof()->dofGlobalProcessIsGhost( k ) )
                        checkGhost->add( k, 1 );
                checkGhost->close();

                for ( size_type k : markedDof )
                    if ( (*checkGhost)(k) > 0)
                        isOnInterProcessBoundary.insert( k );
            }



            std::set< size_type > eltsTodo;
            std::set< size_type > eltsDone;
            std::map< size_type, int > globalClusterDofDone;

            // first, localize the points which are DONE from the given ones
            std::set< size_type > dofDONE;

            // add the points contained in the file is they exist
            listPoints.insert( listPoints.end(), refPoints.begin(), refPoints.end() );

            // just used to know if the point is on this proc
            std::ofstream nodeFile;
            if (exportInnerPoints && (Environment::worldComm().rank() == 0) )
                nodeFile.open( "innerpoints.particles", std::ofstream::out);

            auto ctx = M_spaceP0->context();
            for(auto const& pt : listPoints)
            {
                if (exportInnerPoints && (Environment::worldComm().rank() == 0) )
                    nodeFile << pt[0] << "," << pt[1] << ","<< "0" <<std::endl;
                ctx.add( pt );
            }

            if (exportInnerPoints && (Environment::worldComm().rank() == 0) )
                nodeFile.close();


            for (int i=0; i<ctx.nPoints(); ++i)
            {
                if (Environment::worldComm().localRank() != ctx.processorHavingPoint(i))
                    continue;

                const node_type pt = listPoints[i];

                double minDist = bigdouble;
                size_type minDof = 0;

                for (size_type k : markedDof)
                {
                    const node_type dofCoord = M_spaceP1->dof()->dofPoint( k ).get<0>();
                    const double dist = squareDistToPoint(dofCoord, pt);
                    if (dist < minDist)
                    {
                        minDist = dist;
                        minDof = k;
                    }
                }

                dofDONE.insert( minDof );
            }

            // make sure at least one proc having marked elements has a point setting inner region
            bool okToInitialize = (markedDof.size() != 0 ) && (ctx.nPoints() != 0);
            okToInitialize = mpi::all_reduce( Environment::worldComm(), okToInitialize, std::logical_or<bool>() );
            CHECK( okToInitialize ) << "There is no partition which has at the same time marked elements (around the interface) and points setting the inner region. Consider adding some points.\n";

            std::cout<<"localized elements, dofDONE.size = "<<dofDONE.size()<<std::endl;

            auto doElement = [&](size_type eltId, size_type globIndex)
                {
                    /*
                     find in the element eltId the side of the dof globIndex
                     put all the dof of the element not done to the good side and make them DONE
                     put the elements in which they appear to TODO if not already DONE
                     */
                    if ( eltsDone.count( eltId ) )
                    {
                        if (eltsTodo.count(eltId))
                            eltsTodo.erase(eltId);
                        return;
                    }

                    eltsDone.insert( eltId );

                    side_type sideRef=noSide;
                    // search the side of the element DONE
                    size_type i=0;
                    for (; i<ndofv; ++i)
                        if (M_spaceP1->dof()->localToGlobal( eltId, i, 0).index() == globIndex )
                        {
                            sideRef = dofIsOnSide[ eltId ][i];
                            break;
                        }

                    CHECK( sideRef != noSide ) << "proc " << Environment::worldComm().rank()
                    << ", eltId : "<< eltId << ", has not been able to retrive the sides of the dofs inside.\n";

                    const int signRef = (*phi)(globIndex) > 0 ? 1 : -1;

                    for (int j=0; j<ndofv; ++j)
                    {
                        if (j == i)
                            continue;

                        const size_type indexTodo = M_spaceP1->dof()->localToGlobal( eltId, j, 0).index();

                        if ( dofDONE.count( indexTodo ) )
                            continue;

                        else
                        {
                            const side_type side = dofIsOnSide[ eltId ][j];
                            (*phi)(indexTodo) = std::abs((*phi)(indexTodo)) * ((side == sideRef) ? signRef : -signRef);
                            dofDONE.insert( indexTodo );
                            if (isOnInterProcessBoundary.count(indexTodo))
                                globalClusterDofDone[ processorToCluster( indexTodo ) ] = (*phi)(indexTodo) > 0 ? 1 : -1;

                            // insert the elements having this dof to the list of next element todo (if not already done)
                            for (const size_type& eltCandidate : eltsAtGlobalIndex[indexTodo] )
                                if (! eltsDone.count( eltCandidate ) )
                                    eltsTodo.insert( eltCandidate );
                        }
                    }

                    eltsTodo.erase(eltId);

                }; //doElement


            auto communicateDonePointsOnBoundary = [&]()
                {
                    std::vector< std::map< size_type, int > > all_globalClusterDofDone( Environment::worldComm().size() );
                    mpi::all_gather(Environment::worldComm(), globalClusterDofDone, all_globalClusterDofDone);

                    CHECK( all_globalClusterDofDone.size() == Environment::worldComm().size()) <<"some globlClusterDofDONE might be empty\n";

                    for (int i=0; i<Environment::worldComm().size(); ++i)
                    {
                        if (i==Environment::worldComm().rank())
                            continue;

                        for (std::pair< size_type, int > const& eltSign : all_globalClusterDofDone[i])
                        {
                            size_type index = invalid_size_type_value;
                            if (M_spaceP1->dof()->dofGlobalClusterIsOnProc( eltSign.first ))
                                index = clusterToProcessor( eltSign.first );

                            else if (M_spaceP1->dof()->dofGlobalProcessIsGhost(eltSign.first))
                                index = ghostClusterToProc[ eltSign.first ];

                            if ( index != invalid_size_type_value ) // if the dof is on the proc or ghost
                                if ( !dofDONE.count(index) )
                                {
                                    (*phi)( index ) = eltSign.second * std::abs((*phi)(index));
                                    dofDONE.insert( index );
                                    for (size_type elts : eltsAtGlobalIndex[ index ])
                                        doElement( elts, index );
                                }
                        }
                    }

                    globalClusterDofDone.clear();
                }; //communicateDonePointsOnBoundary

            // in the following loop, dofDONE is modified in doElement and the iterator range would be modified if iterated on it
            auto dofDONEcopy = dofDONE;
            // the set dofDONE contains the dof we are sure are negative. Use them to initialize the loop
            for (size_type dd : dofDONEcopy)
            {
                (*phi)(dd) *= -1;

                if (isOnInterProcessBoundary.count(dd))
                    globalClusterDofDone[ processorToCluster( dd ) ] = -1;

                for (size_type elts : eltsAtGlobalIndex[ dd ])
                    doElement( elts, dd );
            }

            std::cout<<"proc "<<Environment::worldComm().rank()
                     <<" finished doing first dofDONE\n";

            int nbEltTodoGlobal = mpi::all_reduce( Environment::worldComm(), eltsTodo.size(), std::plus<int>() );
            while( nbEltTodoGlobal != 0)
            {
                if (Environment::worldComm().size()>1)
                    communicateDonePointsOnBoundary();

                if ( ! eltsTodo.empty() )
                {
                    // in each element, find the dof DONE and do the other one thanks to it
                    const size_type elt = *eltsTodo.begin();
                    for (int j=0; j<ndofv; ++j)
                    {
                        const size_type index = M_spaceP1->dof()->localToGlobal( elt, j, 0).index();

                        if ( dofDONE.count( index ) )
                        {
                            doElement( elt, index );
                            break;
                        }
                    }
                }

                nbEltTodoGlobal = mpi::all_reduce( Environment::worldComm(), eltsTodo.size(), std::plus<int>() );
            }

            std::cout<<"finished, dofDONE.size = "<<dofDONE.size()<<std::endl;
            std::cout<<"markedDof.size = "<<markedDof.size()<<std::endl;

        }//setInnerRegion

#endif




}; //DistToCurve


template< typename FunctionSpaceP0Type, typename FunctionSpaceP1Type >
boost::shared_ptr< DistToCurve< FunctionSpaceP0Type, FunctionSpaceP1Type > >
distToCurve( boost::shared_ptr<FunctionSpaceP0Type> spaceP0, boost::shared_ptr<FunctionSpaceP1Type> spaceP1 )
{
    auto dtc = DistToCurve< FunctionSpaceP0Type, FunctionSpaceP1Type >::New( spaceP0, spaceP1 );
    return dtc;
}


} //namespace Feel

#endif // DISTTOCURVE_HPP
