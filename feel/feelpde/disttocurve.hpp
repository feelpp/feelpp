/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Vincent Doyeux <vincent.doyeux@ujf-grenoble.fr>
       Date: 2014-02-02

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
            M_mesh( spaceP1->mesh() ),
            M_dt1(0),
            M_dt2(0)
        {
            // make the map ghostClusterToProc
            // TODO : should be in datamap.hpp
            if (Environment::worldComm().size() > 1)
                for (int k=0; k<M_spaceP1->nLocalDof(); ++k)
                    if ( M_spaceP1->dof()->dofGlobalProcessIsGhost( k ) )
                        ghostClusterToProc[ processorToCluster( k ) ] = k;

            eltHavingPoints = M_spaceP0->element();
        }


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
                                               double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2,
                                               bool exportPoints = false,
                                               std::string exportName=""
                                               )
        {
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


        /* 2d version only */
        element_ptrtype fromParametrizedCurve( std::function<double(double)> xexpr,
                                               std::function<double(double)> yexpr,
                                               double t1Start, double t1End, double dt1,
                                               bool broadenCurveForElementDetection = true,
                                               double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2,
                                               bool exportPoints = false,
                                               std::string exportName=""
                                               )
        {
            generatePointsFromParametrization(xexpr, yexpr,
                                                t1Start, t1End, dt1,
                                                exportPoints, exportName);

            locateElementsCrossedByCurve(broadenCurveForElementDetection, broadenessAmplitude );

            auto shape = makeDistanceFunctionSequential();

            clear(); // no need for the maps used to create the distance function

            if ( Environment::worldComm().size() > 1)
                reduceDistanceFunction( shape );

            return shape;
        }



        // same function as the previous one but extract the functions x(t), y(t) from a tuple
        element_ptrtype fromParametrizedCurve( std::tuple< std::function<double(double)>, std::function<double(double)> > paramFct,
                                               double tStart, double tEnd, double dt,
                                               bool broadenCurveForElementDetection = true,
                                               double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2,
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
                                               double broadenessAmplitude = option("gmsh.hsize").as<double>() / 2,
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
        const size_type firstDof;
        const uint16_type ndofv;

        const int dim;
        int periodT2;

        // function spaces
        spaceP0_ptrtype M_spaceP0;
        spaceP1_ptrtype M_spaceP1;
        typename FunctionSpaceP1Type::mesh_ptrtype M_mesh;

        // store a marker on the last elements having been localized
        typename FunctionSpaceP0Type::element_type eltHavingPoints;

        std::map< size_type, size_type > ghostClusterToProc;

        // contains parameter t, node
        std::map< double, node_type > tNodeMap;
        double M_dt1;
        double M_dt2;
        // contains as key the id of an element and value the numbers of the nodes (t) which are crossing it
        std::map<size_type, std::vector<double> > pointsAtIndex;


        // -------- private methods ---------

        inline size_type clusterToProcessor( size_type dof )
        {return M_spaceP1->dof()->mapGlobalClusterToGlobalProcess( dof - firstDof ); }

        inline size_type processorToCluster( size_type dof )
        {return M_spaceP1->dof()->mapGlobalProcessToGlobalCluster( dof ); }

        // todo : readPointsFromFile

        typename FunctionSpaceP0Type::element_type getCrossedElements()
        {return eltHavingPoints;}


        void clear()
        {
            tNodeMap.clear();
            pointsAtIndex.clear();
            M_dt1=0;
            M_dt2=0;
        }


        /* 2d version */
        template < class TExpr >
        void generatePointsFromParametrization( TExpr xexpr, TExpr yexpr,
                                                double t1Start, double t1End, double dt1,
                                                bool exportPoints = false, std::string exportName="" )
        {
            M_dt1 = dt1;

            std::ofstream nodeFile;

            if (exportPoints && (Environment::worldComm().rank() == 0) )
                {
                    std::string expName = exportName.empty() ? "nodes.particles" : exportName+".particles";
                    nodeFile.open(expName, std::ofstream::out);
                }

            int count=0;

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
            M_dt1 = dt1;
            M_dt2 = dt2;

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




        void locateElementsCrossedByCurve(bool randomlyBroadenNodesPositions=false, double randomnessAmplitude = option("gmsh.hsize").as<double>() / 2 )
        {
            // locate the elements crossed by the curve
            // store their ids in a map with the "t" of the nodes being in the element
            // update the marker2 with the elements being crossed

            CHECK( ! tNodeMap.empty() )<<"\n No nodes defining the curve have been loaded.\n";

            // create a P0 elt containing the ids of the elements
            auto ids = M_spaceP0->element();
            for (auto const& it : elements(M_mesh) )
                ids.assign( it.id(), 0, 0, it.id() );

            auto ctx = M_spaceP0->context();

            std::default_random_engine re( (unsigned int)time(0) );
            std::uniform_real_distribution<double> smallRd( -randomnessAmplitude, randomnessAmplitude );

            for (auto const& tnode : tNodeMap)
                {
                    node_type nodeToAdd = tnode.second;

                    if (randomlyBroadenNodesPositions)
                        for (int i=0; i<dim; ++i)
                            nodeToAdd[i] += smallRd(re);

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
                        const double tOfNode = (*tnodeit).first;
                        const size_type index = allIndexes(i);

                        if ( pointsAtIndex.count( index ) )
                            pointsAtIndex[ index ].push_back( tOfNode );
                        else
                            {
                                std::vector<double> v(1, tOfNode);
                                pointsAtIndex[ index ] = v;
                            }

                        eltHavingPoints.assign( index, 0, 0, 1);
                    }

                tnodeit++;

                }

            M_mesh->updateMarker2( eltHavingPoints );

        }




        element_ptrtype makeDistanceFunctionSequential()
        {
            CHECK( M_dt1 >= 0 )<<"\n Intervall between the points of the curve has to be set (either the dt of the parametrized curve or an integer 1 if curve is read from a file\n";
            if (dim==3)
                CHECK( M_dt2 >= 0 )<<"\n Intervall between the points of the curve has to be set (either the dt of the parametrized curve or an integer 1 if curve is read from a file\n";

            const double bigdouble = 1e8;

            auto shape = M_spaceP1->elementPtr();
            *shape = vf::project(M_spaceP1, elements(M_mesh), cst(bigdouble) );


            // squared distance between a point where only its "t" is given, and a node (x,y)
            auto distToPt = [this] (double t, double x, double y, double z) -> double
                {
                    node_type node = this->tNodeMap[ t ];
                    if (dim==2)
                        return (x - node[0])*(x - node[0]) + (y - node[1])*(y - node[1]);
                    else
                        return (x - node[0])*(x - node[0]) + (y - node[1])*(y - node[1]) + (z - node[2])*(z - node[2]);
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
                        const double xdof = M_spaceP1->dof()->dofPoint( indexGlobDof ).template get<0>()[0];
                        const double ydof = M_spaceP1->dof()->dofPoint( indexGlobDof ).template get<0>()[1];
                        const double zdof = dim==3 ? M_spaceP1->dof()->dofPoint( indexGlobDof ).template get<0>()[2] : 0;

                        // find the point in the element having the closest distance with the dof. This distance will be the value of shape at this dof (if a smaller distance on the same dof is not found in an other element).
                        //This method assumes that the distance between the points of the curve is very small compared to the size of the mesh

                        for (auto const& pt : pointsAtIndex[ it_elt->id() ] )
                            {
                                const double dtp = distToPt( pt, xdof, ydof, zdof );
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
                                double t1x, t1y, t1z, t2x, t2y, t2z;
                                // try to get the point next to the closest point (in the particular case where closest point is the last point, get the previous one)

                                try
                                    {
                                        const node_type closestPointPlusDtCoord = tNodeMap.at(closestPoint + 1.);
                                        t1x = closestPointPlusDtCoord[0] - closestPointCoord[0];
                                        t1y = closestPointPlusDtCoord[1] - closestPointCoord[1];
                                        t1z = (dim==3) ? closestPointPlusDtCoord[2] - closestPointCoord[2] : 0;
                                    }

                                catch (const std::out_of_range& oor)
                                    {
                                        const node_type closestPointMinusDtCoord = tNodeMap.at(closestPoint - 1.);
                                        t1x =  closestPointCoord[0] - closestPointMinusDtCoord[0];
                                        t1y =  closestPointCoord[1] - closestPointMinusDtCoord[1];
                                        t1z =  (dim==3) ? closestPointCoord[2] - closestPointMinusDtCoord[2] : 0;
                                    }

                                if (dim==3)
                                    {
                                        try
                                            {
                                                const node_type closestPointPlusDtCoord = tNodeMap.at(closestPoint + periodT2);
                                                t2x = closestPointPlusDtCoord[0] - closestPointCoord[0];
                                                t2y = closestPointPlusDtCoord[1] - closestPointCoord[1];
                                                t2z = closestPointPlusDtCoord[2] - closestPointCoord[2];
                                            }

                                        catch (const std::out_of_range& oor)
                                            {

                                                const node_type closestPointMinusDtCoord = tNodeMap.at(closestPoint - periodT2);
                                                t2x =  closestPointCoord[0] - closestPointMinusDtCoord[0];
                                                t2y =  closestPointCoord[1] - closestPointMinusDtCoord[1];
                                                t2z =  closestPointCoord[2] - closestPointMinusDtCoord[2];
                                            }
                                    }
                                // (vx, vy) = vector pointing from the closest point on curve to the concerned dof
                                const double vx = xdof - closestPointCoord[0];
                                const double vy = ydof - closestPointCoord[1];
                                const double vz = (dim==3) ? zdof - closestPointCoord[2] : 0;

                                // the sign of the distance function is ruled by the vectorial product of the tangent vector and the vector v : sign(v x t)
                                // in 3D, it should be somthing like :  sign( (v x t) . n ) where n is the normal of the param surface pointing outward
                                const double nx = t1y*t2z-t1z*t2y;
                                const double ny = t1z*t2x-t1x*t2z;
                                const double nz = t1x*t2y-t1y*t2x;

                                double signProdVec;
                                if (dim==2)
                                    signProdVec = vx * t1y - vy * t1x > 0 ? 1 : -1;
                                else
                                    signProdVec = vx * nx + vy * ny + vz * nz > 0 ? 1 : -1;

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
