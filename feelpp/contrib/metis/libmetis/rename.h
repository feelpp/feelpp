/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h 13933 2013-03-29 22:20:46Z karypis $
 *
 */


#ifndef _LIBMETIS_RENAME_H_
#define _LIBMETIS_RENAME_H_


/* balance.c */
#define Balance2Way			feel__libmetis__Balance2Way
#define Bnd2WayBalance			feel__libmetis__Bnd2WayBalance
#define General2WayBalance		feel__libmetis__General2WayBalance
#define McGeneral2WayBalance            feel__libmetis__McGeneral2WayBalance

/* bucketsort.c */
#define BucketSortKeysInc		feel__libmetis__BucketSortKeysInc

/* checkgraph.c */
#define CheckGraph                      feel__libmetis__CheckGraph
#define CheckInputGraphWeights          feel__libmetis__CheckInputGraphWeights
#define FixGraph                        feel__libmetis__FixGraph

/* coarsen.c */
#define CoarsenGraph			feel__libmetis__CoarsenGraph
#define Match_RM                        feel__libmetis__Match_RM
#define Match_SHEM                      feel__libmetis__Match_SHEM
#define Match_2Hop                      feel__libmetis__Match_2Hop
#define Match_2HopAny                   feel__libmetis__Match_2HopAny
#define Match_2HopAll                   feel__libmetis__Match_2HopAll
#define PrintCGraphStats                feel__libmetis__PrintCGraphStats
#define CreateCoarseGraph		feel__libmetis__CreateCoarseGraph
#define CreateCoarseGraphNoMask		feel__libmetis__CreateCoarseGraphNoMask
#define CreateCoarseGraphPerm		feel__libmetis__CreateCoarseGraphPerm
#define SetupCoarseGraph		feel__libmetis__SetupCoarseGraph
#define ReAdjustMemory			feel__libmetis__ReAdjustMemory

/* compress.c */
#define CompressGraph			feel__libmetis__CompressGraph
#define PruneGraph			feel__libmetis__PruneGraph

/* contig.c */
#define FindPartitionInducedComponents  feel__libmetis__FindPartitionInducedComponents   
#define IsConnected                     feel__libmetis__IsConnected
#define IsConnectedSubdomain            feel__libmetis__IsConnectedSubdomain
#define FindSepInducedComponents        feel__libmetis__FindSepInducedComponents
#define EliminateComponents             feel__libmetis__EliminateComponents
#define MoveGroupContigForCut           feel__libmetis__MoveGroupContigForCut
#define MoveGroupContigForVol           feel__libmetis__MoveGroupContigForVol

/* debug.c */
#define ComputeCut			feel__libmetis__ComputeCut
#define ComputeVolume			feel__libmetis__ComputeVolume
#define ComputeMaxCut			feel__libmetis__ComputeMaxCut
#define CheckBnd			feel__libmetis__CheckBnd
#define CheckBnd2			feel__libmetis__CheckBnd2
#define CheckNodeBnd			feel__libmetis__CheckNodeBnd
#define CheckRInfo			feel__libmetis__CheckRInfo
#define CheckNodePartitionParams	feel__libmetis__CheckNodePartitionParams
#define IsSeparable			feel__libmetis__IsSeparable
#define CheckKWayVolPartitionParams     feel__libmetis__CheckKWayVolPartitionParams

/* fm.c */
#define FM_2WayRefine                   feel__libmetis__FM_2WayRefine
#define FM_2WayCutRefine                feel__libmetis__FM_2WayCutRefine
#define FM_Mc2WayCutRefine              feel__libmetis__FM_Mc2WayCutRefine
#define SelectQueue                     feel__libmetis__SelectQueue
#define Print2WayRefineStats            feel__libmetis__Print2WayRefineStats

/* fortran.c */
#define Change2CNumbering		feel__libmetis__Change2CNumbering
#define Change2FNumbering		feel__libmetis__Change2FNumbering
#define Change2FNumbering2		feel__libmetis__Change2FNumbering2
#define Change2FNumberingOrder		feel__libmetis__Change2FNumberingOrder
#define ChangeMesh2CNumbering		feel__libmetis__ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		feel__libmetis__ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		feel__libmetis__ChangeMesh2FNumbering2

/* graph.c */
#define SetupGraph			feel__libmetis__SetupGraph
#define SetupGraph_adjrsum              feel__libmetis__SetupGraph_adjrsum
#define SetupGraph_tvwgt                feel__libmetis__SetupGraph_tvwgt
#define SetupGraph_label                feel__libmetis__SetupGraph_label
#define SetupSplitGraph                 feel__libmetis__SetupSplitGraph
#define CreateGraph                     feel__libmetis__CreateGraph
#define InitGraph                       feel__libmetis__InitGraph
#define FreeRData                       feel__libmetis__FreeRData
#define FreeGraph                       feel__libmetis__FreeGraph

/* initpart.c */
#define Init2WayPartition		feel__libmetis__Init2WayPartition
#define InitSeparator			feel__libmetis__InitSeparator
#define RandomBisection			feel__libmetis__RandomBisection
#define GrowBisection			feel__libmetis__GrowBisection
#define McRandomBisection               feel__libmetis__McRandomBisection
#define McGrowBisection                 feel__libmetis__McGrowBisection
#define GrowBisectionNode		feel__libmetis__GrowBisectionNode

/* kmetis.c */
#define MlevelKWayPartitioning		feel__libmetis__MlevelKWayPartitioning
#define InitKWayPartitioning            feel__libmetis__InitKWayPartitioning

/* kwayfm.c */
#define Greedy_KWayOptimize		feel__libmetis__Greedy_KWayOptimize
#define Greedy_KWayCutOptimize		feel__libmetis__Greedy_KWayCutOptimize
#define Greedy_KWayVolOptimize          feel__libmetis__Greedy_KWayVolOptimize
#define Greedy_McKWayCutOptimize        feel__libmetis__Greedy_McKWayCutOptimize
#define Greedy_McKWayVolOptimize        feel__libmetis__Greedy_McKWayVolOptimize
#define IsArticulationNode              feel__libmetis__IsArticulationNode
#define KWayVolUpdate                   feel__libmetis__KWayVolUpdate

/* kwayrefine.c */
#define RefineKWay			feel__libmetis__RefineKWay
#define AllocateKWayPartitionMemory	feel__libmetis__AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	feel__libmetis__ComputeKWayPartitionParams
#define ProjectKWayPartition		feel__libmetis__ProjectKWayPartition
#define ComputeKWayBoundary		feel__libmetis__ComputeKWayBoundary
#define ComputeKWayVolGains             feel__libmetis__ComputeKWayVolGains
#define IsBalanced			feel__libmetis__IsBalanced

/* mcutil */
#define rvecle                          feel__libmetis__rvecle
#define rvecge                          feel__libmetis__rvecge
#define rvecsumle                       feel__libmetis__rvecsumle
#define rvecmaxdiff                     feel__libmetis__rvecmaxdiff
#define ivecle                          feel__libmetis__ivecle
#define ivecge                          feel__libmetis__ivecge
#define ivecaxpylez                     feel__libmetis__ivecaxpylez
#define ivecaxpygez                     feel__libmetis__ivecaxpygez
#define BetterVBalance                  feel__libmetis__BetterVBalance
#define BetterBalance2Way               feel__libmetis__BetterBalance2Way
#define BetterBalanceKWay               feel__libmetis__BetterBalanceKWay
#define ComputeLoadImbalance            feel__libmetis__ComputeLoadImbalance
#define ComputeLoadImbalanceDiff        feel__libmetis__ComputeLoadImbalanceDiff
#define ComputeLoadImbalanceDiffVec     feel__libmetis__ComputeLoadImbalanceDiffVec
#define ComputeLoadImbalanceVec         feel__libmetis__ComputeLoadImbalanceVec

/* mesh.c */
#define CreateGraphDual                 feel__libmetis__CreateGraphDual
#define FindCommonElements              feel__libmetis__FindCommonElements
#define CreateGraphNodal                feel__libmetis__CreateGraphNodal
#define FindCommonNodes                 feel__libmetis__FindCommonNodes
#define CreateMesh                      feel__libmetis__CreateMesh
#define InitMesh                        feel__libmetis__InitMesh
#define FreeMesh                        feel__libmetis__FreeMesh

/* meshpart.c */
#define InduceRowPartFromColumnPart     feel__libmetis__InduceRowPartFromColumnPart

/* minconn.c */
#define ComputeSubDomainGraph           feel__libmetis__ComputeSubDomainGraph
#define UpdateEdgeSubDomainGraph        feel__libmetis__UpdateEdgeSubDomainGraph
#define PrintSubDomainGraph             feel__libmetis__PrintSubDomainGraph
#define EliminateSubDomainEdges         feel__libmetis__EliminateSubDomainEdges
#define MoveGroupMinConnForCut          feel__libmetis__MoveGroupMinConnForCut
#define MoveGroupMinConnForVol          feel__libmetis__MoveGroupMinConnForVol

/* mincover.c */
#define MinCover			feel__libmetis__MinCover
#define MinCover_Augment		feel__libmetis__MinCover_Augment
#define MinCover_Decompose		feel__libmetis__MinCover_Decompose
#define MinCover_ColDFS			feel__libmetis__MinCover_ColDFS
#define MinCover_RowDFS			feel__libmetis__MinCover_RowDFS

/* mmd.c */
#define genmmd				feel__libmetis__genmmd
#define mmdelm				feel__libmetis__mmdelm
#define mmdint				feel__libmetis__mmdint
#define mmdnum				feel__libmetis__mmdnum
#define mmdupd				feel__libmetis__mmdupd


/* ometis.c */
#define MlevelNestedDissection		feel__libmetis__MlevelNestedDissection
#define MlevelNestedDissectionCC	feel__libmetis__MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	feel__libmetis__MlevelNodeBisectionMultiple
#define MlevelNodeBisectionL2		feel__libmetis__MlevelNodeBisectionL2
#define MlevelNodeBisectionL1		feel__libmetis__MlevelNodeBisectionL1
#define SplitGraphOrder			feel__libmetis__SplitGraphOrder
#define SplitGraphOrderCC		feel__libmetis__SplitGraphOrderCC
#define MMDOrder			feel__libmetis__MMDOrder

/* options.c */
#define SetupCtrl                       feel__libmetis__SetupCtrl
#define SetupKWayBalMultipliers         feel__libmetis__SetupKWayBalMultipliers
#define Setup2WayBalMultipliers         feel__libmetis__Setup2WayBalMultipliers
#define PrintCtrl                       feel__libmetis__PrintCtrl
#define FreeCtrl                        feel__libmetis__FreeCtrl
#define CheckParams                     feel__libmetis__CheckParams

/* parmetis.c */
#define MlevelNestedDissectionP		feel__libmetis__MlevelNestedDissectionP
#define FM_2WayNodeRefine1SidedP        feel__libmetis__FM_2WayNodeRefine1SidedP
#define FM_2WayNodeRefine2SidedP        feel__libmetis__FM_2WayNodeRefine2SidedP

/* pmetis.c */
#define MlevelRecursiveBisection	feel__libmetis__MlevelRecursiveBisection
#define MultilevelBisect		feel__libmetis__MultilevelBisect
#define SplitGraphPart			feel__libmetis__SplitGraphPart

/* refine.c */
#define Refine2Way			feel__libmetis__Refine2Way
#define Allocate2WayPartitionMemory	feel__libmetis__Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	feel__libmetis__Compute2WayPartitionParams
#define Project2WayPartition		feel__libmetis__Project2WayPartition

/* separator.c */
#define ConstructSeparator		feel__libmetis__ConstructSeparator
#define ConstructMinCoverSeparator	feel__libmetis__ConstructMinCoverSeparator

/* sfm.c */
#define FM_2WayNodeRefine2Sided         feel__libmetis__FM_2WayNodeRefine2Sided 
#define FM_2WayNodeRefine1Sided         feel__libmetis__FM_2WayNodeRefine1Sided
#define FM_2WayNodeBalance              feel__libmetis__FM_2WayNodeBalance

/* srefine.c */
#define Refine2WayNode			feel__libmetis__Refine2WayNode
#define Allocate2WayNodePartitionMemory	feel__libmetis__Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	feel__libmetis__Compute2WayNodePartitionParams
#define Project2WayNodePartition	feel__libmetis__Project2WayNodePartition

/* stat.c */
#define ComputePartitionInfoBipartite   feel__libmetis__ComputePartitionInfoBipartite
#define ComputePartitionBalance		feel__libmetis__ComputePartitionBalance
#define ComputeElementBalance		feel__libmetis__ComputeElementBalance

/* timing.c */
#define InitTimers			feel__libmetis__InitTimers
#define PrintTimers			feel__libmetis__PrintTimers

/* util.c */
#define iargmax_strd                    feel__libmetis__iargmax_strd 
#define iargmax_nrm                     feel__libmetis__iargmax_nrm
#define iargmax2_nrm                    feel__libmetis__iargmax2_nrm
#define rargmax2                        feel__libmetis__rargmax2
#define InitRandom                      feel__libmetis__InitRandom
#define metis_rcode                     feel__libmetis__metis_rcode

/* wspace.c */
#define AllocateWorkSpace               feel__libmetis__AllocateWorkSpace                  
#define AllocateRefinementWorkSpace     feel__libmetis__AllocateRefinementWorkSpace
#define FreeWorkSpace                   feel__libmetis__FreeWorkSpace
#define wspacemalloc                    feel__libmetis__wspacemalloc
#define wspacepush                      feel__libmetis__wspacepush
#define wspacepop                       feel__libmetis__wspacepop
#define iwspacemalloc                   feel__libmetis__iwspacemalloc
#define rwspacemalloc                   feel__libmetis__rwspacemalloc
#define ikvwspacemalloc                 feel__libmetis__ikvwspacemalloc
#define cnbrpoolReset                   feel__libmetis__cnbrpoolReset
#define cnbrpoolGetNext                 feel__libmetis__cnbrpoolGetNext
#define vnbrpoolReset                   feel__libmetis__vnbrpoolReset
#define vnbrpoolGetNext                 feel__libmetis__vnbrpoolGetNext

#endif


