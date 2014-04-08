/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Date: 2014-02-13

   Copyright (C) 2014 Feel++ Consortium

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
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="partition",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    /*
    //auto mesh = unitSquare();
    auto dgmesh = loadMesh( _mesh=new Mesh<Hypercube<2>>,
                            _h=option(_name="gmsh.hsize2").as<double>(),
                            _partitioner=option(_name="gmsh.partitioner").as<int>() );

    // auto e = exporter( _mesh=dgmesh );
    // e->step(0)->setMesh( dgmesh );
    // e->save();

    auto Vh = Pdh<1>( dgmesh, true );
    auto Xh = Pch<1>( dgmesh );

    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(dgmesh),
                  _expr=id(v));

    auto dgform = form2( _trial=Vh, _test=Vh,
                         _pattern=size_type(Pattern::EXTENDED) );
    dgform = integrate(_range=elements(dgmesh),
                       _expr=gradt(u)*trans(grad(v)) );
    // dgform +=integrate( internalfaces( dgmesh ),
    //                     + trans( jumpt( cst(1.0)/4 ) )*jump( cst( 1.0 )/4 ) / (measFace()) );

    dgform +=integrate( internalfaces( dgmesh ),
                        // - {grad(u)} . [v]
                        -averaget( gradt( u ) )*jump( id( v ) )
                        // - [u] . {grad(v)}
                        -average( grad( v ) )*jumpt( idt( u ) )
                        // penal*[u] . [v]/h_face
                        + 50* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace() );

    dgform += on( _range=boundaryfaces(dgmesh), _element=u, _rhs=l, _expr=cst(0.));

    dgform.solve(_rhs=l,_solution=u);

    // if (Environment::worldComm().rank()==0)
    //     u.printMatlab("dgsol.m");

    //auto dg_backend = backend();
    //auto p = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2, _backend=dg_backend );
    auto p = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2 );
    auto uc = p->project( idv(u) );

    auto er = exporter( _mesh=dgmesh );
    //er->step(0)->setMesh( dgmesh );
    er->step(0)->add( "u", u );
    er->step(0)->add( "uc", uc );
    er->save();

    // std::vector<int> vec = {0,1,2,0};
    // //auto subworldcomm = Environment::worldComm().subWorldComm(vec);
    // auto subworldcomm = Environment::worldComm();
    // subworldcomm.setIsActive(vec);
    // std::cout<<"worldcomm.size()= "<< subworldcomm.size() <<"\n";

    // subworldcomm.showMe();
    */

    //mpi::communicator world;
    //if(world.size() < 2)
    //    throw std::runtime_error("Please run with at least 2 MPI processes!");

    //int group_a_ranks[4] = {1,2,4,6};
    std::vector<int> group_a_ranks = {1,2,4,6};
    mpi::group world_group = Environment::worldComm().group();
    //mpi::group group_a = world_group.include(group_a_ranks,group_a_ranks+4);
    mpi::group group_a = world_group.include(group_a_ranks.begin(),group_a_ranks.end());
    //mpi::communicator comm_a(world,group_a);
    //mpi::communicator comm_b(Environment::worldComm(),group_a);
    mpi::communicator comm_b(Environment::worldComm(),group_a);

    //std::vector<int> newIsActive(Environment::worldComm().godSize(),false);
    //newIsActive[Environment::worldComm().godRank()]=true;
    //std::vector<int> newIsActive = {true,false,true,false};
    std::vector<int> newIsActive = {true,true,true,true};

    std::string value("Hello world!");
    //boost::shared_ptr<WorldComm> commSeq;

    if ( comm_b != MPI_COMM_NULL )
    {
        auto commSeq = WorldComm( comm_b,
                                  comm_b,
                                  comm_b,
                                  //Environment::worldComm().godComm(),
                                  //Environment::worldComm().godRank(), // local color
                                  comm_b.rank(), // local color
                                  newIsActive );


        commSeq.showMe();

        if (commSeq.rank()==0)
        {
            std::cout<<"commSeq.size= "<< commSeq.size() <<"\n";
            std::cout<<"commSeq.localSize= "<< commSeq.localSize() <<"\n";
            std::cout<<"commSeq.globalSize= "<< commSeq.globalSize() <<"\n";
            std::cout<<"commSeq.godSize= "<< commSeq.godSize() <<"\n";
        }

        std::cout<<"PDG: commSeq.godRank()= "<< commSeq.godRank() << " and commSeq.globalRank()= "<< commSeq.globalRank() <<"\n";

        // if(commSeq.rank() == 0)
        // {
        //     value = "Hello group a!";
        // }
        // mpi::broadcast(commSeq, value, 0);

        auto dgmesh = createGMSHMesh( _mesh=new Mesh<Hypercube<2> >(commSeq),
                                      _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                      _desc=domain( _name="dgmesh", _addmidpoint=false, _usenames=false, _shape="hypercube",
                                                    _dim=2, _h=option(_name="gmsh.hsize2").as<double>(),
                                                    _convex="Hypercube",
                                                    _xmin=0.,
                                                    _xmax=1.,
                                                    _ymin=0.,
                                                    _ymax=1.,
                                                    _worldcomm=commSeq
                                                    //_group=2
                                                    ),
                                      //_structured=3,
                                      _partitions=commSeq.globalSize(),
                                      _worldcomm=commSeq
                                      );


        // auto e = exporter( _mesh=dgmesh );
        // e->step(0)->setMesh( dgmesh );
        // e->save();

        std::cout<<"DG MESH WORLDCOMM SIZE= "<< dgmesh->worldComm().globalSize() <<"\n";

        auto Vh = Pdh<1>( dgmesh, true );
        auto Xh = Pch<1>( dgmesh );
        std::cout<<"Xh->nDof= "<< Xh->nDof() <<"\n";
        std::cout<<"Xh->nLocalDof= "<< Xh->nLocalDof() <<"\n";
        std::cout<<"Xh->nLocalDofWithoutGhost= "<< Xh->nLocalDofWithoutGhost() <<"\n";

        auto u = Vh->element();
        auto v = Vh->element();

        auto l = form1( _test=Vh );
        l = integrate(_range=elements(dgmesh),
                      _expr=id(v));

        auto dgform = form2( _trial=Vh, _test=Vh,
                             _pattern=size_type(Pattern::EXTENDED) );

        dgform = integrate(_range=elements(dgmesh),
                           _expr=gradt(u)*trans(grad(v)) );

        dgform +=integrate( internalfaces( dgmesh ),
                            // - {grad(u)} . [v]
                            -averaget( gradt( u ) )*jump( id( v ) )
                            // - [u] . {grad(v)}
                            -average( grad( v ) )*jumpt( idt( u ) )
                            // penal*[u] . [v]/h_face
                            + 50* ( trans( jumpt( idt( u ) ) )*jump( id( v ) ) )/hFace() );

        std::cout<<"apply dirichlet boundary conditions starts\n";
        dgform += on( _range=boundaryfaces(dgmesh), _element=u, _rhs=l, _expr=cst(0.));
        std::cout<<"apply dirichlet boundary conditions done\n";

        std::cout<<"solve starts\n";
        dgform.solve(_rhs=l,_solution=u);
        std::cout<<"solve done\n";

        u.printMatlab("u.m");

        std::cout<<"projection starts\n";
        auto dgbackend = backend(_worldcomm=commSeq, _rebuild=true);
        auto p = opProjection( _domainSpace=Xh, _imageSpace=Xh, _type=L2, _backend=dgbackend );
        auto uc = p->project( idv(u) );
        std::cout<<"projection done\n";

        uc.printMatlab("uc.m");

#if 0
        // auto ddofpt_it = Vh->dof()->dof()->dofPointBegin();
        // auto ddofpt_en = Vh->dof()->dof()->dofPointEnd();

        // for ( ; ddofpt_it != ddofpt_en; ++ddofpt_it )
        // {
        //     auto ddofpt_coord = ddofpt_it->get<0>();
        //     auto ddofpt_id = ddofpt_it->get<1>();
        //     if (commSeq.rank()==0) std::cout<< "INNER proc "<< commSeq.rank() <<": dof_id= "<< ddofpt_id << " and dof_coord= " << ddofpt_coord <<"\n";
        // }
#endif

        for( auto const& dof : Vh->dof()->localDof() )
        {
            if (commSeq.rank()==0)
                std::cout << "PROC "<< commSeq.rank()  << " id:" << dof.first.localDof()
                          << " global dof : " << dof.second.index() << " pts: " << Vh->dof()->dofPoint( dof.second.index() ).template get<0>() <<"\n";
        }
        // commSeq.barrier();
        auto er = exporter( _mesh=dgmesh );
        er->step(0)->add( "u", u );
        er->step(0)->add( "uc", uc );
        er->save();
    }

    // Environment::worldComm().barrier();
    //std::cout << "Process #" << Environment::worldComm().rank() << " says " << value << std::endl;
    MPI::Finalize();
    exit(0);
}
