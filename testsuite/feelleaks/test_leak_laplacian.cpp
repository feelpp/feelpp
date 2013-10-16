// -*- coding: utf-8; mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
#include <feel/feel.hpp>

#if defined(FEELPP_HAS_GPERFTOOLS)
#include <gperftools/heap-checker.h>
#endif /* FEELPP_HAS_GPERFTOOLS */

int main(int argc, char**argv )
{
#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check("checker 1");
    CHECK(check.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */

    //# marker1 #
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="test_leak_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    //# endmarker1 #



    boost::shared_ptr<Mesh<Simplex<2>>> mesh;
    decltype( Pch<2>( mesh ) ) Vh;
    decltype( exporter( _mesh=mesh ) ) e;

#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check0("checker 2");
    CHECK(check0.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */


#if defined(FEELPP_HAS_GPERFTOOLS)
    HeapLeakChecker check3("checker 3");
#endif /* FEELPP_HAS_GPERFTOOLS */
    {
        mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
        Vh = Pch<2>( mesh );
        LOG(INFO) << "Vh.use_count(): " << Vh.use_count();
        auto u = Vh->element();
        auto v = Vh->element();
        LOG(INFO) << "Vh.use_count() after element: " << Vh.use_count();
        //# endmarker2 #

        //# marker3 #
        auto l = form1( _test=Vh );
        l = integrate(_range=elements(mesh),
                      _expr=id(v));
        LOG(INFO) << "Vh.use_count() after form1: " << Vh.use_count();
        auto a = form2( _trial=Vh, _test=Vh);
        a = integrate(_range=elements(mesh),
                      _expr=gradt(u)*trans(grad(v)) );
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
              _expr=expr( option(_name="functions.g").as<std::string>(), symbols<2>() ) );
        a.solve(_rhs=l,_solution=u);
        LOG(INFO) << "Vh.use_count() after form2: " << Vh.use_count();

#if 0
        e = exporter( _mesh=mesh );
        e->add( "u", u );
        e->save();
#endif
        LOG(INFO) << "Destructors being called...";
        LOG(INFO) << "mesh count : " << mesh.use_count();
        LOG(INFO) << "Vh count : " << Vh.use_count();
        LOG(INFO) << "e count : " << e.use_count();
    }
        LOG(INFO) << "after block mesh count : " << mesh.use_count();
        LOG(INFO) << "after block Vh count : " << Vh.use_count();
        LOG(INFO) << "after block e count : " << e.use_count();
        mesh.reset();
        Vh.reset();
        e.reset();
        LOG(INFO) << "Destructors done.";
        CHECK( mesh.use_count() == 0 ) << "Invalid mesh shared_ptr";
        CHECK( Vh.use_count() == 0 ) << "Invalid functionspace shared_ptr";
        CHECK( e.use_count() == 0 ) << "Invalid exporter shared_ptr";

#if defined(FEELPP_HAS_GPERFTOOLS)
    CHECK(check3.NoLeaks()) << "There are leaks";
#endif /* FEELPP_HAS_GPERFTOOLS */

    CHECK( mesh.use_count() == 0 ) << "Invalid mesh shared_ptr";
    CHECK( Vh.use_count() == 0 ) << "Invalid functionspace shared_ptr";
    CHECK( e.use_count() == 0 ) << "Invalid exporter shared_ptr";
}
