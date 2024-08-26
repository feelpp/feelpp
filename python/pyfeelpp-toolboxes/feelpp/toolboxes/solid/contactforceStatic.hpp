#include <feel/feelmodels/modelcore/modelnumerical.hpp>
#include <feel/feelmodels/solid/solidmechanics.hpp>
#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/json.hpp>
#include <fmt/core.h>
#include <Eigen/Geometry> 
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelvf/matvec.hpp>
#include <mpi.h>
#include <feel/feelmesh/meshsupportbase.hpp>
#include <feel/feelmesh/traits.hpp>
#include <feel/feelmesh/filters.hpp>
#include <feel/feeldiscr/localization.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/savegmshmesh.hpp>
#include <feel/feelvf/ginac.hpp>
#include <feel/feelvf/detail/ginacmatrix.hpp>
#include <fmt/ostream.h>
#include <feel/feelmesh/bvh.hpp>

using namespace Feel;
using namespace Feel::FeelModels;
using index_type = uint32_type;


template<typename SolidMechanics>
void
initG(SolidMechanics const& t, std::vector<double> const& direction, std::vector<double> const& initCondition)
{   
    // Init the distance fields 
    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    auto g = Vh->element();

    // Load and export rigid obstacles
    std::string filename = "$cfgdir/wall.geo";
    auto wall = loadMesh(_mesh=new Mesh<Simplex<SolidMechanics::nDim>>("Wall"), _filename = filename);

    auto expWall = exporter(_mesh = wall, _name = fmt::format("Wall"));
    expWall->addRegions();
    expWall->save();

    // Raytracing to compute distance
    using bvh_ray_type = BVHRay<SolidMechanics::nRealDim>;
    Eigen::VectorXd origin(SolidMechanics::nDim);
    Eigen::VectorXd dir(SolidMechanics::nDim);

    if constexpr(SolidMechanics::nDim == 2)
        dir << direction[0], direction[1];
    else if constexpr(SolidMechanics::nDim == 3)
        dir << direction[0], direction[1], direction[2];
    
    for (auto const& theface : boundaryfaces(t.mesh()) )
    {                
        auto & face = boost::unwrap_ref( theface );

        auto &point = face.point(0);
        if (point.isOnBoundary())
        {    

            if constexpr(SolidMechanics::nDim == 2)
                origin << point.node()[0], point.node()[1];
            else if constexpr(SolidMechanics::nDim == 3)
                origin << point.node()[0], point.node()[1], point.node()[2];

            bvh_ray_type ray(origin,dir);
            auto bvh = boundingVolumeHierarchy(_range=boundaryfaces(wall));
            auto rayIntersection = bvh->intersect(ray) ;

            if (!rayIntersection.empty()) 
            {
                for ( auto const& rir : rayIntersection )
                {
                    for (auto const& ldof  : Vh->dof()->faceLocalDof( face.id() ))
                        g[ldof.index()] = rir.distance();
                }
            } 
        }    
    }

    // Export distance field
    auto exp = exporter(_mesh = t.mesh(), _name = fmt::format("Distance"));
    exp->addRegions();
    exp->add("g", g);
    exp->save();

    // Save field
    g.saveHDF5("distanceG.h5");
}

template<typename SolidMechanics>
typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype
getContactRegion(SolidMechanics t, std::vector<double> const& direction, double gamma0, int dispCondition, int raytracing, double distance, double lameFirstExpr, double lameSecondExpr, typename SolidMechanics::element_displacement_type u, int & nbrFaces)
{   
    // Get function spaces
    auto Vh =  Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    auto Xh = t.functionSpaceDisplacement();

    // Get data
    auto nwall = [&direction]() 
    { 
        if constexpr(SolidMechanics::nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(SolidMechanics::nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    auto [ havg, hmin, hmax ] = hMeasures( t.mesh() );

    auto const Id = eye<SolidMechanics::nDim,SolidMechanics::nDim>();
    auto epsv = sym(gradv(u));
    auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

    auto g = Vh->elementPtr();
    g->loadHDF5("/data/scratch/vanlandeghem/feel/ball_newmark/np_1/distanceG.h5");

    auto exp = exporter(_mesh = t.mesh(), _name = fmt::format("DistanceInter"));
    exp->addRegions();
    exp->add("g", idv(g));
    exp->save();

    // Set contact region
    auto contactregion = Vh->element();
    if (dispCondition == 1)
    {
        if (raytracing == 1)
            contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr = trans(nw)*idv(u) - idv(g) );
        else 
            contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr = trans(nw)*idv(u) - abs(cst(distance) - Py()));
    }
        
    else
    {
        if (raytracing == 1)
            contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr = (cst(gamma0)/h())*(trans(nw)*idv(u) - idv(g)) - trans(nw)*sigmav );
        else 
            contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr = (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(distance) - Py())) - trans(nw)*sigmav );
    }
        
    
    typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
        
    nbrFaces = 0;
    auto const& trialDofIdToContainerId =  form2( _test=Vh,_trial=Vh).dofIdToContainerIdTest();
    for (auto const& theface : boundaryfaces(t.mesh()) )
    {                
        auto & face = boost::unwrap_ref( theface );
        int contactDof = 0;
        for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
        {
            index_type thedof = ldof.index();
            thedof = trialDofIdToContainerId[ thedof ];

            if (contactregion[thedof] >= -hmin/10 )
                contactDof++;
                    
            if (SolidMechanics::nOrderDisplacement == 1)
            {
                if (SolidMechanics::nDim == 2)
                {
                    if (contactDof == 2)
                    {
                        nbrFaces++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (SolidMechanics::nDim == 3)
                {
                    std::cout << "TODO" << std::endl;
                }         
            }
            else if (SolidMechanics::nOrderDisplacement == 2)
            {
                if (SolidMechanics::nDim == 2)
                {
                    if (contactDof == 3)
                    {
                        nbrFaces++;
                        myelts->push_back( boost::cref( face ) );
                    }
                }
                else if (SolidMechanics::nDim == 3)
                {
                    std::cout << "TODO" << std::endl;
                }
            }
        }
    }
    myelts->shrink_to_fit();

    // Exports 
    if (SolidMechanics::nOrderDisplacement == 1)
    {
        if (dispCondition == 1)
        {
            if (raytracing == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", trans(nw)*idv(u) - idv(g), boundaryfaces(t.mesh()), "Pch1" ); 
            else 
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", trans(nw)*idv(u) - abs(cst(distance) - Py()), boundaryfaces(t.mesh()), "Pch1" ); 
        } 
        else 
        {
            if (raytracing == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", (cst(gamma0)/h())*(trans(nw)*idv(u) - idv(g)) - trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch1" ); 
            else 
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(distance) - Py())) - trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch1" );  
        } 
    }
    else if (SolidMechanics::nOrderDisplacement == 2)
    {
        if (dispCondition == 1)
        {
            if (raytracing == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", trans(nw)*idv(u) - idv(g), boundaryfaces(t.mesh()), "Pch2" ); 
            else 
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", trans(nw)*idv(u) - abs(cst(distance) - Py()), boundaryfaces(t.mesh()), "Pch2" ); 
        } 
        else 
        {
            if (raytracing == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", (cst(gamma0)/h())*(trans(nw)*idv(u) - idv(g)) - trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch2" ); 
            else 
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(distance) - Py())) - trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch2" );  
        } 
    }

    return myelts;
}

// Unilateral contact 

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralSteady( SolidMechanics t, DataType & data, std::vector<double> const& direction, const double gamma0, const double theta, int setA, std::vector<double> const& pressurePoint, int dispCondition, int withMarker, int raytracing, double distance) 
{
    // Set function spaces
    auto Xh = t.functionSpaceDisplacement();
    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    
    // Set data : deformation field, distance field, nw
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    auto const Id = eye<nDim,nDim>();

    auto g = Vh->elementPtr();
    g->loadHDF5("/data/scratch/vanlandeghem/feel/ball_newmark/np_1/distanceG.h5");

    // Set linear and bilinear forms
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
   
    if (setA == 1)
        A->zero();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            // Get Lame coefficients
            double lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            double lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            // Get deformation and stress tensors
            auto epst = sym(gradt(u));
            auto eps = sym(grad(u));
            auto epsv = sym(gradv(u));

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            if (setA == 1)
                bilinearFormDD = integrate(_range=elements(t.mesh()),_expr= inner((lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst),eps), _geomap=t.geomap());
            
            // Get contact region
            int nbrFaces = 0;
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            if (withMarker == 0)
            {
                myelts = getContactRegion(t, direction, gamma0, dispCondition,raytracing, distance, lameFirstExpr, lameSecondExpr,  u, nbrFaces);
                std::cout << "Number of faces in the contactRegion : " << nbrFaces << std::endl;
            }
            else 
                nbrFaces++;

            // Add contact terms        
            if (nbrFaces > 0)
            {
                if (withMarker == 0)
                {
                    auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
                    
                    // Exports 
                    if (SolidMechanics::nOrderDisplacement == 1)
                    {
                        t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
                        t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch1" );
                    }
                    else if (SolidMechanics::nOrderDisplacement == 2)
                    {
                        t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch2" );
                        t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch2" );
                    }
                    
                    // Add contact terms
                    if (raytracing == 1)
                    {
                        bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                        bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                        linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
                    }
                    else 
                    {
                        bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                        bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                        linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
                    }

                }
                else 
                {
                    if (raytracing == 1)
                    {
                        bilinearFormDD += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                        bilinearFormDD += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                        linearFormDisp += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
                    }
                    else 
                    {
                        bilinearFormDD += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                        bilinearFormDD += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                        linearFormDisp += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
                    }
                }
            }

            // Exports 
            if (SolidMechanics::nOrderDisplacement == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch1" );

            else if (SolidMechanics::nOrderDisplacement == 2)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch2" );
            

            // Contact Pressure at given point
            auto ctx = Vh->context();
            node_type t1(nDim);
            
            if constexpr (nDim == 2)
            {
                t1(0)=pressurePoint[0]; t1(1)=pressurePoint[1];
            }  
            else 
            {
                t1(0)=pressurePoint[0]; t1(1)=pressurePoint[1]; t1(2)=pressurePoint[2];
            }
            ctx.add( t1 );

            auto stress = project(_space=Vh, _range= boundaryfaces(t.mesh()), _expr =  trans(nw)*sigmav );
            auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(stress) );     
            
            std::cout << "Contact pressure at given point : " << evaluateStress(0,0) << std::endl;
        }
    }

    A->close();
    F->close();
}

// Bilateral contact : TODO

// Export solution and mesh

template<typename SolidMechanics>
void
storeData( SolidMechanics const& t ) 
{
    auto const& u = t.fieldDisplacement();
    u.saveHDF5("solution_ref.h5");

    auto const& mesh = t.mesh();
    saveGMSHMesh(_mesh=mesh,_filename= "mesh_ref.msh" );
}

// Compute error 

template<int nDim, typename SolidMechanics>
void
error( SolidMechanics const& t) 
{
    auto meshref = loadMesh(_mesh=new typename SolidMechanics::mesh_type, _filename = "/data/scratch/vanlandeghem/feel/convergence/np_1/mesh_ref.msh");
    auto uref = SolidMechanics::space_displacement_type::New(meshref)->elementPtr() ;
    uref->loadHDF5("/data/scratch/vanlandeghem/feel/benchmark/np_1/solution_ref.h5");

    auto op_inter = opInterpolation(_domainSpace =  SolidMechanics::space_displacement_type::New(t.mesh()), _imageSpace = SolidMechanics::space_displacement_type::New(meshref) );
    auto uinter =  SolidMechanics::space_displacement_type::New(meshref)->element(); 
    op_inter->apply(t.fieldDisplacement(), uinter);

    // Test
    auto exp = exporter(_mesh = meshref, _name = fmt::format("Test_error"));
    exp->addRegions();
    exp->add("uref",*uref);
    exp->add("uinter",uinter);
    exp->save();

    auto h1norm = normH1( _range=elements(meshref), _expr=idv(*uref), _grad_expr=gradv(*uref) );
    auto h1err = normH1( _range=elements(meshref), _expr = idv(uinter) - idv(*uref), _grad_expr = gradv(uinter) - gradv(*uref));
    std::cout << "H1 relative error : " << h1err/h1norm << std::endl;

    auto const Id = eye<nDim,nDim>();

    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
    double theta =  jsonCollisionForce["theta"].get<double>();
    
    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epsvref = sym(gradv(*uref));
    auto epsvinter = sym(gradv(uinter));

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmavref = (lameFirstExpr*trace(epsvref)*Id + 2*lameSecondExpr*epsvref)*N();
            auto sigmavinter = (lameFirstExpr*trace(epsvinter)*Id + 2*lameSecondExpr*epsvinter)*N();

            auto l2norm = normL2( _range=markedfaces(meshref,"Gamma_C"), _expr= cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) );
            auto l2err = normL2( _range=markedfaces(meshref,"Gamma_C"), _expr= sqrt(cst(gamma0)/h())*(cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) - cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(uinter) - trans(nw)*sigmavinter)) );
            std::cout << "Contact stress relative error : " << l2err/l2norm << std::endl;
        }
    }
}

