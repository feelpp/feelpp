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
using size_type_ = index_type;

class Measures
{
    public:
        static int nbr;
        static std::vector<int> it;
        static std::vector<double> Eh;
        static std::vector<double> Eh_1;
        static std::vector<double> Eh_2;
        static std::vector<double> R;
        static std::vector<double> R1;
        static std::vector<double> R2;
        static std::vector<double> disp;
        static std::vector<double> E_s;
        static std::vector<double> E;
        static std::vector<double> Lv;
        static std::vector<double> volume;
        static std::vector<double> evaluateStress;
};

int Measures::nbr = 0;
std::vector<int> Measures::it = {0};
std::vector<double> Measures::Eh = {0.};
std::vector<double> Measures::Eh_1 = {0.};
std::vector<double> Measures::Eh_2 = {0.};
std::vector<double> Measures::R = {0.};
std::vector<double> Measures::R1 = {0.};
std::vector<double> Measures::R2 = {0.};
std::vector<double> Measures::disp = {0.};
std::vector<double> Measures::E_s = {0.};
std::vector<double> Measures::E = {0.};
std::vector<double> Measures::Lv = {0.};
std::vector<double> Measures::volume = {0.};
std::vector<double> Measures::evaluateStress = {0.};

template<int nDim, typename SolidMechanics>
void
setWalls(SolidMechanics const& t)
{  
    std::string filename = "$cfgdir/wall.geo";
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<nDim>>("Wall"), _filename = filename);
    auto exp = exporter(_mesh = mesh, _name = fmt::format("Wall"));
    exp->addRegions();
    exp->save();
}

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralSteady( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, int & it, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());
    
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );

            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                    break;
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
                     
            if (it > 0)
            {
                bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner((cst(gamma0)/h())*trans(nw)*idt(u) - trans(nw)*sigmat, (cst(gamma0)/h())*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner((cst(gamma0)/h())*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }

            // Exports 
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, myfaces, "Pch1" );
        }
    }

    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamic( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, int & it, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());

    
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();


            // Compute faces in contact
            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactRegion", (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch1" ); 

            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                int contactDof = 0;
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                        contactDof++;
                    
                    
                    if (contactDof == 2)
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );

            // Add contact terms      
            if (it > 0)
            {
                bilinearFormDD += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }

            // Exports
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfaces, "Pch1" );
            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, myfaces, "Pch1" );
               
        }
    }
    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}
/*
template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamicIter( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
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
    auto Vh = Pch<1>(t.mesh());
    
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F); 

    auto bilinearFormCopie = form2( _test=Xh,_trial=Xh );
    auto linearFormCopie = form1( _test=Xh);

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto epst0 = sym(gradt(u));
            auto eps0 = sym(grad(u));
            auto epsv0 = sym(gradv(u));

            auto sigmat0 = (lameFirstExpr*trace(epst0)*Id + 2*lameSecondExpr*epst0)*N();
            auto sigma0 = (lameFirstExpr*trace(eps0)*Id + 2*lameSecondExpr*eps0)*N();
            auto sigmav0 = (lameFirstExpr*trace(epsv0)*Id + 2*lameSecondExpr*epsv0)*N();

            // Compute faces in contact
            auto contactregion0 = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav0 );
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts0( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            int it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                int contactDof = 0;
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion0[thedof] >= -1e-8 )
                        contactDof++;
                    
                    
                    if (contactDof == 2)
                    {
                        it++;
                        myelts0->push_back( boost::cref( face ) );
                    }   
                }
            }
            
            myelts0->shrink_to_fit();
            auto myfaces0 = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts0->begin(), myelts0->end(), myelts0 );

            
            // Iterations 
            if (it > 0)
            {
                // Iteration one
                std::cout << "Iteration 1" << std::endl;

                bilinearFormCopie = bilinearFormDD;
                linearFormCopie = linearFormDisp;

                bilinearFormCopie += integrate (_range=myfaces0,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat0, trans(nw)*sigma0) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=myfaces0,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat0, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma0) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=myfaces0,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma0) ,_geomap=t.geomap() );

                auto unew = Xh->element();  
                bilinearFormCopie.solve(_rhs=linearFormCopie,_solution=unew);


                auto epst1 = sym(gradt(unew));
                auto eps1 = sym(grad(unew));
                auto epsv1 = sym(gradv(unew));

                auto sigmat1 = (lameFirstExpr*trace(epst1)*Id + 2*lameSecondExpr*epst1)*N();
                auto sigma1 = (lameFirstExpr*trace(eps1)*Id + 2*lameSecondExpr*eps1)*N();
                auto sigmav1 = (lameFirstExpr*trace(epsv1)*Id + 2*lameSecondExpr*epsv1)*N();

                // Compute faces in contact
                auto contactregion1 = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew) - abs(cst(g) - Py())) - trans(nw)*sigmav1 );
                typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts1( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                it = 0;
                for (auto const& theface : boundaryfaces(t.mesh()) )
                {                
                    auto & face = boost::unwrap_ref( theface );
                    int contactDof = 0;
                    for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                    {
                        index_type thedof = ldof.index();
                        thedof = trialDofIdToContainerId[ thedof ];

                        if (contactregion1[thedof] >= -1e-8 )
                            contactDof++;
                    
                    
                        if (contactDof == 2)
                        {
                            it++;
                            myelts1->push_back( boost::cref( face ) );
                        }   
                    }
                }
            
                myelts1->shrink_to_fit();
                auto myfaces1 = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts1->begin(), myelts1->end(), myelts1 );
            
                // Iteration 2
                std::cout << "Iteration 2" << std::endl;

                bilinearFormCopie = bilinearFormDD;
                linearFormCopie = linearFormDisp;

                bilinearFormCopie += integrate (_range=myfaces1,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat1, trans(nw)*sigma1) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=myfaces1,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew) - trans(nw)*sigmat1, cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigma1) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=myfaces1,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigma1) ,_geomap=t.geomap() );

                auto unew1 = Xh->element();  
                bilinearFormCopie.solve(_rhs=linearFormCopie,_solution=unew1);

                auto epst2 = sym(gradt(unew1));
                auto eps2 = sym(grad(unew1));
                auto epsv2 = sym(gradv(unew1));

                auto sigmat2 = (lameFirstExpr*trace(epst2)*Id + 2*lameSecondExpr*epst2)*N();
                auto sigma2 = (lameFirstExpr*trace(eps2)*Id + 2*lameSecondExpr*eps2)*N();
                auto sigmav2 = (lameFirstExpr*trace(epsv2)*Id + 2*lameSecondExpr*epsv2)*N();

                // Compute faces in contact
                auto contactregion2 = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew1) - abs(cst(g) - Py())) - trans(nw)*sigmav2 );
                typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts2( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                it = 0;
                for (auto const& theface : boundaryfaces(t.mesh()) )
                {                
                    auto & face = boost::unwrap_ref( theface );
                    int contactDof = 0;
                    for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                    {
                        index_type thedof = ldof.index();
                        thedof = trialDofIdToContainerId[ thedof ];

                        if (contactregion2[thedof] >= -1e-8 )
                            contactDof++;
                    
                    
                        if (contactDof == 2)
                        {
                            it++;
                            myelts2->push_back( boost::cref( face ) );
                        }   
                    }
                }
            
                myelts2->shrink_to_fit();
                auto myfaces2 = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts2->begin(), myelts2->end(), myelts2 );


                // Last iteration
                std::cout << "Iteration 3" << std::endl;
                // Add contact terms   
            
                bilinearFormDD += integrate (_range=myfaces2,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat2, trans(nw)*sigma2) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfaces2,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmat2, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigma2) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfaces2,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigma2) ,_geomap=t.geomap() ); 
            
            }   
        }
    }
    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}
*/

template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamicIter( SolidMechanics t, DataType & data, const double g, std::vector<double> const& direction, const double gamma0, const double theta, const double tolerance, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
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
    auto Vh = Pch<1>(t.mesh());
    
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F); 

    auto bilinearFormCopie = form2( _test=Xh,_trial=Xh );
    auto linearFormCopie = form1( _test=Xh);

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto epst = sym(gradt(u));
            auto eps = sym(grad(u));
            auto epsv = sym(gradv(u));

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            // Compute faces in contact
            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            int it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                int contactDof = 0;
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                        contactDof++;
                    
                    
                    if (contactDof == 2)
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );




            // Iterations 
            if (it > 0)
            {
                // Definitions

                auto unew = Xh->element();
                auto unew1 = Xh->element();


                int iteration = 0;
                
                unew = u;
                unew1 = u;

                double error = integrate(_range=elements(t.mesh()), _expr = norm2( idv(unew)-idv(unew1))).evaluate()(0,0); 
                //double norm = integrate(_range=elements(t.mesh()), _expr =norm2(idv(unew))).evaluate()(0,0);
                double norm = 1.;


                while (((error/norm) > tolerance) || (iteration < 1))
                {
                    std::cout << " Relative error : " << error/norm << std::endl;
                    std::cout << "Iteration : " << iteration << std::endl;

                    // Copie current solution
                    unew = unew1; 

                    // Update 
                    auto epstnew = sym(gradt(unew1));
                    auto epsnew = sym(grad(unew1));
                    auto epsvnew = sym(gradv(unew1));

                    auto sigmatnew = (lameFirstExpr*trace(epstnew)*Id + 2*lameSecondExpr*epstnew)*N();
                    auto sigmanew = (lameFirstExpr*trace(epsnew)*Id + 2*lameSecondExpr*epsnew)*N();
                    auto sigmavnew = (lameFirstExpr*trace(epsvnew)*Id + 2*lameSecondExpr*epsvnew)*N();

                    // Compute faces in contact
                    auto contactregionnew = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew1) - abs(cst(g) - Py())) - trans(nw)*sigmavnew );
                    typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myeltsnew( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                    it = 0;
                    for (auto const& theface : boundaryfaces(t.mesh()) )
                    {                
                        auto & face = boost::unwrap_ref( theface );
                        int contactDof = 0;
                        for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                        {
                            index_type thedof = ldof.index();
                            thedof = trialDofIdToContainerId[ thedof ];

                            if (contactregionnew[thedof] >= -1e-8 )
                                contactDof++;
                    
                    
                            if (contactDof == 2)
                            {
                                it++;
                                myeltsnew->push_back( boost::cref( face ) );
                            }   
                        }
                    }
            
                    myeltsnew->shrink_to_fit();
                    auto myfacesnew = boost::make_tuple( mpl::size_t<MESH_FACES>(), myeltsnew->begin(), myeltsnew->end(), myeltsnew );
            
                    // Copie of bilinear and linear form

                    bilinearFormCopie = bilinearFormDD;
                    linearFormCopie = linearFormDisp;


                    // Add contact terms
                    
                    bilinearFormCopie += integrate (_range=myfacesnew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                    bilinearFormCopie += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                    linearFormCopie += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );


                    // Compute unew 
                    bilinearFormCopie.solve(_rhs=linearFormCopie,_solution=unew1);

                    // Update
                    iteration++;
                    error = integrate(_range=elements(t.mesh()), _expr = norm2( idv(unew)-idv(unew1))).evaluate()(0,0); 
                    //norm = integrate(_range=elements(t.mesh()), _expr =norm2(idv(unew))).evaluate()(0,0);
                    norm = 1.;

                    if (iteration == 10)
                        break;



                }
                
                // Last iteration
                std::cout << "Last Iteration" << std::endl;
                std::cout << " Relative error : " << error / norm << std::endl;

                // Update 
                auto epstnew = sym(gradt(unew1));
                auto epsnew = sym(grad(unew1));
                auto epsvnew = sym(gradv(unew1));

                auto sigmatnew = (lameFirstExpr*trace(epstnew)*Id + 2*lameSecondExpr*epstnew)*N();
                auto sigmanew = (lameFirstExpr*trace(epsnew)*Id + 2*lameSecondExpr*epsnew)*N();
                auto sigmavnew = (lameFirstExpr*trace(epsvnew)*Id + 2*lameSecondExpr*epsvnew)*N();

                // Compute faces in contact
                auto contactregionnew = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew1) - abs(cst(g) - Py())) - trans(nw)*sigmavnew );
                typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myeltsnew( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                it = 0;
                for (auto const& theface : boundaryfaces(t.mesh()) )
                {                
                    auto & face = boost::unwrap_ref( theface );
                    int contactDof = 0;
                    for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                    {
                        index_type thedof = ldof.index();
                        thedof = trialDofIdToContainerId[ thedof ];

                        if (contactregionnew[thedof] >= -1e-8 )
                            contactDof++;
                    
                    
                        if (contactDof == 2)
                        {
                            it++;
                            myeltsnew->push_back( boost::cref( face ) );
                        }   
                    }
                }
            
                myeltsnew->shrink_to_fit();
                auto myfacesnew = boost::make_tuple( mpl::size_t<MESH_FACES>(), myeltsnew->begin(), myeltsnew->end(), myeltsnew );
               
                bilinearFormDD += integrate (_range=myfacesnew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(g) - Py()), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
            
            } 
        }
    }
    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}


template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamicIterRaytracing( SolidMechanics t, DataType & data, std::vector<double> const& direction, const double gamma0, const double theta, const double tolerance, double & volume ) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());
    
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F); 

    auto bilinearFormCopie = form2( _test=Xh,_trial=Xh );
    auto linearFormCopie = form1( _test=Xh);


    // Compute distance 
    auto g = Vh->element();

    // Load wall mesh
    std::string filename = "$cfgdir/wall.geo";
    auto wall = loadMesh(_mesh=new Mesh<Simplex<nDim>>("Wall"), _filename = filename);

    // Raytracing
    using bvh_ray_type = BVHRay<SolidMechanics::nRealDim>;
    Eigen::VectorXd origin(nDim);
    Eigen::VectorXd dir(nDim);

    if constexpr(nDim == 2)
        dir << direction[0], direction[1];
    else if constexpr(nDim == 3)
        dir << direction[0], direction[1], direction[2];
    
    for (auto const& theface : boundaryfaces(t.mesh()) )
    {                
        auto & face = boost::unwrap_ref( theface );

        auto &point = face.point(0);
        if (point.isOnBoundary())
        {                
            origin << point.node()[0], point.node()[1];
            
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
    t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "distance", idv(g), boundaryfaces(t.mesh()), "Pch1" );

    // Get outward normal
    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            auto epst = sym(gradt(u));
            auto eps = sym(grad(u));
            auto epsv = sym(gradv(u));

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            // Compute faces in contact
            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(idv(g))) - trans(nw)*sigmav );
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            int it = 0;
            auto const& trialDofIdToContainerId = bilinearFormDD.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                int contactDof = 0;
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                        contactDof++;
                    
                    
                    if (contactDof == 2)
                    {
                        it++;
                        myelts->push_back( boost::cref( face ) );
                    }   
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );

            // Iterations 
            if (it > 0)
            {
                // Definitions

                auto unew = Xh->element();
                auto unew1 = Xh->element();

                int iteration = 0;
                
                unew = u;
                unew1 = u;

                double error = integrate(_range=elements(t.mesh()), _expr = norm2( idv(unew)-idv(unew1))).evaluate()(0,0); 
                double norm = integrate(_range=elements(t.mesh()), _expr =norm2(idv(unew))).evaluate()(0,0);
                //double norm = 1.;


                while (((error/norm) > tolerance) || (iteration < 1))
                {
                    std::cout << " Relative error : " << error/norm << std::endl;
                    std::cout << "Iteration : " << iteration << std::endl;

                    // Copie current solution
                    unew = unew1; 

                    // Update 
                    auto epstnew = sym(gradt(unew1));
                    auto epsnew = sym(grad(unew1));
                    auto epsvnew = sym(gradv(unew1));

                    auto sigmatnew = (lameFirstExpr*trace(epstnew)*Id + 2*lameSecondExpr*epstnew)*N();
                    auto sigmanew = (lameFirstExpr*trace(epsnew)*Id + 2*lameSecondExpr*epsnew)*N();
                    auto sigmavnew = (lameFirstExpr*trace(epsvnew)*Id + 2*lameSecondExpr*epsvnew)*N();

                    // Compute faces in contact
                    auto contactregionnew = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew1) - abs(idv(g))) - trans(nw)*sigmavnew );
                    typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myeltsnew( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                    it = 0;
                    for (auto const& theface : boundaryfaces(t.mesh()) )
                    {                
                        auto & face = boost::unwrap_ref( theface );
                        int contactDof = 0;
                        for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                        {
                            index_type thedof = ldof.index();
                            thedof = trialDofIdToContainerId[ thedof ];

                            if (contactregionnew[thedof] >= -1e-8 )
                                contactDof++;
                    
                    
                            if (contactDof == 2)
                            {
                                it++;
                                myeltsnew->push_back( boost::cref( face ) );
                            }   
                        }
                    }
            
                    myeltsnew->shrink_to_fit();
                    auto myfacesnew = boost::make_tuple( mpl::size_t<MESH_FACES>(), myeltsnew->begin(), myeltsnew->end(), myeltsnew );
            
                    // Copie of bilinear and linear form

                    bilinearFormCopie = bilinearFormDD;
                    linearFormCopie = linearFormDisp;


                    // Add contact terms
                    
                    bilinearFormCopie += integrate (_range=myfacesnew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                    bilinearFormCopie += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                    linearFormCopie += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(idv(g)), cst(gamma0)/h()*trans(nw)*id(unew) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );


                    // Compute unew 
                    bilinearFormCopie.solve(_rhs=linearFormCopie,_solution=unew1);

                    // Update
                    iteration++;
                    error = integrate(_range=elements(t.mesh()), _expr = norm2( idv(unew)-idv(unew1))).evaluate()(0,0); 
                    //norm = integrate(_range=elements(t.mesh()), _expr =norm2(idv(unew))).evaluate()(0,0);
                    //norm = 1.;

                    if (iteration == 10)
                        break;



                }
                
                // Last iteration
                std::cout << "Last Iteration" << std::endl;
                std::cout << " Relative error : " << error / norm << std::endl;

                // Update 
                auto epstnew = sym(gradt(unew1));
                auto epsnew = sym(grad(unew1));
                auto epsvnew = sym(gradv(unew1));

                auto sigmatnew = (lameFirstExpr*trace(epstnew)*Id + 2*lameSecondExpr*epstnew)*N();
                auto sigmanew = (lameFirstExpr*trace(epsnew)*Id + 2*lameSecondExpr*epsnew)*N();
                auto sigmavnew = (lameFirstExpr*trace(epsvnew)*Id + 2*lameSecondExpr*epsvnew)*N();

                // Compute faces in contact
                auto contactregionnew = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(unew1) - abs(idv(g))) - trans(nw)*sigmavnew );
                typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myeltsnew( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
                it = 0;
                for (auto const& theface : boundaryfaces(t.mesh()) )
                {                
                    auto & face = boost::unwrap_ref( theface );
                    int contactDof = 0;
                    for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                    {
                        index_type thedof = ldof.index();
                        thedof = trialDofIdToContainerId[ thedof ];

                        if (contactregionnew[thedof] >= -1e-8 )
                            contactDof++;
                    
                    
                        if (contactDof == 2)
                        {
                            it++;
                            myeltsnew->push_back( boost::cref( face ) );
                        }   
                    }
                }
            
                myeltsnew->shrink_to_fit();
                auto myfacesnew = boost::make_tuple( mpl::size_t<MESH_FACES>(), myeltsnew->begin(), myeltsnew->end(), myeltsnew );
               
                bilinearFormDD += integrate (_range=myfacesnew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                bilinearFormDD += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                linearFormDisp += integrate (_range=myfacesnew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(idv(g)), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
            
            } 
        }
    }
    A->close();
    F->close();

    // Volume
    volume = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
}



template<int nDim, typename SolidMechanics >
void
energy( SolidMechanics t) 
{
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    double g = jsonCollisionForce["walls"]["distances"].get<double>();
    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
    double theta =  jsonCollisionForce["theta"].get<double>();

    auto Xh = t.functionSpaceDisplacement();
    auto const& u = t.fieldDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    
    auto nw = nwall();
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));
    
    auto const Id = eye<nDim,nDim>();
    auto Vh = Pch<1>(t.mesh());

    auto form = form2( _test=Xh,_trial=Xh );

    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            // Get faces in contact
            auto contactregion = project(_space=Vh, _range=boundaryfaces(t.mesh()), _expr =  (cst(gamma0)/h())*(trans(nw)*idv(u) - abs(cst(g) - Py())) - trans(nw)*sigmav );
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            
            auto const& trialDofIdToContainerId = form.dofIdToContainerIdTest();
            for (auto const& theface : boundaryfaces(t.mesh()) )
            {                
                auto & face = boost::unwrap_ref( theface );
                int contactDof = 0;
                for( auto const& ldof : Vh->dof()->faceLocalDof( face.id() ) )
                {
                    index_type thedof = ldof.index();
                    thedof = trialDofIdToContainerId[ thedof ];

                    if (contactregion[thedof] >= -1e-8 )
                        contactDof++;
                    
                    if (contactDof == 2)
                        myelts->push_back( boost::cref( face ) );     
                }
            }
            
            myelts->shrink_to_fit();
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );
           
            // Compute energies
            auto const& v = t.fieldVelocity();
            
            // Strain energy 
            Measures::Eh_1[Measures::nbr] = 0.5*normL2Squared(_range = elements(t.mesh()), _expr = idv(v));
            Measures::Eh_2[Measures::nbr] = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)),gradv(u)) ).evaluate()( 0,0 );
            Measures::Eh[Measures::nbr] = Measures::Eh_1[Measures::nbr] + Measures::Eh_2[Measures::nbr];

            Measures::R1[Measures::nbr] = normL2Squared(_range=myfaces, _expr= sqrt(h()) * trans(nw)*(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)))*N());
            Measures::R2[Measures::nbr] = normL2Squared( _range=myfaces, _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*idv(u)  - abs(cst(g) - Py())) - trans(nw)*(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)))*N() ));
            Measures::disp[Measures::nbr] = integrate( _range= myfaces, _expr= trans(nw)*idv(u)  - abs(cst(g) - Py()) ).evaluate()( 0,0 );

            Measures::R[Measures::nbr] = (Measures::R1[Measures::nbr] - Measures::R2[Measures::nbr])/(2.*gamma0);
            Measures::E_s[Measures::nbr] = Measures::Eh[Measures::nbr] - theta*Measures::R[Measures::nbr];
            
            // Potential energy
            auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
            for ( auto const& bodyForce : physicSolidData->bodyForces() )
            {
                auto f = bodyForce.expr();
                Measures::Lv[Measures::nbr] = integrate( _range=elements(t.mesh()),_expr= inner( f,idv(u) ) ).evaluate()( 0,0 );
            }

            // Energy conservation
            Measures::E[Measures::nbr] = Measures::E_s[Measures::nbr] - Measures::Lv[Measures::nbr];

            auto ctx = Vh->context();
            node_type t1(nDim);
            // Coordiantes of lowest point
            if constexpr (nDim == 2)
            {
                t1(0)=0.0; t1(1)=-20.0;
            }  
            else 
            {
                t1(0)=0.0; t1(1)=-20.0; t1(2)=0.0;
            }
            ctx.add( t1 );
            auto stress = project(_space=Vh, _range= boundaryfaces(t.mesh()), _expr =  trans(nw)*(lameFirstExpr*trace(sym(gradv(u)))*Id + 2*lameSecondExpr*sym(gradv(u)))*N() );
            auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(stress) );     
            Measures::evaluateStress[Measures::nbr] = evaluateStress(0,0);
        }
    }
}

template<typename SolidMechanics>
void
storeData( SolidMechanics const& t ) 
{
    auto const& u = t.fieldDisplacement();
    u.saveHDF5("solution_ref.h5");

    auto const& mesh = t.mesh();
    saveGMSHMesh(_mesh=mesh,_filename= "mesh_ref.msh" );
}


template<int nDim, typename SolidMechanics>
void
error( SolidMechanics const& t) 
{
    auto meshref = loadMesh(_mesh=new typename SolidMechanics::mesh_type, _filename = "/data/scratch/vanlandeghem/feel/convergence/np_1/mesh_ref.msh");
    auto uref = SolidMechanics::space_displacement_type::New(meshref)->elementPtr() ;
    uref->loadHDF5("/data/scratch/vanlandeghem/feel/convergence/np_1/solution_ref.h5");

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

    // Get gamma0
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }
    double gamma0 = jsonCollisionForce["gamma_0"].get<double>();

    // Get nw
    std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    auto const Id = eye<nDim,nDim>();

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

            auto l2norm = normL2( _range=markedfaces(meshref,"contact_wall"), _expr= cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) );
            auto l2err = normL2( _range=markedfaces(meshref,"contact_wall"), _expr= sqrt(cst(gamma0)/h())*(cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(*uref) - trans(nw)*sigmavref) - cst(1.)/(cst(gamma0)/h()) * (cst(gamma0)/h()*trans(nw)*idv(uinter) - trans(nw)*sigmavinter)) );
            std::cout << "Contact stress relative error : " << l2err/l2norm << std::endl;
        }
    }
}

template<int nDim, typename SolidMechanics>
void
raytracing( SolidMechanics const& t, std::vector<double> const& direction) 
{
    auto Vh = Pch<1>(t.mesh());
    // g storing the distance for each boundary node of the body
    auto g = Vh->element();

    // Load wall mesh
    std::string filename = "$cfgdir/wall.geo";
    auto wall = loadMesh(_mesh=new Mesh<Simplex<nDim>>("Wall"), _filename = filename);

    // Raytracing
    using bvh_ray_type = BVHRay<SolidMechanics::nRealDim>;
    Eigen::VectorXd origin(nDim);
    Eigen::VectorXd dir(nDim);

    if constexpr(nDim == 2)
        dir << direction[0], direction[1];
    else if constexpr(nDim == 3)
        dir << direction[0], direction[1], direction[2];
    
    for (auto const& theface : boundaryfaces(t.mesh()) )
    {                
        auto & face = boost::unwrap_ref( theface );

        auto &point = face.point(0);
        if (point.isOnBoundary())
        {                
            origin << point.node()[0], point.node()[1];
            
            bvh_ray_type ray(origin,dir);
            auto bvh = boundingVolumeHierarchy(_range=boundaryfaces(wall));
            auto rayIntersection = bvh->intersect(ray) ;

            if (!rayIntersection.empty()) 
            {
                for ( auto const& rir : rayIntersection )
                {
                    for (auto const& ldof  : Vh->dof()->faceLocalDof( face.id() ))
                        //std::cout << g[ldof.index()] << std::endl;
                        g[ldof.index()] = rir.distance();
                }
            } 
        }    
        
    }

    auto exp = exporter(_mesh = t.mesh(), _name = fmt::format("Distance"));
    exp->addRegions();
    exp->add("g",g);
    exp->save();

               
}

template<typename SolidMechanics>
void
init_Measures(SolidMechanics const& t)
{
    Measures::nbr = 0;
    Measures::it.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh_1.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Eh_2.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::R.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::R1.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::R2.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::disp.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E_s.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Lv.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::volume.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::evaluateStress.resize(t.timeStepBase()->iterationNumber()-1);
}

template<typename SolidMechanics>
void
export_Measures(SolidMechanics const& t)
{
    std::ofstream myFile("measures.csv");
    
    myFile << "iter,Eh,Eh_1,Eh_2,R,R1,R2,disp,Estrain,E,Lv,it,volume,evaluateStress\n";
    
    for(int i = 0; i < Measures::nbr; ++i)
    {
        myFile << i << "," << Measures::Eh[i] << "," << Measures::Eh_1[i] << "," << Measures::Eh_2[i] << "," << Measures::R[i] << "," << Measures::R1[i] << "," << Measures::R2[i] << ","  << Measures::disp[i] << "," <<  Measures::E_s[i] << "," << Measures::E[i] << "," << Measures::Lv[i] << "," << Measures::it[i] << "," << Measures::volume[i] << "," << Measures::evaluateStress[i] << "\n";
    }
    
    myFile.close();
    
    // Reset
    Measures::nbr = 0;
    Measures::it.resize(1);
    Measures::Eh.resize(1);
    Measures::Eh_1.resize(1);
    Measures::Eh_2.resize(1);
    Measures::R.resize(1);
    Measures::R1.resize(1);
    Measures::R2.resize(1);
    Measures::disp.resize(1);
    Measures::E_s.resize(1);
    Measures::E.resize(1);
    Measures::Lv.resize(1);
    Measures::volume.resize(1);
    Measures::evaluateStress.resize(1);
}


/*
Read JSON and execute function
*/

template<int nDim, std::size_t residualType, typename SolidMechanics, typename DataType>
void
contactForceModels(SolidMechanics const& t, DataType & data)
{
    bool buildCstPart = data.buildCstPart();
    if(buildCstPart)
        return;
    
    fs::path path (Environment::expand(soption(_name="solid.filename")));
    json jsonCollisionForce;

    if (fs::exists(path))
    {
        std::ifstream i(path.string().c_str());
        json j = json::parse(i);
        jsonCollisionForce = j["CollisionForce"]["body"]["setup"];
    }

    std::string model = jsonCollisionForce["model"].get<std::string>(); // steady or dynamic
    std::string type = jsonCollisionForce["type"].get<std::string>(); // unilateral or bilateral
    
    // Model : steady
    if (model.compare("steady") == 0)
    {
        double volume = 0.;
        int it = 0;

        // Type case : unilateral
        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>(); // 1 or more

            if (nbr_walls == 1)
            {
                double g = jsonCollisionForce["walls"]["distances"].get<double>();
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                
                SingleUnilateralSteady<nDim, residualType>(t, data, g, direction, gamma0, theta, it, volume);

                // Store measures
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<double> g = jsonCollisionForce["walls"]["distances"].get<std::vector<double>>();
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }        
        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }    

    // Model : dynamic
    if (model.compare("dynamic") == 0)
    {
        int it = 0; 
        double volume = 0.;

        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

            if (nbr_walls == 1)
            {
                double g = jsonCollisionForce["walls"]["distances"].get<double>();
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                
                SingleUnilateralDynamic<nDim, residualType>(t, data, g, direction, gamma0, theta,it,volume);

                // Store measures
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<double> g = jsonCollisionForce["walls"]["distances"].get<std::vector<double>>();
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }

        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }

    if (model.compare("dynamicIter") == 0)
    {
        int it = 0; 
        double volume = 0.;

        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

            if (nbr_walls == 1)
            {
                double g = jsonCollisionForce["walls"]["distances"].get<double>();
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                double tolerance = jsonCollisionForce["tolerance"].get<double>();
                
                SingleUnilateralDynamicIter<nDim, residualType>(t, data, g, direction, gamma0, theta, tolerance,volume);

                // Store measures
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<double> g = jsonCollisionForce["walls"]["distances"].get<std::vector<double>>();
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }

        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }

    if (model.compare("dynamicIterRay") == 0)
    {
        int it = 0; 
        double volume = 0.;

        if (type.compare("unilateral") == 0)
        {
            int nbr_walls = jsonCollisionForce["walls"]["nbr"].get<int>();

            if (nbr_walls == 1)
            {
                std::vector<double> direction = jsonCollisionForce["walls"]["directions"].get<std::vector<double>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
                double tolerance = jsonCollisionForce["tolerance"].get<double>();
                
                SingleUnilateralDynamicIterRaytracing<nDim, residualType>(t, data, direction, gamma0, theta, tolerance,volume);

                // Store measures
                Measures::it[Measures::nbr] = it;
                Measures::volume[Measures::nbr] = volume;
                Measures::nbr += 1;

            }   
            else if (nbr_walls > 1)
            {
                std::vector<std::vector<double>> directions = jsonCollisionForce["walls"]["directions"].get<std::vector<std::vector<double>>>();
                double gamma0 = jsonCollisionForce["gamma_0"].get<double>();
                double theta =  jsonCollisionForce["theta"].get<double>();
            }

        }
        else if (type.compare("bilateral") == 0)
        {

        }
    }

}

