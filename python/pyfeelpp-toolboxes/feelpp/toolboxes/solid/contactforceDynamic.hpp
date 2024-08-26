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
        static std::vector<double> R2tmp;
        static std::vector<double> disp;
        static std::vector<double> E_s;
        static std::vector<double> E;
        static std::vector<double> Lv;
        static std::vector<double> EnoC;
        static std::vector<double> volume;
        static std::vector<double> evaluateStress;
        static std::vector<double> evaluateDisp;
        static std::vector<double> time;
};


int Measures::nbr = 0;
std::vector<int> Measures::it = {0};
std::vector<double> Measures::Eh = {0.};
std::vector<double> Measures::Eh_1 = {0.};
std::vector<double> Measures::Eh_2 = {0.};
std::vector<double> Measures::R = {0.};
std::vector<double> Measures::R1 = {0.};
std::vector<double> Measures::R2 = {0.};
std::vector<double> Measures::R2tmp = {0.};
std::vector<double> Measures::disp = {0.};
std::vector<double> Measures::E_s = {0.};
std::vector<double> Measures::E = {0.};
std::vector<double> Measures::Lv = {0.};
std::vector<double> Measures::EnoC = {0.};
std::vector<double> Measures::volume = {0.};
std::vector<double> Measures::evaluateStress = {0.};
std::vector<double> Measures::evaluateDisp = {0.};
std::vector<double> Measures::time = {0.};


template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamic( SolidMechanics t, DataType & data, std::vector<double> const& direction, const double gamma0, const double theta, int setA, int dispCondition, int withMarker, int raytracing, double distance ) 
{
    // Get function sapces
    auto Xh = t.functionSpaceDisplacement();
    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    
    // Get data
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
            {
                auto density = t.materialsProperties()->density( matName ).exprScalar().evaluate()(0,0);
                bilinearFormDD = integrate(_range=elements(t.mesh()),_expr= inner((lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst),eps), _geomap=t.geomap());
                bilinearFormDD += integrate( _range=elements(t.mesh()), _expr= t.timeStepNewmark()->polySecondDerivCoefficient()*density*inner(idt(u),id(u)),_geomap=t.geomap() );
            }

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
            
        }
    }
    A->close();
    F->close();
}


template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamicPenalty( SolidMechanics t, DataType & data, std::vector<double> const& direction, const double epsilon, double distance ) 
{
    // Get function sapces
    auto Xh = t.functionSpaceDisplacement();
    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    
    // Get data
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

    
    // Set linear and bilinear forms
    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();
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

            auto density = t.materialsProperties()->density( matName ).exprScalar().evaluate()(0,0);
            bilinearFormDD = integrate(_range=elements(t.mesh()),_expr= inner((lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst),eps), _geomap=t.geomap());
            bilinearFormDD += integrate( _range=elements(t.mesh()), _expr= t.timeStepNewmark()->polySecondDerivCoefficient()*density*inner(idt(u),id(u)),_geomap=t.geomap() );
            
            // Get contact region 
            int nbrFaces = 0;
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );

            myelts = getContactRegion(t, direction, 0, 1,0, distance, lameFirstExpr, lameSecondExpr,  u, nbrFaces);
            std::cout << "Number of faces in the contactRegion : " << nbrFaces << std::endl;

            // Add contact terms      
            if (nbrFaces > 0)
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

                bilinearFormDD += integrate (_range=myfaces,_expr= cst(1.)/cst(epsilon) * inner(trans(nw)*idt(u),trans(nw)*id(u)));
                linearFormDisp += integrate (_range=myfaces,_expr= cst(1.)/cst(epsilon) * inner(abs(cst(distance) - Py()),trans(nw)*id(u)));  
            }  

            // Exports 
            if (SolidMechanics::nOrderDisplacement == 1)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch1" );
            else if (SolidMechanics::nOrderDisplacement == 2)
                t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactPressure", trans(nw)*sigmav, boundaryfaces(t.mesh()), "Pch2" );
            
        }
    }
    A->close();
    F->close();

}


template<typename SolidMechanics, typename DataType>
typename SolidMechanics::element_displacement_type
oneIteration( SolidMechanics t, DataType & data, double lameFirstExpr, double lameSecondExpr, double density, int setA, std::vector<double> const& direction, const double gamma0, int dispCondition, int withMarker, double theta, int raytracing, double distance, typename SolidMechanics::element_displacement_type u, std::vector<int> & liste_faces, std::vector<double> & liste_error, std::vector<double> & liste_energy) 
{
    auto Xh = t.functionSpaceDisplacement();
    auto const& uTrue = t.fieldDisplacement();
    auto unew = t.functionSpaceDisplacement()->element();
    auto const Id = eye<SolidMechanics::nDim,SolidMechanics::nDim>();

    // Compute deformation and stress tensors
    auto epst = sym(gradt(u));
    auto eps = sym(grad(u));
    auto epsv = sym(gradv(u));

    auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
    auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
    auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

    auto g = Pch<SolidMechanics::nOrderDisplacement>(t.mesh())->elementPtr();
    g->loadHDF5("/data/scratch/vanlandeghem/feel/ball_newmark/np_1/distanceG.h5");

    auto nwall = [&direction]() 
    { 
        if constexpr(SolidMechanics::nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(SolidMechanics::nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    


    // Define linear and bilinear forms
    auto bilinearFormCopie = form2( _test=Xh,_trial=Xh );
    auto linearFormCopie = form1( _test=Xh);

    sparse_matrix_ptrtype& A = data.matrix();
    vector_ptrtype& F = data.rhs();

    if (setA == 1)
        A->zero();
    
    auto bilinearFormDD = form2( _test=Xh,_trial=Xh,_matrix=A );
    auto linearFormDisp = form1( _test=Xh, _vector=F);

    if (setA == 1)
    {
        bilinearFormDD = integrate(_range=elements(t.mesh()),_expr= inner((lameFirstExpr*trace(sym(gradt(uTrue)))*Id + 2*lameSecondExpr*sym(gradt(uTrue))),sym(grad(uTrue))), _geomap=t.geomap());
        bilinearFormDD += integrate( _range=elements(t.mesh()), _expr= t.timeStepNewmark()->polySecondDerivCoefficient()*density*inner(idt(uTrue),idv(uTrue)),_geomap=t.geomap() );
    }
    
    bilinearFormCopie = bilinearFormDD;
    linearFormCopie = linearFormDisp;

    // Get contact region 
    int nbrFaces = 0;
    typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
    if (withMarker == 0)
    {
        myelts = getContactRegion(t, direction, gamma0, dispCondition,raytracing, distance, lameFirstExpr, lameSecondExpr,  u, nbrFaces);
        std::cout << "Number of faces in the contactRegion : " << nbrFaces << std::endl;
        liste_faces.push_back(nbrFaces);
    }
    else 
        nbrFaces++;
    
    auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );

    // Add contact terms      
    if (nbrFaces > 0)
    {
        if (withMarker == 0)
        {           
            if (raytracing == 1)
            {
                bilinearFormCopie += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }
            else 
            {
                bilinearFormCopie += integrate (_range=myfaces,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=myfaces,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }
        }
        else 
        {
            if (raytracing == 1)
            {
                bilinearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }
            else 
            {
                bilinearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmat, trans(nw)*sigma) ,_geomap=t.geomap() ); 
                bilinearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(u) - trans(nw)*sigmat, cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() );
                linearFormCopie += integrate (_range=markedfaces(t.mesh(), "Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(u) - cst(theta)*trans(nw)*sigma) ,_geomap=t.geomap() ); 
            }
        }
    }  
                   
    // Compute unew 
    bilinearFormCopie.solve(_rhs=linearFormCopie,_solution=unew);

    // Compute error                
    double error = integrate(_range=elements(t.mesh()), _expr = norm2( idv(u)-idv(unew))).evaluate()(0,0); 
    liste_error.push_back(error);

    // Compute energy
    
    double E1 = 0.5*normL2Squared(_range = elements(t.mesh()), _expr = t.timeStepNewmark()->polyFirstDerivCoefficient()*idv(unew)-idv(t.timeStepNewmark()->polyFirstDeriv()));
    double E2 = 0.;
    if (setA == 1)
        E2 = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)),sym(gradv(unew))) ).evaluate()( 0,0 );
    else 
        E2 = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)),gradv(unew)) ).evaluate()( 0,0 );
            
    double R1 = 0; 
    if (withMarker == 0)
        R1 = normL2Squared(_range=myfaces, _expr= sqrt(h()) * trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N());
    else 
        R1 = normL2Squared(_range=markedfaces(t.mesh(), "Gamma_C"), _expr= sqrt(h()) * trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N());

    double R2 = 0;
    if (withMarker == 0)
    {
        if (raytracing == 1)
            R2 = normL2Squared( _range=myfaces, _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(unew))  - idv(g)) - trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N() ));
        else 
            R2 = normL2Squared( _range=myfaces, _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*idv(unew)  - abs(cst(distance) - Py())) - trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N() ));
    }
    else 
    {
        if (raytracing == 1)
            R2 = normL2Squared( _range=markedfaces(t.mesh(),"Gamma_C"), _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(unew))  - idv(g)) - trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N() ));
        else 
            R2 = normL2Squared( _range=markedfaces(t.mesh(),"Gamma_C"), _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*idv(unew)  - abs(cst(distance) - Py())) - trans(nw)*(lameFirstExpr*trace(sym(gradv(unew)))*Id + 2*lameSecondExpr*sym(gradv(unew)))*N() ));
    }
    
    double Lv = 0;
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<SolidMechanics::nDim>>(physicData);
        for ( auto const& bodyForce : physicSolidData->bodyForces() )
        {
            auto f = bodyForce.expr();
            Lv = integrate( _range=elements(t.mesh()),_expr= inner( f,idv(unew) ) ).evaluate()( 0,0 );
        }
    }

    double E = E1 + E2 - theta*(R1 - R2)/(2.*gamma0) - Lv;   
    liste_energy.push_back(E);

    return unew;
}


template<int nDim, std::size_t residualType,typename SolidMechanics,typename DataType>
void
SingleUnilateralDynamicFixedPoint( SolidMechanics t, DataType & data, std::vector<double> const& direction, const double gamma0, const double theta, const double tolerance, int setA, int dispCondition, int withMarker, int raytracing, double distance ) 
{
    // Get function spaces
    auto Xh = t.functionSpaceDisplacement();
    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    
    // Get data
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
            
            // Get Lame oefficients
            double lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            double lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);

            // Get deformation and stress tensors
            auto epst = sym(gradt(u));
            auto eps = sym(grad(u));
            auto epsv = sym(gradv(u));

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();

            double density = 0;
            if (setA == 1)
            {
                density = t.materialsProperties()->density( matName ).exprScalar().evaluate()(0,0);
                bilinearFormDD = integrate(_range=elements(t.mesh()),_expr= inner((lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst),eps), _geomap=t.geomap());
                bilinearFormDD += integrate( _range=elements(t.mesh()), _expr= t.timeStepNewmark()->polySecondDerivCoefficient()*density*inner(idt(u),idv(u)),_geomap=t.geomap() );
            }

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
            

            if (nbrFaces > 0)
            {
                // Definie data : deformation templates, listes
                auto unew = Xh->element();
                auto unew1 = Xh->element();

                unew = u;
                unew1 = u;

                int iteration = 0;
                double error = 0.; 
                std::vector<double> liste_error;
                std::vector<double> liste_energy;
                std::vector<int> liste_faces;

                // Start iterations

                while ((error > tolerance) || (iteration < 1))
                {
                    std::cout << " Iteration : " << iteration << " with error : " << error << std::endl;

                    // Copie current solution
                    unew = unew1; 

                    // Compute new solution
                    unew1 = oneIteration(t, data, lameFirstExpr, lameSecondExpr, density, setA, direction, gamma0, dispCondition, withMarker, theta, raytracing, distance, unew, liste_faces, liste_error, liste_energy); 

                    double Ecurrent = liste_energy[iteration];
                    std::cout << "Difference : " << abs(Ecurrent - Measures::E[Measures::nbr-1]) << std::endl;
                    std::cout << "Last : " << Measures::E[Measures::nbr] << std::endl;

                    // Update
                    error = liste_error[iteration];
                    iteration++;

                    if (abs(Ecurrent - Measures::E[Measures::nbr-1]) < 0.0001)
                        break;
                    
                    if (iteration == 20)
                        break;

                }

                unew1 = unew;
                
                // Export listes
                std::cout << "List error : " << liste_error << std::endl;
                std::cout << "List energy : " << liste_energy << std::endl;
                std::cout << "List faces : " << liste_faces << std::endl;

                // Update 
                auto epstnew = sym(gradt(unew1));
                auto epsnew = sym(grad(unew1));
                auto epsvnew = sym(gradv(unew1));

                auto sigmatnew = (lameFirstExpr*trace(epstnew)*Id + 2*lameSecondExpr*epstnew)*N();
                auto sigmanew = (lameFirstExpr*trace(epsnew)*Id + 2*lameSecondExpr*epsnew)*N();
                auto sigmavnew = (lameFirstExpr*trace(epsvnew)*Id + 2*lameSecondExpr*epsvnew)*N();

                // Get contact region 
                nbrFaces = 0;
                typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myeltsNew( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
                if (withMarker == 0)
                {
                    myeltsNew = getContactRegion(t, direction, gamma0, dispCondition,raytracing, distance, lameFirstExpr, lameSecondExpr,  unew1, nbrFaces);
                    std::cout << "Number of faces in the contactRegion : " << nbrFaces << std::endl;
                }
                else 
                    nbrFaces++;

                if (nbrFaces > 0)
                {
                    if (withMarker == 0)
                    {
                        auto myfacesNew = boost::make_tuple( mpl::size_t<MESH_FACES>(), myeltsNew->begin(), myeltsNew->end(), myeltsNew );
                    
                        // Exports 
                        if (SolidMechanics::nOrderDisplacement == 1)
                        {
                            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch1" );
                            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfacesNew, "Pch1" );
                        }
                        else if (SolidMechanics::nOrderDisplacement == 2)
                        {
                            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(0.), elements(t.mesh()), "Pch2" );
                            t.modelMesh().template updateField<typename SolidMechanics::mesh_type>( "contactFaces", cst(1.), myfacesNew, "Pch2" );
                        }

                        if (raytracing == 1)
                        {
                            bilinearFormDD += integrate (_range=myfacesNew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                            bilinearFormDD += integrate (_range=myfacesNew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                            linearFormDisp += integrate (_range=myfacesNew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                        }
                        else 
                        {
                            bilinearFormDD += integrate (_range=myfacesNew,_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                            bilinearFormDD += integrate (_range=myfacesNew,_expr=  cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew)  ,_geomap=t.geomap() );
                            linearFormDisp += integrate (_range=myfacesNew,_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew)  ,_geomap=t.geomap() ); 
                        }
                    }
                    else 
                    {
                        if (raytracing == 1)
                        {
                            bilinearFormDD += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr= - cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                            bilinearFormDD += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() );
                            linearFormDisp += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*id(g), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                        }
                        else 
                        {
                            bilinearFormDD += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr=- cst(theta)/(cst(gamma0)/h())*inner(trans(nw)*sigmatnew, trans(nw)*sigmanew),_geomap=t.geomap() ); 
                            bilinearFormDD += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr=  cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*trans(nw)*idt(unew1) - trans(nw)*sigmatnew, cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew)  ,_geomap=t.geomap() );
                            linearFormDisp += integrate (_range=markedfaces(t.mesh(),"Gamma_C"),_expr= cst(1.)/(cst(gamma0)/h())*inner(cst(gamma0)/h()*abs(cst(distance) - Py()), cst(gamma0)/h()*trans(nw)*id(unew1) - cst(theta)*trans(nw)*sigmanew) ,_geomap=t.geomap() ); 
                        }
                    }
                }   
            }
        }
    }
    A->close();
    F->close();
}

template<int nDim, typename SolidMechanics >
void
energy( SolidMechanics t) 
{
    // Get parameters from json
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
    std::vector<double> pressurePoint = jsonCollisionForce["pressurePoint"].get<std::vector<double>>();
    int setA = jsonCollisionForce["updateA"].get<int>();
    int dispCondition = jsonCollisionForce["dispCondition"].get<int>();
    int withMarker = jsonCollisionForce["withMarker"].get<int>();
    int raytracing = jsonCollisionForce["raytracing"].get<int>();
    double distance = 0;
    distance = jsonCollisionForce["walls"]["distances"].get<double>(); 
    
    // Get data 
    auto const& u = t.fieldDisplacement();
    auto const& v = t.fieldVelocity();

    auto Vh = Pch<SolidMechanics::nOrderDisplacement>(t.mesh());
    auto Xh = t.functionSpaceDisplacement();

    auto nwall = [&direction]() 
    { 
        if constexpr(nDim == 2) 
            return vec(cst(direction[0]),cst(direction[1]));
        else if constexpr(nDim == 3)
            return vec(cst(direction[0]),cst(direction[1]),cst(direction[2]));  
    };
    auto nw = nwall();
    auto const Id = eye<nDim,nDim>();
    
    auto g = Vh->elementPtr() ;
    g->loadHDF5("/data/scratch/vanlandeghem/feel/ball_newmark/np_1/distanceG.h5");
    
    for ( auto const& [physicName,physicData] : t.physicsFromCurrentType() )
    {
        for ( std::string const& matName : t.materialsProperties()->physicToMaterials( physicName ) )
        {
            auto const& matProperties = t.materialsProperties()->materialProperties( matName );
            
            // Get Lame coefficients
            auto lameFirstExpr = matProperties.property( "Lame-first-parameter" ).exprScalar().evaluate()(0,0);
            auto lameSecondExpr = matProperties.property( "Lame-second-parameter" ).exprScalar().evaluate()(0,0);
            
            // Get deformation and stress tensors
            auto epst = sym(gradt(u));
            auto eps = sym(grad(u));
            auto epsv = sym(gradv(u));

            auto sigmat = (lameFirstExpr*trace(epst)*Id + 2*lameSecondExpr*epst)*N();
            auto sigma = (lameFirstExpr*trace(eps)*Id + 2*lameSecondExpr*eps)*N();
            auto sigmav = (lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv)*N();
            auto sig = lameFirstExpr*trace(epsv)*Id + 2*lameSecondExpr*epsv;


            // Get contact region 
            int nbrFaces = 0;
            typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_ptrtype myelts( new typename MeshTraits<typename SolidMechanics::mesh_type>::faces_reference_wrapper_type );
            if (withMarker == 0)
                myelts = getContactRegion(t, direction, gamma0, dispCondition,raytracing, distance, lameFirstExpr, lameSecondExpr,  u, nbrFaces);
            else 
                nbrFaces++;
            
            auto myfaces = boost::make_tuple( mpl::size_t<MESH_FACES>(), myelts->begin(), myelts->end(), myelts );

            // Compute energies 
            //Measures::Eh_1[Measures::nbr] = 0.5*normL2Squared(_range = elements(t.mesh()), _expr = idv(v));
            Measures::Eh_1[Measures::nbr] = 0.5*normL2Squared(_range=elements(t.mesh()),_expr=t.timeStepNewmark()->polyFirstDerivCoefficient()*idv(u)-idv(t.timeStepNewmark()->polyFirstDeriv()));
            if (setA == 1)
                Measures::Eh_2[Measures::nbr] = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(sig,epsv) ).evaluate()( 0,0 );
            else 
                Measures::Eh_2[Measures::nbr] = 0.5*integrate( _range= elements( t.mesh() ), _expr= inner(sig,gradv(u)) ).evaluate()( 0,0 );
            
            Measures::Eh[Measures::nbr] = Measures::Eh_1[Measures::nbr] + Measures::Eh_2[Measures::nbr];

            if (withMarker == 0)
            {
                Measures::R1[Measures::nbr] = normL2Squared(_range=myfaces, _expr= sqrt(h()) * trans(nw)*sigmav);
                
                if (raytracing == 1)
                {
                    Measures::disp[Measures::nbr] = integrate( _range= myfaces, _expr= trans(nw)*(idv(u))  - idv(g) ).evaluate()( 0,0 );
                    Measures::R2[Measures::nbr] = normL2Squared( _range=myfaces, _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(u))  - idv(g)) - trans(nw)*sigmav ));
                }
                    
                else 
                {
                    Measures::disp[Measures::nbr] = integrate( _range= myfaces, _expr= trans(nw)*(idv(u))  - abs(cst(distance) - Py()) ).evaluate()( 0,0 );
                    Measures::R2[Measures::nbr] = normL2Squared( _range=myfaces, _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(u))  - abs(cst(distance) - Py())) - trans(nw)*sigmav ));
                }
            }
            else 
            {
                Measures::R1[Measures::nbr] = normL2Squared(_range=markedfaces(t.mesh(),"Gamma_C"), _expr= sqrt(h()) * trans(nw)*sigmav);
                if (raytracing == 1)
                {
                    Measures::disp[Measures::nbr] = integrate( _range= markedfaces(t.mesh(),"Gamma_C"), _expr= trans(nw)*(idv(u))  - idv(g) ).evaluate()( 0,0 );
                    Measures::R2[Measures::nbr] = normL2Squared( _range=markedfaces(t.mesh(),"Gamma_C"), _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(u))  - idv(g)) - trans(nw)*sigmav ));
                }
                    
                else 
                {
                    Measures::disp[Measures::nbr] = integrate( _range= markedfaces(t.mesh(),"Gamma_C"), _expr= trans(nw)*(idv(u))  - abs(cst(distance) - Py()) ).evaluate()( 0,0 );
                    Measures::R2[Measures::nbr] = normL2Squared( _range=markedfaces(t.mesh(),"Gamma_C"), _expr= sqrt(h()) * (cst(gamma0)/h() * (trans(nw)*(idv(u))  - abs(cst(distance) - Py())) - trans(nw)*sigmav ));
                }
            }

            Measures::R[Measures::nbr] = (Measures::R1[Measures::nbr] - Measures::R2[Measures::nbr])/(2.*gamma0);
            Measures::E_s[Measures::nbr] = Measures::Eh[Measures::nbr] - theta*Measures::R[Measures::nbr];

            auto physicSolidData = std::static_pointer_cast<ModelPhysicSolid<nDim>>(physicData);
            for ( auto const& bodyForce : physicSolidData->bodyForces() )
            {
                auto f = bodyForce.expr();
                Measures::Lv[Measures::nbr] = integrate( _range=elements(t.mesh()),_expr= inner( f,idv(u) ) ).evaluate()( 0,0 );
            }

            // Energy conservation
            Measures::E[Measures::nbr] = Measures::E_s[Measures::nbr] - Measures::Lv[Measures::nbr];
            Measures::EnoC[Measures::nbr] = Measures::Eh[Measures::nbr] - Measures::Lv[Measures::nbr];

            auto ctx = Vh->context();
            node_type t1(nDim);
            
            // Coordiantes of given point
            if constexpr (nDim == 2)
            {
                t1(0)=pressurePoint[0]; t1(1)=pressurePoint[1];
            }  
            else 
            {
                t1(0)=pressurePoint[0]; t1(1)=pressurePoint[1]; t1(2)=pressurePoint[2];
            }

            ctx.add( t1 );
   
            auto stress = project(_space=Vh, _range= boundaryfaces(t.mesh()), _expr =  trans(nw)*sig*nw );
            auto Disp = project(_space=Vh, _range= boundaryfaces(t.mesh()), _expr = trans(nw)*(idv(u)) );

            auto evaluateStress = evaluateFromContext( _context=ctx, _expr= idv(stress) ); 
            auto evaluateDisp = evaluateFromContext( _context=ctx, _expr= -idv(Disp) );     
            
            Measures::evaluateStress[Measures::nbr] = evaluateStress(0,0);
            Measures::evaluateDisp[Measures::nbr] = evaluateDisp(0,0);
        }
    }

    Measures::volume[Measures::nbr] = integrate (_range = elements(t.mesh()), _expr = cst(1.0)).evaluate()( 0, 0 );
    Measures::time[Measures::nbr] = t.timeStepBase()->time();            
    Measures::nbr += 1;
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
    Measures::R2tmp.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::disp.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E_s.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::E.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::Lv.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::EnoC.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::volume.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::evaluateStress.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::evaluateDisp.resize(t.timeStepBase()->iterationNumber()-1);
    Measures::time.resize(t.timeStepBase()->iterationNumber()-1);
}

template<typename SolidMechanics>
void
export_Measures(SolidMechanics const& t)
{
    std::ofstream myFile("measures.csv");
    
    myFile << "iter,Eh,Eh_1,Eh_2,R,R1,R2,R2tmp,disp,Estrain,E,Lv,EnoC,it,volume,evaluateStress,evaluateDisp,time\n";
    
    for(int i = 0; i < Measures::nbr; ++i)
    {
        myFile << i << "," << Measures::Eh[i] << "," << Measures::Eh_1[i] << "," << Measures::Eh_2[i] << "," << Measures::R[i] << "," << Measures::R1[i] << "," << Measures::R2[i] << "," << Measures::R2tmp[i] << "," << Measures::disp[i] << "," <<  Measures::E_s[i] << "," << Measures::E[i] << "," << Measures::Lv[i] << "," << Measures::EnoC[i] << "," << Measures::it[i] << "," << Measures::volume[i] << "," << Measures::evaluateStress[i] << "," << Measures::evaluateDisp[i] << "," << Measures::time[i] <<  "\n";
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
    Measures::R2tmp.resize(1);
    Measures::disp.resize(1);
    Measures::E_s.resize(1);
    Measures::E.resize(1);
    Measures::Lv.resize(1);
    Measures::EnoC.resize(1);
    Measures::volume.resize(1);
    Measures::evaluateStress.resize(1);
    Measures::evaluateDisp.resize(1);
    Measures::time.resize(1);
}