/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2012-02-01

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file nirb.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2012-02-01
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

//#include <feel/feelalg/solvereigen.hpp>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <vector>
#include <algorithm>

#include<boost/range/algorithm/max_element.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline po::options_description makeOptions()
{
    po::options_description NIRBoptions("Nirb options");
    NIRBoptions.add_options()
        // meshes parameters
        ("hfinsize", po::value<double>()->default_value( 0.05 ), "fine mesh size")
        ("hcoarsesize", po::value<double>()->default_value( 0.1 ), "coarse mesh size")
        //("ReadingCoarseMesh",po::value<int>()->default_value(0),"Reading COARSE MESH in file if set to 1");
        //("ReadingFineMesh",po::value<int>()->default_value(0),"Reading FINE MESH in file if set to 1");
    
        //("FineMeshFilename",po::value<string>()->default_value("FineMesh.msh"),"Fine Mesh filename");
        //("CoarseMeshFilename",po::value<string>()->default_value("CoarseMesh.msh"),"Fine Mesh filename");
        // Reduced basis parameters
        ("NbSnapshot",po::value<int>()->default_value(100),"numbers of snapshot computed")
        ("sizeRB",po::value<int>()->default_value(15),"size of reduced basis")

        ("muMin", po::value<double>()->default_value( 0 ), "angle in [0,pi/2]")
        ("muMax", po::value<double>()->default_value( M_PI/2.),"angle in [0,pi/2]")
    
        ("mu",po::value<double>()->default_value(1.),"angle in [0,pi/2]")
    
        ("Sampling",po::value<int>()->default_value(1),"Does not compute sampling if set to 0, (set to 1 by default)")

        ("Offline",po::value<int>()->default_value(1),"integer equal to 0 if the offline has not  to be done")
    
        ("ComputeError",po::value<int>()->default_value(1),"integer equal to 0 if the error computation has not to be done")

        ("polynomialOrder",po::value<int>()->default_value(3),"polynomial order");
        ;
    return NIRBoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline AboutData makeAbout()
{
    AboutData about( "nirb-test" ,
                     "nirb-test" ,
                     "0.2",
                     "Non intrusive reduced basis",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

//inline
//bool SortFunction( double i, double j) { return std::abs(i)<std::abs(j) ;}

template<int PolynomialOrder> class NIRBTEST
    :
    public Simget
{
    typedef Simget super;
public:

    //static const uint16_type Order = 3;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension 2 of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //! the basis type of our approximation space
    typedef bases<Lagrange<PolynomialOrder,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef std::vector<element_type> un_type;

    /**
     * Constructor
     */
    NIRBTEST( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        FineMeshSize( this->vm()["hfinsize"].template as<double>() ),
        CoarseMeshSize( this->vm()["hcoarsesize"].template as<double>() ),
        //ReadingCoarseMesh( this->vm()["ReadingCoarseMesh"].template as<int>() ),
        //ReadingFineMesh( this->vm()["ReadingFineMesh"].template as<int>() ),
        NbSnapshot( this->vm()["NbSnapshot"].template as<int>() ),
        sizeRB(this->vm()["sizeRB"].template as<int>()),
        muMin( this->vm()["muMin"].template as<double>() ),
        muMax( this->vm()["muMax"].template as<double>() ),
        mu(this->vm()["mu"].template as<double>() ),
        Sampling(this->vm()["Sampling"].template as<int>()),
        Offline( this->vm()["Offline"].template as<int>() ),
        ComputeError( this->vm()["ComputeError"].template as<int>() )
    {
    }

    void run();
    
    void run( const double* X, unsigned long P, double* Y, unsigned long N ); 
    element_type blackbox( space_ptrtype Xh, double param);
    element_type BuildNirbSolution(space_ptrtype XhFine,space_ptrtype XhCoarse,double param);
    void ConstructNIRB (space_ptrtype Xh);
    void ComputeSnapshot(space_ptrtype Xh);
    void ChooseRBFunction(space_ptrtype Xh);
    void OrthogonalisationRBFunction(space_ptrtype Xh);
    void OrthogonalisationRBFunctionL2GrammSchmidt(space_ptrtype Xh,Eigen::MatrixXd &M );
    
    double OrthogonalisationRBFunctionL2GrammSchmidt(space_ptrtype Xh, Eigen::MatrixXd &MU,int n);
 

    Eigen::MatrixXd  ConstructStiffMatrixSnapshot(space_ptrtype Xh);
    Eigen::MatrixXd  ConstructStiffMatrixRB(space_ptrtype Xh,std::string filename);
    Eigen::MatrixXd  ConstructMassMatrixRB(space_ptrtype Xh,std::string filename);
    
    
    gmsh_ptrtype createGeo(double hsize,std::string MeshFileName);

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double CoarseMeshSize;
    double FineMeshSize;
    //int ReadingCoarseMesh; // By Default = 0
    //int ReadingFineMesh;   // By Default = 0

    //! Reduced basis parameter
    int NbSnapshot,sizeRB;
    double muMin,muMax,mu;

    // Paramater to set OFF or ON the "offline" procedure
    
    int Sampling;
    // if sampling = 1 ==> compute sampling (by default)
    // if sampling = 0 ==> read sampling from a file

    int Offline;
    // if Offline = 1 then ConstructNIRB is called
    // if Offline = 0 then ConstructNIRB is not called, the Nirb Basis function are read from a file

    // Paramater to set OFF or ON the computation of the error between the F.E and the NIRB Methods
    int ComputeError;
    //if ComputeError = 1  == ON
    //if ComputeError = 0 == OFF
    
    int polynomialOrder;
    }; // NIRBTEST
//const uint16_type NIRBTEST::Order;
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder>::run()
{

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute NIRBTEST<" << PolynomialOrder << ">\n";
    std::vector<double> X(10);
    X[0] = FineMeshSize;
    X[1] = CoarseMeshSize;
    X[2] = NbSnapshot;
    X[3] = sizeRB;
    X[4] = muMin;
    X[5] = muMax;
    X[6] = mu;
    X[7] = Sampling;
    X[8] = Offline;
    X[9] = ComputeError;
    std::vector<double> Y(1);
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    /*
    std::cout << "fine mesh size : " << FineMeshSize << std::endl;
    std::cout << "coarse mesh size : " << CoarseMeshSize << std::endl;
    */
    std:: cout << "Size reduced basis : " << sizeRB << endl;
    std:: cout << "Polynomial order : P" << PolynomialOrder << endl;

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "%1%/P%2%/" )
                                       % this->about().appName()
                                       % PolynomialOrder );
 
     
    mesh_ptrtype meshExtraCoarse, meshCoarse, meshFine ;
    
    std::ifstream MeshFiles ("/home/chakir/Mesh/Mesh_H0.msh");
    if (MeshFiles){
        std::cout << "Meshes read in file " << endl;
        meshExtraCoarse = loadGMSHMesh(_mesh=new mesh_type,
                                       _filename="/home/chakir/Mesh/Mesh_H0.msh"); 
        
        meshCoarse  =  loadGMSHMesh( _mesh=new mesh_type,
                                    _filename="/home/chakir/Mesh/Mesh_H1.msh");  

        
        meshFine  =  loadGMSHMesh( _mesh=new mesh_type,
                                  _filename="/home/chakir/Mesh/Mesh_H2.msh");  
    }
    else {
        std::cout << "Meshes build using 'createGMSHMesh' " << endl;
        meshExtraCoarse = createGMSHMesh( _mesh = new mesh_type,
                                         _desc = createGeo(X[1],"Mesh_ExtraCoarse"),
                                         _update=MESH_UPDATE_FACES | MESH_UPDATE_EDGES );
        meshCoarse = createGMSHMesh(_mesh= new mesh_type,
                                    _desc = createGeo(X[1]/2,"Mesh_Coarse"),
                                    _refine=1);
        if(X[1]/4 > X[0]){
            meshFine = createGMSHMesh(_mesh= new mesh_type,
                                      _desc = createGeo(X[0],"Mesh_Fine"),
                                      _refine=1);
        }
        else{
            meshFine = createGMSHMesh(_mesh= new mesh_type,
                                      _desc = createGeo(X[1]/4,"Mesh_Fine"),
                                      _refine=1);
        }

    }
    
    
    
    space_ptrtype XhFine = space_type::New( meshFine );
    //Fine F.E space used to build the NIRB basis
    
    space_ptrtype XhCoarse = space_type::New( meshCoarse );
    space_ptrtype XhExtraCoarse = space_type::New( meshExtraCoarse );



    //STEP ONE : Construction of the "non intruisive reduced basis (nirb) functions"
    //Computation of the X[3] snapshots solution on Xhfine
    //Selection of X[2] fonctions to build the "reduced basis" using a POD technique
    //Orthogonalisation in L2 and H1 norm of "reduced basis function"
    //Save this final functions
    
    
    //STEP ONE : Construction of the "non intruisive reduced basis (nirb) functions"
    //Computation of the X[3] snapshots solution on Xhfine
    //Selection of X[2] fonctions to build the "reduced basis" using a POD technique
    //Orthogonalisation in L2 and H1 norm of "reduced basis function"
    //Save this final functions
    
     
    if (!Offline){
        std :: string path;
        //Checking if the NIRB basis files exist
        //WARNING : We also must check if the basis have been computed on the right FE space (Not done yet)
        for (int i = 0 ; i < sizeRB; i++){
            path = (boost::format("./RB%1%NIRB_BasisFile_%2%")%sizeRB%i).str();
            auto ui = XhFine->element();
            ui.load(_path=path);
            int NdofNirb;
            NdofNirb = ui.size();
            if (NdofNirb != XhFine->nLocalDof()){
                std::cout <<"WARNING Error in  NIRB_BasisFile_" << i << " Ndof not the same as in Xhfine => OFFline parameter set to 1 " << endl;
                Offline = 1;
                break;
            }
        }
    }
     
    if (Offline){
        ConstructNIRB(XhFine) ;
    }
    
    //STEP TWO : Approximation of the solution using the "nirb" functions for a choosen mu
    
     
    
    
    double p = mu;
    
    

    
    
    auto uNirbCoarse = XhFine->element();   // Fine/Coarse Grid NIRB solution
     
    double TimeCoarse,TimeFine;
    boost::timer ti;   
    ti.restart();
    uNirbCoarse = BuildNirbSolution(XhFine,XhCoarse,p);
    TimeCoarse =  ti.elapsed();
    std :: cout << "Construction of NIRB solution (uNirbCoarse) - Fine/Coarse Grid (saved in nirb2GridCoarse):" << endl;
    std::cout << "Time to build solution " << TimeCoarse << " sec" << endl;    
    
    export_ptrtype exporter2GridCoarse(export_type::New( this->vm(), "nirb2GridCoarse"));
    exporter2GridCoarse->step(0)->setMesh( meshFine );
    exporter2GridCoarse->step(0)->add("uNirbCoarse",uNirbCoarse);
    exporter2GridCoarse->save();
    
    if (ComputeError){
        std :: cout << "Error calculation " << endl;
        
        
        auto u1Grid = XhFine->element();
        ti.restart();
        u1Grid = BuildNirbSolution(XhFine,XhFine,p);
        TimeFine =  ti.elapsed();
        std::cout << "Construction of uNirbFine - Fine/Fine Grid  (saved in nirb1Grid): done "<< endl;
        std::cout << "Time to build  " << TimeFine << " sec" << endl;
        
        
        
        mesh_ptrtype meshRef;
        
        if (MeshFiles){
            meshRef  =  loadGMSHMesh( _mesh=new mesh_type,
                                    _filename="/home/chakir/Mesh/Mesh_H3.msh"); 
        }
        else{
            if(X[1]/4 > X[0]){
                meshRef = createGMSHMesh(_mesh= new mesh_type,
                                          _desc = createGeo(X[0]/2,"Mesh_Fine"),
                                          _refine=1);
            }
            else{
                meshRef = createGMSHMesh(_mesh= new mesh_type,
                                          _desc = createGeo(X[1]/8,"Mesh_Fine"),
                                          _refine=1);
            }
        }
       
        
       
        space_ptrtype XhRef = space_type::New(meshRef);
        
        
        
        //Reference F.E space used to compute the "reference" solution
        //which supposed to be accurate enough
        
        
        auto uRef = XhRef->element();
        
        uRef = blackbox( XhRef, p );
        std :: cout << "Computation of FE solution (uRef) - Ref Grid (not saved): done" << endl;
        
        
        //Computation of the H1 norm of uRef
        double L2NormUref = std :: sqrt(integrate(_range=elements(XhRef->mesh()),
                                                  _expr=( idv(uRef)*idv(uRef) )).evaluate()(0,0));
        double SemiH1NormUref = std :: sqrt(integrate(_range=elements(XhRef->mesh()),
                                                      _expr=(gradv(uRef)*trans(gradv(uRef)))).evaluate()(0,0));
        double H1NormUref  = std :: sqrt(SemiH1NormUref *SemiH1NormUref + L2NormUref*L2NormUref);
        
        //std::cout << "L2NormURef = " << L2NormUref <<  " -- H1NormUref = " << H1NormUref << endl;
        //std::cout << endl;
        
        
        
        // //Computation of relative error measured in  H1 norm
        
        
        // // auto DFine = M_backend->newMatrix( _test=XhRef, _trial=XhFine  );//Sparse F.E mass matrix D
        // // form2( _test=XhRef, _trial=XhFine, _matrix=DFine ) =
        // //     integrate( _range=elements(XhFine->mesh()), _expr=idt(uRef)*id(uFine));
        // // auto DCoarse = M_backend->newMatrix( _test=XhRef, _trial=XhCoarse  );//Sparse F.E mass matrix D
        // // form2( _test=XhRef, _trial=XhCoarse, _matrix=DCoarse ) =
        // //     integrate( _range=elements(XhFine->mesh()), _expr=idt(uRef)*id(uCoarse));
        
        
        
  
        double  ErrH1uNirbCoarse,ErrH1u1Grid;
       
        
        
        // uRef - uNirbCoarse
        ErrH1uNirbCoarse = integrate(_range=elements(XhRef->mesh()),
                                     _expr=( (gradv(uRef)-gradv(uNirbCoarse))*trans(gradv(uRef)-gradv(uNirbCoarse))
                                            + (idv(uRef)-idv(uNirbCoarse))*(idv(uRef)-idv(uNirbCoarse)) )).evaluate()(0,0);
        ErrH1uNirbCoarse = sqrt(ErrH1uNirbCoarse)/H1NormUref;
        
        
        //uRef - uNirbFine (= uRef - u1Grid)
        ErrH1u1Grid = integrate(_range=elements(XhRef->mesh()),
                                _expr=( (gradv(uRef)-gradv(u1Grid))*trans(gradv(uRef)-gradv(u1Grid))
                                       + (idv(uRef)-idv(u1Grid))*(idv(uRef)-idv(u1Grid)) )).evaluate()(0,0);
        ErrH1u1Grid = sqrt(ErrH1u1Grid)/H1NormUref;
          
        std :: cout << "H1-norm NIRB Error " << endl;
        std :: cout << "||u_ref - u_NirbCoarse ||_{H1}  = " << ErrH1uNirbCoarse << endl;
        std :: cout << "||u_ref - u_NirbFine ||_{H1}  = " << ErrH1u1Grid<< endl;
        std :: cout <<  endl;
        
        
        
        
        
        export_ptrtype exporter1Grid(export_type::New( this->vm(), "nirbFine"));
        exporter1Grid->step(0)->setMesh( meshFine );
        exporter1Grid->step(0)->add("uNirbFine",u1Grid);    
        exporter1Grid->save();
        
        
         
        
        auto uFine = XhFine->element();
        auto uCoarse = XhCoarse->element();
        uCoarse = blackbox( XhCoarse, p);
        std :: cout << "Computation of FE solution (uCoarse) -  Coarse Grid (not saved):"<< endl;
        
        uFine = blackbox( XhFine, p );
        std :: cout << "Computation of FE solution (uFine) - Fine Grid (saved in EFFine): done" << endl;
        
        
         
        // uref - uFine 
        double ErrSemiH1uFine,ErrL2uFine,ErrH1uFine;
        
        
        ErrL2uFine = std::sqrt(integrate(_range=elements(XhRef->mesh()),
                                         _expr=(  (idv(uRef)-idv(uFine))*(idv(uRef)-idv(uFine)) )).evaluate()(0,0));
        
        
        ErrSemiH1uFine = std::sqrt(integrate(_range=elements(XhRef->mesh()),
                                             _expr=( (gradv(uRef)-gradv(uFine))*trans(gradv(uRef)-gradv(uFine)))).evaluate()(0,0));
        
        
        ErrH1uFine = std::sqrt(ErrSemiH1uFine *ErrSemiH1uFine + ErrL2uFine*ErrL2uFine);
        ErrH1uFine = ErrH1uFine/H1NormUref;
        ErrL2uFine = ErrL2uFine/L2NormUref;
        
        
        
        // uref - uCoarse      
        double ErrL2uCoarse ,ErrSemiH1uCoarse ,ErrH1uCoarse;
        ErrL2uCoarse = std::sqrt(integrate(_range=elements(XhRef->mesh()),
                                           _expr=(  (idv(uRef)-idv(uCoarse))*(idv(uRef)-idv(uCoarse)) )).evaluate()(0,0));
        
        
        ErrSemiH1uCoarse = std::sqrt(integrate(_range=elements(XhRef->mesh()),
                                               _expr=( (gradv(uRef)-gradv(uCoarse))*trans(gradv(uRef)-gradv(uCoarse)))).evaluate()(0,0));
        
        ErrH1uCoarse = std::sqrt(ErrSemiH1uCoarse * ErrSemiH1uCoarse + ErrL2uCoarse*ErrL2uCoarse);
        ErrH1uCoarse = ErrH1uCoarse/H1NormUref;
        ErrL2uCoarse = ErrL2uCoarse/L2NormUref;
        
        
        std :: cout << "H1-Norm F.E Error " << endl;
        std :: cout << "||u_ref - u_coarse ||_{H1}  = " << ErrH1uCoarse << endl;
        std :: cout << "||u_ref - u_fine ||_{H1}  = " << ErrH1uFine << endl;
        std :: cout << endl;
        
        
        export_ptrtype exporterFine(export_type::New( this->vm(), "EF_Fine" ) );
        export_ptrtype exporterCoarse(export_type::New( this->vm(), "EF_Coarse" ) );
        exporterFine->step(0)->setMesh( meshFine );
        exporterCoarse->step(0)->setMesh( meshCoarse );  
        
        exporterFine->step(0)->add("uFine", uFine );
        exporterCoarse->step(0)->add("uCoarse", uCoarse );
        
        exporterFine->save();
        exporterCoarse->save();
         
     }
    MeshFiles.close();
    
    
    



}
//-----------------------------------------
//-----------------------------------------


//-----------------------------------------
//-----------------------------------------

template< int PolynomialOrder>
gmsh_ptrtype NIRBTEST<PolynomialOrder>::createGeo(double hsize, std::string MeshFileName){
    double hCornersize1 = hsize/15.;
    double hCornersize2 = hsize/30.;
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = 2.2;" << "\n";
    ostr << "Mesh.CharacteristicLengthExtendFromBoundary=1;" << endl;
    ostr << "Mesh.CharacteristicLengthFromPoints=1;" << endl;
    ostr << "Mesh.ElementOrder=1;" << endl;
    ostr << "Mesh.SecondOrderIncomplete = 0;" << endl;
    ostr << "Mesh.Algorithm = 6;" << endl;
    ostr << "Mesh.OptimizeNetgen=1;" << endl;
    ostr << "// partitioning data" << endl;
    ostr << "Mesh.Partitioner=1;" << endl;
    ostr << "Mesh.NbPartitions=1;" << endl;
    ostr << "Mesh.MshFilePartitioned=0;" << endl;
    ostr << "Point(1) = {0,0,0.0," << hsize << "};" << endl;
    ostr << "Point(2) = {1,0,0.0," << hCornersize2 << "};" << endl;
    ostr << "Point(3) = {1,1,0.0," << hCornersize1 << "};"<< endl;
    ostr << "Point(4) = {0,1,0.0," << hCornersize2 << "};"<< endl;
    ostr << "Line(1) = {4,1};" << endl;
    ostr << "Line(2) = {1,2};" << endl;
    ostr << "Line(3) = {2,3};" << endl;
    ostr << "Line(4) = {3,4};" << endl;
    ostr << "Line Loop(5) = {1,2,3,4};" << endl;
    ostr << "Plane Surface(6) = {5};" << endl;
    ostr << "Physical Line(\"Dirichlet1\") = {1,2};" << endl;
    ostr << "Physical Line(\"Dirichlet2\") = {3};" << endl;
    ostr << "Physical Line(\"Dirichlet3\") = {4};" << endl;
    ostr << "Physical Surface(\"Mat1\") = {6};" << endl;
    std::ostringstream nameStr;
    nameStr.precision( 3 );
    nameStr << MeshFileName;
    //nameStr << MeshFileName.c_str();
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( ostr.str() );
    return gmshp;
    
}

//-----------------------------------------
//-----------------------------------------

template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder>::ComputeSnapshot(space_ptrtype Xh){

    int XhNdof = Xh->nLocalDof();
    double pi = 3.14159265;

    auto ui = Xh->element();
    for (int i =0; i<NbSnapshot;i++)
    {
        double theta = (i+1)*(muMax-muMin)/(113.);
        ui = blackbox(Xh,theta);
        std::string path = (boost::format("./Sol_%1%") %i).str() ;
        ui.save(_path=path);
   }

}
//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructStiffMatrixSnapshot(space_ptrtype Xh){

    auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );//Sparse F.E stiffness matrix D
 	form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=gradt(ui)*trans(grad(uj)));

	Eigen::MatrixXd S (NbSnapshot,NbSnapshot);//Dense RB stiffness matrix S

    for (int i = 0; i< NbSnapshot;i++)
    {
        ui.zero();
        std::string path = (boost::format("./Sol_%1%") %i).str() ;
        ui.load(_path=path);
        if (D->energy(ui,ui)==0.){
            std::cerr << "In 'ConstructStiffMatrixSnapshot': ERROR IN LOADING FILE " << path << endl;
            exit(0);
        }

        for(int j=0; j< i;j++)
        {
            uj.zero();
            path = (boost::format("./Sol_%1%") %j).str() ;
            uj.load(_path=path); 
            if (D->energy(uj,uj)==0.){
                std::cerr << "In 'ConstructStiffMatrixSnapshot': ERROR IN LOADING FILE " << path << endl;
                exit(0);
            }
            double Sij = D->energy(ui,uj);
            S(i,j) = Sij;
            S(j,i) = Sij;
        }
        double Sii = D->energy(ui,ui);
        S(i,i) = Sii;
    }
    
    std::string filename = (boost::format("StiffMatrixP%1%")%PolynomialOrder).str();
    std :: ofstream fileMat(filename);
    fileMat << NbSnapshot << endl;
    for (int i = 0;i<NbSnapshot;i++){
        for (int j =0;j<NbSnapshot;j++){
            fileMat << S(i,j) << endl;
        }
    }
    fileMat.close();

 	return S;
}
 

//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructStiffMatrixRB(space_ptrtype Xh, std::string filename){
	auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE stiffness matrix
	form2( _test=Xh, _trial=Xh, _matrix=D ) =
    integrate( _range=elements(Xh->mesh()), _expr=gradt(ui)*trans(grad(uj)));

	Eigen::MatrixXd S (sizeRB,sizeRB); //Dense RB stiffness matrix

	for (int i = 0; i< sizeRB;i++)
    {
        ui.zero();
        std::string path = (boost::format(filename+"%1%") %i).str() ;
		ui.load(_path=path);
        if (D->energy(ui,ui)==0.){
            std::cerr << "In 'ConstructStiffMatrixRB': ERROR IN LOADING FILE " << path << endl;
            exit(0);
        }

		for(int j=0; j< i;j++)
        {
            uj.zero();
            path = (boost::format(filename+"%1%") %j).str() ;
			uj.load(_path=path);
            if (D->energy(uj,uj)==0.){
                std::cerr << "In 'ConstructStiffMatrixRB': ERROR IN LOADING FILE " << path << endl;
                exit(0);
            }
			double Sij = D->energy(ui,uj);
			S(i,j) = Sij;
			S(j,i) = Sij;
		}
 		double Sii = D->energy(ui,ui);
		S(i,i) = Sii;
	}

	return S;
}
  
 
//--------------------------------------------------------------
//--------------------------------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructMassMatrixRB(space_ptrtype Xh, std::string filename){
	auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE mass matrix
	form2( _test=Xh, _trial=Xh, _matrix=D ) =
	integrate( _range=elements(Xh->mesh()), _expr=idt(ui)*id(uj));

	Eigen::MatrixXd S (sizeRB,sizeRB); //Dense RB mass matrix

	for (int i = 0; i< sizeRB;i++)
    {
        ui.zero();
        std::string path = (boost::format(filename+"%1%") %i).str() ;
		ui.load(_path=path);
        if (D->energy(ui,ui)==0.){
            std::cerr << "In 'ConstructMassMatrixRB': ERROR IN LOADING FILE " << path << endl;
            exit(0);
        }

		for(int j=0; j< i;j++)
        {
            uj.zero();
            path = (boost::format(filename+"%1%") %j).str() ;
			uj.load(_path=path);
            if (D->energy(uj,uj)==0.){
                std::cerr << "In 'ConstructMassMatrixRB': ERROR IN LOADING FILE " << path << endl;
                exit(0);
            }
			double Sij = D->energy(ui,uj);
			S(i,j) = Sij;
			S(j,i) = Sij;
		}
		double Sii = D->energy(ui,ui);
		S(i,i) = Sii;
	}
	return S;

}



 

//-----------------------------------------
//-----------------------------------------
 template< int PolynomialOrder>
 void NIRBTEST<PolynomialOrder> ::ChooseRBFunction(space_ptrtype Xh){
    
     //std :: cout << "Inside ChooseRBFunction" << endl;
     Eigen::MatrixXd S (NbSnapshot,NbSnapshot);//Dense Stiffness Matrix
     if (!Sampling){
         std::ifstream fileMat((boost::format("StiffMatrixP%1%")%PolynomialOrder).str());
         if(!fileMat){
             std::cerr << "Error in reading file " << (boost::format("StiffMatrixP%1%")%PolynomialOrder).str() <<  "does not exist " <<endl;
             std::cout << "Need to re-build the Stiffness Matrix -> SAMPLING SET TO 1 " << endl;
             Sampling = 1;
         }
         else
         {
             int Itemp;
             fileMat >> Itemp;
             if (Itemp != NbSnapshot){
                 std :: cerr << "Error in reading file" << (boost::format("StiffMatrixP%1%")%PolynomialOrder).str() <<  " --  " << Itemp << " != NbSnapshot(" << NbSnapshot << ")" <<  endl;
                 std::cout << "Need to re-build the Stiffness Matrix -> SAMPLING SET TO 1 " << endl;
                 Sampling = 1;
             }
             else{
                 for (int i=0;i<NbSnapshot;i++){
                     for (int j =0;j<NbSnapshot;j++){
                         double Dtemp;
                         fileMat >> Dtemp;
                         S(i,j) = Dtemp; 
                     }
                 }
             }
         }
         
     }
     if (Sampling)
     {
         S = ConstructStiffMatrixSnapshot(Xh);
     }
     Eigen::EigenSolver <Eigen::MatrixXd> eigen_solver(S);
     
     //int nb_EigenValue = eigen_solver.eigenvalues().size();
     //eigen_solver.eigenvectors().col(i)  =  eigenvector #i
     
     std :: cout << "Computation of the eigenvalues of the stiffness matrix S (NbSnapshot,NbSnapshot) : done" << endl;
     
     //Sorting N = "SizeBR" F.E solutions uh(mu_k) the more represented in the eigenvectors
     //associated to the N largest eigenvalues
     
     
     
     
     // saving the index's number to identify the NIRB basis functions
     std::string path = (boost::format("./IndBR%1%") %sizeRB).str() ;
     std::ofstream find(path);
     if (!find){
         std :: cerr <<" 'ChooseRBFunction routine' - ERROR IN OPENING FILE:" << path << endl;
         exit(0);
     }
     
     std::vector<int> Uind(NbSnapshot);
     for(int i=0;i<NbSnapshot;i++)
     {
         Uind[i] = 0; //Set to zero if u(Uind[i]) is not NIRB basis function.
     }
     
     
     int Ndof = Xh->nLocalDof(); 
     auto ui = Xh->element();
     auto uj = Xh->element();
     
     
     auto D = M_backend->newMatrix( _test=Xh , _trial=Xh  ); //Sparse FE mass matrix
     form2( _test=Xh , _trial=Xh, _matrix=D ) = integrate( _range=elements(Xh->mesh()), _expr=id(ui)*idt(uj));
     
     
     
     Eigen::MatrixXd Mu (Ndof,sizeRB);
     double Tol = 1e-30;
     std::ofstream fout("fout");
     //Sorting eigenvector #i
     std::vector<double>Vi(NbSnapshot);
     
     for (int i=0;i<sizeRB;i++)
     {
         for(int ri = 0; ri < NbSnapshot;ri++)
         {
             //If u(Uind[ri]) is a NIRB basis function, the "Uind[ri]" component of Vi is set to zero
             //to avoid it.
             if (Uind[ri] == 0){
                 Vi[ri] = std::abs(real(eigen_solver.eigenvectors().col(i)[ri]));
             }
             else {
                 Vi[ri] = 0;
             }
         }
         //Using max_element to find the position of the maximum component of Vi
         
         int IndMax = std::distance(Vi.begin(),boost::max_element(Vi)); 
         
         //copy u[IndMax] in Mu
         ui.zero();
         path = (boost::format("./Sol_%1%") %IndMax).str();
         ui.load(_path=path);
         //cout << "LOAD DONE " << endl;
         if (D->energy(ui,ui)==0.){
             std::cerr <<"'ChooseRBFunction routine': ERROR IN LOADING FILE  " << path << endl;
             exit(0);
         }
         for (int rj = 0; rj < Ndof; rj++){
             Mu(rj,i) = ui(rj);
         }
         if (i>0){
             double normL2 = OrthogonalisationRBFunctionL2GrammSchmidt(Xh,Mu,i);
             //std :: cout << "i = " << i << " --- NormU(" << IndMax << ") = " << normL2 << endl;
             fout << "i = " << i << " --- NormU(" << IndMax << ") = " << normL2 << endl;
             if (normL2 > Tol){
                 Uind[IndMax] = i;
                 find << IndMax << endl;
             }
             else{
                 Uind[IndMax] = -1;
                 i--;
             }
         }
         else
         {
             Uind[IndMax] = i;
             find << IndMax << endl;
         }
         
     }
     find.close();
     fout.close();
     
     
     
     //OrthogonalisationRBFunctionL2GrammSchmidt(Xh,Mu);
     
     
     //save in a file the matrix A
     for (int i=0;i<sizeRB;i++)
     {
         ui.zero();
         for (int j =0;j<Ndof;j++)
         {
             ui(j) = Mu(j,i);
         }
         double L2norm = std::sqrt(D->energy(ui,ui));
         for (int j =0;j<Ndof;j++)
         {
             ui(j) = Mu(j,i)/L2norm;
         }
         ui.save(_path=(boost::format("./Sol_OR%1%") %i).str());
     }
       
     
 }

//-----------------------------------------
template< int PolynomialOrder>
double NIRBTEST <PolynomialOrder>  :: OrthogonalisationRBFunctionL2GrammSchmidt(space_ptrtype Xh, Eigen::MatrixXd & MU, int n){
    
    //We suppose that the n-1 first function has been already orthogonalized
    auto un = Xh->element();
    auto uk = Xh->element();
    
    auto D = M_backend->newMatrix( _test=Xh , _trial=Xh  ); //Sparse FE mass matrix
    form2( _test=Xh , _trial=Xh, _matrix=D ) = integrate( _range=elements(Xh->mesh()), _expr=id(un)*idt(uk));
    
    int Ndof = Xh->nLocalDof();
    
    double Dtemp1,Dtemp2,Dtemp3,Dtemp4;
    for (int j=0;j<Ndof;j++){
        un(j) = MU(j,n);
    }
    for (int k=0;k<n;k++){
        for (int j=0;j<Ndof;j++){
            uk(j) = MU(j,k);
        }
        Dtemp1 = D->energy(un,uk); 
        Dtemp2 = D->energy(uk,uk);
        Dtemp3 = -Dtemp1/Dtemp2; 
        un.add(Dtemp3,uk);
        
    }
    for (int j = 0; j<Ndof;j++){ 
        MU(j,n) = un(j);
    }
    Dtemp4 = D->energy(un,un);
    return Dtemp4;
}


//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> :: OrthogonalisationRBFunctionL2GrammSchmidt
	      (space_ptrtype Xh,Eigen::MatrixXd & M ){

   double Ndof = Xh->nLocalDof();
   double Dtemp,Dtemp1,Dtemp2,Dtemp3, L2norm;
   auto ui = Xh->element();
   auto uk = Xh->element();

   auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE mass matrix
     form2( _test=Xh, _trial=Xh, _matrix=D ) =
	integrate( _range=elements(Xh->mesh()), _expr=idt(ui)*id(uk));

//    std::cout << "In 'OrthogonalisationRBFunctionL2GrammSchmidt' Warning : PRE-ORTHOGONALISATION NOT DONE !" << endl;
//
     Eigen::MatrixXd MTemp(Ndof,sizeRB);
     MTemp = M;

     // In MTemp non orthogonalized vector
     // In M orthogonalized vector
     for (int i=1;i<sizeRB;i++)
     {
         ui.zero();
         for (int j =0;j<Ndof;j++)
         {
             ui(j) = MTemp(j,i) ;
         }
         for (int k=0;k<i;k++)
         {
             for (int j =0;j<Ndof;j++)
             {
                 uk(j) = M(j,k);
             }
             Dtemp1 = D->energy(ui,uk);
             Dtemp2 = D->energy(uk,uk);
             Dtemp3 = - Dtemp1/Dtemp2;
             ui.add(Dtemp3,uk);
         }
         for (int j=0;j<Ndof;j++){
             M(j,i) = ui(j);
         }
     }

    for (int i=0;i<sizeRB;i++)
    {
        ui.zero();
        for (int j =0;j<Ndof;j++)
        {
            ui(j) = M(j,i);
        }
        L2norm = D->energy(ui,ui);
        L2norm = std::sqrt(L2norm);
        //	std :: cout << "i = " << i << "L2norm = " << L2norm << endl;
        for (int j =0;j<Ndof;j++){
            M(j,i) = ui(j)/L2norm;
        }
    }

    //  std::cout << "In 'OrthogonalisationRBFunctionL2GrammSchmidt' : L2-Normalisation done " << endl;

}
//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> ::OrthogonalisationRBFunction(space_ptrtype Xh){
    //First step : L2 -Pre-orthogonalization and normalization using a Gramm-schmidt method
    
   
    auto ui = Xh->element();
    auto uj = Xh->element();
    
    auto D = M_backend->newMatrix( _test=Xh , _trial=Xh ); //Sparse FE mass matrix
    form2( _test=Xh , _trial=Xh , _matrix=D ) =
    integrate( _range=elements(Xh ->mesh()), _expr=id(uj)*idt(ui));
    
    int Ndof = Xh->nLocalDof();
    Eigen :: MatrixXd M(Ndof,sizeRB);
    std::string path;
    int TindBR;
    
    for (int i=0;i<sizeRB;i++)
    {
        ui.zero();
        std::string path = (boost::format("./Sol_OR%1%") %i).str()  ;
        ui.load(_path=path);
        if( D->energy(ui,ui) == 0.){
            std::cerr << " In OrthogonalisationRBFunction : Sol_OR" << i << " is equal to zero " << endl;
            exit(0);
        }
        for (int j=0;j<Ndof;j++)
        {
            M(j,i) = ui(j);
        }
    }
    
    
    
    //Second step : Solving a generalized eigenvalue problem
    //A*Epsilon_i = \lambda_i*B*Epsilon_i
    //A : Dense stiffness matrix in the pre-orthogonalized nirb basis
    //B : Dense mass matrix in the pre-orthogonalized nirb basis  = Idendity matrix
    
    //Construction of the matrix A and B
    
    //     M.resize(sizeRB,sizeRB);
    Eigen :: MatrixXd A (sizeRB,sizeRB);
    Eigen :: MatrixXd B (sizeRB,sizeRB);
    Eigen :: MatrixXd C (sizeRB,sizeRB);
    
    A = ConstructStiffMatrixRB(Xh,"./Sol_OR");
    B = ConstructMassMatrixRB(Xh,"./Sol_OR");
    //     Eigen:: GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
    C = B.inverse()*A;
    
    Eigen:: EigenSolver<Eigen::MatrixXd> es;
    
    es.compute(C);
    
    A.resize(Ndof,sizeRB);//to upload the final NIRB function after H1-L2 orthogonalisation
    A.setZero(Ndof,sizeRB);
    for (int i = 0; i < sizeRB;i++)
    {
        for (int j =0;j< Ndof;j++)
        {
            for (int k =0;k< sizeRB;k++){
                A(j,i)=A(j,i)+ real(es.eigenvectors().col(i)[k])*M(j,k);
            }
        }
    }
    
    
    OrthogonalisationRBFunctionL2GrammSchmidt(Xh,A);
    
    //save in a file the final NIRB basis functions
    for (int i=0;i<sizeRB;i++)
    {
        ui.zero();
        for (int j =0;j<Ndof;j++)
        {
            ui(j) = A(j,i);
        }
        path = (boost::format("./RB%1%NIRB_BasisFile_%2%") %sizeRB%i).str() ;
        ui.save(_path=path);
        
    }
    
}
////-----------------------------------------
////-----------------------------------------
//template< int PolynomialOrder>
//void NIRBTEST<PolynomialOrder> ::OrthogonalisationRBFunction(space_ptrtype Xh){
//     //First step : L2 -Pre-orthogonalization and normalization using a Gramm-schmidt method
//
//    int Ndof = Xh->nLocalDof();
//    Eigen :: MatrixXd M(Ndof,sizeRB);
//    
//    int TindBR;
//    auto ui = Xh->element();
//    auto uj = Xh->element();
//    
//    auto D = M_backend->newMatrix( _test=Xh , _trial=Xh ); //Sparse FE mass matrix
//    form2( _test=Xh , _trial=Xh , _matrix=D ) =
//    integrate( _range=elements(Xh ->mesh()), _expr=id(uj)*idt(ui));
//    
//    //reading the index's number to identify the NIRB basis functions
//    std::string path = (boost::format("./IndBR%1%") %sizeRB).str() ;
//    std :: ifstream f_in(path);
//    if (!f_in)
//    {
//      std :: cerr << "In 'OrthogonalisationRBFunction' subroutine : Error in opening file -> " << path << endl;
//    }
//    for (int i=0;i<sizeRB;i++)
//    {
//        f_in >> TindBR;
//        ui.zero();
//        std::string path = (boost::format("./Sol_%1%") %TindBR).str() ;
//        ui.load(_path=path);
//        if( D->energy(ui,ui) == 0.){
//            std::cerr << " In OrthogonalisationRBFunction : Sol_" << i << " is equal to zero " << endl;
//            exit(0);
//        }
//        for (int j=0;j<Ndof;j++)
//        {
//            M(j,i) = ui(j);
//        }
//    }
//
//    double Dtemp;
//    OrthogonalisationRBFunctionL2GrammSchmidt(Xh,M);
//
//    //save in a file the matrix A
//    for (int i=0;i<sizeRB;i++)
//    {
//        ui.zero();
//        for (int j =0;j<Ndof;j++)
//        {
//            ui(j) = M(j,i);
//        }
//        std::string path = (boost::format("./Sol_OR%1%") %i).str() ;
//        ui.save(_path=path);
//    }
//
//    //Second step : Solving a generalized eigenvalue problem
//    //A*Epsilon_i = \lambda_i*B*Epsilon_i
//    //A : Dense stiffness matrix in the pre-orthogonalized nirb basis
//    //B : Dense mass matrix in the pre-orthogonalized nirb basis  = Idendity matrix
//
//    //Construction of the matrix A and B
//
//    //     M.resize(sizeRB,sizeRB);
//    Eigen :: MatrixXd A (sizeRB,sizeRB);
//    Eigen :: MatrixXd B (sizeRB,sizeRB);
//    Eigen :: MatrixXd C (sizeRB,sizeRB);
//
//    A = ConstructStiffMatrixRB(Xh,"./Sol_OR");
//    B = ConstructMassMatrixRB(Xh,"./Sol_OR");
//    //     Eigen:: GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
//    C = B.inverse()*A;
//
//    Eigen:: EigenSolver<Eigen::MatrixXd> es;
//
//    es.compute(C);
//
//    A.resize(Ndof,sizeRB);//to upload the final NIRB function after H1-L2 orthogonalisation
//    A.setZero(Ndof,sizeRB);
//    for (int i = 0; i < sizeRB;i++)
//    {
//       for (int j =0;j< Ndof;j++)
//       {
//           for (int k =0;k< sizeRB;k++){
//               A(j,i)=A(j,i)+ real(es.eigenvectors().col(i)[k])*M(j,k);
//           }
//       }
//    }
//
//
//    OrthogonalisationRBFunctionL2GrammSchmidt(Xh,A);
//
//    //save in a file the final NIRB basis functions
//    for (int i=0;i<sizeRB;i++)
//    {
//        ui.zero();
//        for (int j =0;j<Ndof;j++)
//        {
//            ui(j) = A(j,i);
//        }
//        path = (boost::format("./RB%1%NIRB_BasisFile_%2%") %sizeRB%i).str() ;
//        ui.save(_path=path);
//
//    }
//
//}
//---------------------------------------------------
//---------------------------------------------------
template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder> ::BuildNirbSolution(space_ptrtype XhFine,space_ptrtype XhCoarse, double param){

    //std :: cout << "In BuildNirbSolution" << endl;


    auto uNirb = XhFine->element();
    auto uCoarse = XhCoarse->element();
    auto uCoarseInterpolate = XhFine->element();
    auto ui = XhFine->element();




    auto D = M_backend->newMatrix( _test=XhFine, _trial=XhFine  ); //Sparse FE mass matrix
    form2( _test=XhFine, _trial=XhFine, _matrix=D ) =
    integrate( _range=elements(XhFine->mesh()), _expr=id(uCoarseInterpolate)*idt(ui));

//    auto D = M_backend->newMatrix( _test=XhCoarse, _trial=XhFine  ); //Sparse FE mass matrix
//    form2( _test=XhCoarse, _trial=XhFine, _matrix=D ) =
//    integrate( _range=elements(XhFine->mesh()), _expr=id(uCoarse)*idt(ui));

    uCoarse = blackbox(XhCoarse,param);//Computation of the coarse solution
    interpolate (XhFine,uCoarse,uCoarseInterpolate);
    //Interpolation of the coarse solution on the fine Mesh to compute the coefficiant BetaiH

    //Computation of the coefficiant \BetaiH = \int uCoarse*\Epsilon_i
    //with \Epsilon_i being the final Nirb basis functions

    //Reading Nirb basis function #i and building the Nirb solution
    Eigen :: VectorXd BetaiH(sizeRB);
    uNirb.zero();
    std::string path;

    for (int i =0;i<sizeRB;i++)
    {

        ui.zero();
        path =(boost::format("./RB%1%NIRB_BasisFile_%2%") %sizeRB%i).str()  ;
        ui.load(_path=path);
        if( D->energy(ui,ui) == 0.){
            std::cerr << "ERROR: In BuilNirbSolution : NIRB Basis function (" << i << ") equal to zero " << endl;
            exit(0);
        }
        //BetaiH[i] = integrate( _range=elements(XhFine->mesh()),_expr=(idv(uCoarse)*idv(ui)) ).evaluate()(0,0);
        //BetaiH[i] = D->energy(uCoarse,ui);

        BetaiH[i] = D->energy(uCoarseInterpolate,ui);
        uNirb.add( BetaiH[i], ui );
    }


    return uNirb;


}

//-----------------------------------------

//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> :: ConstructNIRB(space_ptrtype Xh){
   
    std :: cout << "OFFLINE PROCEDURE :  Construction of the 'non intruisive' reduced basis (nirb) " << endl;
    /*
    if(!Sampling){
        //check if the sampling function exist
        std::string filename;
        for (int i=0;i<NbSnapshot;i++){
            filename = (boost::format("./Sol_%1%/u.fdb") %i).str();
            auto ui = Xh->element();
            ui.load(filename);
            int NdofNirb;
            NdofNirb = ui.size();
            if (NdofNirb != Xh->nLocalDof()){
                std::cout <<"WARNING Error in  Sol_" << i << " Ndof not the same as in Xh=> Sampling parameter set to 1 " << endl;
                Sampling = 1;
                break;
            }
        }
    }
    */
    boost::timer ti;
    double Time_snapshot;
    if(Sampling){
        std::cout << "Sampling Procedure :" << endl;
        ComputeSnapshot(Xh);
        std::cout << "Computation of the " << NbSnapshot << " snapshots : done  -- ";
        std::cout << "Time per snapshot: " << Time_snapshot/NbSnapshot << endl;
    }
    Time_snapshot = ti.elapsed();
    ti.restart();
    
	ChooseRBFunction(Xh);
    double TimeChooseRB =  ti.elapsed();
	std::cout << "Choice of " << sizeRB << " reduced basis functions : done  -- ";
    std::cout << "Time: " << TimeChooseRB << endl;
    ti.restart();
    
	//Orthogonalisation de Gram-schmidt
	OrthogonalisationRBFunction(Xh);
    double Time_buildNIRBbasis = ti.elapsed();
	std::cout << "H1 and L2 orthogonalisation of the reduced basis functions : done  -- ";
    std::cout << "Time: " << Time_buildNIRBbasis << endl;
    double TotalTime = Time_buildNIRBbasis +TimeChooseRB + Time_snapshot;
    std::cout <<"Total time for the OFFLINE procedure: " << TotalTime << endl;

}


//---------------------------------------------------
template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder>::blackbox( space_ptrtype Xh, double param )
{

    auto u = Xh->element();
    auto v = Xh->element();
    // std::cout << "Xh ndof=" << Xh->nLocalDof() << "\n";
    value_type pi = M_PI;
    double penaltyTerm = 30.;

    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto velocity = vec(cst(cos(param)),cst(sin(param)));
    /*
    auto g = ( chi(abs(Px()-1) < 1e-10 )*Px()*Px()+
               chi(abs(Py()-1) < 1e-10 )*Py()*Py() );
    */
    
    auto g = (Px()*Px()*Py()*Py());
    auto F = M_backend->newVector( Xh );
    /*
    form1( _test=Xh, _vector=F ) =
        integrate( _range=boundaryfaces(Xh->mesh()),
                   _expr=g*(-0.01*dn(v)+penaltyTerm*id(v)/hFace())  );
     */
    auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );
    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=0.01*gradt(u)*trans(grad(v))+ (gradt(u)*velocity)*id(v) );

    /*
    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate( boundaryfaces(Xh->mesh()),
                   -0.01*dnt(u)*id(v)-0.01*dn(v)*idt(u)+penaltyTerm*id(v)*idt(u)/hFace());
    */
    
    form2( _test=Xh, _trial=Xh, _matrix=D ) += on(boundaryfaces(Xh->mesh()), u, F,g );
    backend_type::build()->solve( _matrix=D, _solution=u, _rhs=F );

    return u;
} // NIRBTEST::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    Environment env( argc, argv );
    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }
    /**
     * register the simgets
     */

    const int polynomialOrder =  app.vm()["polynomialOrder"].as<int>();


    if(polynomialOrder == 1)
    {
        app.add( new NIRBTEST<1>( app.vm(), app.about() ) );
    }
    else if(polynomialOrder == 2)
    {
        app.add( new NIRBTEST<2>( app.vm(), app.about() ) );
    }
    else if(polynomialOrder == 3)
    {
        app.add( new NIRBTEST<3>( app.vm(), app.about() ) );
    }
    else
    {
        throw std::logic_error( "Error with polynomialOrder variable, this application allows only P1 P2 and P3" );
    }

    /**
     * run the application
     */
    app.run();
}






