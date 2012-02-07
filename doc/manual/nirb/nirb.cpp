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
inline
po::options_description
makeOptions()
{
    po::options_description laplacianoptions("Laplacian options");
    laplacianoptions.add_options()
	// meshs parameters
        ("hfinsize", po::value<double>()->default_value( 0.1 ), "fine mesh size")
        ("hcoarsesize", po::value<double>()->default_value( 0.5 ), "coarse mesh size")
	
	// Reduced basis parameters        
	("NbSnapshot",po::value<int>()->default_value(100),"numbers of snapshot computed")
        ("sizeBR",po::value<int>()->default_value(15),"size of reduced basis")
	
	("muMin", po::value<double>()->default_value( 0 ), "angle in [0,pi/2]")
	("muMax", po::value<double>()->default_value( M_PI/2),"angle in [0,pi/2]")
        ;
    return laplacianoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
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


class NIRBTEST
    :
    public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 3;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 1
    typedef Simplex<2> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar> > basis_type;

    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;


    /**
     * Constructor
     */
    NIRBTEST( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        FineMeshSize( this->vm()["hfinsize"].as<double>() ),
        CoarseMeshSize( this->vm()["hcoarsesize"].as<double>() ),
	NbSnapshot( this->vm()["NbSnapshot"].as<int>() ),
	sizeBR(this->vm()["sizeBR"].as<int>()),
        muMin( this->vm()["muMin"].as<double>() ),
        muMax( this->vm()["muMax"].as<double>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );
    element_type blackbox( space_ptrtype Xh, double param);
    element_type BuildNirbSolution(space_ptrtype XhFine,space_ptrtype XhCoarse,double param);
    void ConstructNIRB (space_ptrtype Xh); 
    void CalculSnapshot(space_ptrtype Xh);
    void ChooseRBFunction(space_ptrtype Xh);
    void OrthogonalisationRBFunction(space_ptrtype Xh);
    void OrthogonalisationRBFunctionL2GrammSchmidt(space_ptrtype Xh,Eigen::MatrixXd &M );

    Eigen::MatrixXd  ConstructStiffMatrixSnapshot(space_ptrtype Xh);
    Eigen::MatrixXd  ConstructStiffMatrixRB(space_ptrtype Xh);
    Eigen::MatrixXd  ConstructMassMatrixRB(space_ptrtype Xh); 
    
private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double FineMeshSize;
    double CoarseMeshSize;

    // Reduced basis parameter
    int NbSnapshot,sizeBR;
    double muMin,muMax;
}; // NIRBTEST
const uint16_type NIRBTEST::Order;
void
NIRBTEST::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute NIRBTEST<" << 2 << ">\n";
    std::vector<double> X(6);
    X[0] = FineMeshSize;
    X[1] = CoarseMeshSize;
    X[2] = NbSnapshot;
    X[3] = sizeBR;
    X[4] = muMin;
    X[5] = muMax;
    std::vector<double> Y( 1 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
void
NIRBTEST::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "%1%/P%2%/" )
                                       % this->about().appName()
                                       % Order );


    mesh_ptrtype meshFine = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name="Omega1",
                                                       _usenames=true,
                                                       _shape="hypercube",
                                                       _dim=2,
                                                       _h=X[0] ) );
    mesh_ptrtype meshCoarse = createGMSHMesh( _mesh=new mesh_type,
                                         _desc=domain( _name="Omega2",
                                                       _usenames=true,
                                                       _shape="hypercube",
                                                       _dim=2,
                                                       _h=X[1] ));

    space_ptrtype XhFine = space_type::New( meshFine );
    space_ptrtype XhCoarse = space_type::New( meshCoarse );

    export_ptrtype exporterFine( export_type::New( this->vm(), "nirbInFine" ) );
    export_ptrtype exporterCoarse( export_type::New( this->vm(), "nirbInCoarse" ) );
    export_ptrtype exporter2Grid(export_type::New( this->vm(), "nirbOut"));

    exporterFine->step(0)->setMesh( meshFine );
    exporterCoarse->step(0)->setMesh( meshCoarse );
    exporter2Grid->step(0)->setMesh( meshFine );


    //STEP ONE : Construction of the "non intruisive reduced basis (nirb) functions"
    //Computation of the X[3] snapshots solution on Xhfine
    //Selection of X[2] fonctions to build the "reduced basis" using a POD technique
    //Orthogonalisation in L2 and H1 norm of "reduced basis function"
    //Save this final functions
    std:: cout << "STEP ONE  : Construction of the 'non intruisive' reduced basis (nirb) functions " << endl;
    ConstructNIRB(XhFine) ;	
     

   //STEP TWO : Approximation of the solution using the "nirb" functions for a choosen mu 
    	
    std:: cout << "STEP TWO " << endl; 
 
   
    
    double p = 5*M_PI/(2*9); 
    auto uFine = XhFine->element();
    auto uCoarse = XhCoarse->element();
    auto uNirb = XhFine->element();
    uFine = blackbox( XhFine, p );
    uCoarse = blackbox( XhCoarse, p);
    uNirb = BuildNirbSolution(XhFine,XhCoarse,p);
    std :: cout << "After BuildNirbSolution " << endl;
   
    exporterFine->step(0)->add("uFine", uFine ); 
    exporterCoarse->step(0)->add("uCoarse", uCoarse ); 
    exporter2Grid->step(0)->add( "uNirb", uNirb ); 
 
    exporterFine->save(); 
    exporterCoarse->save();  
    exporter2Grid->save(); 

}
//-----------------------------------------
//-----------------------------------------
void NIRBTEST::CalculSnapshot(space_ptrtype Xh){
  
	int XhNdof = Xh->nLocalDof(); 
	double pi = 3.14159265;
	char path [50], buf[3]; 
	auto ui = Xh->element();
        for (int i =0; i<NbSnapshot;i++){
		double theta = i*(muMax-muMax)/(113.);
		ui = blackbox(Xh,theta);
		strcpy(path,"./Sol_");
		sprintf(buf,"%d",i);
		strcat(path,buf);
		ui.save(_path=path);  
	}
}
//-----------------------------------------
//-----------------------------------------
Eigen::MatrixXd NIRBTEST :: ConstructStiffMatrixSnapshot(space_ptrtype Xh){ 
// 	std::cout << "Entering ConstructStiffMatrix routine "<<endl;
	char path [50], buf[3]; 
	auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );//Sparse F.E stiffness matrix D 
 	form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=gradt(ui)*trans(grad(uj)));
	
	Eigen::MatrixXd S (NbSnapshot,NbSnapshot);//Dense RB stiffness matrix S
	 	 
        for (int i = 0; i< NbSnapshot;i++){  
		strcpy(path,"./Sol_");
		sprintf(buf,"%d",i);
		strcat(path,buf);
		ui.load(_path=path);
	 
 		for(int j=0; j< i;j++){   
			strcpy(path,"./Sol_");
			sprintf(buf,"%d",j);
			strcat(path,buf);
			uj.load(_path=path); 
			double Sij = D->energy(ui,uj); 
			S(i,j) = Sij;
			S(j,i) = Sij;
		}	
 		double Sii = D->energy(ui,ui); 
		S(i,i) = Sii;	
	}  
	
// 	std::cout << "Quitting ConstructStiffMatrix routine " << endl;
	
 	return S;
}
 
//-----------------------------------------
//-----------------------------------------
Eigen::MatrixXd NIRBTEST :: ConstructStiffMatrixRB(space_ptrtype Xh){
	char path [50], buf[3]; 
	auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE stiffness matrix
	form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=gradt(ui)*trans(grad(uj)));
	
	Eigen::MatrixXd S (sizeBR,sizeBR); //Dense RB stiffness matrix
	
	for (int i = 0; i< sizeBR;i++){  
		strcpy(path,"./Sol_OR");
		sprintf(buf,"%d",i);
		strcat(path,buf);
		ui.load(_path=path);
	
		for(int j=0; j< i;j++){   
			strcpy(path,"./Sol_OR");
			sprintf(buf,"%d",j);
			strcat(path,buf);
			uj.load(_path=path); 
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
Eigen::MatrixXd NIRBTEST :: ConstructMassMatrixRB(space_ptrtype Xh){
	char path [50], buf[3]; 
	auto ui = Xh->element();
	auto uj = Xh->element();
	auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE mass matrix
	form2( _test=Xh, _trial=Xh, _matrix=D ) =
	integrate( _range=elements(Xh->mesh()), _expr=idt(ui)*id(uj));
	
	Eigen::MatrixXd S (sizeBR,sizeBR); //Dense RB mass matrix
	
	for (int i = 0; i< sizeBR;i++){  
		strcpy(path,"./Sol_OR");
		sprintf(buf,"%d",i);
		strcat(path,buf);
		ui.load(_path=path);
	
		for(int j=0; j< i;j++){   
			strcpy(path,"./Sol_OR");
			sprintf(buf,"%d",j);
			strcat(path,buf);
			uj.load(_path=path); 
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
void NIRBTEST ::ChooseRBFunction(space_ptrtype Xh){
   
	char path [50], buf[3];
	
	Eigen::MatrixXd S (NbSnapshot,NbSnapshot);//Dense Stiffness Matrix 
	S = ConstructStiffMatrixSnapshot(Xh);
	Eigen::EigenSolver <Eigen::MatrixXd> eigen_solver(S); 
	//int nb_EigenValue = eigen_solver.eigenvalues().size();
	//eigen_solver.eigenvectors().col(i)  =  eigenvector #i
	std :: cout << "Computation of the eigenvalues of the stiffness matrix S (NbSnapshot,NbSnapshot) s: OK" << endl;
	//Sorting N = "SizeBR" F.E solutions uh(mu_k) the more represented in the eigenvectors
	//associated to the N largest eigenvalues
	 
 	
	
	int *Lind = new int [NbSnapshot];
	int *Uind = new int [NbSnapshot];
	
	for (int i =0;i<NbSnapshot;i++){
	    Lind[i] = i;
	    Uind[i] = -1; 
	} 
	double memory;
	int compt,tempi;
	bool marqueur;
	
	
	
	strcpy(path,"./IndBR"); 
	// saving the index's number to identify the NIRB basis functions 
	sprintf(buf,"%d",sizeBR);
	strcat(path,buf);
	std::ofstream find(path); 
	if (!find){
	  std :: cerr <<" 'ChooseRBFunction routine' - Error in opening file :" << path << endl;
	} 
	//Sorting eigenvector #i 
	Eigen::VectorXd Vi(NbSnapshot); 
	for (int i=0;i<sizeBR;i++){ 
	  for(int ri = 0; ri < NbSnapshot;ri++){
	     Vi[ri] =  real(eigen_solver.eigenvectors().col(i)[ri]);
	     Lind[ri] = ri;
	  }
	 
	  for (int j = 1; j< NbSnapshot;j++){  
	    memory = Vi[j];
	    tempi = Lind[j];
	    compt = j-1;
	    marqueur = std::abs (Vi[compt])>std::abs(memory);
	    while(marqueur){
	      if (std::abs(Vi[compt])>std::abs(memory)){
		Vi[compt+1] = Vi[compt];
		Lind[compt+1]= Lind[compt,1];
		compt--;
		marqueur = true;
	      }
	      else marqueur = false;
	      if (compt<0)
		marqueur = false;
	    }
	    Vi[compt+1] = memory;
	    Lind[compt+1] = tempi;
	  } 
	  //Extracting reduced basis function #i + update of Uind
	  bool pasoki = true;
	  int rc = NbSnapshot-1;
	  while (pasoki && rc >0){
	    if (Uind [Lind[rc]] ==-1){
	      pasoki= false; 
	      Uind [Lind[rc]] = i;    
	      find << Lind[rc]<<endl; 
	    }
	    else rc--; 
	  }
	  
	}
 	find.close();  
    delete [] Lind;
    delete [] Uind;
}

//-----------------------------------------
//-----------------------------------------
void NIRBTEST :: OrthogonalisationRBFunctionL2GrammSchmidt
	      (space_ptrtype Xh,Eigen::MatrixXd & M ){
//      auto ui = Xh->element(); 
//      auto uj = Xh->element();
//      auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  ); //Sparse FE mass matrix
//      form2( _test=Xh, _trial=Xh, _matrix=D ) =
// 	integrate( _range=elements(Xh->mesh()), _expr=idt(ui)*id(uj));
   std::cout << " In OrthogonalisationRBFunctionL2GrammSchmidt :ATTENTION : sORTHOGONALISATION NOT DONE !" << endl;
 
   
}

//-----------------------------------------
//-----------------------------------------
void NIRBTEST ::OrthogonalisationRBFunction(space_ptrtype Xh){
     //First step : L2 -Pre-orthogonalization and normalization using a Gramm-schmidt method
     
    int Ndof = Xh->nLocalDof();
    char path [50], buf[3];
    int TindBR;
    auto ui = Xh->element(); 
    Eigen :: MatrixXd M(Ndof,sizeBR);
    strcpy(path,"./IndBR"); 
//     reading the index's number to identify the NIRB basis functions 
    sprintf(buf,"%d",sizeBR);
    strcat(path,buf);
    std :: ifstream f_in(path);
    if (!f_in){
      std :: cerr << "In 'OrthogonalisationRBFunction' subroutine : Error in opening file -> " << path << endl;
    }
    for (int i=0;i<sizeBR;i++){
      f_in >> TindBR; 
      strcpy(path,"./Sol_");
      sprintf(buf,"%d",TindBR);
      strcat(path,buf);
      ui.load(_path=path);
      for (int j=0;j<Ndof;j++){
	M(j,i) = ui(j);
      }
    } 
     
//      element_ptrtype Pu (new element_type(u.functionSpace()));
//      Pu->zero();
     double Dtemp;
     OrthogonalisationRBFunctionL2GrammSchmidt(Xh,M);
     
     //save in a file the matrix A
     for (int i=0;i<sizeBR;i++){
	ui.zero();
	for (int j =0;j<Ndof;j++){
	  Dtemp = M(j,i);
	  ui(j)=Dtemp;
	}
	strcpy(path,"./Sol_OR_");
	sprintf(buf,"%d",i);
	strcat(path,buf);
	ui.save(_path=path);
    }
     
     
     
     
     //Second step : Solving a generalized eigenvalue problem
     //A*Epsilon_i = \lambda_i*B*Epsilon_i
     //A : Dense stiffness matrix in the pre-orthogonalized nirb basis
     //B : Dense mass matrix in the pre-orthogonalized nirb basis  = Idendity matrix
     
     //Construction of the matrix A and B
     
    M.resize(sizeBR,sizeBR);  
    Eigen :: MatrixXd A (sizeBR,sizeBR);
    Eigen :: MatrixXd B (sizeBR,sizeBR); 
    Eigen :: MatrixXd C (sizeBR,sizeBR);
    A = ConstructStiffMatrixRB(Xh);
    B = ConstructMassMatrixRB(Xh);
    
//     Eigen:: GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
    M = B.inverse()*A;
    Eigen:: EigenSolver<Eigen::MatrixXd> es (M);
    M.resize(Ndof,sizeBR);//to download the NIRB function after the Gramm-schmidt orthogonalisation
    
    for (int i=0;i<sizeBR;i++){
	ui.zero();
	strcpy(path,"./Sol_OR_");
	sprintf(buf,"%d",i);
	strcat(path,buf);
	ui.load(_path=path);
	for (int j =0;j<Ndof;j++){
	   M(j,i) = ui(j);;
	}
    }
     
    
    A.resize(Ndof,sizeBR);//to upload the final NIRB function after H1-L2 orthogonalisation
    A.setZero(Ndof,sizeBR);
    for (int i = 0; i < sizeBR;i++){
      for (int j =0;j< Ndof;j++){
	for (int k =0;k< sizeBR;k++){
	    A(j,i)=A(j,i)+ real(es.eigenvectors().col(i)[k])*M(j,k);
	}
      }
    }
    
     OrthogonalisationRBFunctionL2GrammSchmidt(Xh,A);
    
    
    //save in a file the final NIRB basis functions
     for (int i=0;i<sizeBR;i++){
	ui.zero();
	for (int j =0;j<Ndof;j++){
	  Dtemp = A(j,i);
	  ui(j)=Dtemp;
	}
	strcpy(path,"./NIRB_BasisF_");
	sprintf(buf,"%d",i);
	strcat(path,buf);
	ui.save(_path=path);
    }
     
    

    
}
//---------------------------------------------------
//---------------------------------------------------
NIRBTEST::element_type NIRBTEST ::BuildNirbSolution(space_ptrtype XhFine,space_ptrtype XhCoarse, double param){
  
  std :: cout << "In BuildNirbSolution" << endl;
  char path [50], buf[3];
  
  auto uNirb = XhFine->element();
  auto uCoarse = XhCoarse->element();
  auto uCoarseInterpolate = XhFine->element();
  auto ui = XhFine->element();
 
  auto D = M_backend->newMatrix( _test=XhFine, _trial=XhFine  ); //Sparse FE mass matrix
	form2( _test=XhFine, _trial=XhFine, _matrix=D ) =
	integrate( _range=elements(XhFine->mesh()), _expr=idt(uNirb)*id(ui));
   
  uCoarse = blackbox(XhCoarse,param);//Computation of the coarse solution 
  interpolate (XhFine,uCoarse,uCoarseInterpolate); 
  //Interpolation of the coarse solution on the fine Mesh to compute the coefficiant BetaiH
  
  //Computation of the coefficiant \BetaiH = \int uCoarse*\Epsilon_i 
  //with \Epsilon_i being the final Nirb basis functions
   
  double BetaiH;
  uNirb.zero();
  for (int i =0;i<sizeBR;i++){
    ui.zero();
     //Reading Nirb basis function #i and building the Nirb solution  
    
    strcpy(path,"./NIRB_BasisF_");
    sprintf(buf,"%d",i);
    strcat(path,buf);
    ui.load(_path=path);
    
    BetaiH  = D->energy(uCoarseInterpolate,ui);
    for (int j=0;j<XhFine->nLocalDof();j++){
      BetaiH = ui(j)*BetaiH;
      ui(j) = BetaiH;
    }
    uNirb = uNirb + ui; 
  }
  return uNirb;
  
  
}

//-----------------------------------------

//-----------------------------------------
//-----------------------------------------
void NIRBTEST :: ConstructNIRB(space_ptrtype Xh){
	
	CalculSnapshot(Xh);
	std :: cout << "Computation of the " << NbSnapshot << " snapshots : OK " << endl;
	ChooseRBFunction(Xh); 
	std:: cout << "Choice of " << sizeBR << " Reduced basis functions :OK " << endl; 
	//Orthogonalisation de Gram-schmidt
	OrthogonalisationRBFunction(Xh);	
}


//---------------------------------------------------
NIRBTEST::element_type NIRBTEST::blackbox( space_ptrtype Xh, double param )
{

    auto u = Xh->element();
    auto v = Xh->element();
   // std::cout << "Xh ndof=" << Xh->nLocalDof() << "\n";
    value_type pi = M_PI;

    //! deduce from expression the type of g (thanks to keyword 'auto')
    auto velocity = vec(cst(cos(param)),cst(sin(param)));
    auto g = ( chi(abs(Px()-1) < 1e-10 )*Px()*Px()+
               chi(abs(Py()-1) < 1e-10 )*Py()*Py() );
    auto F = M_backend->newVector( Xh );

    form1( _test=Xh, _vector=F ) =
        integrate( _range=boundaryfaces(Xh->mesh()),
                   _expr=g*(-0.01*dn(v)+30*id(v)/hFace())  );

    auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );
    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements(Xh->mesh()), _expr=0.01*gradt(u)*trans(grad(v))+ (gradt(u)*velocity)*id(v) );

    form2( _test=Xh, _trial=Xh, _matrix=D ) +=
        integrate( boundaryfaces(Xh->mesh()),
                   -0.01*dnt(u)*id(v)-0.01*dn(v)*idt(u)+30*id(v)*idt(u)/hFace());

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
    app.add( new NIRBTEST( app.vm(), app.about() ) );

    /**
     * run the application
     */
    app.run();
}






