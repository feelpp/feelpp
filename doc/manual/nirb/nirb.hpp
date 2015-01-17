//
//  nirb.hpp
//


#ifndef _nirb_hpp
#define _nirb_hpp



#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/loadmesh.hpp>
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
#include<boost/filesystem.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

using namespace std;

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
    //! vector type associated with backend
    typedef typename backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;


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
    typedef std::vector<element_type> vector_of_element_type;

    /**
     * Constructor
     */
    NIRBTEST()
        :
        super(),
        M_backend( backend_type::build( soption("backend")) ),
        CoarseMeshSize( doption("hcoarsesize"   ) ),
        FineMeshSize(   doption("hfinsize"      ) ),
        ReadingMeshes(  ioption("ReadingMeshes" ) ),
        NbSnapshot(     ioption("NbSnapshot"    ) ),
        sizeRB(         ioption("sizeRB"        ) ),
        muMin(          doption("muMin"         ) ),
        muMax(          doption("muMax"         ) ),
        mu(             doption("mu"            ) ),
        Sampling(       ioption("Sampling"      ) ),
        SamplingCoarse( ioption("SamplingCoarse") ),
        Offline(        ioption("Offline"       ) ),
        ComputeError(   ioption("ComputeError"  ) )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );


    // computation of the "blackbox" solution
    element_type blackbox( space_ptrtype Xh, double param );


    //Construction of  NIRB basis

    void ComputeSnapshot( space_ptrtype Xh, std::string filename );

    void ChooseRBFunction( space_ptrtype Xh,
                           vector_of_element_type &VuBasis,
                           sparse_matrix_ptrtype const &MassMatrix,
                           sparse_matrix_ptrtype const &StiffMatrix );

    void OrthogonalisationRBFunction( space_ptrtype Xh,
                                      vector_of_element_type &Vu,
                                      vector_of_element_type &M_Vu,
                                      sparse_matrix_ptrtype const &StiffMatrix,
                                      sparse_matrix_ptrtype const &MassMatrix );

    void OrthogonalisationRBFunctionL2GrammSchmidt( space_ptrtype Xh,
            vector_of_element_type &Vu,
            sparse_matrix_ptrtype const &MassMatrix );

    double OrthogonalisationRBFunctionL2GrammSchmidt( space_ptrtype Xh,
            vector_of_element_type &Vu,
            int n,
            sparse_matrix_ptrtype const &MassMatrix );


    Eigen::MatrixXd  ConstructStiffMatrixSnapshot( space_ptrtype Xh
            ,sparse_matrix_ptrtype  const &StiffMatrix );

    Eigen::MatrixXd  ConstructStiffMatrixRB( space_ptrtype Xh,
            vector_of_element_type &Vu,
            sparse_matrix_ptrtype const &StiffMatrix );

    Eigen::MatrixXd  ConstructMassMatrixRB( space_ptrtype Xh,
                                            vector_of_element_type &Vu,
                                            sparse_matrix_ptrtype const &MassMatrix );

    void ConstructNIRB ( space_ptrtype Xh,
                         vector_of_element_type  &M_VNirbBasis,
                         vector_of_element_type  &VNirbBasis );


    //Construction of NIRB solution
    element_type  BuildNirbSolutionWithoutPostProcess( space_ptrtype XhFine,
            element_type uCoarseInterpolate,
            vector_of_element_type const &M_VNirbBasis,
            vector_of_element_type const &VNirbBasis,
            Eigen::VectorXd & BetaiH );
    element_type BuildCoarseInterpolation( space_ptrtype XhFine,
                                           space_ptrtype XhCoarse,
                                           double param );

    element_type BuildNirbSolution( space_ptrtype XhFine,
                                    space_ptrtype XhCoarse,
                                    vector_of_element_type const &M_VNirbBasis,
                                    vector_of_element_type const &VNirbBasis,
                                    double param,
                                    Eigen::VectorXd BetaiH );

    //Post-processing of NIRB solution
    Eigen::MatrixXd BuildBetaH( space_ptrtype XhFine,
                                space_ptrtype XhCoarse,
                                vector_of_element_type  &M_VNirbBasis );

    Eigen::MatrixXd BuildBetah( space_ptrtype XhFine,
                                vector_of_element_type  &M_VNirbBasis );

    Eigen::MatrixXd BuildPostProcessMatrix( space_ptrtype XhFine,
                                            space_ptrtype XhCoarse,
                                            vector_of_element_type  &M_VNirbBasis );

    element_type BuildNirbSolutionWithPostProcess( space_ptrtype XhFine,
            space_ptrtype XhCoarse,
            vector_of_element_type  &M_VNirbBasis,
            vector_of_element_type  &VNirbBasis,
            Eigen::VectorXd &BetaiuH,
            element_type &uCoarseInterpolate );


    //Create geometry
    gmsh_ptrtype createGeo( double hsize,std::string MeshFileName );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double CoarseMeshSize;
    double FineMeshSize;
    int ReadingMeshes; // By Default = 0

    //! Reduced basis parameter
    int NbSnapshot,sizeRB;
    double muMin,muMax,mu;

    // Paramater to set OFF or ON the "offline" procedure

    int Sampling;
    int SamplingCoarse;
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


//-----------------------------------------
//-----------------------------------------
//const uint16_type NIRBTEST::Order;
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder>::run()
{

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute NIRBTEST<" << PolynomialOrder << ">\n";
    std::vector<double> X( 12 );
    X[0] = FineMeshSize;
    X[1] = CoarseMeshSize;
    X[2] = ReadingMeshes;
    X[3] = NbSnapshot;
    X[4] = sizeRB;
    X[5] = muMin;
    X[6] = muMax;
    X[7] = mu;
    X[8] = Sampling;
    X[9] = SamplingCoarse;
    X[10] = Offline;
    X[11] = ComputeError;
    std::vector<double> Y( 1 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{

    std::cout<< "Size reduced basis : " << sizeRB <<std::endl;
    std::cout<< "Polynomial order : P" << PolynomialOrder <<std::endl;

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "%1%/P%2%/" )
                                       % this->about().appName()
                                       % PolynomialOrder );


    mesh_ptrtype meshExtraCoarse, meshCoarse, meshFine ;

    if ( ReadingMeshes )
    {
        meshCoarse  =  loadGMSHMesh( _mesh=new mesh_type,
                                     _filename="/home/chakir/Mesh/Mesh_H1.msh",
                                    _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
        meshFine  =  loadGMSHMesh( _mesh=new mesh_type,
                                  _filename="/home/chakir/Mesh/Mesh_H2.msh",
                                  _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK  );
        std::cout << "Meshes read in file : done " <<std::endl;
    }

    else
    {
        std::cout << "Meshes build using 'createGMSHMesh' " <<std::endl;
        std::cout << "Coarse mesh size : " << X[1] << std::endl;
        meshCoarse = createGMSHMesh( _mesh= new mesh_type,
                                     _desc = createGeo( X[1],"Mesh_Coarse" ),
                                     _refine=1 );

        if ( X[1]> X[0] )
        {

            std::cout << "Fine mesh size : " << X[0] << std::endl;
            meshFine = createGMSHMesh( _mesh= new mesh_type,
                                       _desc = createGeo( X[0],"Mesh_Fine" ),
                                       _refine=1 );
        }

        else
        {

            std::cout << "WARNING : 'hfinesize > 'hcoarsesize' --> Fine mesh size set to hcoarsesize/2 : " << X[1]/2 << std::endl;
            meshFine = createGMSHMesh( _mesh= new mesh_type,
                                       _desc = createGeo( X[1]/2,"Mesh_Fine" ),
                                       _refine=1 );
        }

    }



    space_ptrtype XhFine = space_type::New( meshFine );
    //Fine F.E space used to build the NIRB basis

    space_ptrtype XhCoarse = space_type::New( meshCoarse );



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


    auto ui = XhFine->element();
    ui.zero();
    vector_of_element_type VNirbBasis( sizeRB,ui );
    vector_of_element_type MassMat_x_VuNirb( sizeRB,ui );

    if ( !Offline )
    {
        std :: string path;

        //Checking if the NIRB basis files exist
        //WARNING : We also must check if the basis have been computed on the right FE space (Not done yet)
        for ( int i = 0 ; i < sizeRB; i++ )
        {
            path = ( boost::format( "./RB%1%NIRB_BasisFile_%2%" )%sizeRB%i ).str();
            ui.load( _path=path );

            if ( ui.l2Norm()== 0. )
            {
                std::cout <<"WARNING Error : RB "  << sizeRB << "NIRB_BasisFile_" << i << " is equal to zero => OFFline parameter set to 1 " <<std::endl;
                Offline = 1;
                break;
            }

            VNirbBasis[i] = ui;
            path =( boost::format( "./RB%1%D_Nirb_Basis%2%File" ) %sizeRB%i ).str()  ;
            ui.load( _path=path );

            if ( ui.l2Norm()== 0. )
            {
                std::cout <<"WARNING Error : RB" << sizeRB << "D_NIRB_Basis" << i << "File is equal to zero => Re-calculating M*Ui" <<std::endl;
                auto Mui = M_backend->newVector( XhFine );
                auto Ui = M_backend->newVector( XhFine );

                auto uj = XhFine->element();

                sparse_matrix_ptrtype MassMatrix = M_backend->newMatrix( _test=XhFine, _trial=XhFine ); //Sparse FE mass matrix
                form2( _test=XhFine , _trial=XhFine, _matrix= MassMatrix ) =
                    integrate( _range=elements( XhFine->mesh() ), _expr=id( uj )*idt( ui ) );

                *Ui = VNirbBasis[i];
                Ui->close();
                MassMatrix->multVector( Ui,Mui );
                MassMat_x_VuNirb[i] = *Mui;
            }

            else
            {
                MassMat_x_VuNirb[i] = ui;
            }
        }
    }

    if ( Offline )
    {
        ConstructNIRB( XhFine,MassMat_x_VuNirb,VNirbBasis );

        if ( SamplingCoarse )
        {
            std::cout << "Computation of Coarse snapshot " <<std::endl;
            boost::timer ti;
            ComputeSnapshot( XhCoarse,"_Coarse_" );

            double Time_snapshot_Coarse = ti.elapsed();
            std::cout << "Computation of the " << NbSnapshot << " snapshots : done  -- ";
            std::cout << "Time per snapshot: " << Time_snapshot_Coarse/NbSnapshot << " sec " <<std::endl;
            SamplingCoarse = 0;
        }
    }

    //STEP TWO : Approximation of the solution using the "nirb" functions for a choosen mu




    double p = mu;


    std :: cout << "OFFLINE Procedure has been skipped" <<std::endl;


    auto uCoarseInterpolate = XhFine->element();
    auto uNirbCoarse = XhFine->element();   // Fine/Coarse Grid NIRB solution
    Eigen::VectorXd BetaiH( sizeRB );
    double TimeCoarse,TimeCoarsePostProcess,TimeFine;
    boost::timer ti;
    ti.restart();
    uCoarseInterpolate = BuildCoarseInterpolation( XhFine,XhCoarse,p );
    TimeCoarse =  ti.elapsed();
    std :: cout << "Calculation of uCoarse and it's interpolation on XhFine :"   << TimeCoarse << " sec" <<std::endl;
    uNirbCoarse = BuildNirbSolutionWithoutPostProcess( XhFine,uCoarseInterpolate,MassMat_x_VuNirb,VNirbBasis,BetaiH );
    //uNirbCoarse = BuildNirbSolution(XhFine,XhCoarse,MassMat_x_VuNirb,VNirbBasis,p,BetaiH);
    TimeCoarse =  ti.elapsed();
    std :: cout << "Construction of NIRB solution (uNirbCoarse) - Fine/Coarse Grid (saved in nirb2GridCoarse):" <<std::endl;
    std::cout << "Time to build solution " << TimeCoarse << " sec" <<std::endl;

    export_ptrtype exporter2GridCoarse( export_type::New( this->vm(), "nirb2GridCoarse" ) );
    exporter2GridCoarse->step( 0 )->setMesh( meshFine );
    exporter2GridCoarse->step( 0 )->add( "uNirbCoarse",uNirbCoarse );
    exporter2GridCoarse->save();

    auto uNirbCoarsePostProcess = XhFine->element();   // Fine/Coarse Grid NIRB solution
    ti.restart();
    if (SamplingCoarse)
    {
        std::cout << "Computation of Coarse snapshot " <<std::endl;
        boost::timer ti;
        ComputeSnapshot( XhCoarse,"_Coarse_" );

        double Time_snapshot_Coarse = ti.elapsed();
        std::cout << "Computation of the " << NbSnapshot << " snapshots : done  -- ";
        std::cout << "Time per snapshot: " << Time_snapshot_Coarse/NbSnapshot << " sec " <<std::endl;
    }
    uNirbCoarsePostProcess = BuildNirbSolutionWithPostProcess( XhFine, XhCoarse, MassMat_x_VuNirb,
                             VNirbBasis,BetaiH,uCoarseInterpolate );

    TimeCoarsePostProcess =  ti.elapsed();
    std :: cout << "Post-processing of NIRB solution (uNirbCoarsePostProcess) - Fine/Coarse Grid (saved in nirb2GridCoarse):" <<std::endl;
    std::cout << "Time  " << TimeCoarsePostProcess << " sec" <<std::endl;

    export_ptrtype exporter2GridCoarsePostProcess( export_type::New( this->vm(), "nirb2GridCoarsePostProcess" ) );
    exporter2GridCoarsePostProcess->step( 0 )->setMesh( meshFine );
    exporter2GridCoarsePostProcess->step( 0 )->add( "uNirbCoarsePostProcess",uNirbCoarsePostProcess );
    exporter2GridCoarsePostProcess->save();

    if ( ComputeError )
    {
        std :: cout << "Error calculation " <<std::endl;
        Eigen::VectorXd Betaih( sizeRB );

        auto u1Grid = XhFine->element();
        ti.restart();
        u1Grid = BuildNirbSolution( XhFine,XhFine,MassMat_x_VuNirb,VNirbBasis,p,Betaih );
        TimeFine =  ti.elapsed();
        std::cout << "Construction of uNirbFine - Fine/Fine Grid  (saved in nirb1Grid): done "<<std::endl;
        std::cout << "Time to build  " << TimeFine << " sec" <<std::endl;



        mesh_ptrtype meshRef;

        if ( ReadingMeshes )
        {
            meshRef  =  loadGMSHMesh( _mesh=new mesh_type,
                                      _filename="/home/chakir/Mesh/Mesh_H3.msh",
                                     _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK  );
        }

        else
        {
            if ( X[1] > X[0] )
            {
                meshRef = createGMSHMesh( _mesh= new mesh_type,
                                          _desc = createGeo( X[0]/2,"Mesh_Fine" ),
                                          _refine=1 );
            }

            else
            {
                meshRef = createGMSHMesh( _mesh= new mesh_type,
                                          _desc = createGeo( X[1]/4,"Mesh_Fine" ),
                                          _refine=1 );
            }
        }



        space_ptrtype XhRef = space_type::New( meshRef );



        //Reference F.E space used to compute the "reference" solution
        //which supposed to be accurate enough


        auto uRef = XhRef->element();

        uRef = blackbox( XhRef, p );
        std :: cout << "Computation of FE solution (uRef) - Ref Grid (not saved): done" <<std::endl;


        //Computation of the H1 norm of uRef
        double L2NormUref = std :: sqrt( integrate( _range=elements( XhRef->mesh() ),
                                         _expr=( idv( uRef )*idv( uRef ) ) ).evaluate()( 0,0 ) );
        double SemiH1NormUref = std :: sqrt( integrate( _range=elements( XhRef->mesh() ),
                                             _expr=( gradv( uRef )*trans( gradv( uRef ) ) ) ).evaluate()( 0,0 ) );
        double H1NormUref  = std :: sqrt( SemiH1NormUref *SemiH1NormUref + L2NormUref*L2NormUref );

        //std::cout << "L2NormURef = " << L2NormUref <<  " -- H1NormUref = " << H1NormUref <<std::endl;
        //std::cout <<std::endl;



        // //Computation of relative error measured in  H1 norm


        // // auto DFine = M_backend->newMatrix( _test=XhRef, _trial=XhFine  );//Sparse F.E mass matrix D
        // // form2( _test=XhRef, _trial=XhFine, _matrix=DFine ) =
        // //     integrate( _range=elements(XhFine->mesh()), _expr=idt(uRef)*id(uFine));
        // // auto DCoarse = M_backend->newMatrix( _test=XhRef, _trial=XhCoarse  );//Sparse F.E mass matrix D
        // // form2( _test=XhRef, _trial=XhCoarse, _matrix=DCoarse ) =
        // //     integrate( _range=elements(XhFine->mesh()), _expr=idt(uRef)*id(uCoarse));




        double  ErrH1uNirbCoarse,ErrH1uNirbCoarsePostProcess,ErrH1u1Grid;



        // uRef - uNirbCoarse
        ErrH1uNirbCoarsePostProcess = integrate( _range=elements( XhRef->mesh() ),
                                      _expr=( ( gradv( uRef )-gradv( uNirbCoarsePostProcess ) )*trans( gradv( uRef )-gradv( uNirbCoarsePostProcess ) )
                                              + ( idv( uRef )-idv( uNirbCoarsePostProcess ) )*( idv( uRef )-idv( uNirbCoarsePostProcess ) ) ) ).evaluate()( 0,0 );
        ErrH1uNirbCoarsePostProcess = sqrt( ErrH1uNirbCoarsePostProcess )/H1NormUref;


        // uRef - uNirbCoarse
        ErrH1uNirbCoarse = integrate( _range=elements( XhRef->mesh() ),
                                      _expr=( ( gradv( uRef )-gradv( uNirbCoarse ) )*trans( gradv( uRef )-gradv( uNirbCoarse ) )
                                              + ( idv( uRef )-idv( uNirbCoarse ) )*( idv( uRef )-idv( uNirbCoarse ) ) ) ).evaluate()( 0,0 );
        ErrH1uNirbCoarse = sqrt( ErrH1uNirbCoarse )/H1NormUref;

        //uRef - uNirbFine (= uRef - u1Grid)
        ErrH1u1Grid = integrate( _range=elements( XhRef->mesh() ),
                                 _expr=( ( gradv( uRef )-gradv( u1Grid ) )*trans( gradv( uRef )-gradv( u1Grid ) )
                                         + ( idv( uRef )-idv( u1Grid ) )*( idv( uRef )-idv( u1Grid ) ) ) ).evaluate()( 0,0 );
        ErrH1u1Grid = sqrt( ErrH1u1Grid )/H1NormUref;

        std :: cout << "H1-norm NIRB Error " <<std::endl;
        std :: cout << "||u_ref - u_NirbCoarse ||_{H1}  = " << ErrH1uNirbCoarse <<std::endl;
        std :: cout << "||u_ref - u_NirbCoarsePostProcess ||_{H1}  = " << ErrH1uNirbCoarsePostProcess <<std::endl;
        std :: cout << "||u_ref - u_NirbFine ||_{H1}  = " << ErrH1u1Grid<<std::endl;
        std :: cout << std::endl;





        export_ptrtype exporter1Grid( export_type::New( this->vm(), "nirbFine" ) );
        exporter1Grid->step( 0 )->setMesh( meshFine );
        exporter1Grid->step( 0 )->add( "uNirbFine",u1Grid );
        exporter1Grid->save();



        auto uFine = XhFine->element();
        auto uCoarse = XhCoarse->element();

        uCoarse = blackbox( XhCoarse, p );
        std :: cout << "Computation of FE solution (uCoarse) -  Coarse Grid (not saved):"<<std::endl;

        uFine = blackbox( XhFine, p );
        std :: cout << "Computation of FE solution (uFine) - Fine Grid (saved in EFFine): done" <<std::endl;

        /*

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


        std :: cout << "H1-Norm F.E Error " <<std::endl;
        std :: cout << "||u_ref - u_coarse ||_{H1}  = " << ErrH1uCoarse <<std::endl;
        std :: cout << "||u_ref - u_fine ||_{H1}  = " << ErrH1uFine <<std::endl;
        std :: cout <<std::endl;

        */
        export_ptrtype exporterFine( export_type::New( this->vm(), "EF_Fine" ) );
        export_ptrtype exporterCoarse( export_type::New( this->vm(), "EF_Coarse" ) );
        exporterFine->step( 0 )->setMesh( meshFine );
        exporterCoarse->step( 0 )->setMesh( meshCoarse );

        exporterFine->step( 0 )->add( "uFine", uFine );
        exporterCoarse->step( 0 )->add( "uCoarse", uCoarse );

        exporterFine->save();
        exporterCoarse->save();

    }







}


//-----------------------------------------
//-----------------------------------------

template< int PolynomialOrder>
gmsh_ptrtype NIRBTEST<PolynomialOrder>::createGeo( double hsize, std::string MeshFileName )
{
    double hCornersize1 = hsize/15.;
    double hCornersize2 = hsize/30.;
    std::ostringstream ostr;
    ostr << "Mesh.MshFileVersion = 2.2;" << "\n";
    ostr << "Mesh.CharacteristicLengthExtendFromBoundary=1;" <<std::endl;
    ostr << "Mesh.CharacteristicLengthFromPoints=1;" <<std::endl;
    ostr << "Mesh.ElementOrder=1;" <<std::endl;
    ostr << "Mesh.SecondOrderIncomplete = 0;" <<std::endl;
    ostr << "Mesh.Algorithm = 6;" <<std::endl;
    ostr << "Mesh.OptimizeNetgen=1;" <<std::endl;
    ostr << "// partitioning data" <<std::endl;
    ostr << "Mesh.Partitioner=1;" <<std::endl;
    ostr << "Mesh.NbPartitions=1;" <<std::endl;
    ostr << "Mesh.MshFilePartitioned=0;" <<std::endl;
    ostr << "Point(1) = {0,0,0.0," << hsize << "};" <<std::endl;
    ostr << "Point(2) = {1,0,0.0," << hCornersize2 << "};" <<std::endl;
    ostr << "Point(3) = {1,1,0.0," << hCornersize1 << "};"<<std::endl;
    ostr << "Point(4) = {0,1,0.0," << hCornersize2 << "};"<<std::endl;
    ostr << "Line(1) = {4,1};" <<std::endl;
    ostr << "Line(2) = {1,2};" <<std::endl;
    ostr << "Line(3) = {2,3};" <<std::endl;
    ostr << "Line(4) = {3,4};" <<std::endl;
    ostr << "Line Loop(5) = {1,2,3,4};" <<std::endl;
    ostr << "Plane Surface(6) = {5};" <<std::endl;
    ostr << "Physical Line(\"Dirichlet1\") = {1,2};" <<std::endl;
    ostr << "Physical Line(\"Dirichlet2\") = {3};" <<std::endl;
    ostr << "Physical Line(\"Dirichlet3\") = {4};" <<std::endl;
    ostr << "Physical Surface(\"Mat1\") = {6};" <<std::endl;
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
void NIRBTEST<PolynomialOrder>::ComputeSnapshot( space_ptrtype Xh,std::string filename )
{

    int XhNdof = Xh->nLocalDof();
    double pi = 3.14159265;

    auto ui = Xh->element();

    for ( int i =0; i<NbSnapshot; i++ )
    {
        double theta = ( i+1 )*( muMax-muMin )/( 113. );
        ui = blackbox( Xh,theta );
        if (ui.l2Norm()==0.)
        {
            std::cerr << "In ComputeSnapshot : ERROR IN USING BLACKBOX -> ui = 0" << endl;
            exit(0);
        }
        std::string path = "./Sol" + filename + ( boost::format( "%1%" ) %i ).str() ;
        ui.save( _path=path );
        //boost::filesystem::path full_path_Sol( boost::filesystem::current_path() );
        //cout << "Path where the sol are save " << full_path_Sol << endl;

    }

}

//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructStiffMatrixSnapshot( space_ptrtype Xh,sparse_matrix_ptrtype const & StiffMatrix )
{


    auto ui = Xh->element();
    auto uj = Xh->element();
    /*
     auto StiffMatrix = M_backend->newMatrix( _test=Xh, _trial=Xh  );//Sparse F.E stiffness matrix D
     form2( _test=Xh, _trial=Xh, _matrix=D ) =
     integrate( _range=elements(Xh->mesh()), _expr=gradt(ui)*trans(grad(uj)));
    */
    Eigen::MatrixXd S ( NbSnapshot,NbSnapshot ); //Dense RB stiffness matrix S
    cout << " Dans ConstructStiffMatrixSnapshot " << endl;
    //boost::filesystem::path full_path_Sol( boost::filesystem::current_path() );
    //cout << "Path where the sol are load " << full_path_Sol << endl;

    for ( int i = 0; i< NbSnapshot; i++ )
    {
        ui.zero();
        std::string path = ( boost::format( "./Sol_%1%" ) %i ).str() ;
        //std::string suffixe = "-";
        //ui.load( _path=path,_suffix=suffixe );

        ui.load( _path=path);

        if ( ui.l2Norm()==0. )
        {
            std::cerr << "In 'ConstructStiffMatrixSnapshot': ERROR IN LOADING FILE " << path <<std::endl;
            auto t = Xh->element();
            t.load(_path=path);
            cout << "Reloading of ui done -> L2norm(ui)= " << t.l2Norm() << endl;
            exit( 0 );
        }

        for ( int j=0; j< i; j++ )
        {
            uj.zero();
            path = ( boost::format( "./Sol_%1%" ) %j ).str() ;
            uj.load( _path=path );

            if ( uj.l2Norm()==0. )
            {
                std::cerr << "In 'ConstructStiffMatrixSnapshot': ERROR IN LOADING FILE " << path <<std::endl;
                exit( 0 );
            }

            double Sij = StiffMatrix->energy( ui,uj );
            S( i,j ) = Sij;
            S( j,i ) = Sij;
        }

        double Sii = StiffMatrix->energy( ui,ui );
        S( i,i ) = Sii;
    }

    std::string filename = ( boost::format( "StiffMatrixP%1%" )%PolynomialOrder ).str();
    std::ofstream fileMat( filename.c_str() );
    fileMat << NbSnapshot <<std::endl;

    for ( int i = 0; i<NbSnapshot; i++ )
    {
        for ( int j =0; j<NbSnapshot; j++ )
        {
            fileMat << S( i,j ) <<std::endl;
        }
    }

    fileMat.close();

    return S;
}


//--------------------------------------------------------------
//--------------------------------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructStiffMatrixRB( space_ptrtype Xh,vector_of_element_type &Vu,sparse_matrix_ptrtype const & StiffMatrix )
{

    Eigen::MatrixXd S ( sizeRB,sizeRB ); //Dense RB stiffness matrix

    for ( int i = 0; i< sizeRB; i++ )
    {
        for ( int j=0; j< i; j++ )
        {
            double Sij = StiffMatrix->energy( Vu[i],Vu[j] );
            S( i,j ) = Sij;
            S( j,i ) = Sij;
        }

        double Sii = StiffMatrix->energy( Vu[i],Vu[i] );
        S( i,i ) = Sii;
    }

    return S;
}
//--------------------------------------------------------------
//--------------------------------------------------------------
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: ConstructMassMatrixRB( space_ptrtype Xh,vector_of_element_type &Vu,sparse_matrix_ptrtype const &MassMatrix )
{

    Eigen::MatrixXd S ( sizeRB,sizeRB ); //Dense RB mass matrix

    for ( int i = 0; i< sizeRB; i++ )
    {
        for ( int j=0; j< i; j++ )
        {
            double Sij = MassMatrix->energy( Vu[i],Vu[j] );
            S( i,j ) = Sij;
            S( j,i ) = Sij;
        }

        double Sii = MassMatrix->energy( Vu[i],Vu[i] );
        S( i,i ) = Sii;
    }

    return S;

}

//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> ::ChooseRBFunction( space_ptrtype Xh,vector_of_element_type &VuBasis,sparse_matrix_ptrtype const &MassMatrix,sparse_matrix_ptrtype const &StiffMatrix )
{

    //std :: cout << "Inside ChooseRBFunction" <<std::endl;
    Eigen::MatrixXd S ( NbSnapshot,NbSnapshot ); //Dense Stiffness Matrix

    if ( !Sampling )
    {
        std::string str = ( boost::format( "StiffMatrixP%1%" )%PolynomialOrder ).str();
        std::ifstream fileMat( str.c_str() );

        if ( !fileMat )
        {
            std::cerr << "Error in reading file " << ( boost::format( "StiffMatrixP%1%" )%PolynomialOrder ).str() <<  "does not exist " << std::endl;
            std::cout << "Need to re-build the Stiffness Matrix -> SAMPLING SET TO 1 " << std::endl;
            Sampling = 1;
        }

        else
        {
            int Itemp;
            fileMat >> Itemp;

            if ( Itemp != NbSnapshot )
            {
                std :: cerr << "Error in reading file" << ( boost::format( "StiffMatrixP%1%" )%PolynomialOrder ).str() <<  " --  " << Itemp << " != NbSnapshot(" << NbSnapshot << ")" << std::endl;
                std::cout << "Need to re-build the Stiffness Matrix -> SAMPLING SET TO 1 " << std::endl;
                Sampling = 1;
            }

            else
            {
                for ( int i=0; i<NbSnapshot; i++ )
                {
                    for ( int j =0; j<NbSnapshot; j++ )
                    {
                        double Dtemp;
                        fileMat >> Dtemp;
                        S( i,j ) = Dtemp;
                    }
                }
            }
        }

    }

    if ( Sampling )
    {
        S = ConstructStiffMatrixSnapshot( Xh,StiffMatrix );
    }


    Eigen::EigenSolver <Eigen::MatrixXd> eigen_solver( S );

    //int nb_EigenValue = eigen_solver.eigenvalues().size();
    //eigen_solver.eigenvectors().col(i)  =  eigenvector #i

    std :: cout << "Computation of the eigenvalues of the stiffness matrix S (NbSnapshot,NbSnapshot) : done" <<std::endl;

    //Sorting N = "SizeBR" F.E solutions uh(mu_k) the more represented in the eigenvectors
    //associated to the N largest eigenvalues




    // saving the index's number to identify the NIRB basis functions
    std::string path = ( boost::format( "./IndBR%1%" ) %sizeRB ).str() ;
    std::ofstream find( path.c_str() );

    if ( !find )
    {
        std :: cerr <<" 'ChooseRBFunction routine' - ERROR IN OPENING FILE:" << path <<std::endl;
        exit( 0 );
    }

    std::vector<int> Uind( NbSnapshot );

    for ( int i=0; i<NbSnapshot; i++ )
    {
        Uind[i] = 0; //Set to zero if u(Uind[i]) is not NIRB basis function.
    }


    auto ui = Xh->element();


    double Tol = 1e-30;
    //std::ofstream fout("fout");
    //Sorting eigenvector #i
    std::vector<double>Vi( NbSnapshot );

    for ( int i=0; i<sizeRB; i++ )
    {
        for ( int ri = 0; ri < NbSnapshot; ri++ )
        {
            //If u(Uind[ri]) is a NIRB basis function, the "Uind[ri]" component of Vi is set to zero
            //to avoid it.
            if ( Uind[ri] == 0 )
            {
                Vi[ri] = std::abs( real( eigen_solver.eigenvectors().col( i )[ri] ) );
            }

            else
            {
                Vi[ri] = 0;
            }
        }

        //Using max_element to find the position of the maximum component of Vi

        int IndMax = std::distance( Vi.begin(),boost::max_element( Vi ) );

        //copy u[IndMax] in Mu
        ui.zero();
        path = ( boost::format( "./Sol_%1%" ) %IndMax ).str();
        ui.load( _path=path );

        //cout << "LOAD DONE " <<std::endl;
        if ( ui.l2Norm()==0. )
        {
            std::cerr <<"'ChooseRBFunction routine': ERROR IN LOADING FILE  " << path <<std::endl;
            exit( 0 );
        }

        VuBasis[i]= ui;

        if ( i>0 )
        {
            double normL2 = OrthogonalisationRBFunctionL2GrammSchmidt( Xh,VuBasis,i,MassMatrix );

            //std :: cout << "i = " << i << " --- NormU(" << IndMax << ") = " << normL2 <<std::endl;
            //fout << "i = " << i << " --- NormU(" << IndMax << ") = " << normL2 <<std::endl;
            if ( normL2 > Tol )
            {
                Uind[IndMax] = i;
                find << IndMax <<std::endl;
            }

            else
            {
                Uind[IndMax] = -1;
                i--;
            }
        }

        else
        {
            Uind[IndMax] = i;
            find << IndMax <<std::endl;
        }

    }

    find.close();
    //fout.close();


    double L2normVui;

    //save in a file the matrix A
    for ( int i=0; i<sizeRB; i++ )
    {
        L2normVui = std::sqrt( MassMatrix->energy( VuBasis[i],VuBasis[i] ) );
        VuBasis[i].scale( 1./L2normVui );
        VuBasis[i].save( _path=( boost::format( "./Sol_OR%1%" ) %i ).str() );
    }


}


//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
double NIRBTEST <PolynomialOrder>  :: OrthogonalisationRBFunctionL2GrammSchmidt( space_ptrtype Xh, vector_of_element_type &Vu,int n,sparse_matrix_ptrtype const &MassMatrix )
{

    if ( MassMatrix->linftyNorm()==0. )
    {
        std::cout << "ERROR in OrthogonalisationRBFunctionL2GrammSchmidt : MassMatrix is null " <<std::endl;
        exit( 0. );
    }

    double Dtemp1,Dtemp2,Dtemp3,Dtemp4;
    //We suppose that the n-1 first function has been already orthogonalized

    for ( int k=0; k<n; k++ )
    {

        Dtemp1 = MassMatrix->energy( Vu[n],Vu[k] );
        Dtemp2 = MassMatrix->energy( Vu[k],Vu[k] );
        Dtemp3 = -Dtemp1/Dtemp2;
        Vu[n].add( Dtemp3,Vu[k] );

    }

    Dtemp4 = MassMatrix->energy( Vu[n],Vu[n] );
    return Dtemp4;
}

//-----------------------------------------
//-----------------------------------------

//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> :: OrthogonalisationRBFunctionL2GrammSchmidt
( space_ptrtype Xh,vector_of_element_type &Vu, sparse_matrix_ptrtype const &MassMatrix  )
{

    double Dtemp,Dtemp1,Dtemp2,Dtemp3, L2norm;

    for ( int i=1; i<sizeRB; i++ )
    {
        for ( int k=0; k<i; k++ )
        {
            Dtemp1 = MassMatrix->energy( Vu[i],Vu[k] );
            Dtemp2 = MassMatrix->energy( Vu[k],Vu[k] );
            Dtemp3 = - Dtemp1/Dtemp2;
            Vu[i].add( Dtemp3,Vu[k] );
        }

    }

    for ( int i=0; i<sizeRB; i++ )
    {
        L2norm = std::sqrt( MassMatrix->energy( Vu[i],Vu[i] ) );
        Vu[i].scale( 1./L2norm ) ;
    }

}

//-----------------------------------------
//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> ::OrthogonalisationRBFunction( space_ptrtype Xh,vector_of_element_type &Vu,vector_of_element_type &M_Vu,sparse_matrix_ptrtype const &StiffMatrix,sparse_matrix_ptrtype const &MassMatrix )
{
    //First step : L2 -Pre-orthogonalization and normalization using a Gramm-schmidt method

    std::string path;
    auto ui = Xh->element();
    ui.zero();

    //Second step : Solving a generalized eigenvalue problem
    //A*Epsilon_i = \lambda_i*B*Epsilon_i
    //A : Dense stiffness matrix in the pre-orthogonalized nirb basis
    //B : Dense mass matrix in the pre-orthogonalized nirb basis  = Idendity matrix

    //Construction of the matrix A and B

    //     M.resize(sizeRB,sizeRB);
    Eigen :: MatrixXd A ( sizeRB,sizeRB );
    Eigen :: MatrixXd B ( sizeRB,sizeRB );
    Eigen :: MatrixXd C ( sizeRB,sizeRB );

    A = ConstructStiffMatrixRB( Xh,Vu,StiffMatrix );
    B = ConstructMassMatrixRB( Xh,Vu,MassMatrix );
    //     Eigen:: GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
    C = B.inverse()*A;

    Eigen:: EigenSolver<Eigen::MatrixXd> es;

    es.compute( C );
    vector_of_element_type VuNirb( sizeRB,ui );


    for ( int i =0; i<sizeRB; i++ )
    {
        for ( int k=0; k<sizeRB; k++ )
        {
            VuNirb[i].add( real( es.eigenvectors().col( i )[k] ),Vu[k] );
        }
    }


    OrthogonalisationRBFunctionL2GrammSchmidt( Xh,VuNirb,MassMatrix );

    auto Mui = M_backend->newVector( Xh );
    auto Ui = M_backend->newVector( Xh );

    //save in a file the final NIRB basis functions
    for ( int i=0; i<sizeRB; i++ )
    {
        path = ( boost::format( "./RB%1%NIRB_BasisFile_%2%" ) %sizeRB%i ).str() ;
        VuNirb[i].save( _path=path );
        Vu[i] = VuNirb[i];
        *Ui = VuNirb[i];
        Ui->close();
        MassMatrix->multVector( Ui,Mui );
        M_Vu[i] = *Mui;

        if ( M_Vu[i].l2Norm()>1e-20 )
        {
            path = ( boost::format( "./RB%1%D_Nirb_Basis%2%File" )%sizeRB%i ).str();
            M_Vu[i].save( _path=path );
        }

        else
        {
            std::cerr << "ERROR: in computation of D*NIRB Basis function (" << i << ") : Result equal to zero " <<std::endl;
            exit( 0 );
        }

    }
}

//---------------------------------------------------
//---------------------------------------------------

template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder> :: BuildNirbSolutionWithoutPostProcess ( space_ptrtype XhFine, element_type uCoarseInterpolate,vector_of_element_type const &M_VNirbBasis, vector_of_element_type const &VNirbBasis, Eigen::VectorXd & BetaiH )
{



    auto uNirb = XhFine->element();

    auto UcoarseInterpol = M_backend->newVector( XhFine );
    *UcoarseInterpol =  uCoarseInterpolate;

    //Computation of the coefficiant \BetaiH = \int uCoarse*\Epsilon_i
    //with \Epsilon_i being the final Nirb basis functions

    uNirb.zero();
    auto Mui = M_backend->newVector( XhFine );

    for ( int i =0; i<sizeRB; i++ )
    {
        *Mui=M_VNirbBasis[i];
        BetaiH[i] = inner_product( UcoarseInterpol,Mui );
        uNirb.add( BetaiH[i], VNirbBasis[i] );
    }


    return uNirb;


}
//---------------------------------------------------
//---------------------------------------------------

template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder> ::BuildCoarseInterpolation( space_ptrtype XhFine,space_ptrtype XhCoarse, double param )
{
    //boost::timer ti;
    //ti.restart();
    auto uCoarse = XhCoarse->element();
    auto uCoarseInterpolate = XhFine->element();
    auto ui = XhFine->element();
    uCoarse = blackbox( XhCoarse,param );
    //std :: cout << "Calculation of uCoarse "   << ti.elapsed() << " sec" <<std::endl;
    //ti.restart();
    ////Interpolation of the coarse solution on the fine Mesh to compute the coefficiant BetaiH

    //auto opI = opInterpolation(_domainSpace=XhCoarse,_imageSpace=XhFine);
    //std :: cout << "Construction of the operator"   << ti.elapsed() << " sec" <<std::endl;
    //ti.restart();
    //opI->apply(uCoarse,uCoarseInterpolate);
    //std :: cout << "Application of the operator"   << ti.elapsed() << " sec" <<std::endl;
    interpolate ( XhFine,uCoarse,uCoarseInterpolate );
    return uCoarseInterpolate;

}


// //---------------------------------------------------
////---------------------------------------------------

template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder> :: BuildNirbSolution( space_ptrtype XhFine,space_ptrtype XhCoarse, vector_of_element_type const &M_VNirbBasis, vector_of_element_type const &VNirbBasis, double param,Eigen::VectorXd BetaiH )
{

    //Eigen::VectorXd BetaiH(sizeRB);
    auto uCoarseInterpolate = XhFine->element();
    auto uNirb = XhFine->element();
    uCoarseInterpolate = BuildCoarseInterpolation( XhFine,XhCoarse,param );
    uNirb = BuildNirbSolutionWithoutPostProcess( XhFine,uCoarseInterpolate,M_VNirbBasis,VNirbBasis,BetaiH );
    return uNirb;

}

//-----------------------------------------
template< int PolynomialOrder>
void NIRBTEST<PolynomialOrder> :: ConstructNIRB( space_ptrtype Xh, vector_of_element_type  &M_VNirbBasis, vector_of_element_type  &VNirbBasis )
{

    std :: cout << "OFFLINE PROCEDURE :  Construction of the 'non intruisive' reduced basis (nirb) " <<std::endl;
    boost::timer ti;
    double Time_snapshot = 0.;

    if ( Sampling )
    {
        std::cout << "Sampling Procedure :" <<std::endl;
        ComputeSnapshot( Xh,"_" );
        Time_snapshot = ti.elapsed();
        std::cout << "Computation of the " << NbSnapshot << " snapshots : done  -- ";
        std::cout << "Time per snapshot: " << Time_snapshot/NbSnapshot << " sec " <<std::endl;
    }

    ti.restart();

    //construction of elementary mass and stiffness matrix
    auto ui = Xh->element();
    auto uj = Xh->element();

    sparse_matrix_ptrtype MassMatrix = M_backend->newMatrix( _test=Xh, _trial=Xh ); //Sparse FE mass matrix
    form2( _test=Xh , _trial=Xh, _matrix= MassMatrix ) =
        integrate( _range=elements( Xh->mesh() ), _expr=id( uj )*idt( ui ) );

    sparse_matrix_ptrtype StiffMatrix = M_backend->newMatrix( _test=Xh, _trial=Xh ); //Sparse FE mass matrix
    form2( _test=Xh , _trial=Xh, _matrix= StiffMatrix ) =
        integrate( _range=elements( Xh->mesh() ), _expr=gradt( ui )*trans( grad( uj ) ) );
    double Time_Construct_Elementary_Matrix = ti.elapsed();
    std::cout << "Time to build elementary F.E matrix = " << Time_Construct_Elementary_Matrix << " sec " <<std::endl;

    ti.restart();
    ChooseRBFunction( Xh,VNirbBasis,MassMatrix,StiffMatrix );
    double TimeChooseRB =  ti.elapsed();
    std::cout << "Choice of " << sizeRB << " reduced basis functions : done  -- ";
    std::cout << "Time: " << TimeChooseRB << " sec "<<std::endl;
    ti.restart();

    //Orthogonalisation de Gram-schmidt

    OrthogonalisationRBFunction( Xh,VNirbBasis,M_VNirbBasis,StiffMatrix,MassMatrix );
    //OrthogonalisationRBFunction(Xh);
    double Time_buildNIRBbasis = ti.elapsed();
    std::cout << "H1 and L2 orthogonalisation of the reduced basis functions : done  -- ";
    std::cout << "Time: " << Time_buildNIRBbasis << " sec " <<std::endl;
    double TotalTime = Time_buildNIRBbasis +TimeChooseRB + Time_Construct_Elementary_Matrix + Time_snapshot;
    std::cout <<"Total time for the OFFLINE procedure: " << TotalTime << " sec "<<std::endl;



}



//---------------------------------------------------
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
    auto velocity = vec( cst( cos( param ) ),cst( sin( param ) ) );
    /*
     auto g = ( chi(abs(Px()-1) < 1e-10 )*Px()*Px()+
     chi(abs(Py()-1) < 1e-10 )*Py()*Py() );
     */

    auto g = ( Px()*Px()*Py()*Py() );
    auto F = M_backend->newVector( Xh );
    /*
     form1( _test=Xh, _vector=F ) =
     integrate( _range=boundaryfaces(Xh->mesh()),
     _expr=g*(-0.01*dn(v)+penaltyTerm*id(v)/hFace())  );
     */
    auto D = M_backend->newMatrix( _test=Xh, _trial=Xh  );
    form2( _test=Xh, _trial=Xh, _matrix=D ) =
        integrate( _range=elements( Xh->mesh() ), _expr=0.01*gradt( u )*trans( grad( v ) )+ ( gradt( u )*velocity )*id( v ) );

    /*
     form2( _test=Xh, _trial=Xh, _matrix=D ) +=
     integrate( boundaryfaces(Xh->mesh()),
     -0.01*dnt(u)*id(v)-0.01*dn(v)*idt(u)+penaltyTerm*id(v)*idt(u)/hFace());
     */

    form2( _test=Xh, _trial=Xh, _matrix=D ) += on( boundaryfaces( Xh->mesh() ), u, F,g );
    backend_type::build(soption("backend"))->solve( _matrix=D, _solution=u, _rhs=F );

    return u;
}


//---------------------------------------------------
//---------------------------------------------------
//Post-processing of NIRB solution
template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> :: BuildBetaH( space_ptrtype XhFine,
        space_ptrtype XhCoarse,vector_of_element_type  &M_VNirbBasis )
{

    //boost::timer ti;
    //ti.restart();
    //auto opI = opInterpolation(_domainSpace=XhCoarse,_imageSpace=XhFine);
    //std :: cout << "Construction of the operator"   << ti.elapsed() << " sec" <<std::endl;
    //ti.restart();
    auto uCoarse = XhCoarse->element();
    auto uCoarseInterpolation = XhFine->element();

    std::vector<int> IndTab( sizeRB );
    Eigen::MatrixXd BetaH( sizeRB,sizeRB );
    std::string path = ( boost::format( "./IndBR%1%" ) %sizeRB ).str() ;
    std::ifstream find( path.c_str() );

    if ( !find )
    {
        std::cerr << "ERROR : Problem in opening 'IndBR'" << sizeRB << " file " <<std::endl;
        exit( 0 );


    }

    for ( int i =0 ; i < sizeRB; i++ )
    {
        int Itemp;
        find >> Itemp;
        IndTab[i] = Itemp;
    }

    path =  ( boost::format( "./AlphaH_%1%" ) %sizeRB ).str();
    std::ofstream fB( path.c_str() );
    fB << sizeRB << " " << sizeRB <<std::endl;

    for ( int i =0; i<sizeRB; i++ )
    {
        uCoarse.zero();
        uCoarseInterpolation.zero();
        path = ( boost::format( "./Sol_Coarse_%1%" ) %IndTab[i] ).str();
        uCoarse.load( _path=path );

        if ( uCoarse.l2Norm()==0. )
        {
            std::cerr <<"'BuildBetaH routine': ERROR IN LOADING FILE  " << path <<std::endl;
            std::cerr <<"Coarse sampling has to be done " << std::endl;
            ComputeSnapshot( XhCoarse,"_Coarse_" );
            i--;

        }

        //opI->apply(uCoarse,uCoarseInterpolation);
        interpolate( XhFine,uCoarse,uCoarseInterpolation );
        auto UcoarseInterpol_i = M_backend->newVector( XhFine );
        *UcoarseInterpol_i = uCoarseInterpolation;
        auto Muj = M_backend->newVector( XhFine );

        for ( int j =0; j<sizeRB; j++ )
        {
            *Muj=M_VNirbBasis[j];
            BetaH( j,i ) = inner_product( UcoarseInterpol_i,Muj ); //beta_j(\mu_i)
            fB << BetaH( j,i ) <<std::endl;
        }
    }

    fB.close();
    return BetaH;

}
//---------------------------------------------------
//---------------------------------------------------

template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> ::  BuildBetah( space_ptrtype XhFine,
        vector_of_element_type  &M_VNirbBasis )
{

    auto uFine = XhFine->element();

    std::vector<int> IndTab( sizeRB );
    Eigen::MatrixXd BetaH( sizeRB,sizeRB );
    std::string path = ( boost::format( "./IndBR%1%" ) %sizeRB ).str() ;
    std::ifstream find( path.c_str() );

    if ( !find )
    {
        std::cerr << "ERROR : Problem in opening 'IndBR'" << sizeRB << " file " <<std::endl;
        exit( 0 );
    }

    for ( int i =0 ; i < sizeRB; i++ )
    {
        int Itemp;
        find >> Itemp;
        IndTab[i] = Itemp;
    }

    path =  ( boost::format( "./Betah_%1%" ) %sizeRB ).str();
    std::ofstream fB( path.c_str() );
    fB << sizeRB << " " << sizeRB <<std::endl;

    for ( int i =0; i<sizeRB; i++ )
    {
        uFine.zero();
        path = ( boost::format( "./Sol_%1%" ) %IndTab[i] ).str();
        uFine.load( _path=path );

        if ( uFine.l2Norm()==0. )
        {
            std::cerr <<"'ChooseRBFunction routine': ERROR IN LOADING FILE  " << path <<std::endl;
            exit( 0 );
        }

        auto Ufine_i = M_backend->newVector( XhFine );
        *Ufine_i = uFine;
        auto Muj = M_backend->newVector( XhFine );

        for ( int j =0; j<sizeRB; j++ )
        {
            *Muj=M_VNirbBasis[j];
            BetaH( j,i ) = inner_product( Ufine_i,Muj ); // beta_j(\mu_i)
            fB << BetaH( j,i ) <<std::endl;
        }
    }

    fB.close();

    return BetaH;

}
//---------------------------------------------------
//---------------------------------------------------

template< int PolynomialOrder>
Eigen::MatrixXd NIRBTEST<PolynomialOrder> ::BuildPostProcessMatrix( space_ptrtype XhFine, space_ptrtype XhCoarse,vector_of_element_type  &M_VNirbBasis )
{

    Eigen::MatrixXd BetaH( sizeRB,sizeRB );
    Eigen::MatrixXd Betah( sizeRB,sizeRB );
    Eigen::MatrixXd T( sizeRB,sizeRB );

    int FileOkA = 1,FileOkB = 1 ;

    if ( !Offline )
    {
        std::string path =  ( boost::format( "./AlphaH_%1%" ) %sizeRB ).str();
        std::ifstream fA ( path.c_str() );

        if ( !fA )
        {
            FileOkA = 0;
            BetaH = BuildBetaH( XhFine,XhCoarse,M_VNirbBasis );
        }

        path =  ( boost::format( "./betah_%1%" ) %sizeRB ).str();
        std::ifstream fB ( path.c_str() );

        if ( !fB )
        {
            FileOkB = 0;
            Betah = BuildBetah( XhFine,M_VNirbBasis );
        }

        int Itemp1,Itemp2;
        double dTempA,dTempB;

        if ( FileOkA && FileOkB )
        {
            fA >> Itemp1 >> Itemp2;
            fB >> Itemp1 >> Itemp2;

            for ( int i = 0; i< sizeRB; i++ )
            {
                for ( int j=0; j<sizeRB; j++ )
                {
                    fA >> dTempA;
                    fB >> dTempB;
                    BetaH( j,i ) = dTempA;
                    Betah( j,i ) = dTempB;
                }
            }
        }

        else
        {
            if ( FileOkA )
            {
                fA >> Itemp1 >> Itemp2;

                for ( int i = 0; i< sizeRB; i++ )
                {
                    for ( int j=0; j<sizeRB; j++ )
                    {
                        fA >> dTempA;
                        BetaH( j,i ) = dTempA;
                    }
                }
            }

            if ( FileOkB )
            {
                fB >> Itemp1 >> Itemp2;

                for ( int i = 0; i< sizeRB; i++ )
                {
                    for ( int j=0; j<sizeRB; j++ )
                    {
                        fB >> dTempB;
                        Betah( j,i ) = dTempB;
                    }
                }
            }
        }

        fA.close();
        fB.close();

    }

    if ( Offline )
    {
        BetaH = BuildBetaH( XhFine,XhCoarse,M_VNirbBasis );
        Betah = BuildBetah( XhFine,M_VNirbBasis );
    }

    T = Betah* BetaH.inverse() ;

    return T;


}


// //---------------------------------------------------
////---------------------------------------------------

template< int PolynomialOrder>
typename NIRBTEST<PolynomialOrder>::element_type NIRBTEST<PolynomialOrder> ::
BuildNirbSolutionWithPostProcess( space_ptrtype XhFine,
                                  space_ptrtype XhCoarse,
                                  vector_of_element_type  &M_VNirbBasis,
                                  vector_of_element_type  &VNirbBasis,
                                  Eigen::VectorXd &BetaiuH,
                                  element_type &uCoarseInterpolate )
{

    Eigen::MatrixXd T ( sizeRB,sizeRB );
    T = BuildPostProcessMatrix( XhFine,XhCoarse,M_VNirbBasis );

    //Eigen::VectorXd AlphaiH (sizeRB);
    //AlphaiH = T*BetaiuH;
    double AlphaiH;
    auto uNirb = XhFine->element();
    uNirb.zero();

    for ( int i =0; i<sizeRB; i++ )
    {
        AlphaiH =0.;

        for ( int k=0; k<sizeRB; k++ )
        {
            AlphaiH += T( i,k )*BetaiuH[k];
        }

        uNirb.add( AlphaiH, VNirbBasis[i] );
    }


    return uNirb;
}



#endif
