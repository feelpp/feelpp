/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Stephane Veys <stephane.veys@imag.fr>
       Date: 2009-09-30

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
   \file pod.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
           Stephane Veys <stephane.veys@imag.fr>
   \date 2009-09-30
 */

#ifndef __POD_H
#define __POD_H 1


#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include "boost/tuple/tuple_io.hpp"
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/support/lambda.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>


#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>

#include <vector>

#include <feel/feelalg/solvereigen.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcrb/parameterspace.hpp>
#include <feel/feelcrb/crbdb.hpp>
#include <feel/feelcrb/crbscm.hpp>
#include <feel/feelcore/serialization.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelts/bdf.hpp>

namespace Feel
{



/**
 * \class POD
 * \brief POD class
 *
 * This class implements POD method useful to treat transient problems with the
 * certified reduced basis method.
 *
 * @author Christophe Prud'homme
 * @author Stephane Veys
 */
template<typename TruthModelType>
class POD
{
public :

    typedef TruthModelType truth_model_type;
    typedef truth_model_type model_type;
    typedef boost::shared_ptr<truth_model_type> truth_model_ptrtype;


    typedef typename model_type::value_type value_type;

    //! element of the functionspace type
    typedef typename model_type::element_type element_type;
    typedef typename model_type::element_ptrtype element_ptrtype;



    //! mesh type
    typedef typename model_type::mesh_type mesh_type;

    //! mesh shared_ptr
    typedef typename model_type::mesh_ptrtype mesh_ptrtype;

    //! function space type
    typedef typename model_type::functionspace_type functionspace_type;
    typedef typename model_type::functionspace_ptrtype functionspace_ptrtype;

    typedef typename model_type::space_type space_type;

    //! time discretization
    typedef Bdf<space_type>  bdf_type;
    typedef boost::shared_ptr<bdf_type> bdf_ptrtype;


    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    typedef std::vector<element_type> wn_type;

    typedef std::vector<element_type> mode_set_type;


    typedef Eigen::VectorXd vectorN_type;
    typedef Eigen::MatrixXd matrixN_type;

    typedef typename model_type::vector_ptrtype vector_ptrtype;

    typedef typename model_type::backend_type backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    POD()
        :
        M_store_pod_matrix ( false ),
        M_store_pod_matrix_format_octave ( false ),
        M_Nm( 1 ),
        M_pod_matrix(),
        M_model(),
        M_use_solutions( true )
    {}

    /**
     * copy constructor
     */
    POD( POD const & o )
        :
        M_store_pod_matrix( o.M_store_pod_matrix ),
        M_store_pod_matrix_format_octave( o.M_store_pod_matrix_format_octave ),
        M_Nm( o.M_Nm ),
        M_pod_matrix( o.M_matrix ),
        M_model( o.M_model ),
        M_use_solutions( o.M_use_solutions )
    {}

    //! destructor
    ~POD()
    {}



    //! return a bool to indicate if we store the pod matrix
    bool storePodMatrix() const
    {
        return M_store_pod_matrix;
    }

    //! return a bool to indicate if we store the pod matrix with octave format
    bool storePodMatrixFormatOctave() const
    {
        return M_store_pod_matrix_format_octave;
    }

    //! return the pod matrix
    const matrixN_type& podMatrix()
    {
        return M_pod_matrix;
    }

    //! return number of mode used per mu
    const int nm()
    {
        return M_Nm;
    }

    //! return model used
    const truth_model_ptrtype & model()
    {
        return M_model;
    }

    const double timeInitial()
    {
        return M_time_initial;
    }

    void setBdf( bdf_ptrtype& bdf )
    {
        M_bdf = bdf;
    }

    void setNm ( const int Nm )
    {
        M_Nm=Nm;
    }

    void setModel ( truth_model_ptrtype Model )
    {
        M_model=Model;
    }

    void setTimeInitial( double Ti )
    {
        M_time_initial = Ti;
    }

    //! fill the matrix which will be used to perform the POD
    void fillPodMatrix( const wn_type&  elements_set );

    /**
     * input/output : MpdeSet (set of modes to add in the reduced basis)
     * input : is_primal ( bool which indicates if the problem is the primal one or not )
     * optional input : elements set to perform POD
     */
    int pod( mode_set_type& ModeSet, bool is_primal, const wn_type& elements_set=std::vector<element_type>(), bool use_solutions=true );

    void exportMode( double time, element_ptrtype& mode );

    void projectionOnPodSpace();

private :

    bool M_store_pod_matrix;

    bool M_store_pod_matrix_format_octave;

    //number of modes
    int M_Nm;

    matrixN_type M_pod_matrix;

    truth_model_ptrtype M_model;

    backend_ptrtype M_backend;

    export_ptrtype exporter;

    bdf_ptrtype M_bdf;

    double M_time_initial;

    bool M_use_solutions;
};//class POD

po::options_description podOptions( std::string const& prefix = "" );

template<typename TruthModelType>
void POD<TruthModelType>::exportMode( double time, element_ptrtype& mode )
{

    LOG(INFO) << "exportResults starts\n";

    functionspace_ptrtype function_space = M_model->functionSpace();
    mesh_ptrtype mesh = function_space->mesh();

    //1D
    std::ofstream mode_file;
    mode_file.open( "mode.dat",std::ios::out );
    std::map<double,double> data;

    for ( auto it = mesh->beginElement( ),
            en = mesh->endElement( );
            it!=en; ++it )
    {
        for ( size_type i = 0; i < space_type::basis_type::nLocalDof; ++i )
        {
            value_type a = it->point( 0 ).node()[0];
            value_type b = it->point( 1 ).node()[0];
            value_type x = 0;

            if ( i == 0 )
                x=a;

            else if ( i == 1 )
                x=b;

            else
                x= a + ( i-1 )*( b-a )/( space_type::basis_type::nLocalDof-1 );

            data[x] = mode->localToGlobal( it->id(), i, 0 );
        }

    }

    BOOST_FOREACH( auto d, data )
    {
        mode_file << d.first << " " << d.second << "\n";
    }
    mode_file.close();

    //2D - 3D
    //exporter->step(time)->setMesh( mesh );
    //exporter->step(time)->add( "mode", *mode );
    //exporter->save();

}

template<typename TruthModelType>
void POD<TruthModelType>::fillPodMatrix( const wn_type& elements_set)
{
    if( M_use_solutions )
    {
        boost::mpi::timer timer;

        //M_bdf->setRestart( true );
        int K = M_bdf->timeValues().size()-1;
        M_pod_matrix.resize( K,K );

        auto bdfi = M_bdf->deepCopy();
        auto bdfj = M_bdf->deepCopy();

        bdfi->setRestart( true );
        bdfj->setRestart( true );
        bdfi->setTimeInitial( M_time_initial );
        bdfj->setTimeInitial( M_time_initial );
        bdfi->setRestartAtLastSave(false);
        bdfj->setRestartAtLastSave(false);
        for ( bdfi->restart(); !bdfi->isFinished(); bdfi->next() )
        {
            int i = bdfi->iteration()-1;
            bdfi->loadCurrent();

            //here is a barrier. For now, if this barrier is removed,
            //during the bdfj->restart() the program will wait and the memory increase.
            //during the init step in bdf, the metadata are read,
            //and it is during this step that the memory could increase dangerously.
            Environment::worldComm().barrier();

            for ( bdfj->restart(); !bdfj->isFinished() && ( bdfj->iteration() < bdfi->iteration() ); bdfj->next() )
            {
                int j = bdfj->iteration()-1;
                bdfj->loadCurrent();

                M_pod_matrix( i,j ) = M_model->scalarProductForPod( bdfj->unknown( 0 ), bdfi->unknown( 0 ) );
                M_pod_matrix( j,i ) = M_pod_matrix( i,j );
            }

            M_pod_matrix( i,i ) = M_model->scalarProductForPod( bdfi->unknown( 0 ), bdfi->unknown( 0 ) );
        }
        double time=timer.elapsed();
        if( Environment::worldComm().isMasterRank() )
        {
            std::cout<<"POD matrix filled in  "<<time<<" s"<<std::endl;
        }
    }//fill pod matrix with solutions
    else
    {
        //use elements given by CRB
        int size = elements_set.size();
        M_pod_matrix.resize( size , size );
        CHECK( size > 0 )<<" elements set doesn't conatin any element\n";
        for(int i=0; i<size; i++)
        {
            for(int j=i+1; j<size; j++)
            {
                M_pod_matrix( i,j ) = M_model->scalarProductForPod( elements_set[i], elements_set[j] );
                M_pod_matrix( j,i ) = M_pod_matrix( i,j );
            }
            M_pod_matrix( i,i ) = M_model->scalarProductForPod( elements_set[i] , elements_set[i] );
        }
    }
}//fillPodMatrix


//the bool is_primal is used to determine the max number of modes ( M_Nm ) whithout taking eigenvectors
//associated with too small eigenvalues
//Since we always use M_Nm = 1 or 2 (so we never have to modify M_Nm) we only correct M_Nm in the primal case
template<typename TruthModelType>
int POD<TruthModelType>::pod( mode_set_type& ModeSet, bool is_primal, const wn_type& elements_set, bool use_solutions )
{
    //M_backend = backend_type::build( BACKEND_PETSC, Environment::worldComm() );
    M_backend = backend();

    Eigen::SelfAdjointEigenSolver< matrixN_type > eigen_solver;
    M_use_solutions = use_solutions;

    fillPodMatrix( elements_set );
    int size = M_pod_matrix.cols();

    //store the matrix
    if ( M_store_pod_matrix )
    {
        std::ofstream matrix_file;
        LOG(INFO)<<"saving Pod matrix in a file \n";
        matrix_file.open( "PodMatrix",std::ios::out );
        matrix_file<<M_pod_matrix.rows();
        matrix_file<<"\n";
        matrix_file<<M_pod_matrix.cols();
        matrix_file<<"\n";

        for ( int i=0; i<M_pod_matrix.rows(); i++ )
        {
            for ( int j=0; j<M_pod_matrix.cols(); j++ )
            {
                matrix_file<< std::setprecision( 16 ) <<M_pod_matrix( i,j )<<" ";
            }

            matrix_file<<"\n";
        }

        LOG(INFO)<<" matrix wrote in file named PodMatrix \n";
    }

    else if ( M_store_pod_matrix_format_octave )
    {
        std::ofstream matrix_file;
        LOG(INFO)<<"saving Pod matrix in a file \n";
        matrix_file.open( "PodMatrixOctave.mat",std::ios::out );
        matrix_file<<"# name: A\n";
        matrix_file<<"# type: matrix\n";
        matrix_file<<"# rows: "<<M_pod_matrix.rows()<<"\n";
        matrix_file<<"# columns: "<<M_pod_matrix.cols()<<"\n";

        for ( int i=0; i<M_pod_matrix.rows(); i++ )
        {
            for ( int j=0; j<M_pod_matrix.cols(); j++ )
            {
                matrix_file<< std::setprecision( 16 ) <<M_pod_matrix( i,j )<<" ";
            }

            matrix_file<<"\n";
        }

        LOG(INFO)<<" matrix wrote in file named PodMatrix \n";
        matrix_file.close();
    }


    eigen_solver.compute( M_pod_matrix ); // solve M_pod_matrix psi = lambda psi

    int number_of_eigenvalues =  eigen_solver.eigenvalues().size();
    LOG(INFO)<<"Number of eigenvalues  : "<<number_of_eigenvalues<<"\n";
    //we copy eigenvalues in a std::vector beacause it's easier to manipulate it
    std::vector<double> eigen_values( number_of_eigenvalues );

    int too_small_index = 0;

    for ( int i=0; i<number_of_eigenvalues; i++ )
    {
        if ( imag( eigen_solver.eigenvalues()[i] )>1e-12 )
        {
            throw std::logic_error( "[POD::pod] ERROR : complex eigenvalues were found" );
        }

        eigen_values[i]=real( eigen_solver.eigenvalues()[i] );
    }


    //tab will contains index of maximum eigenvalues
    std::vector<int> max_idx;
    std::vector<double> copy_eigen_values ( eigen_values );
    int number_of_good_eigenvectors=0;
    bool go = true;
    double min_eigenvalue = doption(_name="pod.minimum-eigenvalue");
    //for(int maxmode=0; maxmode<M_Nm; maxmode++)
    auto max = std::max_element( copy_eigen_values.begin(), copy_eigen_values.end() );
    double idx = std::distance( copy_eigen_values.begin(), max );
    if( *max < min_eigenvalue )
        go=false;
    while( go )
    {
        number_of_good_eigenvectors++;
        max_idx.push_back( idx );
        copy_eigen_values[idx]=0;
        max = std::max_element( copy_eigen_values.begin(), copy_eigen_values.end() );
        idx = std::distance( copy_eigen_values.begin(), max );
        if( *max < min_eigenvalue )
            go = false;
        if( number_of_good_eigenvectors == M_Nm )
            go = false;
        if( number_of_good_eigenvectors == number_of_eigenvalues )
            go = false;
    }

    CHECK( number_of_good_eigenvectors > 0 )<<"The max eigenvalue "<<*max<<" is under the minimum eigenvalue set by the user "<<min_eigenvalue<<" and we have zero eigenvectors !\n";

    int position_of_largest_eigenvalue=max_idx[0];

    if( M_Nm == -1 )
    {
        //in this case the user doesn't know a priori
        //what value give to M_Nm
        //this case appears for when apply-POD-to-WN=true
        M_Nm = number_of_good_eigenvectors;
    }
    if ( M_Nm > number_of_good_eigenvectors && number_of_good_eigenvectors>0 && is_primal )
    {
        //in this case, the user set M_Nm but with respect to
        //the minimum eigenvalue set, we don't find enough eigenvectors
        M_Nm = number_of_good_eigenvectors;
    }

    for ( int i=0; i<M_Nm; i++ )
    {
        element_ptrtype mode ( new element_type( M_model->functionSpace() ) );
        mode->zero();

        int index=0;

        if( M_use_solutions )
        {
            M_bdf->setRestart( true );

            for ( M_bdf->restart(); !M_bdf->isFinished(); M_bdf->next() )
            {
                M_bdf->loadCurrent();
                double psi_k = real( eigen_solver.eigenvectors().col( position_of_largest_eigenvalue )[index] );
                M_bdf->unknown( 0 ).scale( psi_k );
                mode->add( 1 , M_bdf->unknown( 0 ) );
                index++;
            }
        }//if use solutions
        else
        {
            double eigenvalue = eigen_values[ position_of_largest_eigenvalue ];
            for(int j=0; j<size; j++)
            {
                double psi_k = real( eigen_solver.eigenvectors().col( position_of_largest_eigenvalue )[index] );
                mode->add( psi_k , elements_set[j] );
                index++;
            }
        }//if use elements set


        //check
        auto eigenvector = eigen_solver.eigenvectors().col( position_of_largest_eigenvalue );
        double eigenvalue = eigen_values[position_of_largest_eigenvalue];
        auto Aw = M_pod_matrix*eigenvector;
        auto lambdaw = eigenvalue*eigenvector;
        double Awnorm = Aw.norm();
        double lambdawnorm = lambdaw.norm();
        double eigenvectornorm = eigenvector.norm();
        double check_tol = doption(_name="pod.check-tol");
        CHECK( math::abs(Awnorm - lambdawnorm) < check_tol )<<" A w : "<<Awnorm<<" and lambda w : "<<lambdawnorm<<" so math::abs(A w - lambda w) : "<<math::abs(Awnorm - lambdawnorm)<<" and pod.check-tol : "<<check_tol<<" -- eigenvalue : "<<eigenvalue<<"\n";

        if( (i+1) < max_idx.size() )
            position_of_largest_eigenvalue=max_idx[i+1];
        ModeSet.push_back( *mode );
    }// loop on number of modes


    //check orthogonality
    if( boption(_name="pod.check-orthogonality") )
    {
        position_of_largest_eigenvalue=max_idx[0];
        for(int i=1; i<M_Nm; i++)
        {
            auto modei = ModeSet[i-1];
            for(int j=i+1; j<M_Nm; j++)
            {
                auto modej = ModeSet[j-1];
                double prod = M_model->scalarProductForPod( modei , modej );
                CHECK( prod < 1e-10 )<<"scalar product between mode "<<i<<" and mode "<<j<<" is not null and is "<<prod<<"\n";
            }
            double eigenvalue = eigen_values[ position_of_largest_eigenvalue ];
            double prod = M_model->scalarProductForPod( modei , modei );
            {
                CHECK( (prod - eigenvalue) < 1e-12 )<<"scalar product between mode "<<i<<" and mode "<<i<<" is "<<prod<<" and associated eigenvalue : "<<eigenvalue<<" ( i = "<<i<<") \n";
            }
            if( i < max_idx.size() )
                position_of_largest_eigenvalue=max_idx[i];
        }
    }

    if ( M_store_pod_matrix_format_octave )
    {
        std::cout<<"M_store_pod_matrix_format_octave = "<<M_store_pod_matrix_format_octave<<std::endl;
        std::ofstream eigenvalue_file;
        eigenvalue_file.open( "eigen_values.mat",std::ios::out );
        eigenvalue_file<<"# name: E\n";
        eigenvalue_file<<"# type: matrix\n";
        eigenvalue_file<<"# rows: "<<number_of_eigenvalues<<"\n";
        eigenvalue_file<<"# columns: 1\n";

        for ( int i=0; i<number_of_eigenvalues; i++ )
        {
            eigenvalue_file<<std::setprecision( 16 )<<eigen_values[i]<<"\n";
        }

        eigenvalue_file.close();


        std::ofstream eigenvector_file( ( boost::format( "eigen_vectors" ) ).str().c_str() );

        for ( int j=0; j<M_pod_matrix.cols(); j++ )
        {
            for ( int i=0; i<number_of_eigenvalues; i++ )
            {
                eigenvector_file<<std::setprecision( 16 )<<eigen_solver.eigenvectors().col( i )[j] <<" ";
            }

            eigenvector_file<<"\n";
        }

        eigenvector_file.close();


        for ( int i=0; i<M_Nm; i++ )
        {
            std::ofstream mode_file( ( boost::format( "mode_%1%" ) %i ).str().c_str() );
            element_type e = ModeSet[i];

            for ( size_type j=0; j<e.size(); j++ ) mode_file<<std::setprecision( 16 )<<e( j )<<"\n";

            mode_file.close();
        }

    }

    return M_Nm;

}


}//namespace Feel

#endif /* __POD_H */


