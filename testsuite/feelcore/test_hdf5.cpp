/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
       Date: 2015-04-09

  Copyright (C) 2015 UJF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_hdf5.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2015-04-09
 */

#define BOOST_TEST_MODULE test_hdf5
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/hdf5.hpp>


namespace test_hdf5
{
// example take from https://www.hdfgroup.org/HDF5/doc/Intro/IntroExamples.html#CreateExample
void run1()
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[2];              /* dataset dimensions */
    herr_t      status;
    int NX=5, NY=6;
    std::string FILE = "SDS.h5";
    int RANK = 2;
    std::string DATASETNAME = "IntArray";
    int **      data;          /* data to write */
    int         i, j;

    /* allocate data */
    data = new int*[NX];
    for(i = 0; i < NX; i++)
    {
        data[i] = new int[NY];
    }

    /*
     * Data  and output buffer initialization.
     */
    for (j = 0; j < NX; j++) {
	for (i = 0; i < NY; i++)
	    data[j][i] = i + j;
    }
    /*
     * 0 1 2 3 4 5
     * 1 2 3 4 5 6
     * 2 3 4 5 6 7
     * 3 4 5 6 7 8
     * 4 5 6 7 8 9
     */

    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    file = H5Fcreate(FILE.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Describe the size of the array and create the data space for fixed
     * size dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY;
    dataspace = H5Screate_simple(RANK, dimsf, NULL);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    datatype = H5Tcopy(H5T_NATIVE_INT);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
#ifdef H5_USE_16_API
    dataset = H5Dcreate(file, DATASETNAME.c_str(), datatype, dataspace,
			H5P_DEFAULT);
#else
    dataset = H5Dcreate(file, DATASETNAME.c_str(), datatype, dataspace,
                        H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT );
#endif
    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, data);

    /*
     * Close/release resources.
     */
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);

    /*
     * Release memory
     */
    for(i = 0; i < NX; i++)
    {
        delete[] data[i];
    }
    delete[] data;
}

void run2()
{
    using namespace Feel;

    auto mycomm = Environment::worldComm();
    HDF5 hdf5;
    std::string filename="file.h5";

    //create file
    std::vector<double> eltValues = { 5.2, 4.8 };
    std::vector<double> eltValues2 = { 3.6, 2.7, 8.2 };
    std::vector<uint> sizeValues(1, ((mycomm.rank()%2)==0)?eltValues.size():eltValues2.size() );

    hdf5.openFile( filename, Environment::worldComm(), false );
    int rank=2;
    hsize_t * dims;
    hsize_t * dims2;
    hsize_t * offset;
    hsize_t * dimsElt;
    hsize_t * dimsElt2;
    hsize_t * offsetElt;

    /* memory allocation */
    dims = new hsize_t[rank];
    dims2 = new hsize_t[rank];
    offset = new hsize_t[rank];
    dimsElt = new hsize_t[rank];
    dimsElt2 = new hsize_t[rank];
    offsetElt = new hsize_t[rank];

    dims[0] = mycomm.size();dims[1] = 1;
    dims2[0] = sizeValues.size();dims2[1] = 1;
    offset[0] = mycomm.rank(); offset[1] = 0;

    // create size tab
    hdf5.createTable( "size", H5T_NATIVE_UINT, dims );
    hdf5.write( "size", H5T_NATIVE_UINT, dims2, offset, sizeValues.data() );
    hdf5.closeTable( "size" );

    if ( (mycomm.size()%2) == 0 )
        dimsElt[0] = (mycomm.size()/2)*(eltValues.size()+eltValues2.size());
    else
        dimsElt[0] = ((mycomm.size()-1)/2)*(eltValues.size()+eltValues2.size())+eltValues.size();
    dimsElt[1] = 1;

    if  ( (mycomm.rank()%2) == 0 )
    {
        offsetElt[0] = (mycomm.rank()/2)*(eltValues.size()+eltValues2.size());
        dimsElt2[0] = eltValues.size();
    }
    else
    {
        offsetElt[0] = ((mycomm.rank()-1)/2)*(eltValues.size()+eltValues2.size())+eltValues.size();
        dimsElt2[0] = eltValues2.size();
    }
    dimsElt2[1] = 1;
    offsetElt[1] = 0;

    // create double tab
    hdf5.createTable( "element", H5T_NATIVE_DOUBLE, dimsElt );
    hdf5.write( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, ((mycomm.rank()%2)==0)?eltValues.data() : eltValues2.data() );
    hdf5.closeTable( "element" );

    hdf5.closeFile();

    // reload file
    HDF5 hdf5b;
    std::vector<uint> sizeValuesReload( 1 );
    std::vector<double> eltValuesReload;

    hdf5b.openFile( filename, Environment::worldComm(), true );

    hdf5b.openTable( "size",dims );
    hdf5b.read( "size", H5T_NATIVE_UINT, dims2, offset, sizeValuesReload.data() );
    hdf5b.closeTable( "size" );

    eltValuesReload.resize( sizeValuesReload[0] );

    hdf5b.openTable( "element",dimsElt );
    hdf5b.read( "element", H5T_NATIVE_DOUBLE, dimsElt2, offsetElt, eltValuesReload.data() );
    hdf5b.closeTable( "element" );

    hdf5b.closeFile();

    // check reload
    BOOST_CHECK( sizeValues.size() == sizeValuesReload.size() );
    for ( int k = 0 ; k < sizeValues.size() ; ++k)
        BOOST_CHECK( sizeValues[k] == sizeValuesReload[k] );

    if  ( (mycomm.rank()%2) == 0 )
    {
        BOOST_CHECK( eltValues.size() == eltValuesReload.size() );
        for ( int k = 0 ; k < eltValues.size() ; ++k)
            BOOST_CHECK_CLOSE( eltValues[k], eltValuesReload[k], 1e-9 );
    }
    else
    {
        BOOST_CHECK( eltValues2.size() == eltValuesReload.size() );
        for ( int k = 0 ; k < eltValues2.size() ; ++k)
            BOOST_CHECK_CLOSE( eltValues2[k], eltValuesReload[k], 1e-9 );
    }

    /*
     * Release memory
     */
    delete[] dims;
    delete[] dims2;
    delete[] offset;
    delete[] dimsElt;
    delete[] dimsElt2;
    delete[] offsetElt;
}

// Test groups.
void run3()
{
    using namespace Feel;

    auto comm = Environment::worldComm();
    HDF5 h5;

    bool isLoad = false;
    h5.openFile( "groups.h5", comm, isLoad );
    isLoad = true;

    h5.openGroups( "class_1/sub_1" );
    h5.closeGroups( "class_1/sub_1" );

    h5.openGroups( "class_2" );
    h5.closeGroups( "class_2" );

    // Check slash.
    h5.openGroups( "/class_3" );
    h5.closeGroups( "/class_3" );

    h5.openGroups( "/class_4/" );
    h5.closeGroups( "/class_4/" );

    h5.openGroups( "class_5/" );
    h5.closeGroups( "class_5/" );

    // Check slash with sub.
    h5.openGroups( "class_6/sub_6" );
    h5.closeGroups( "class_6/sub_6" );

    h5.openGroups( "/class_7/sub_7" );
    h5.closeGroups( "/class_7/sub_7" );

    h5.openGroups( "/class_8/sub_8/" );
    h5.closeGroups( "/class_8/sub_8/" );

    h5.openGroups( "class_9/sub_9/" );
    h5.closeGroups( "class_9/sub_9/" );

    // Check reopen group.
    h5.openGroups( "class_1" );
    h5.closeGroups( "class_1" );

    // Check create when top group exist.
    h5.openGroups( "class_1/sub_2/subsub_1" );
    h5.closeGroups( "class_1/sub_2/subsub_1" );

    h5.closeFile();


    // Check that groups exist in the file.
    h5.openFile( "groups.h5", comm, isLoad );

    BOOST_CHECK( h5.groupExist( "/class_1" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_2" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_3" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_4" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_5" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_6" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_7" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_8" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_6/sub_6" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_7/sub_7" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_8/sub_8" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_9/sub_9" ) == true );
    BOOST_CHECK( h5.groupExist( "/class_1/sub_2/subsub_1" ) == true );

    h5.closeFile();
}


} // namespace test_hdf5

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( hdf5 )

BOOST_AUTO_TEST_CASE( hdf5 )
{
    test_hdf5::run1();
    test_hdf5::run2();
    test_hdf5::run3();
}

BOOST_AUTO_TEST_SUITE_END()

