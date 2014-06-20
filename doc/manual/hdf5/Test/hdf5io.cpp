#include <iostream>
#include <string>

#include "hdf5.h"
#define FILE "test.h5"

int main (void)
{
    const std::string DATASET_NAME ("dset") ;
    std::cout << "HDF5 I/O TEST" << std::endl ;

    hid_t       file_id;   /* file identifier */
    herr_t      status;

    int data [4][6] ;

    for (int i = 0 ; i < 4 ; i ++)
        for (int j = 0 ; j < 6 ; j ++)
            data[i][j] = i * 6 + j ;

    int data1 [4][6] ;

    for (int i = 0 ; i < 4 ; i ++)
        for (int j = 0 ; j < 6 ; j ++)
            data1[i][j] = i * 6 + j + 24 ;

    /* Create a new file */
    file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dims [2] ;
    dims[0] = 4 ;
    dims[1] = 6 ;

    /* Create a dataspace */
    hid_t dataspace_id = H5Screate_simple (2, dims, dims) ;

   
    hid_t group_id = H5Gcreate (file_id, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
    hid_t group_id1 = H5Gcreate (file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
    
    /*Create a dataset */
    hid_t dataset_id = H5Dcreate (group_id, "./Nodes", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
    hid_t dataset1_id = H5Dcreate (group_id, "./Elements", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ; 
    hid_t dataset2_id = H5Dcreate (group_id1, "./Nodes", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ; 
    hid_t dataset3_id = H5Dcreate (group_id1, "./Elements", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ; 

    H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) ;
    H5Dwrite (dataset1_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1) ;
    /* Close dataspace properly */
    status = H5Sclose (dataspace_id) ;
    /* Close dataset properly */
    status = H5Dclose (dataset_id) ;
    status = H5Dclose (dataset1_id) ;
    status = H5Dclose (dataset2_id) ;
    status = H5Dclose (dataset3_id) ;
    /* Close group properly */
    status = H5Gclose (group_id) ;
    status = H5Gclose (group_id1) ;
    /* Terminate access to the file. */
    status = H5Fclose(file_id); 
    return 0;  // successfully terminated

}
