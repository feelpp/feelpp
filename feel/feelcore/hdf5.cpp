/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-16

  Copyright (C) 2013 Université de Strasbourg

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
   \file hdf5.cpp
   \author Radu Popescu <radu.popescu@epfl.ch>
   \author Christophe Prud'homme <prudhomme@unistra.fr> (adaptation from LifeV to Feel++)
   \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
   \date 2013-10-16
 */

#include <feel/feelcore/hdf5.hpp>

#ifdef FEELPP_HAS_HDF5

// ===================================================
// Constructor
// ===================================================

Feel::HDF5::HDF5 (const std::string& fileName, const comm_type& comm,
                       const bool& existing)
{
    openFile (fileName, comm, existing);
}

// ===================================================
// Public Methods
// ===================================================

void Feel::HDF5::openFile (const std::string& fileName,
                           const comm_type& comm,
                           const bool& existing)
{
    hid_t plistId;
    MPI_Comm mpiComm = comm.comm();
    MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plistId = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (plistId, mpiComm, info);

    // Create/open a file collectively and release property list identifier.
    if (existing)
    {
        M_fileId = H5Fopen (fileName.c_str(), H5F_ACC_RDONLY, plistId);
    }
    else
    {
        M_fileId = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC,
                              H5P_DEFAULT, plistId);
    }
    H5Pclose (plistId);
}

void Feel::HDF5::createGroup (const std::string& GroupName)
{
#ifdef H5_USE_16_API
        M_groupList [GroupName] = H5Gcreate (M_fileId, GroupName.c_str(), H5P_DEFAULT) ;
#else
        M_groupList [GroupName] = H5Gcreate (M_fileId, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
}

void Feel::HDF5::createTable (const std::string& GroupName, const std::string& tableName, hid_t& fileDataType, 
                                     hsize_t tableDimensions[], const bool& existing)
{
    if (!existing)
    {   
#ifdef H5_USE_16_API
        M_groupList [GroupName] = H5Gcreate (M_fileId, GroupName.c_str(), H5P_DEFAULT) ;
#else
        M_groupList [GroupName] = H5Gcreate (M_fileId, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
    }

    tableHandle& currentTable = M_tableList[GroupName+tableName] ;
    hid_t group_id = M_groupList[GroupName] ;

    currentTable.filespace = H5Screate_simple (2, tableDimensions,
                                               tableDimensions);
#ifdef H5_USE_16_API
    currentTable.dataset = H5Dcreate (group_id, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT);
#else
    currentTable.dataset = H5Dcreate (group_id, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
#endif
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
}

void Feel::HDF5::createTable (const std::string& tableName,
                                 hid_t& fileDataType,
                                 hsize_t tableDimensions[])
{
    tableHandle& currentTable = M_tableList[tableName];

    currentTable.filespace = H5Screate_simple (2, tableDimensions,
                                               tableDimensions);
#ifdef H5_USE_16_API
    currentTable.dataset = H5Dcreate (M_fileId, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT);
#else
    currentTable.dataset = H5Dcreate (M_fileId, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
#endif
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_INDEPENDENT);
}

void Feel::HDF5::openTable (const std::string& tableName,
                               hsize_t tableDimensions[])
{
    tableHandle& currentTable = M_tableList[tableName];

#ifdef H5_USE_16_API
    currentTable.dataset = H5Dopen (M_fileId, tableName.c_str());
#else
    currentTable.dataset = H5Dopen (M_fileId, tableName.c_str(), H5P_DEFAULT);
#endif
    currentTable.filespace = H5Dget_space (currentTable.dataset);
    H5Sget_simple_extent_dims (currentTable.filespace, tableDimensions, NULL);
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
}

void Feel::HDF5::write (const std::string& tableName,
                           hid_t& memDataType, hsize_t currentCount[],
                           hsize_t currentOffset[], void* buffer)
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (2, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dwrite (currentTable.dataset, memDataType, memspace,
              currentTable.filespace, currentTable.plist, buffer);

    H5Sclose (memspace);
}

void Feel::HDF5::read (const std::string& tableName,
                          hid_t& memDataType, hsize_t currentCount[],
                          hsize_t currentOffset[], void* buffer)
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (2, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dread (currentTable.dataset, memDataType, memspace,
             currentTable.filespace, currentTable.plist,
             buffer);

    H5Sclose (memspace);
}

void Feel::HDF5::closeTable (const std::string& tableName)
{
    tableHandle& currentTable = M_tableList[tableName];
    H5Dclose (currentTable.dataset);
    H5Sclose (currentTable.filespace);
    H5Pclose (currentTable.plist);
    M_tableList.erase (tableName);
}

void Feel::HDF5::closeGroup (const std::string& groupName)
{
    H5Gclose (M_groupList [groupName]) ;
    M_groupList.erase (groupName) ;
}

void Feel::HDF5::closeFile()
{
    for (std::map<std::string, hid_t>::iterator it = M_groupList.begin() ; it != M_groupList.end() ; it++)
        H5Gclose (it->second) ;
    H5Fclose (M_fileId);
}

#endif /* FEELPP_HAS_HDF5 */
