/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org> (adaptation from LifeV to Feel++)
   \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
   \author Guillaume Dollé <gdolle@unistra.fr>
   \date 2015-10-01
 */

#include <feel/feelcore/hdf5.hpp>

#ifdef FEELPP_HAS_HDF5


// ===================================================
// Constructor
// ===================================================

Feel::HDF5::HDF5()
    : 
    useCollectiveMetadataOps_(false),
    useCollMetadataWrite_(false)
{

}

Feel::HDF5::HDF5( const std::string& fileName, const comm_type& comm,
                  const bool& existing, bool useCollectiveMetadataOps, bool useCollMetadataWrite )
    :
    useCollectiveMetadataOps_(useCollectiveMetadataOps),
    useCollMetadataWrite_(useCollMetadataWrite)
{
    openFile (fileName, comm, existing);
}

// ===================================================
// Public Methods
// ===================================================

bool Feel::HDF5::groupExist( const std::string& groupName )
{
    herr_t status;
    // Turn off error print.
    H5Eset_auto( H5E_DEFAULT, NULL, NULL );
    //status = H5Oget_info_by_name( M_fileId, groupName.c_str(), 0, H5P_DEFAULT );
    status = H5Gget_objinfo( M_fileId, groupName.c_str(), 0, NULL );
    bool exist = ( status >= 0 );
    return exist;
}

void Feel::HDF5::openFile( const std::string& fileName,
                           const comm_type& comm,
                           const bool& existing, const bool& rdwr )
{
    hid_t plistId;
    MPI_Comm mpiComm = comm;
    MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plistId = H5Pcreate (H5P_FILE_ACCESS);
#if defined( H5_HAVE_PARALLEL )
    H5Pset_fapl_mpio (plistId, mpiComm, info);
    // Set properties based on member variables
    if (useCollectiveMetadataOps_) 
    {
        H5Pset_all_coll_metadata_ops(plistId, true);
    }
    if (useCollMetadataWrite_) 
    {
        H5Pset_coll_metadata_write(plistId, true);
    }
#endif

    // Create/open a file collectively and release property list identifier.
    if (existing)
    {
        /* if the file does not already exists
         * this is an error case: The user marked the file as existing
         * and it does not exists. Create the file so we don't get this error
         */
        if( !fs::exists( fileName ) )
        { M_fileId = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId); }
        /* Case where the file exists */
        else
        {
            if(rdwr)
            { M_fileId = H5Fopen (fileName.c_str(), H5F_ACC_RDWR, plistId); }
            else
            { M_fileId = H5Fopen (fileName.c_str(), H5F_ACC_RDONLY, plistId); }
        }
    }
    else
    {
        M_fileId = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC,
                              H5P_DEFAULT, plistId);
    }
    H5Pclose (plistId);
}

void Feel::HDF5::createGroup( const std::string& groupName )
{
#ifdef H5_USE_16_API
        M_groupList [groupName] = H5Gcreate (M_fileId, groupName.c_str(), H5P_DEFAULT) ;
#else
        M_groupList [groupName] = H5Gcreate (M_fileId, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) ;
#endif
        LOG(INFO) << "Create HDF5 group: " << groupName << "\n";
}

void Feel::HDF5::openGroup( const std::string& groupName,
                            const bool& createIfNotExist )
{
    bool exist = groupExist( groupName );
    if( exist )
    {
        M_groupList [groupName] = H5Gopen (M_fileId, groupName.c_str(), H5P_DEFAULT);
        LOG(INFO) << "Open HDF5 group: " << groupName;
    }
    else
    {
        if( createIfNotExist )
        {
            createGroup( groupName );
        }
    }
}

void Feel::HDF5::openGroups( const std::string& groupName,
                             const bool& createIfNotExist )
{
    const boost::char_separator<char> sep("/");
    boost::tokenizer<boost::char_separator<char>> tokens(groupName, sep);

    // Begin at root.
    std::string group="/";
    for (const auto& t : tokens)
    {
        group+=t;
        openGroup( group, createIfNotExist );
        group+="/";
    }
}

void Feel::HDF5::createTable ( const std::string& groupName,
                               const std::string& tableName,
                               hid_t& fileDataType,
                               hsize_t tableDimensions[],
                               const bool& existing,
                               unsigned int nbDims )
{
    if (!existing)
    {
       openGroups(groupName);
    }

    tableHandle& currentTable = M_tableList[groupName+tableName] ;
    hid_t group_id = M_groupList[groupName] ;

    currentTable.filespace = H5Screate_simple (nbDims, tableDimensions,
                                               tableDimensions);
#ifdef H5_USE_16_API
    currentTable.dataset = H5Dcreate (group_id, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT);
#else
    currentTable.dataset = H5Dcreate (group_id, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
#endif
#if defined( H5_HAVE_PARALLEL )
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
#endif    
}

void Feel::HDF5::createTable( const std::string& tableName,
                              hid_t& fileDataType,
                              hsize_t tableDimensions[],
                              unsigned int nbDims )
{
    tableHandle& currentTable = M_tableList[tableName];

    currentTable.filespace = H5Screate_simple (nbDims, tableDimensions,
                                               tableDimensions);
#ifdef H5_USE_16_API
    currentTable.dataset = H5Dcreate (M_fileId, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT);
#else
    currentTable.dataset = H5Dcreate (M_fileId, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
#endif
#if defined( H5_HAVE_PARALLEL )
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_INDEPENDENT);
#endif    
}

void Feel::HDF5::openTable( const std::string& tableName,
                            hsize_t tableDimensions[] )
{
    tableHandle& currentTable = M_tableList[tableName];

#ifdef H5_USE_16_API
    currentTable.dataset = H5Dopen (M_fileId, tableName.c_str());
#else
    currentTable.dataset = H5Dopen (M_fileId, tableName.c_str(), H5P_DEFAULT);
#endif
    currentTable.filespace = H5Dget_space (currentTable.dataset);
    H5Sget_simple_extent_dims (currentTable.filespace, tableDimensions, NULL);
#if defined( H5_HAVE_PARALLEL )    
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
#endif    
}

void Feel::HDF5::write( const std::string& tableName,
                        hid_t& memDataType, hsize_t currentCount[],
                        hsize_t currentOffset[], const void* buffer, unsigned int nbDims )
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (nbDims, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dwrite (currentTable.dataset, memDataType, memspace,
              currentTable.filespace, currentTable.plist, buffer);

    H5Sclose (memspace);
}

void Feel::HDF5::read( const std::string& tableName,
                       hid_t& memDataType, hsize_t currentCount[],
                       hsize_t currentOffset[], void* buffer, int nbDims )
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (nbDims, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dread (currentTable.dataset, memDataType, memspace,
             currentTable.filespace, currentTable.plist,
             buffer);

    H5Sclose (memspace);
}

void Feel::HDF5::read_elements( const std::string& tableName,
                                hid_t& memDataType, hsize_t currentCount[], hsize_t num_elem, const hsize_t * coord,
                                void* buffer, int nbDims )
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (nbDims, currentCount, currentCount);

    H5Sselect_elements( currentTable.filespace, H5S_SELECT_SET, num_elem, coord );

    H5Dread (currentTable.dataset, memDataType, memspace,
             currentTable.filespace, currentTable.plist,
             buffer);

    H5Sclose (memspace);
}

void Feel::HDF5::closeTable( const std::string& tableName )
{
    tableHandle& currentTable = M_tableList[tableName];
    H5Dclose (currentTable.dataset);
    H5Sclose (currentTable.filespace);
    H5Pclose (currentTable.plist);
    M_tableList.erase (tableName);
}

void Feel::HDF5::closeGroup( const std::string& groupName )
{
    H5Gclose (M_groupList [groupName]);
    M_groupList.erase (groupName);
}

void Feel::HDF5::closeGroups( const std::string& groupName )
{
    const boost::char_separator<char> sep("/");
    boost::tokenizer<boost::char_separator<char>> tokens(groupName, sep);
    std::vector<std::string> vtokens(tokens.begin(),tokens.end());

    std::string group = groupName;
    int gnsize = group.size();

    // Transform to begin and not end with a slash.
    if(group[0] != '/') group.insert( group.begin(), '/' );
    if(group[gnsize-1] == '/') group.pop_back();

    for( auto it=vtokens.rbegin(); it!=vtokens.rend(); ++it )
    {
        if( M_groupList.count(group) != 0 )
        {
            closeGroup(group);
            LOG(INFO) << "HDF5 close group: " << group << "\n";
        }
        if( group.size() > it->size() )
        {
            group.erase( group.end()-it->size()-1,group.end() );
        }
    }
}

void Feel::HDF5::closeFile()
{
    for (std::map<std::string, hid_t>::iterator it = M_groupList.begin() ; it != M_groupList.end() ; it++)
        H5Gclose (it->second) ;
    H5Fclose (M_fileId);
}

#endif /* FEELPP_HAS_HDF5 */
