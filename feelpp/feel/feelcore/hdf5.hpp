/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
   \file hdf5.hpp
   \author Radu Popescu <radu.popescu@epfl.ch> (LifeV)
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org> (adaptation from LifeV to Feel++)
   \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
   \author Guillaume Dollé <gdolle@unistra.fr>
   \date 2015-10-01
 */
#ifndef FEELPP_HDF5_HPP
#define FEELPP_HDF5_HPP

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/worldcomm.hpp>

#include <map>
#include <string>

#ifdef FEELPP_HAS_HDF5
#undef OMPI_SKIP_MPICXX
#include <hdf5.h>

//Tell the compiler to restore the warning previously silented
//#pragma GCC diagnostic warning "-Wunused-variable"
//#pragma GCC diagnostic warning "-Wunused-parameter"

namespace Feel
{
/*!
  @brief Convenience wrapper for the C interface of the HDF5 library
  @author Radu Popescu <radu.popescu@epfl.ch>
  @author Christophe Prud'homme <christophe.prudhomme@feelpp.org> (adaptation from LifeV to Feel++)

  This class provides an easy way to write and read data from an HDF5 container.
  It is designed to handle a single open file at one time, with multiple open
  datasets simultaneously.
  It does not make any assumptions about data types, or dataset dimensions.
  It is a very thin wrapper on top of the HDF5 library and it is designed for
  high performance parallel operations.

  Usage:
      - open (or create) a file with HDF5IO::openFile
      - create or open existing data table with HDF5IO::createTable or
        HDF5IO::openTable (you can have multiple open tables at one time)
      - read or write blocks of data with HDF5IO::read or HDF5IO::write
      - close the tables after you are finished, with HDF5::closeTable
      - close the file: HDF5::closeFile
*/
class FEELPP_EXPORT HDF5
{
public:
    //! @name Public Types
    //@{
    typedef /*WorldComm*/boost::mpi::communicator comm_type;
    //@}

    //! @name Constructors and Destructor
    //@{
    //! Default empty constructor
    HDF5();

    //! Constructor
    /*!
     * Constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     * \param useCollectiveMetadataOps enable or disable collective metadata opts
     * \param useCollMetadataWrite enable or disable collective metadata write
     */
    HDF5( const std::string& fileName, const comm_type& comm,
          const bool& existing = false, bool useCollectiveMetadataOps = true, bool useCollMetadataWrite = true );

    // Copy constructor and assignment operator are disabled
    HDF5 (const HDF5&) = delete;
    HDF5& operator= (const HDF5&) = delete;

    //! Empty destructor
    virtual ~HDF5() {}
    //@}

    //! @name Public Methods
    //@{

    void enableCollectiveMetadataOps(bool enable = true) 
    {
        useCollectiveMetadataOps_ = enable;
    }

    void enableCollMetadataWrite(bool enable = true) 
    {
        useCollMetadataWrite_ = enable;
    }

    //! Check if group exist.
    /*!
     * \param groupName full path to the group under root.
     */
    bool groupExist( const std::string& groupName );

    //! Create or open a file.
    /*!
     * Create a file or open an existing file
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     * \param rdwr when opening an existing file, specifies whether to open it 
     *        in read-only or read/write mode.
     */
    void openFile( const std::string& fileName, const comm_type& comm,
                   const bool& existing, const bool& rdwr = false );

    //! Create a new group
    /*!
     * Create a new group in the open file
     */
    void createGroup( const std::string& tableName );

    //! Open group
    /*!
     * Open an existing group.
     * \param groupName Name of the group.
     * \param createIfNotExist Create all groups which does not exist in
     *        the given group path (default: true).
     */
    void openGroup( const std::string& groupName,
                    const bool& createIfNotExist = true );

    //! Open groups
    /*!
     * Open recursively all existing group from the given group name.
     * \param groupName Name of the group.
     * \param createIfNotExist Create all groups which does not exist in
     *        the given group path (default: true).
     */
    void openGroups( const std::string& groupName,
                     const bool& createIfNotExist = true );

    //! Create a new table
    /*!
     * Create a new table in the open file
     * \param tableName a string containing the table name
     * \param fileDataType data type that is to be used in the HDF5 container
     *        should be a standard HDF5 type, not a machine native type;
     *        consult HDF5 documentation for more information
     * \param tableDimensions array of hsize_t of size nbDims which holds the
     *        dimensions of the table
     * \param nbDims the number of dimensions of the array to store 
     */
    void createTable( const std::string& tableName,
                      hid_t& fileDataType,
                      hsize_t tableDimensions[],
                      unsigned int nbDims = 2 );

    //! Create a new table
    /*!
     * Create a new table in the open file under the given group.
     * \param groupName A string containing the group name.
     * \param tableName A string containing the table name.
     * \param fileDataType Data type that is to be used in the HDF5 container
     *        should be a standard HDF5 type, not a machine native type;
     *        Consult HDF5 documentation for more information.
     * \param tableDimensions Array of hsize_t of size nbDims which holds the
     *        dimensions of the table.
     * \param nbDims the number of dimensions of the array to store 
     */
    void createTable( const std::string& GroupName,
                      const std::string& tableName,
                      hid_t& fileDataType,
                      hsize_t tableDimensions[],
                      const bool& existing,
                      unsigned int nbDims = 2 );

    //! Open a new table
    /*!
     * Open a new table in the open file
     * \param tableName a string containing the table name
     * \param tableDimensions array of hsize_t of size 2 which will hold the
     *        dimensions of the table (output parameter)
     */
    void openTable( const std::string& tableName, hsize_t tableDimensions[] );

    //! Write
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the buffer that is to be written
     * \param currentCount an array of hsize_t of size nbDims describing the shape
     *        of the block to be written (see HDF5 documentation)
     * \param currentOffset an array of hsize_t of size nbDims describing the
     *        stride of the block to be written (see HDF5 documentation)
     * \param buffer pointer to a memory region containing the data to be
     *        written
     * \param nbDims the number of dimensions of the array to store 
     */
    void write( const std::string& tableName,
                hid_t& memDataType,
                hsize_t currentCount[],
                hsize_t currentOffset[],
                const void* buffer,
                unsigned int nbDims = 2 );

    //! Read
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the destination buffer
     * \param currentCount an array of hsize_t of size two describing the shape
     *        of the block to be read (see HDF5 documentation)
     * \param currentOffset an array of hsize_t of size two describing the
     *        stride of the block to be read (see HDF5 documentation)
     * \param buffer pointer to a memory region that represents the destination
     *        of the read operation
     * \param nbDims the number of dimensions of the array to store 
     */
    void read( const std::string& tableName,
               hid_t& memDataType,
               hsize_t currentCount[],
               hsize_t currentOffset[],
               void* buffer,
               int nbDims = 2 );

    //! Read_element
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the destination buffer
     * \param currentCount an array of hsize_t of size two describing the shape
     *        of the block to be read (see HDF5 documentation)
     * \param num_elem number of elements to be selected
     * \param coord a pointer to a buffer containing a serialized copy of a 2-dimensional array of zero-based values specifying the coordinates of the elements in the point selection
     * \param buffer pointer to a memory region that represents the destination
     *        of the read operation
     * \param nbDims the number of dimensions of the array to store
     */
     void read_elements( const std::string& tableName,
                         hid_t& memDataType, hsize_t currentCount[], hsize_t num_elem, const hsize_t * coord,
                         void* buffer, int nbDims = 2 );
    //! Close an open group.
    /*
     * \param groupName A string containing the group name.
     */
    void closeGroup( const std::string& groupName );

    //! Close recursively open groups.
    /*
     * \param groupName A string containing the group name.
     */
    void closeGroups( const std::string& groupName );

    //! Close open table.
    /*!
     * \param tableName A string containing the table name.
     */
    void closeTable( const std::string& tableName );

    //! Close an open file
    /*!
     * Call this when finished operating with a file.
     */
    void closeFile();

    //@}

private:
    
    bool useCollectiveMetadataOps_;
    bool useCollMetadataWrite_;

    // typedef for internal use
    typedef struct
    {
        hid_t filespace;
        hid_t dataset;
        hid_t plist;
    } tableHandle;

    //! Private Data Members
    //@{
    // HDF5 handles
    std::map<std::string, hid_t> M_groupList ;
    std::map<std::string, tableHandle> M_tableList;
    hid_t M_fileId;
    //@}
}; // class HDF5

} /* namespace Feel */

#endif /* FEELPP_HAS_HDF5 */

#endif /* FEELPP_HDF5_HPP */
