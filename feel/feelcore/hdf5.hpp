/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

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
   \file hdf5.hpp
  @author Radu Popescu <radu.popescu@epfl.ch> (LifeV)
   \author Christophe Prud'homme <prudhomme@unistra.fr> (adaptation from LifeV to Feel++)
   \author Benjamin Vanthong <benjamin.vanthong@gmail.com>
   \date 2013-10-16
 */
#ifndef FEELPP_HDF5_HPP
#define FEELPP_HDF5_HPP

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/worldcomm.hpp>

#include <map>
#include <string>

#ifdef FEELPP_HAS_HDF5

#include <hdf5.h>

//Tell the compiler to restore the warning previously silented
//#pragma GCC diagnostic warning "-Wunused-variable"
//#pragma GCC diagnostic warning "-Wunused-parameter"

namespace Feel
{

/*!
  @brief Convenience wrapper for the C interface of the HDF5 library
  @author Radu Popescu <radu.popescu@epfl.ch>
  @author Christophe Prud'homme <prudhomme@unistra.fr> (adaptation from LifeV to Feel++)

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
class HDF5
{
public:
    //! @name Public Types
    //@{
    typedef WorldComm comm_type;
    //@}

    //! @name Constructors and Destructor
    //@{
    //! Default empty constructor
    HDF5() {}

    //! Constructor
    /*!
     * Constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     */
    HDF5 (const std::string& fileName, const comm_type& comm,
            const bool& existing = false);

    //! Empty destructor
    virtual ~HDF5() {}
    //@}

    //! @name Public Methods
    //@{
    //! Open
    /*!
     * Create a file or open an existing file
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     */
    void openFile (const std::string& fileName, const comm_type& comm,
                   const bool& existing);

    //! Create a new group
    /*!
     * Create a new group in the open file
     */
    void createGroup (const std::string& tableName);
    //! Create a new table
    /*!
     * Create a new table in the open file
     * \param tableName a string containing the table name
     * \param fileDataType data type that is to be used in the HDF5 container
     *        should be a standard HDF5 type, not a machine native type;
     *        consult HDF5 documentation for more information
     * \param tableDimensions array of hsize_t of size 2 which holds the
     *        dimensions of the table
     */
    void createTable (const std::string& tableName, hid_t& fileDataType,
                      hsize_t tableDimensions[]);

    void createTable (const std::string& GroupName, const std::string& tableName, hid_t& fileDataType,
                             hsize_t tableDimensions[], const bool& existing);

    //! Open a new table
    /*!
     * Open a new table in the open file
     * \param tableName a string containing the table name
     * \param tableDimensions array of hsize_t of size 2 which will hold the
     *        dimensions of the table (output parameter)
     */
    void openTable (const std::string& tableName, hsize_t tableDimensions[]);
    //! Write
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the buffer that is to be written
     * \param currentCount an array of hsize_t of size two describing the shape
     *        of the block to be written (see HDF5 documentation)
     * \param currentOffset an array of hsize_t of size two describing the
     *        stride of the block to be written (see HDF5 documentation)
     * \param buffer pointer to a memory region containing the data to be
     *        written
     */
    void write (const std::string& tableName,
                hid_t& memDataType, hsize_t currentCount[],
                hsize_t currentOffset[], void* buffer);
    //! Write
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
     */
    void read (const std::string& tableName,
               hid_t& memDataType, hsize_t currentCount[],
               hsize_t currentOffset[], void* buffer);
    //! Write
    /*!
     * \param tableName a string containing the table name
     */
    void closeGroup (const std::string& groupName);
    void closeTable (const std::string& tableName);
    //! Close an open file
    /*!
     * Call this when finished operating with a file
     */
    void closeFile();
    //@}

private:
    // typedef for internal use
    typedef struct
    {
        hid_t filespace;
        hid_t dataset;
        hid_t plist;
    } tableHandle;

    // Copy constructor and assignment operator are disabled
    HDF5 (const HDF5&);
    HDF5& operator= (const HDF5&);

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
