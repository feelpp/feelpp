/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-03-19

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file wrapper.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-19
 */
#ifndef __Wrapper_H
#define __Wrapper_H 1

#define FUNC_EXEC_BODY_IN_TEMPDIR( classname, dirname )            \
{                                                                       \
char* currentWorkingDirectory = getCurrentWorkingDirectory (0) ;        \
char* temporaryDirectory=createTemporaryDirectory(#dirname,p_exchangedData,0); \
                                                                        \
int rc = 0;                                                             \
try                                                                     \
{                                                                       \
    CAST(classname*,p_state)->run( INPOINT_ARRAY, INPOINT_SIZE, OUTPOINT_ARRAY, OUTPOINT_SIZE ); \
}                                                                       \
catch( ... )                                                            \
{                                                                       \
    rc = 1;                                                             \
}                                                                       \
if (rc) {                                                               \
    PRINT( "Error in class "#classname );                               \
    return WRAPPER_EXECUTION_ERROR;                                     \
}                                                                       \
deleteTemporaryDirectory ( temporaryDirectory , rc, 0 ) ;               \
free ( currentWorkingDirectory ) ;                                      \
}

#endif /* __Wrapper_H */
