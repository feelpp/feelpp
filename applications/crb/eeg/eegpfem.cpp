/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-03-22

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
   \file eegpfem.cpp
   \author Sylvain Vallaghé <sylvain.vallaghe@gmail.com>
   \date 2011-06-27
 */
#include <eeg.hpp>
#include <feel/feelcrb/pfemapp.hpp>



int main( int argc, char** argv )
{
    Feel::PFemApp<Feel::EEG> app( argc, argv, Feel::makeEEGAbout( "eegpfem" ), Feel::makeEEGOptions()  );
    app.run();
}

