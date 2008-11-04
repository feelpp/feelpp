/* -*- c++ -*- (enables emacs c++ mode)                                    */
/* *********************************************************************** */
/*                                                                         */
/* Library :  Generic Matrix Methods  (gmm)                                */
/* File    :  gmm_inoutput.h : input and output of matrices.               */
/*     									   */
/* Date : July 8, 2003.                                                    */
/* Author : Yves Renard, Yves.Renard@gmm.insa-tlse.fr,                     */
/*          Julien Pommier, Julien.Pommier@gmm.insa-tlse.fr                */
/*                                                                         */
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2003  Yves Renard.                                        */
/*                                                                         */
/* This file is a part of GMM++                                            */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#ifndef GMM_INOUTPUT_H
#define GMM_INOUTPUT_H

#include <stdio.h>

namespace gmm {

  /*************************************************************************/
  /*                                                                       */
  /*  Functions to read and write Harwell Boeing format.                   */
  /*                                                                       */
  /*************************************************************************/

  // Fri Aug 15 16:29:47 EDT 1997
  // 
  //                      Harwell-Boeing File I/O in C
  //                               V. 1.0
  // 
  //          National Institute of Standards and Technology, MD.
  //                            K.A. Remington
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                                NOTICE
  //
  // Permission to use, copy, modify, and distribute this software and
  // its documentation for any purpose and without fee is hereby granted
  // provided that the above copyright notice appear in all copies and
  // that both the copyright notice and this permission notice appear in
  // supporting documentation.
  //
  // Neither the Author nor the Institution (National Institute of Standards
  // and Technology) make any representations about the suitability of this 
  // software for any purpose. This software is provided "as is" without 
  // expressed or implied warranty.
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  

  inline void IOHBTerminate(const char *a) { DAL_THROW(dal::failure_error, a); }

  inline bool is_complex_double__(std::complex<double>) { return true; }
  inline bool is_complex_double__(double) { return false; }

  inline int ParseIfmt(const char *fmt, int* perline, int* width) {
    if (sscanf(fmt, " (%dI%d)", perline, width) != 2) 
      DAL_THROW(dal::failure_error, "invalid HB I-format : " << fmt);
    return *width;
  }
  
  inline int ParseRfmt(const char *fmt, int* perline, int* width,
		       int* prec, int* flag) {
    char p;
    *perline = *width = *flag = *prec = 0;
    if (sscanf(fmt, " (%d%c%d.%d)", perline, &p, width, prec) < 3 || 
	!strchr("PEDF", p)) 
      DAL_THROW(dal::failure_error, "invalid HB REAL format : " << fmt);
    *flag = p;
    return *width;
  }

  /** matrix input/output for Harwell-Boeing format */
  struct HarwellBoeing_IO {
    int nrows() const { return Nrow; }
    int ncols() const { return Ncol; }
    int nnz() const { return Nnzero; }
    int is_complex() const { return Type[0] == 'C'; }
    int is_symmetric() const { return Type[1] == 'S'; }
    int is_hermitian() const { return Type[1] == 'H'; }
    HarwellBoeing_IO() { clear(); }
    HarwellBoeing_IO(const char *filename) { clear(); open(filename); }
    ~HarwellBoeing_IO() { close(); }
    /* open filename and reads header */
    void open(const char *filename);
    /* read the opened file */
    template <typename T, int shift> void read(csc_matrix<T, shift>& A);
    template <typename MAT> void read(MAT &M);
    /* save the matrix */
    template <typename T, int shift> static void write(const char *filename, const csc_matrix<T, shift>& A);
    template <typename MAT> static void write(const char *filename, const MAT& A);
  private:
    FILE *f;
    char Title[73], Key[9], Rhstype[4], Type[4];
    int Nrow, Ncol, Nnzero, Nrhs;
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd; 
    int lcount;


    void close() { if (f) fclose(f); clear(); }
    void clear() { 
      Nrow = Ncol = Nnzero = Nrhs = 0; f = 0; lcount = 0;
      memset(Type, 0, sizeof Type); 
      memset(Key, 0, sizeof Key); 
      memset(Title, 0, sizeof Title); 
    }
    char *getline(char *buf) { 
      fgets(buf, BUFSIZ, f); ++lcount;
      if (sscanf(buf,"%*s") < 0) 
	DAL_THROW(dal::failure_error, "blank line in HB file at line " << lcount);
      return buf;
    }

    int substrtoi(const char *p, size_type len) {
      char s[100]; len = std::min(len, sizeof s - 1);
      strncpy(s,p,len); s[len] = 0; return atoi(s);
    }
    double substrtod(const char *p, size_type len, int Valflag) {
      char s[100]; len = std::min(len, sizeof s - 1);
      strncpy(s,p,len); s[len] = 0;
      if ( Valflag != 'F' && !strchr(s,'E')) {
	/* insert a char prefix for exp */
	int last = strlen(s);
	for (int j=last+1;j>=0;j--) {
	  s[j] = s[j-1];
	  if ( s[j] == '+' || s[j] == '-' ) {
	    s[j-1] = Valflag;                    
	    break;
	  }
	}
      }
      return atof(s);
    }
    template <typename IND_TYPE>   
    int readHB_data(IND_TYPE colptr[], IND_TYPE rowind[], 
		    double val[]) {
      /************************************************************************/
      /*  This function opens and reads the specified file, interpreting its  */
      /*  contents as a sparse matrix stored in the Harwell/Boeing standard   */
      /*  format and creating compressed column storage scheme vectors to hold*/
      /*  the index and nonzero value information.                            */
      /*                                                                      */
      /*    ----------                                                        */
      /*    **CAVEAT**                                                        */
      /*    ----------                                                        */
      /*  Parsing real formats from Fortran is tricky, and this file reader   */
      /*  does not claim to be foolproof.   It has been tested for cases when */
      /*  the real values are printed consistently and evenly spaced on each  */
      /*  line, with Fixed (F), and Exponential (E or D) formats.             */
      /*                                                                      */
      /*  **  If the input file does not adhere to the H/B format, the  **    */
      /*  **             results will be unpredictable.                 **    */
      /*                                                                      */
      /************************************************************************/
      int i,ind,col,offset,count;
      int Ptrperline, Ptrwidth, Indperline, Indwidth;
      int Valperline, Valwidth, Valprec, Nentries;
      int Valflag;           /* Indicates 'E','D', or 'F' float format */
      char line[BUFSIZ];

      /*  Parse the array input formats from Line 3 of HB file  */
      ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
      ParseIfmt(Indfmt,&Indperline,&Indwidth);
      if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
	ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
      }
    
      /*  Read column pointer array:   */
      offset = 0;         /* if base 0 storage is declared (via macro def),  */
      /* then storage entries are offset by 1            */
    
      for (count = 0, i=0;i<Ptrcrd;i++) {
	getline(line);
	for (col = 0, ind = 0;ind<Ptrperline;ind++) {
	  if (count > Ncol) break;
	  colptr[count] = substrtoi(line+col,Ptrwidth)-offset;
	  count++; col += Ptrwidth;
	}
      }
    
      /*  Read row index array:  */    
      for (count = 0, i=0;i<Indcrd;i++) {
	getline(line);
	for (col = 0, ind = 0;ind<Indperline;ind++) {
	  if (count == Nnzero) break;
	  rowind[count] = substrtoi(line+col,Indwidth)-offset;
	  count++; col += Indwidth;
	}
      }
    
      /*  Read array of values:  */
      if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
	if ( Type[0] == 'C' ) Nentries = 2*Nnzero;
	else Nentries = Nnzero;
      
	count = 0;
	for (i=0;i<Valcrd;i++) {
	  getline(line);
	  if (Valflag == 'D')  {
            // const_cast Due to aCC excentricity
	    char *p; while( (p = const_cast<char *>(strchr(line,'D')) )) *p = 'E';
	  }
	  for (col = 0, ind = 0;ind<Valperline;ind++) {
	    if (count == Nentries) break;
	    val[count] = substrtod(line+col, Valwidth, Valflag);
	    count++; col += Valwidth;
	  }
	}
      }
      return 1;
    }
  };
  
  inline void HarwellBoeing_IO::open(const char *filename) {
    int Totcrd,Neltvl,Nrhsix;
    char line[BUFSIZ];
    close();
    f = fopen(filename, "r");
    if (!f) { DAL_THROW(dal::failure_error, "could not open " << filename); }
    /* First line: */
    sscanf(getline(line), "%72c%8s", Title, Key);
    Key[8] = Title[72] = 0;
    /* Second line: */
    Totcrd = Ptrcrd = Indcrd = Valcrd = Rhscrd = 0;
    sscanf(getline(line), "%d%d%d%d%d", &Totcrd, &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd);
    
    /* Third line: */
    Nrow = Ncol = Nnzero = Neltvl = 0;
    if (sscanf(getline(line), "%3c%d%d%d%d", Type, &Nrow, &Ncol, &Nnzero, &Neltvl) < 1)
      IOHBTerminate("Invalid Type info, line 3 of Harwell-Boeing file.\n");
    std::for_each(Type, Type+3, toupper);
    
      /*  Fourth line:  */
    if ( sscanf(getline(line), "%16c%16c%20c%20c",Ptrfmt,Indfmt,Valfmt,Rhsfmt) < 3)
      IOHBTerminate("Invalid format info, line 4 of Harwell-Boeing file.\n"); 
    Ptrfmt[16] = Indfmt[16] = Valfmt[20] = Rhsfmt[20] = 0;
    
    /*  (Optional) Fifth line: */
    if (Rhscrd != 0 ) { 
      Nrhs = Nrhsix = 0;
      if ( sscanf(getline(line), "%3c%d%d", Rhstype, &Nrhs, &Nrhsix) != 1) 
	IOHBTerminate("Invalid RHS type information, line 5 of"
		      " Harwell-Boeing file.\n");
    }
  }

  /* only valid for double and complex<double> csc matrices */
  template <typename T, int shift> void
  HarwellBoeing_IO::read(csc_matrix<T, shift>& A) {

    typedef typename csc_matrix<T, shift>::IND_TYPE IND_TYPE;

    if (!f) DAL_THROW(dal::failure_error, "no file opened!");    
    if (Type[0] == 'P')
      DAL_THROW(dal::failure_error, "Bad HB matrix format (pattern matrices not supported)");
    if (is_complex_double__(T()) && Type[0] == 'R') 
      DAL_THROW(dal::failure_error, "Bad HB matrix format (file contains a REAL matrix)");
    if (!is_complex_double__(T()) && Type[0] == 'C') 
      DAL_THROW(dal::failure_error, "Bad HB matrix format (file contains a COMPLEX matrix)");
    if (A.pr) { delete[] A.pr; delete[] A.ir; delete[] A.jc; }
    A.nc = ncols(); A.nr = nrows();
    A.pr = 0;
    A.jc = new IND_TYPE[ncols()+1];
    A.ir = new IND_TYPE[nnz()];
    A.pr = new T[nnz()];
    readHB_data(A.jc, A.ir, (double*)A.pr);
    for (int i = 0; i <= ncols(); ++i) { A.jc[i] += shift; A.jc[i] -= 1; }
    for (int i = 0; i < nnz(); ++i)    { A.ir[i] += shift; A.ir[i] -= 1; }
  }

  template <typename MAT> void 
  HarwellBoeing_IO::read(MAT &M) {
    csc_matrix<typename gmm::linalg_traits<MAT>::value_type> csc;
    read(csc); 
    resize(M, mat_nrows(csc), mat_ncols(csc));
    copy(csc, M);
  }
  
  template <typename IND_TYPE> 
  inline int writeHB_mat_double(const char* filename, int M, int N, int nz,
				const IND_TYPE colptr[],
				const IND_TYPE rowind[], 
				const double val[], int Nrhs,
				const double /*rhs*/[], const double /*guess*/[],
				const double /*exact*/[], const char* Title,
				const char* Key, const char* Type, 
				const char* Ptrfmt, const char* Indfmt,
				const char* Valfmt, const char* Rhsfmt,
				const char* Rhstype, int shift) {
    /************************************************************************/
      /*  The writeHB function opens the named file and writes the specified  */
      /*  matrix and optional right-hand-side(s) to that file in              */
      /*  Harwell-Boeing format.                                              */
      /*                                                                      */
      /*  For a description of the Harwell Boeing standard, see:              */
      /*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989          */
      /*                                                                      */
      /************************************************************************/
      FILE *out_file;
      int i,entry,offset/* , j, acount, linemod */;
      int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
      int nvalentries, nrhsentries;
      int Ptrperline, Ptrwidth, Indperline, Indwidth;
      int Rhsperline, Rhswidth, Rhsprec, Rhsflag;
      int Valperline, Valwidth, Valprec;
      int Valflag;           /* Indicates 'E','D', or 'F' float format */
      char pformat[16],iformat[16],vformat[19],rformat[19];
    
      if ( Type[0] == 'C' )
	{ nvalentries = 2*nz; nrhsentries = 2*M; }
      else
	{ nvalentries = nz; nrhsentries = M; }
    
      if ( filename != NULL ) {
	if ( (out_file = fopen( filename, "w")) == NULL )
	  DAL_THROW(gmm::failure_error,"Error: Cannot open file: " << filename);
      } else out_file = stdout;
    
      if ( Ptrfmt == NULL ) Ptrfmt = "(8I10)";
      ParseIfmt(Ptrfmt, &Ptrperline, &Ptrwidth);
      sprintf(pformat,"%%%dd",Ptrwidth);
      ptrcrd = (N+1)/Ptrperline;
      if ( (N+1)%Ptrperline != 0) ptrcrd++;
    
      if ( Indfmt == NULL ) Indfmt =  Ptrfmt;
      ParseIfmt(Indfmt, &Indperline, &Indwidth);
      sprintf(iformat,"%%%dd",Indwidth);
      indcrd = nz/Indperline;
      if ( nz%Indperline != 0) indcrd++;
    
      if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
	if ( Valfmt == NULL ) Valfmt = "(4E20.13)";
	ParseRfmt(Valfmt, &Valperline, &Valwidth, &Valprec, &Valflag);
	if (Valflag == 'D') *strchr(Valfmt,'D') = 'E';
	if (Valflag == 'F')
	  sprintf(vformat, "%% %d.%df", Valwidth, Valprec);
	else
	  sprintf(vformat, "%% %d.%dE", Valwidth, Valprec);
	valcrd = nvalentries/Valperline;
	if ( nvalentries%Valperline != 0) valcrd++;
      } else valcrd = 0;
    
      if ( Nrhs > 0 ) {
	if ( Rhsfmt == NULL ) Rhsfmt = Valfmt;
	ParseRfmt(Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec, &Rhsflag);
	if (Rhsflag == 'F')
	  sprintf(rformat,"%% %d.%df",Rhswidth,Rhsprec);
	else
	  sprintf(rformat,"%% %d.%dE",Rhswidth,Rhsprec);
	if (Rhsflag == 'D') *strchr(Rhsfmt,'D') = 'E';
	rhscrd = nrhsentries/Rhsperline; 
	if ( nrhsentries%Rhsperline != 0) rhscrd++;
	if ( Rhstype[1] == 'G' ) rhscrd+=rhscrd;
	if ( Rhstype[2] == 'X' ) rhscrd+=rhscrd;
	rhscrd*=Nrhs;
      } else rhscrd = 0;
    
      totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;
    
    
      /*  Print header information:  */
    
      fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
	      ptrcrd, indcrd, valcrd, rhscrd);
      fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
      fprintf(out_file,"%-16s%-16s%-20s", Ptrfmt, Indfmt, Valfmt);
      //     if ( Nrhs != 0 ) {
      //       /*    Print Rhsfmt on fourth line and                                 */
      //       /*      optional fifth header line for auxillary vector information:  */
      //       fprintf(out_file,"%-20s\n%-14s%d\n",Rhsfmt,Rhstype,Nrhs);
      //     } else
      fprintf(out_file,"\n");
    
      offset = 1 - shift;  /* if base 0 storage is declared (via macro def), */
      /* then storage entries are offset by 1           */
    
      /*  Print column pointers:   */
      for (i = 0; i < N+1; i++) {
	entry = colptr[i]+offset;
	fprintf(out_file,pformat,entry);
	if ( (i+1)%Ptrperline == 0 ) fprintf(out_file,"\n");
      }
    
      if ( (N+1) % Ptrperline != 0 ) fprintf(out_file,"\n");
    
      /*  Print row indices:       */
      for (i=0;i<nz;i++) {
	entry = rowind[i]+offset;
	fprintf(out_file,iformat,entry);
	if ( (i+1)%Indperline == 0 ) fprintf(out_file,"\n");
      }
    
      if ( nz % Indperline != 0 ) fprintf(out_file,"\n");
    
      /*  Print values:            */
    
      if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
	for (i=0;i<nvalentries;i++) {
	  fprintf(out_file,vformat,val[i]);
	  if ( (i+1)%Valperline == 0 ) fprintf(out_file,"\n");
	}
	if ( nvalentries % Valperline != 0 ) fprintf(out_file,"\n");
      }
    
      if ( fclose(out_file) != 0) {
	DAL_THROW(gmm::failure_error,"Error closing file in writeHB_mat_double().");
      } else return 1;
    }

  template <typename T, int shift> void
  HarwellBoeing_IO::write(const char *filename, const csc_matrix<T, shift>& A) {
    const char *t = 0;    
    if (is_complex_double__(T()))
      if (mat_nrows(A) == mat_ncols(A)) t = "CUA"; else t = "CRA";
    else
      if (mat_nrows(A) == mat_ncols(A)) t = "RUA"; else t = "RRA";
    writeHB_mat_double(filename, mat_nrows(A), mat_ncols(A),
		       A.jc[mat_ncols(A)], A.jc, A.ir,
		       (double *)A.pr,
		       0, 0, 0, 0, "GETFEM++ CSC MATRIX", "CSCMAT",
		       t, 0, 0, 0, 0, "F", shift);
  }

  template <typename MAT> void
  HarwellBoeing_IO::write(const char *filename, const MAT& A) {
    gmm::csc_matrix<typename gmm::linalg_traits<MAT>::value_type> 
      tmp(gmm::mat_nrows(A), gmm::mat_ncols(A));
    gmm::copy(A,tmp); 
    HarwellBoeing_IO::write(filename, tmp);
  }
  

  /** save a "double" or "std::complex<double>" matrix into a HarwellBoeing file */
  template <typename T, int shift> inline void
  Harwell_Boeing_save(const char *filename, const csc_matrix<T, shift>& A) {
    HarwellBoeing_IO h; h.write(filename, A);
  }

  /** load a "double" or "std::complex<double>" matrix from a HarwellBoeing file */
  template <typename T, int shift> void
  Harwell_Boeing_load(const char *filename, csc_matrix<T, shift>& A) {
    HarwellBoeing_IO h(filename); h.read(A);
  }


  /*************************************************************************/
  /*                                                                       */
  /*  Functions to read and write MatrixMarket format.                     */
  /*                                                                       */
  /*************************************************************************/

  /* 
   *   Matrix Market I/O library for ANSI C
   *
   *   See http://math.nist.gov/MatrixMarket for details.
   *
   *
   */

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

  typedef char MM_typecode[4];


  /******************* MM_typecode query functions *************************/

#define mm_is_matrix(typecode)	        ((typecode)[0]=='M')
  
#define mm_is_sparse(typecode)	        ((typecode)[1]=='C')
#define mm_is_coordinate(typecode)      ((typecode)[1]=='C')
#define mm_is_dense(typecode)	        ((typecode)[1]=='A')
#define mm_is_array(typecode)	        ((typecode)[1]=='A')
  
#define mm_is_complex(typecode)	        ((typecode)[2]=='C')
#define mm_is_real(typecode)	        ((typecode)[2]=='R')
#define mm_is_pattern(typecode)	        ((typecode)[2]=='P')
#define mm_is_integer(typecode)         ((typecode)[2]=='I')
  
#define mm_is_symmetric(typecode)       ((typecode)[3]=='S')
#define mm_is_general(typecode)	        ((typecode)[3]=='G')
#define mm_is_skew(typecode)	        ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)       ((typecode)[3]=='H')
  
  /******************* MM_typecode modify fucntions ************************/

#define mm_set_matrix(typecode)	        ((*typecode)[0]='M')
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define mm_set_array(typecode)	        ((*typecode)[1]='A')
#define mm_set_dense(typecode)	        mm_set_array(typecode)
#define mm_set_sparse(typecode)	        mm_set_coordinate(typecode)

#define mm_set_complex(typecode)        ((*typecode)[2]='C')
#define mm_set_real(typecode)	        ((*typecode)[2]='R')
#define mm_set_pattern(typecode)        ((*typecode)[2]='P')
#define mm_set_integer(typecode)        ((*typecode)[2]='I')


#define mm_set_symmetric(typecode)      ((*typecode)[3]='S')
#define mm_set_general(typecode)        ((*typecode)[3]='G')
#define mm_set_skew(typecode)	        ((*typecode)[3]='K')
#define mm_set_hermitian(typecode)      ((*typecode)[3]='H')

#define mm_clear_typecode(typecode)     ((*typecode)[0]=(*typecode)[1]= \
			       	        (*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


  /******************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17


  /******************** Matrix Market internal definitions *****************

   MM_matrix_typecode: 4-character sequence

	                object 	    sparse/   	data        storage 
	                            dense     	type        scheme

   string position:	 [0]        [1]		[2]         [3]

   Matrix typecode:     M(atrix)    C(oord)	R(eal)      G(eneral)
		                    A(array)    C(omplex)   H(ermitian)
	                                        P(attern)   S(ymmetric)
                                                I(nteger)   K(kew)

  ***********************************************************************/

#define MM_MTX_STR	   "matrix"
#define MM_ARRAY_STR	   "array"
#define MM_DENSE_STR	   "array"
#define MM_COORDINATE_STR  "coordinate" 
#define MM_SPARSE_STR	   "coordinate"
#define MM_COMPLEX_STR	   "complex"
#define MM_REAL_STR	   "real"
#define MM_INT_STR	   "integer"
#define MM_GENERAL_STR     "general"
#define MM_SYMM_STR	   "symmetric"
#define MM_HERM_STR	   "hermitian"
#define MM_SKEW_STR	   "skew-symmetric"
#define MM_PATTERN_STR     "pattern"

  inline char  *mm_typecode_to_str(MM_typecode matcode) {
    char buffer[MM_MAX_LINE_LENGTH];
    const char *types[4];
    /*    int error =0; */
    /*   int i; */
    
    /* check for MTX type */
    if (mm_is_matrix(matcode)) 
      types[0] = MM_MTX_STR;
    /*
      else
      error=1;
    */
    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
      types[1] = MM_SPARSE_STR;
    else
      if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
      else
        return NULL;
    
    /* check for element data type */
    if (mm_is_real(matcode))
      types[2] = MM_REAL_STR;
    else
      if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
      else
	if (mm_is_pattern(matcode))
	  types[2] = MM_PATTERN_STR;
	else
	  if (mm_is_integer(matcode))
	    types[2] = MM_INT_STR;
	  else
	    return NULL;
    
    
    /* check for symmetry type */
    if (mm_is_general(matcode))
      types[3] = MM_GENERAL_STR;
    else if (mm_is_symmetric(matcode))
      types[3] = MM_SYMM_STR;
    else if (mm_is_hermitian(matcode))
      types[3] = MM_HERM_STR;
    else  if (mm_is_skew(matcode))
      types[3] = MM_SKEW_STR;
    else
      return NULL;
    
    sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
    return strdup(buffer);
    
  }
  
  inline int mm_read_banner(FILE *f, MM_typecode *matcode) {
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH]; 
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;
    /*    int ret_code; */
    
    mm_clear_typecode(matcode);  
    
    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL) 
      return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, 
	       storage_scheme) != 5)
      return MM_PREMATURE_EOF;

    for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
    for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
    for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
      return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
      return  MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
       storgae) or a dense array */


    if (strcmp(crd, MM_SPARSE_STR) == 0)
      mm_set_sparse(matcode);
    else
      if (strcmp(crd, MM_DENSE_STR) == 0)
	mm_set_dense(matcode);
      else
        return MM_UNSUPPORTED_TYPE;
    

    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
      mm_set_real(matcode);
    else
      if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
      else
	if (strcmp(data_type, MM_PATTERN_STR) == 0)
	  mm_set_pattern(matcode);
	else
	  if (strcmp(data_type, MM_INT_STR) == 0)
	    mm_set_integer(matcode);
	  else
	    return MM_UNSUPPORTED_TYPE;
    

    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
      mm_set_general(matcode);
    else
      if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
      else
	if (strcmp(storage_scheme, MM_HERM_STR) == 0)
	  mm_set_hermitian(matcode);
	else
	  if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
	    mm_set_skew(matcode);
	  else
	    return MM_UNSUPPORTED_TYPE;
        
    return 0;
  }

  inline int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz ) {
    char line[MM_MAX_LINE_LENGTH];
    /* int ret_code;*/
    int num_items_read;
    
    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;
    
    /* now continue scanning until you reach the end-of-comments */
    do {
      if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
	return MM_PREMATURE_EOF;
    } while (line[0] == '%');
    
    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3) return 0;
    else
      do { 
	num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
	if (num_items_read == EOF) return MM_PREMATURE_EOF;
      }
      while (num_items_read != 3);
    
    return 0;
  }


  inline int mm_read_mtx_crd_data(FILE *f, int, int, int nz, int I[],
				  int J[], double val[], MM_typecode matcode) {
    int i;
    if (mm_is_complex(matcode)) {
      for (i=0; i<nz; i++)
	if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
	    != 4) return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode)) {
      for (i=0; i<nz; i++) {
	if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
	    != 3) return MM_PREMATURE_EOF;
	
      }
    }
    else if (mm_is_pattern(matcode)) {
      for (i=0; i<nz; i++)
	if (fscanf(f, "%d %d", &I[i], &J[i])
	    != 2) return MM_PREMATURE_EOF;
    }
    else return MM_UNSUPPORTED_TYPE;

    return 0;
  }

  inline int mm_write_mtx_crd(const char *fname, int M, int N, int nz, int I[],
			      int J[], double val[], MM_typecode matcode) {
    FILE *f;
    int i;
    
    if (strcmp(fname, "stdout") == 0) 
      f = stdout;
    else
      if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;
    
    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    char *str = mm_typecode_to_str(matcode);
    fprintf(f, "%s\n", str);
    free(str);
    
    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);
    
    /* print values */
    if (mm_is_pattern(matcode))
      for (i=0; i<nz; i++)
	fprintf(f, "%d %d\n", I[i], J[i]);
    else
      if (mm_is_real(matcode))
        for (i=0; i<nz; i++)
	  fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
      else
	if (mm_is_complex(matcode))
	  for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2*i], 
		    val[2*i+1]);
	else {
	  if (f != stdout) fclose(f);
	  return MM_UNSUPPORTED_TYPE;
	}
    
    if (f !=stdout) fclose(f); 
    return 0;
  }
  

  /** matrix input/output for MatrixMarket storage */
  class MatrixMarket_IO {
    FILE *f;
    bool isComplex, isSymmetric, isHermitian;
    int row, col, nz;
    MM_typecode matcode;
  public:
    MatrixMarket_IO() : f(0) {}
    MatrixMarket_IO(const char *filename) : f(0) { open(filename); }
    template <typename Matrix> MatrixMarket_IO(const char *filename) : f(0) { open(filename); }
    ~MatrixMarket_IO() { if (f) fclose(f); f = 0; }

    int nrows() const { return row; }
    int ncols() const { return col; }
    int nnz() const { return nz; }
    int is_complex() const { return isComplex; }
    int is_symmetric() const { return isSymmetric; }
    int is_hermitian() const { return isHermitian; }

    /* open filename and reads header */
    void open(const char *filename);
    /* read opened file */
    template <typename Matrix> void read(Matrix &A);
    /* write a matrix */
    template <typename T, int shift> static void 
    write(const char *filename, const csc_matrix<T, shift>& A);  
    template <typename MAT> static void 
    write(const char *filename, const MAT& A);  
  };

  /** load a matrix-market file */
  template <typename Matrix> inline void
  MatrixMarket_load(const char *filename, Matrix& A) {
    MatrixMarket_IO mm; mm.open(filename);
    mm.read(A);
  }
  /** write a matrix-market file */
  template <typename T, int shift> void
  MatrixMarket_save(const char *filename, const csc_matrix<T, shift>& A) {
    MatrixMarket_IO mm; mm.write(filename, A);
  }

  inline void MatrixMarket_IO::open(const char *filename) {
    if (f) { fclose(f); }
    f = fopen(filename, "r"); if (!f) DAL_THROW(failure_error, "Sorry, we can not open " << filename);
    if (mm_read_banner(f, &matcode) != 0) {
      DAL_THROW(failure_error,
		"Sorry, we cannnot find the matrix market banner in " << filename);
    }
    if (mm_is_coordinate(matcode) == 0 || mm_is_matrix(matcode) == 0) {
      DAL_THROW(failure_error,
		"file is not coordinate storage or is not a matrix");
    }
    if (mm_is_pattern(matcode)) {
      DAL_THROW(failure_error, "the file does only contain the pattern of a sparse matrix");
    }
    if (mm_is_skew(matcode)) {
      DAL_THROW(failure_error, "not currently supporting skew symmetric");
    }
    isSymmetric = mm_is_symmetric(matcode) || mm_is_hermitian(matcode); 
    isHermitian = mm_is_hermitian(matcode); 
    isComplex =   mm_is_complex(matcode);
    mm_read_mtx_crd_size(f, &row, &col, &nz);
  }

  template <typename Matrix> void MatrixMarket_IO::read(Matrix &A) {
    if (!f) DAL_THROW(dal::failure_error, "no file opened!");
    typedef typename linalg_traits<Matrix>::value_type T;
    
    if (is_complex_double__(T()) && !isComplex)
      DAL_THROW(dal::failure_error, "Bad MM matrix format (complex matrix expected)");
    if (!is_complex_double__(T()) && isComplex)
      DAL_THROW(dal::failure_error, "Bad MM matrix format (real matrix expected)");
    
    A = Matrix(row, col);
    gmm::clear(A);
    
    std::vector<int> I(nz), J(nz);
    std::vector<typename Matrix::value_type> PR(nz);
    mm_read_mtx_crd_data(f, row, col, nz, &I[0], &J[0], (double*)&PR[0], matcode);
    
    for (size_type i = 0; i < size_type(nz); ++i) A(I[i]-1, J[i]-1) = PR[i];
  }

  template <typename T, int shift> void
  MatrixMarket_IO::write(const char *filename, const csc_matrix<T, shift>& A) {
    static MM_typecode t1 = {'M', 'C', 'R', 'G'};
    static MM_typecode t2 = {'M', 'C', 'C', 'G'};
    MM_typecode t;
    
    if (is_complex_double__(T())) std::copy(&(t2[0]), &(t2[0])+4, &(t[0]));
    else std::copy(&(t1[0]), &(t1[0])+4, &(t[0]));
    size_type nz = A.jc[mat_ncols(A)];
    std::vector<int> I(nz), J(nz);
    for (size_type j=0; j < mat_ncols(A); ++j) {      
      for (size_type i = A.jc[j]; i < A.jc[j+1]; ++i) {
	I[i] = A.ir[i] + 1 - shift;
	J[i] = j + 1;
      }
    }
    mm_write_mtx_crd(filename, mat_nrows(A), mat_ncols(A),
		     nz, &I[0], &J[0], (double *)A.pr, t);
  }
  template <typename MAT> void
  MatrixMarket_IO::write(const char *filename, const MAT& A) {
    gmm::csc_matrix<typename gmm::linalg_traits<MAT>::value_type> 
      tmp(gmm::mat_nrows(A), gmm::mat_ncols(A));
    gmm::copy(A,tmp); 
    MatrixMarket_IO::write(filename, tmp);
  }
  
}


#endif //  GMM_INOUTPUT_H
