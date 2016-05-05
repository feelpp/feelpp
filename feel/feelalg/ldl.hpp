/* ========================================================================== */
/* === ldl.c: sparse LDL' factorization and solve package =================== */
/* ========================================================================== */

/* LDL:  a simple set of routines for sparse LDL' factorization.  These routines
 * are not terrifically fast (they do not use dense matrix kernels), but the
 * code is very short.  The purpose is to illustrate the algorithms in a very
 * concise manner, primarily for educational purposes.  Although the code is
 * very concise, this package is slightly faster than the built-in sparse
 * Cholesky factorization in MATLAB 7.0 (chol), when using the same input
 * permutation.
 *
 * The routines compute the LDL' factorization of a real sparse symmetric
 * matrix A (or PAP' if a permutation P is supplied), and solve upper
 * and lower triangular systems with the resulting L and D factors.  If A is
 * positive definite then the factorization will be accurate.  A can be
 * indefinite (with negative values on the diagonal D), but in this case no
 * guarantee of accuracy is provided, since no numeric pivoting is performed.
 *
 * The n-by-n sparse matrix A is in compressed-column form.  The nonzero values
 * in column j are stored in Ax [Ap [j] ... Ap [j+1]-1], with corresponding row
 * indices in Ai [Ap [j] ... Ap [j+1]-1].  Ap [0] = 0 is required, and thus
 * nz = Ap [n] is the number of nonzeros in A.  Ap is an int array of size n+1.
 * The int array Ai and the double array Ax are of size nz.  This data structure
 * is identical to the one used by MATLAB, except for the following
 * generalizations.  The row indices in each column of A need not be in any
 * particular order, although they must be in the range 0 to n-1.  Duplicate
 * entries can be present; any duplicates are summed.  That is, if row index i
 * appears twice in a column j, then the value of A (i,j) is the sum of the two
 * entries.  The data structure used here for the input matrix A is more
 * flexible than MATLAB's, which requires sorted columns with no duplicate
 * entries.
 *
 * Only the diagonal and upper triangular part of A (or PAP' if a permutation
 * P is provided) is accessed.  The lower triangular parts of the matrix A or
 * PAP' can be present, but they are ignored.
 *
 * The optional input permutation is provided as an array P of length n.  If
 * P [k] = j, the row and column j of A is the kth row and column of PAP'.
 * If P is present then the factorization is LDL' = PAP' or L*D*L' = A(P,P) in
 * 0-based MATLAB notation.  If P is not present (a null pointer) then no
 * permutation is performed, and the factorization is LDL' = A.
 *
 * The lower triangular matrix L is stored in the same compressed-column
 * form (an int Lp array of size n+1, an int Li array of size Lp [n], and a
 * double array Lx of the same size as Li).  It has a unit diagonal, which is
 * not stored.  The row indices in each column of L are always returned in
 * ascending order, with no duplicate entries.  This format is compatible with
 * MATLAB, except that it would be more convenient for MATLAB to include the
 * unit diagonal of L.  Doing so here would add additional complexity to the
 * code, and is thus omitted in the interest of keeping this code short and
 * readable.
 *
 * The elimination tree is held in the Parent [0..n-1] array.  It is normally
 * not required by the user, but it is required by ldl::numeric.  The diagonal
 * matrix D is held as an array D [0..n-1] of size n.
 *
 * --------------------
 * C-callable routines:
 * --------------------
 *
 *	ldl::symbolic:  Given the pattern of A, computes the Lp and Parent arrays
 *	    required by ldl::numeric.  Takes time proportional to the number of
 *	    nonzeros in L.  Computes the inverse Pinv of P if P is provided.
 *	    Also returns Lnz, the count of nonzeros in each column of L below
 *	    the diagonal (this is not required by ldl::numeric).
 *	ldl::numeric:  Given the pattern and numerical values of A, the Lp array,
 *	    the Parent array, and P and Pinv if applicable, computes the
 *	    pattern and numerical values of L and D.
 *	ldl::lsolve:  Solves Lx=b for a dense vector b.
 *	ldl::dsolve:  Solves Dx=b for a dense vector b.
 *	ldl::ltsolve: Solves L'x=b for a dense vector b.
 *	ldl::perm:    Computes x=Pb for a dense vector b.
 *	ldl::permt:   Computes x=P'b for a dense vector b.
 *	ldl::valid_perm:  checks the validity of a permutation vector
 *	ldl::valid_matrix:  checks the validity of the sparse matrix A
 *
 * ----------------------------
 * Limitations of this package:
 * ----------------------------
 *
 * In the interest of keeping this code simple and readable, ldl::symbolic and
 * ldl::numeric assume their inputs are valid.  You can check your own inputs
 * prior to calling these routines with the ldl::valid_perm and ldl::valid_matrix
 * routines.  Except for the two ldl::valid_* routines, no routine checks to see
 * if the array arguments are present (non-NULL).  Like all C routines, no
 * routine can determine if the arrays are long enough and don't overlap.
 *
 * The ldl::numeric does check the numerical factorization, however.  It returns
 * n if the factorization is successful.  If D (k,k) is zero, then k is
 * returned, and L is only partially computed.
 *
 * No pivoting to control fill-in is performed, which is often critical for
 * obtaining good performance.  I recommend that you compute the permutation P
 * using AMD or SYMAMD (approximate minimum degree ordering routines), or an
 * appropriate graph-partitioning based ordering.  See the ldldemo.m routine for
 * an example in MATLAB, and the ldlmain.c stand-alone C program for examples of
 * how to find P.  Routines for manipulating compressed-column matrices are
 * available in UMFPACK.  AMD, SYMAMD, UMFPACK, and this LDL package are all
 * available at http://www.cise.ufl.edu/research/sparse.
 *
 * -------------------------
 * Possible simplifications:
 * -------------------------
 *
 * These routines could be made even simpler with a few additional assumptions.
 * If no input permutation were performed, the caller would have to permute the
 * matrix first, but the computation of Pinv, and the use of P and Pinv could be
 * removed.  If only the diagonal and upper triangular part of A or PAP' are
 * present, then the tests in the "if (i < k)" statement in ldl::symbolic and
 * "if (i <= k)" in ldl::numeric, are always true, and could be removed (i can
 * equal k in ldl::symbolic, but then the body of the if statement would
 * correctly do no work since Flag [k] == k).  If we could assume that no
 * duplicate entries are present, then the statement Y [i] += Ax [p] could be
 * replaced with Y [i] = Ax [p] in ldl::numeric.
 *
 * --------------------------
 * Description of the method:
 * --------------------------
 *
 * LDL computes the symbolic factorization by finding the pattern of L one row
 * at a time.  It does this based on the following theory.  Consider a sparse
 * system Lx=b, where L, x, and b, are all sparse, and where L comes from a
 * Cholesky (or LDL') factorization.  The elimination tree (etree) of L is
 * defined as follows.  The parent of node j is the smallest k > j such that
 * L (k,j) is nonzero.  Node j has no parent if column j of L is completely zero
 * below the diagonal (j is a root of the etree in this case).  The nonzero
 * pattern of x is the union of the paths from each node i to the root, for
 * each nonzero b (i).  To compute the numerical solution to Lx=b, we can
 * traverse the columns of L corresponding to nonzero values of x.  This
 * traversal does not need to be done in the order 0 to n-1.  It can be done in
 * any "topological" order, such that x (i) is computed before x (j) if i is a
 * descendant of j in the elimination tree.
 *
 * The row-form of the LDL' factorization is shown in the MATLAB function
 * ldlrow.m in this LDL package.  Note that row k of L is found via a sparse
 * triangular solve of L (1:k-1, 1:k-1) \ A (1:k-1, k), to use 1-based MATLAB
 * notation.  Thus, we can start with the nonzero pattern of the kth column of
 * A (above the diagonal), follow the paths up to the root of the etree of the
 * (k-1)-by-(k-1) leading submatrix of L, and obtain the pattern of the kth row
 * of L.  Note that we only need the leading (k-1)-by-(k-1) submatrix of L to
 * do this.  The elimination tree can be constructed as we go.
 *
 * The symbolic factorization does the same thing, except that it discards the
 * pattern of L as it is computed.  It simply counts the number of nonzeros in
 * each column of L and then constructs the Lp index array when it's done.  The
 * symbolic factorization does not need to do this in topological order.
 * Compare ldl::symbolic with the first part of ldl::numeric, and note that the
 * while (len > 0) loop is not present in ldl::symbolic.
 *
 * LDL Version 1.1 (Apr. 22, 2005), Copyright (c) 2005 by Timothy A Davis,
 * University of Florida.  All Rights Reserved.  Developed while on sabbatical
 * at Stanford University and Lawrence Berkeley National Laboratory.  Refer to
 * the README file for the License.  Available at
 * http://www.cise.ufl.edu/research/sparse.
 */
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/traits.hpp>

namespace Feel
{
namespace glas
{
namespace ublas = boost::numeric::ublas;

template <typename T = double>
class ldl
{
  public:
    typedef T value_type;
    typedef ublas::compressed_matrix<value_type, ublas::column_major> matrix_type;

    ldl( matrix_type& A )
        : n( A.size1() ),
          Ap( A.index1_data().begin() ),
          Ai( A.index2_data().begin() ),
          Ax( A.value_data().begin() ),
          Lp( n + 1 ),
          Parent( n ),
          Lnz( n ),
          Li(),
          Lx(),
          D( n ),
          Y( n ),
          Flag( n ),
          Pattern( n ),
          P( 0 ),
          Pinv( 0 )
    {
    }
    void symbolic();

    size_type numeric();

    void lsolve( ublas::vector<value_type>& X );

    void dsolve( ublas::vector<value_type>& X );

    void ltsolve( ublas::vector<value_type>& X );

    /**
     * shortcut for solving LDL^T X = B
     */
    void solve( ublas::vector<value_type>& X )
    {
        lsolve( X );
        dsolve( X );
        ltsolve( X );
    }

    void perm( ublas::vector<value_type>& X, ublas::vector<value_type> const& B );
    void permt( ublas::vector<value_type>& X, ublas::vector<value_type> const& B );

    bool valid_perm();
    bool valid_matrix();

  private:
    size_type n;

    size_type* Ap;
    size_type* Ai;
    value_type* Ax;

    ublas::vector<size_type> Lp;

    ublas::vector<int64_type> Parent;

    /* output of size n, not defined on input */
    ublas::vector<size_type> Lnz;
    ublas::vector<size_type> Li;  /* output of size lnz=Lp[n], not defined on input */
    ublas::vector<value_type> Lx; /* output of size lnz=Lp[n], not defined on input */
    ublas::vector<value_type> D;  /* output of size n, not defined on input */
    ublas::vector<value_type> Y;  /* workspace of size n, not defn. on input or output */

    /* workspace of size n, not defn. on input or output */
    ublas::vector<size_type> Flag;
    ublas::vector<size_type> Pattern;

    /* optional input of size n */
    size_type* P;
    /* optional output of size n (used if P is not NULL) */
    size_type* Pinv;
};

/* ========================================================================== */
/* === ldl::symbolic ========================================================= */
/* ========================================================================== */

/* The input to this routine is a sparse matrix A, stored in column form, and
 * an optional permutation P.  The output is the elimination tree
 * and the number of nonzeros in each column of L.  Parent [i] = k if k is the
 * parent of i in the tree.  The Parent array is required by ldl::numeric.
 * Lnz [k] gives the number of nonzeros in the kth column of L, excluding the
 * diagonal.
 *
 * One workspace vector (Flag) of size n is required.
 *
 * If P is NULL, then it is ignored.  The factorization will be LDL' = A.
 * Pinv is not computed.  In this case, neither P nor Pinv are required by
 * ldl::numeric.
 *
 * If P is not NULL, then it is assumed to be a valid permutation.  If
 * row and column j of A is the kth pivot, the P [k] = j.  The factorization
 * will be LDL' = PAP', or A (p,p) in MATLAB notation.  The inverse permutation
 * Pinv is computed, where Pinv [j] = k if P [k] = j.  In this case, both P
 * and Pinv are required as inputs to ldl::numeric.
 *
 * The floating-point operation count of the subsequent call to ldl::numeric
 * is not returned, but could be computed after ldl::symbolic is done.  It is
 * the sum of (Lnz [k]) * (Lnz [k] + 2) for k = 0 to n-1.
 */
template <typename T>
void ldl<T>::symbolic()
{
    size_type i, k, p, kk, p2;

    if ( P )
    {
        /* If P is present then compute Pinv, the inverse of P */
        for ( k = 0; k < n; k++ )
        {
            Pinv[P[k]] = k;
        }
    }

    for ( k = 0; k < n; k++ )
    {
        /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
        Parent[k] = -1;                /* parent of k is not yet known */
        Flag[k] = k;                   /* mark node k as visited */
        Lnz[k] = 0;                    /* count of nonzeros in column k of L */
        kk = ( P ) ? ( P[k] ) : ( k ); /* kth original, or permuted, column */
        p2 = Ap[kk + 1];

        for ( p = Ap[kk]; p < p2; p++ )
        {
            /* A (i,k) is nonzero (original or permuted A) */
            i = ( Pinv ) ? ( Pinv[Ai[p]] ) : ( Ai[p] );

            if ( i < k )
            {
                /* follow path from i to root of etree, stop at flagged node */
                for ( ; Flag[i] != k; i = Parent[i] )
                {
                    /* find parent of i if not yet determined */
                    if ( Parent[i] == -1 ) Parent[i] = k;

                    Lnz[i]++;    /* L (k,i) is nonzero */
                    Flag[i] = k; /* mark i as visited */
                }
            }
        }
    }

    /* construct Lp index array from Lnz column counts */
    Lp[0] = 0;

    for ( k = 0; k < n; k++ )
    {
        Lp[k + 1] = Lp[k] + Lnz[k];
    }

    Li.resize( Lp[n] );
    Lx.resize( Lp[n] );
}

/* ========================================================================== */
/* === ldl::numeric ========================================================== */
/* ========================================================================== */

/* Given a sparse matrix A (the arguments n, Ap, Ai, and Ax) and its symbolic
 * analysis (Lp and Parent, and optionally P and Pinv), compute the numeric LDL'
 * factorization of A or PAP'.  The outputs of this routine are arguments Li,
 * Lx, and D.  It also requires three size-n workspaces (Y, Pattern, and Flag).
 */

/* returns n if successful, k if D (k,k) is zero */
template <typename T>
size_type ldl<T>::numeric()
{
    value_type yi, l_ki;
    size_type i, k, p, kk, p2;
    size_type len, top;

    for ( k = 0; k < n; k++ )
    {
        /* compute nonzero Pattern of kth row of L, in topological order */
        Y[k] = 0.0;                    /* Y(0:k) is now all zero */
        top = n;                       /* stack for pattern is empty */
        Flag[k] = k;                   /* mark node k as visited */
        Lnz[k] = 0;                    /* count of nonzeros in column k of L */
        kk = ( P ) ? ( P[k] ) : ( k ); /* kth original, or permuted, column */
        p2 = Ap[kk + 1];

        for ( p = Ap[kk]; p < p2; p++ )
        {
            i = ( Pinv ) ? ( Pinv[Ai[p]] ) : ( Ai[p] ); /* get A(i,k) */

            if ( i <= k )
            {
                Y[i] += Ax[p]; /* scatter A(i,k) into Y (sum duplicates) */

                for ( len = 0; Flag[i] != k; i = Parent[i] )
                {
                    Pattern[len++] = i; /* L(k,i) is nonzero */
                    Flag[i] = k;        /* mark i as visited */
                }

                while ( len > 0 )
                    Pattern[--top] = Pattern[--len];
            }
        }

        /* compute numerical values kth row of L (a sparse triangular solve) */
        D[k] = Y[k]; /* get D(k,k) and clear Y(k) */
        Y[k] = 0.0;

        for ( ; top < n; top++ )
        {
            i = Pattern[top]; /* Pattern [top:n-1] is pattern of L(:,k) */
            yi = Y[i];        /* get and clear Y(i) */
            Y[i] = 0.0;
            p2 = Lp[i] + Lnz[i];

            for ( p = Lp[i]; p < p2; p++ )
            {
                Y[Li[p]] -= Lx[p] * yi;
            }

            l_ki = yi / D[i]; /* the nonzero entry L(k,i) */
            D[k] -= l_ki * yi;
            Li[p] = k; /* store L(k,i) in column form of L */
            Lx[p] = l_ki;
            Lnz[i]++; /* increment count of nonzeros in col i */
        }

        if ( D[k] == 0.0 ) return ( k ); /* failure, D(k,k) is zero */
    }

    return ( n ); /* success, diagonal of D is all nonzero */
}

/* ========================================================================== */
/* === ldl::lsolve:  solve Lx=b ============================================== */
/* ========================================================================== */
template <typename T>
void ldl<T>::lsolve( ublas::vector<value_type>& X )
{
    size_type j, p, p2;

    for ( j = 0; j < n; j++ )
    {
        p2 = Lp[j + 1];

        for ( p = Lp[j]; p < p2; p++ )
        {
            X[Li[p]] -= Lx[p] * X[j];
        }
    }
}

/* ========================================================================== */
/* === ldl::dsolve:  solve Dx=b ============================================== */
/* ========================================================================== */
template <typename T>
void ldl<T>::dsolve( ublas::vector<value_type>& X )
{
    size_type j;

    for ( j = 0; j < n; j++ )
    {
        X[j] /= D[j];
    }
}

/* ========================================================================== */
/* === ldl::ltsolve: solve L'x=b  ============================================ */
/* ========================================================================== */
template <typename T>
void ldl<T>::ltsolve( ublas::vector<value_type>& X )
{
    size_type p, p2;

    for ( int j = n - 1; j >= 0; j-- )
    {
        p2 = Lp[j + 1];

        for ( p = Lp[j]; p < p2; p++ )
        {
            X[j] -= Lx[p] * X[Li[p]];
        }
    }
}

/* ========================================================================== */
/* === ldl::perm: permute a vector, x=Pb ===================================== */
/* ========================================================================== */
template <typename T>
void ldl<T>::perm( ublas::vector<value_type>& X,
                   ublas::vector<value_type> const& B )
{
    size_type j;

    for ( j = 0; j < n; j++ )
    {
        X[j] = B[P[j]];
    }
}

/* ========================================================================== */
/* === ldl::permt: permute a vector, x=P'b =================================== */
/* ========================================================================== */
template <typename T>
void ldl<T>::permt( ublas::vector<value_type>& X,
                    ublas::vector<value_type> const& B )
{
    size_type j;

    for ( j = 0; j < n; j++ )
    {
        X[P[j]] = B[j];
    }
}

/* ========================================================================== */
/* === ldl::valid_perm: check if a permutation vector is valid =============== */
/* ========================================================================== */

/* returns 1 if valid, 0 otherwise */
template <typename T>
bool ldl<T>::valid_perm()
{
    size_type j, k;

    if ( Flag.size() == 0 )
    {
        return false; /* n must be >= 0, and Flag must be present */
    }

    if ( !P )
    {
        return true; /* If NULL, P is assumed to be the identity perm. */
    }

    for ( j = 0; j < n; j++ )
    {
        Flag[j] = 0; /* clear the Flag array */
    }

    for ( k = 0; k < n; k++ )
    {
        j = P[k];

        if ( j >= n || Flag[j] != 0 )
        {
            return false; /* P is not valid */
        }

        Flag[j] = 1;
    }

    return true; /* P is valid */
}

/* ========================================================================== */
/* === ldl::valid_matrix: check if a sparse matrix is valid ================== */
/* ========================================================================== */

/* This routine checks to see if a sparse matrix A is valid for input to
 * ldl::symbolic and ldl::numeric.  It returns 1 if the matrix is valid, 0
 * otherwise.  A is in sparse column form.  The numerical values in column j
 * are stored in Ax [Ap [j] ... Ap [j+1]-1], with row indices in
 * Ai [Ap [j] ... Ap [j+1]-1].  The Ax array is not checked.
 */
template <typename T>
bool ldl<T>::valid_matrix()
{
    size_type j, p;

    /*
    if (n < 0 || !Ap || !Ai || Ap [0] != 0)
    {
    return false ;	    // n must be >= 0, and Ap and Ai must be present
    }*/
    for ( j = 0; j < n; j++ )
    {
        if ( Ap[j] > Ap[j + 1] )
        {
            return false; /* Ap must be monotonically nondecreasing */
        }
    }

    for ( p = 0; p < Ap[n]; p++ )
    {
        if ( Ai[p] >= n )
        {
            return false; /* row indices must be in the range 0 to n-1 */
        }
    }

    return true; /* matrix is valid */
}

} // namespace glas
} // namespace Feel
