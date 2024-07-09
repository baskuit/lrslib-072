/* lrslib.h (vertex enumeration using lexicographic reverse search) */

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */
/*Ver 6.1   major change is new lrsnash driver and library coded by Terje
 * Lensberg */
/*Ver 6.0   major change is mplrs wrapper for multithreading coded by Skip
 * Jordan  */
/*Ver 5.0   major change is plrs wrapper for multithreading coded by Gary
 * Roumanis */
/*Ver 4.2*  library version                                      */
/******************************************************************************/
/*  See http://cgm.cs.mcgill.ca/~avis/C/lrs.html for usage instructions */
/******************************************************************************/

#ifdef LRSLONG
#define ARITH "lrslong.h" /* lrs long integer arithmetic package */
#else
#if defined(GMP)
#define ARITH                                                                  \
  "lrsgmp.h" /* lrs wrapper for gmp multiple precsion arithmetic    */
#else
#define ARITH "lrsmp.h" /* lrs multiple precsion arithmetic    */
#endif
#endif

#include ARITH

/*********************/
/*global constants   */
/*********************/
#define MAX_LRS_GLOBALS 10000L /* number of allocated dictionaries */
#define MAXIMIZE 1L            /* maximize the lp  */
#define MINIMIZE 0L            /* maximize the lp  */
#define GE 1L                  /* constraint is >= */
#define EQ 0L                  /* constraint is linearity */

/*************/
/* typedefs  */
/*************/

/******************************************************************************/
/*                   Indexing after initialization                            */
/*               Basis                                    Cobasis             */
/*   ---------------------------------------    ----------------------------- */
/*  |  i  |0|1| .... |lastdv|lastdv+1|...|m|   | j  | 0 | 1 | ... |d-1|  d  | */
/*  |-----|+|+|++++++|++++++|--------|---|-|   |----|---|---|-----|---|+++++| */
/*  |B[i] |0|1| .... |lastdv|lastdv+1|...|m|   |C[j]|m+1|m+2| ... |m+d|m+d+1| */
/*   -----|+|+|++++++|++++++|????????|???|?|    ----|???|???|-----|???|+++++| */
/*                                                                            */
/* Row[i] is row location for B[i]         Col[j] is column location for C[j] */
/*  -----------------------------              -----------------------------  */
/* |   i   |0|1| ..........|m-1|m|            | j    | 0 | 1 | ... |d-1| d  | */
/* |-------|+|-|-----------|---|-|            |------|---|---|---  |---|++++| */
/* |Row[i] |0|1|...........|m-1|m|            |Col[j]| 1 | 2 | ... | d |  0 | */
/* --------|+|*|***********|***|*|             ------|***|***|*****|***|++++| */
/*                                                                            */
/*  + = remains invariant   * = indices may be permuted ? = swapped by pivot  */
/*                                                                            */
/*  m = number of input rows   n= number of input columns                     */
/*  input dimension inputd = n-1 (H-rep) or n (V-rep)                         */
/*  lastdv = inputd-nredundcol  (each redundant column removes a dec. var)    */
/*  working dimension d=lastdv-nlinearity (an input linearity removes a slack)
 */
/*  obj function in row 0, index 0=B[0]  col 0 has index m+d+1=C[d]           */
/*  H-rep: b-vector in col 0, A matrix in columns 1..n-1                      */
/*  V-rep: col 0 all zero, b-vector in col 1, A matrix in columns 1..n        */
/******************************************************************************/

struct lrs_dic /* dynamic dictionary data */
{
  lrs_mp_matrix A;
  long m;        /* A has m+1 rows, row 0 is cost row            */
  long m_A;      /* =m or m-d if nonnegative flag set            */
  long d;        /* A has d+1 columns, col 0 is b-vector         */
  long d_orig;   /* value of d as A was allocated  (E.G.)        */
  long lexflag;  /* true if lexmin basis for this vertex         */
  long depth;    /* depth of basis/vertex in reverse search tree */
  long i, j;     /* last pivot row and column pivot indices      */
  lrs_mp det;    /* current determinant of basis                 */
  lrs_mp objnum; /* objective numerator value                    */
  lrs_mp objden; /* objective denominator value                  */
  long *B, *Row; /* basis, row location indices                  */
  long *C, *Col; /* cobasis, column location indices             */
  lrs_dic *prev;
  lrs_dic *next;
};

struct lrs_dat /* global problem data   */
{
  long *redineq;        /* holds indices of rows for redundancy check  */
  lrs_mp_matrix Ain;    /* used only in redund to hold input matrix    */
  lrs_mp_vector Gcd;    /* Gcd of each row of numerators               */
  lrs_mp_vector Lcm;    /* Lcm for each row of input denominators      */
  lrs_mp_vector output; /* One line of output dimensioned to n         */

  lrs_mp sumdet;  /* sum of determinants                          */
  lrs_mp Nvolume; /* volume numerator                             */
  lrs_mp Dvolume; /* volume denominator                           */
  lrs_mp boundn;  /* objective bound numerator                    */
  lrs_mp boundd;  /* objective bound denominator                  */
  long unbounded; /* lp unbounded */

  /* initially holds order used to find starting  */
  /* basis, default: m,m-1,...,2,1                */
  long *facet;         /* cobasic indices for restart if needed        */
  long *redundcol;     /* holds columns which are redundant            */
  long *inequality;    /* indices of inequalities corr. to cobasic ind */
  long *linearity;     /* holds cobasic indices of input linearities   */
  long *vars;          /* ordered list of vars for extract,project etc */
  long *minratio;      /* used for lexicographic ratio test            */
  long *temparray;     /* for sorting indices, dimensioned to d        */
  long *isave, *jsave; /* arrays for estimator, malloc'ed at start     */
  long inputd;         /* input dimension: n-1 for H-rep, n for V-rep  */

  long m;         /* number of rows in input file                 */
  long n;         /* number of columns in input file              */
  long lastdv;    /* index of last dec. variable after preproc    */
                  /* given by inputd-nredundcol                   */
  long count[10]; /* count[0]=rays(facets)[1]=vertices(linearities)*/
                  /*  [2]=cobases [3]=pivots [4]=integer vertices  */
                  /*  [5-7]=mplrs R [8]=max vertex/facet depth     */

  long deepest;    /* max depth ever reached in search             */
  long nredundcol; /* number of redundant columns                  */
  long nlinearity; /* number of input linearities                  */
  long totalnodes; /* count total number of tree nodes evaluated   */
  long runs;       /* probes for estimate function                 */
  long seed;       /* seed for random number generator             */
  double cest[10]; /* ests: 0=rays,1=vert,2=bases,3=vol,4=int vert */
  long nextineq;   /* start checking redundancy from this row:def=1*/

  /**** flags  **********                         */
  long dualdeg;       /* TRUE if start dictionary is dual degenerate  */
  long homogeneous;   /* TRUE if all entries in column one are zero   */
  long lponly;        /* true if only lp solution wanted              */
  long long maxdepth; /* max depth to search to in treee              */
  long maximize;      /* flag for LP maximization                     */
  long minimize;   /* flag for LP minimization                     */
  long long mindepth; /* do not backtrack above mindepth              */
  long polytope;      /* TRUE for facet computation of a polytope     */
  long redund;        /* TRUE for computing nash equilibria           */

  /* Variables for saving/restoring cobasis,  db */

  long id;    /* numbered sequentially */
  char *name; /* passed by user */

  /* Variables for cacheing dictionaries, db */
  lrs_dic *Qhead, *Qtail, *olddic;
};

/*******************************/
/* functions  for external use */
/*******************************/

lrs_dat *
lrs_alloc_dat(const char *name); /* allocate for lrs_dat structure "name" */
lrs_dic *
lrs_alloc_dic(lrs_dat *Q); /* allocate for lrs_dic structure corr. to Q   */
long lrs_estimate(lrs_dic *P, lrs_dat *Q); /* get estimates only and returns est
                                              number of cobases in subtree */
long lrs_checkbound(
    lrs_dic *P,
    lrs_dat *Q); /* TRUE if current objective value exceeds specified bound */
long lrs_getfirstbasis(
    lrs_dic **P_p, lrs_dat *Q, lrs_mp_matrix *Lin,
    long no_output); /* gets first basis, FALSE if none,P may get changed if
                        lin. space Lin found  no_output is TRUE supresses output
                        headers P may get changed if lin. space Lin found    */
long lrs_getnextbasis(lrs_dic **dict_p, lrs_dat *Q,
                      long prune); /* gets next lrs tree basis, FALSE if none
                                      backtrack if prune is TRUE */
long lrs_getsolution(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output, long col);
long lrs_getray(lrs_dic *P, lrs_dat *Q, long col, long comment,
                lrs_mp_vector output);
long lrs_getvertex(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output);
long lrs_init();
long lrs_solvelp(
    lrs_dic *P, lrs_dat *Q,
    long maximize); /* solve primal feas LP:TRUE bounded else FALSE */

/*******************************/
/* functions  for internal use */
/*******************************/

/*******************************/
/* basic dictionary functions  */
/*******************************/
long getabasis(lrs_dic *P, lrs_dat *Q,
               long order[]); /* Try to find a starting basis  */
void getnextoutput(lrs_dic *P, lrs_dat *Q, long i, long col,
                   lrs_mp out); /* get A[B[i][col] and copy to out */
long ismin(lrs_dic *P, lrs_dat *Q, long r,
           long s); /* test if A[r][s] is a min ratio for col s */
long lexmin(lrs_dic *P, lrs_dat *Q,
            long col); /* test A to see if current basis is lexmin       */
void pivot(lrs_dic *P, lrs_dat *Q, long bas,
           long cob); /* Qpivot routine for array A  */
long primalfeasible(lrs_dic *P,
                    lrs_dat *Q); /* Do dual pivots to get primal feasibility */
long lrs_ratio(lrs_dic *P, lrs_dat *Q, long col); /* find lex min. ratio  */
long removecobasicindex(lrs_dic *P, lrs_dat *Q,
                        long k); /* remove C[k] from problem  */
long linextractcols(lrs_dic *P,
                    lrs_dat *Q); /* preprocess to remove input linearities */
long restartpivots(lrs_dic *P,
                   lrs_dat *Q); /* restart problem from given cobasis   */
long reverse(lrs_dic *P, lrs_dat *Q, long *r,
             long s); /* TRUE if B[*r] C[s] is a reverse lex-pos pivot  */
long selectpivot(lrs_dic *P, lrs_dat *Q, long *r,
                 long *s); /* select pivot indices using lexicographic rule  */
long dan_selectpivot(lrs_dic *P, lrs_dat *Q, long *r,
                     long *s); /* select pivot indices using dantzig-lex rule */
long ran_selectpivot(
    lrs_dic *P, lrs_dat *Q, long *r,
    long *s); /* select pivot indices using randomedge-lex rule    */
void update(lrs_dic *P, lrs_dat *Q, long *i,
            long *j); /* update the B,C, LOC arrays after a pivot       */
void updatevolume(lrs_dic *P,
                  lrs_dat *Q); /* rescale determinant and update the volume */

/*******************************/
/* other functions using P,Q   */
/*******************************/
long lrs_degenerate(
    lrs_dic *P, lrs_dat *Q); /* TRUE if the dictionary is primal degenerate */
long extractcols(lrs_dic *P,
                 lrs_dat *Q); /* preprocess to just extract given cols */
void rescaledet(lrs_dic *P, lrs_dat *Q, lrs_mp Vnum,
                lrs_mp Vden); /* rescale determinant to get its volume */
void rescalevolume(lrs_dic *P, lrs_dat *Q, lrs_mp Vnum,
                   lrs_mp Vden); /* adjust volume for dimension          */
long lrs_leaf(
    lrs_dic *P,
    lrs_dat *Q); /* true if current dictionary is leaf of reverse search tree */

/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/
long checkredund(lrs_dic *P,
                 lrs_dat *Q); /* solve primal lp to check redund of obj fun.
                                 returns TRUE if redundant, else FALSE */
long checkcobasic(
    lrs_dic *P, lrs_dat *Q,
    long index); /* TRUE if index is cobasic and nondegenerate  FALSE if basic,
                    or degen. cobasic, where it will get pivoted out  */
long checkindex(lrs_dic *P, lrs_dat *Q,
                long index); /* index=0 non-red.,1 red., 2 input linearity NOTE:
                                row is returned all zero if redundant!!  */

/***************************************************/
/* Routines for caching and restoring dictionaries */
/***************************************************/
void lrs_free_dic(lrs_dic *P, lrs_dat *Q);
void lrs_free_dat(lrs_dat *Q);
void copy_dict(lrs_dat *global, lrs_dic *dest, lrs_dic *src);
lrs_dic *alloc_memory(lrs_dat *Q);
lrs_dic *lrs_getdic(lrs_dat *Q);
lrs_dic *resize(lrs_dic *P, lrs_dat *Q);

/*******************************/
/* utilities                   */
/*******************************/
void reorder(long a[], long range); /* reorder array in increasing order with
                                       one misplaced element   */
void reorder1(
    long a[], long b[], long newone,
    long range); /* reorder array a in increasing order with misplaced element
                    newone elements of b go along for the ride */

/***************************/
/* lp_solve like functions */
/***************************/
void lrs_set_row(
    lrs_dic *P, lrs_dat *Q, long row, long num[], long den[],
    long ineq); /* load row i of dictionary from num[]/den[] ineq=GE       */
void lrs_set_row_mp(
    lrs_dic *P, lrs_dat *Q, long row, lrs_mp_vector num, lrs_mp_vector den,
    long ineq); /* same as lrs_set_row except num/den is lrs_mp type       */

/**************************************************************************/
