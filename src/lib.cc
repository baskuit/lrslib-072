/*******************************************************/
/* lrsnashlib is a library of routines for computing   */
/* computing all nash equilibria for two person games  */
/* given by mxn payoff matrices A,B                    */
/*                                                     */
/*                                                     */
/* Main user callable function is                      */
/*         lrs_solve_nash(game *g)                     */
/*                                                     */
/* Requires lrsnashlib.h lrslib.h lrslib.c             */
/*                                                     */
/* Sample driver: lrsnash.c                            */
/* Derived from nash.c in lrslib-060                   */
/* by Terje Lensberg, October 26, 2015:                */
/*******************************************************/

#include "lrslib.h"
// this must be last
#include "lib.h"

// hack
long nash2_main(lrs_dic *P1, lrs_dat *Q1, lrs_dic *P2orig, lrs_dat *Q2,
                long *numequilib, lrs_mp_vector output, FloatOneSumOutput *gg,
                long linindex[]);
/* lrs driver, argv[2]= 2nd input file for nash equilibria */

long lrs_getfirstbasis2(lrs_dic **D_p, lrs_dat *Q, lrs_dic *P2orig,
                        lrs_mp_matrix *Lin, long no_output, long linindex[]);

long getabasis2(lrs_dic *P, lrs_dat *Q, lrs_dic *P2orig, long order[],
                long linindex[]);

long lrs_nashoutput(lrs_dat *Q, lrs_mp_vector output, FloatOneSumOutput *gg,
                    long player);
/* returns TRUE and prints output if not the origin */

void BuildRepP1Is0(lrs_dic *P, lrs_dat *Q, const FastInput *g);
void BuildRepP1Is1(lrs_dic *P, lrs_dat *Q, const FastInput *g);
void FillConstraintRowsP1Is0(lrs_dic *P, lrs_dat *Q, const FastInput *g,
                             int firstRow);
void FillConstraintRowsP1Is1(lrs_dic *P, lrs_dat *Q, const FastInput *g,
                             int firstRow);
void FillFirstRow(lrs_dic *P, lrs_dat *Q, int n);
void FillLinearityRow(lrs_dic *P, lrs_dat *Q, int m, int n);
void FillNonnegativityRows(lrs_dic *P, lrs_dat *Q, int firstRow, int lastRow,
                           int n);

#define MAXSTRAT 200
#define ROW 0
#define COL 1

#include <stdio.h>
#include <string.h>

//========================================================================
// Standard solver. Modified version of main() from lrsNash
//========================================================================
static long FirstTime; /* set this to true for every new game to be solved */

int solve_fast(const FastInput *g, FloatOneSumOutput *gg) {

  lrs_init_no_header();

  lrs_dic *P1;      /* structure for holding current dictionary and indices */
  lrs_dat *Q1, *Q2; /* structure for holding static problem data            */

  lrs_mp_vector
      output1; /* holds one line of output; ray,vertex,facet,linearity */
  lrs_mp_vector
      output2;       /* holds one line of output; ray,vertex,facet,linearity */
  lrs_mp_matrix Lin; /* holds input linearities if any are found             */
  lrs_mp_matrix A2orig;
  lrs_dic *P2orig; /* we will save player 2's dictionary in getabasis      */

  long *linindex; /* for faster restart of player 2                       */

  long col; /* output column index for dictionary                   */
  long startcol = 0;
  long prune = FALSE; /* if TRUE, getnextbasis will prune tree and backtrack  */
  long numequilib = 0; /* number of nash equilibria found */
  long oldnum = 0;

  /*********************************************************************************/
  /* Step 1: Allocate lrs_dat, lrs_dic and set up the problem */
  /*********************************************************************************/
  FirstTime = TRUE; /* This is done for each new game */

  Q1 = lrs_alloc_dat(
      "LRS globals"); /* allocate and init structure for static problem data */
  if (Q1 == NULL) {
    return 0;
  }

  Q1->n = g->rows + 2;
  Q1->m = g->rows + g->cols + 1;

  P1 = lrs_alloc_dic(Q1); /* allocate and initialize lrs_dic */
  if (P1 == NULL) {
    return 0;
  }

  BuildRepP1Is1(P1, Q1, g);

  output1 = lrs_alloc_mp_vector(
      Q1->n + Q1->m); /* output holds one line of output from dictionary     */

  /* allocate and init structure for player 2's problem data */
  Q2 = lrs_alloc_dat("LRS globals");
  if (Q2 == NULL) {
    return 0;
  }

  Q2->n = g->cols + 2;
  Q2->m = g->rows + g->cols + 1;

  P2orig = lrs_alloc_dic(Q2); /* allocate and initialize lrs_dic */
  if (P2orig == NULL) {
    return 0;
  }
  BuildRepP1Is0(P2orig, Q2, g);
  A2orig = P2orig->A;

  output2 = lrs_alloc_mp_vector(
      Q1->n + Q1->m); /* output holds one line of output from dictionary     */

  linindex = static_cast<long *>(
      calloc((P2orig->m + P2orig->d + 2), sizeof(long))); /* for next time */

  /*********************************************************************************/
  /* Step 2: Find a starting cobasis from default of specified order */
  /*         P1 is created to hold  active dictionary data and may be cached */
  /*         Lin is created if necessary to hold linearity space */
  /*         Print linearity space if any, and retrieve output from first dict.
   */
  /*********************************************************************************/

  if (!lrs_getfirstbasis(&P1, Q1, &Lin, TRUE))
    return 1;

  /* Pivot to a starting dictionary                      */
  /* There may have been column redundancy               */
  /* If so the linearity space is obtained and redundant */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */

  if (Q1->homogeneous && Q1->hull)
    startcol++; /* col zero not treated as redundant   */

  col = Q1->nredundcol;

  /*********************************************************************************/
  /* Step 3: Terminate if lponly option set, otherwise initiate a reverse */
  /*         search from the starting dictionary. Get output for each new dict.
   */
  /*********************************************************************************/

  /* We initiate reverse search from this dictionary       */
  /* getting new dictionaries until the search is complete */
  /* User can access each output line from output which is */
  /* vertex/ray/facet from the lrs_mp_vector output         */
  /* prune is TRUE if tree should be pruned at current node */
  do {
    prune = lrs_checkbound(P1, Q1);
    if (!prune && lrs_getsolution(P1, Q1, output1, col)) {
      oldnum = numequilib;
      nash2_main(P1, Q1, P2orig, Q2, &numequilib, output2, gg, linindex);
      if (numequilib > oldnum) {
        lrs_nashoutput(Q1, output1, gg, 1L);
        break;
      }
    }
  } while (lrs_getnextbasis(&P1, Q1, prune));

  lrs_clear_mp_vector(output1, Q1->m + Q1->n);
  lrs_clear_mp_vector(output2, Q1->m + Q1->n);

  lrs_free_dic(P1, Q1); /* deallocate lrs_dic */
  lrs_free_dat(Q1);     /* deallocate lrs_dat */

  /* 2015.10.10  new code to clear P2orig */
  Q2->Qhead = P2orig; /* reset this or you crash free_dic */
  P2orig->A = A2orig; /* reset this or you crash free_dic */

  lrs_free_dic(P2orig, Q2); /* deallocate lrs_dic */
  lrs_free_dat(Q2);         /* deallocate lrs_dat */

  free(linindex);

  gg->value /= g->den;

  return 0;
}

/*********************************************/
/* end of nash driver                        */
/*********************************************/

/**********************************************************/
/* nash2_main is a second driver used in computing nash   */
/* equilibria on a second polytope interleaved with first */
/**********************************************************/

long nash2_main(lrs_dic *P1, lrs_dat *Q1, lrs_dic *P2orig, lrs_dat *Q2,
                long *numequilib, lrs_mp_vector output, FloatOneSumOutput *gg,
                long linindex[]) {

  lrs_dic *P2;       /* This can get resized, cached etc. Loaded from P2orig */
  lrs_mp_matrix Lin; /* holds input linearities if any are found             */
  long col;          /* output column index for dictionary                   */
  long startcol = 0;
  long prune = FALSE;
  long nlinearity;
  long *linearity;
  static long firstwarning = TRUE;
  static long firstunbounded = TRUE;

  long i, j;

  P2 = lrs_getdic(Q2);
  copy_dict(Q2, P2, P2orig);

  linearity = Q2->linearity;
  nlinearity = 0;
  for (i = Q1->lastdv + 1; i <= P1->m; i++) {
    if (!zero(P1->A[P1->Row[i]][0])) {
      j = Q1->inequality[P1->B[i] - Q1->lastdv];
      if (Q1->nlinearity == 0 || j < Q1->linearity[0])
        linearity[nlinearity++] = j;
    }
  }

  if (Q1->nlinearity > 0)
    linearity[nlinearity++] = Q1->linearity[0];

  for (i = 1; i < nlinearity; i++)
    reorder(linearity, nlinearity);

  Q2->nlinearity = nlinearity;
  Q2->polytope = FALSE;

  if (!lrs_getfirstbasis2(&P2, Q2, P2orig, &Lin, TRUE, linindex))
    goto sayonara;

  if (Q2->homogeneous && Q2->hull)
    startcol++; /* col zero not treated as redundant   */

  do {
    prune = lrs_checkbound(P2, Q2);
    col = 0;
    if (!prune && lrs_getsolution(P2, Q2, output, col)) {
      if (lrs_nashoutput(Q2, output, gg, 2L)) {
        (*numequilib)++;
        break;
      }
    }
  } while (lrs_getnextbasis(&P2, Q2, prune));

sayonara:
  lrs_free_dic(P2, Q2);
  return 0;
}

/*********************************************/
/* end of nash2_main                          */
/*********************************************/

/* In lrs_getfirstbasis and lrs_getnextbasis we use D instead of P */
/* since the dictionary P may change, ie. &P in calling routine    */

#define D (*D_p)

long lrs_getfirstbasis2(lrs_dic **D_p, lrs_dat *Q, lrs_dic *P2orig,
                        lrs_mp_matrix *Lin, long no_output, long linindex[])
/* gets first basis, FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is TRUE supresses output headers   */
{
  long i, j, k;

  /* assign local variables to structures */

  lrs_mp_matrix A;
  long *B, *C, *Col;
  long *inequality;
  long *linearity;
  long hull = Q->hull;
  long m, d, nlinearity, nredundcol;

  m = D->m;
  d = D->d;

  nredundcol = 0L;            /* will be set after getabasis        */
  nlinearity = Q->nlinearity; /* may be reset if new linearity read */
  linearity = Q->linearity;

  A = D->A;
  B = D->B;
  C = D->C;
  Col = D->Col;
  inequality = Q->inequality;

  /* default is to look for starting cobasis using linearies first, then     */
  /* filling in from last rows of input as necessary                         */
  /* linearity array is assumed sorted here                                  */
  /* note if restart/given start inequality indices already in place         */
  /* from nlinearity..d-1                                                    */

  for (i = 0; i < nlinearity; i++) /* put linearities first in the order */
    inequality[i] = linearity[i];

  k = 0; /* index for linearity array   */

  if (Q->givenstart)
    k = d;
  else
    k = nlinearity;
  for (i = m; i >= 1; i--) {
    j = 0;
    while (j < k && inequality[j] != i)
      j++; /* see if i is in inequality  */
    if (j == k)
      inequality[k++] = i;
  }

  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      itomp(ZERO, A[0][j]);

  /* Now we pivot to standard form, and then find a primal feasible basis */
  /* Note these steps MUST be done, even if restarting, in order to get */
  /* the same index/inequality correspondance we had for the original prob. */
  /* The inequality array is used to give the insertion order */
  /* and is defaulted to the last d rows when givenstart=FALSE */

  if (!getabasis2(D, Q, P2orig, inequality, linindex)) {
    return FALSE;
  }

  nredundcol = Q->nredundcol;
  d = D->d;

  /* Reset up the inequality array to remember which index is which input
   * inequality */
  /* inequality[B[i]-lastdv] is row number of the inequality with index B[i] */
  /* inequality[C[i]-lastdv] is row number of the inequality with index C[i] */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0) {              /* some cobasic indices will be removed */
    for (i = 0; i < nlinearity; i++) /* remove input linearity indices */
      inequality[linearity[i]] = 0;
    k = 1; /* counter for linearities         */
    for (i = 1; i <= m - nlinearity; i++) {
      while (k <= m && inequality[k] == 0)
        k++; /* skip zeroes in corr. to linearity */
      inequality[i] = inequality[k++];
    }
  } /* end if linearity */
  if (nredundcol > 0) {
    const unsigned int Qn = Q->n;
    *Lin = lrs_alloc_mp_matrix(nredundcol, Qn);

    for (i = 0; i < nredundcol; i++) {
      if (!(Q->homogeneous && Q->hull &&
            i == 0)) { /* skip redund col 1 for homog. hull */
        lrs_getray(D, Q, Col[0], D->C[0] + i - hull,
                   (*Lin)[i]); /* adjust index for deletions */
      }

      if (!removecobasicindex(D, Q, 0L)) {
        lrs_clear_mp_matrix(*Lin, nredundcol, Qn);
        return FALSE;
      }
    }
  } /* end if nredundcol > 0 */

  /* Do dual pivots to get primal feasibility */
  if (!primalfeasible(D, Q)) {
    return FALSE;
  }

  /* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize) {
    Q->unbounded = !lrs_solvelp(D, Q, Q->maximize);

    /* check to see if objective is dual degenerate */
    j = 1;
    while (j <= d && !zero(A[0][j]))
      j++;
    if (j <= d)
      Q->dualdeg = TRUE;
  } else
  /* re-initialize cost row to -det */
  {
    for (j = 1; j <= d; j++) {
      copy(A[0][j], D->det);
      storesign(A[0][j], NEG);
    }

    itomp(ZERO, A[0][0]); /* zero optimum objective value */
  }

  /* reindex basis to 0..m if necessary */
  /* we use the fact that cobases are sorted by index value */
  while (C[0] <= m) {
    i = C[0];
    // j = inequality[B[i] - lastdv];
    // inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
    // inequality[C[0] - lastdv] = j;
    C[0] = B[i];
    B[i] = i;
    reorder1(C, Col, ZERO, d);
  }

  /* Check to see if necessary to resize */
  if (Q->inputd > D->d)
    *D_p = resize(D, Q);

  return TRUE;
}

/********* end of lrs_getfirstbasis  ***************/
long getabasis2(lrs_dic *P, lrs_dat *Q, lrs_dic *P2orig, long order[],
                long linindex[])

/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */
{
  /* 2015.10.10 linindex now preallocated and received as parameter so we can
   * free it */

  //  static long firsttime = TRUE; /* stays true until first valid dictionary
  //  built */

  long i, j, k;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *linearity = Q->linearity;
  long *redundcol = Q->redundcol;
  long m, d, nlinearity;
  long nredundcol = 0L; /* will be calculated here */

  m = P->m;
  d = P->d;
  nlinearity = Q->nlinearity;
  // 2015.9.15
  /* after first time we update the change in linearities from the last time,
   * saving many pivots */
  if (!FirstTime) {
    for (i = 1; i <= m + d; i++)
      linindex[i] = FALSE;
    for (i = 0; i < nlinearity; i++) {
      linindex[d + linearity[i]] = TRUE;
    }

    for (i = 1; i <= m; i++) {
      if (linindex[B[i]]) { /* pivot out unwanted linearities */
        k = 0;
        while (k < d && (linindex[C[k]] || zero(A[Row[i]][Col[k]])))
          k++;

        if (k < d) {
          j = i; /* note this index changes in update, cannot use i!) */

          if (C[k] > B[j]) /* decrease i or we may skip a linearity */
            i--;
          pivot(P, Q, j, k);
          update(P, Q, &j, &k);
        } else {
          /* this is not necessarily an error, eg. two identical rows/cols in
           * payoff matrix */
          if (!zero(A[Row[i]][0])) { /* error condition */
            return (FALSE);
          }
        }

      } /* if linindex */
    } /* for i   .. */
  } else { /* we have not had a successful dictionary built from the given
              linearities */

    /* standard lrs processing is done on only the first call to getabasis2 */

    for (j = 0; j < m; j++) {
      i = 0;
      while (i <= m && B[i] != d + order[j])
        i++;                         /* find leaving basis index i */
      if (j < nlinearity && i > m) { /* cannot pivot linearity to cobasis */
        return FALSE;
      }
      if (i <= m) { /* try to do a pivot */
        k = 0;
        while (C[k] <= d && zero(A[Row[i]][Col[k]]))
          k++;

        if (C[k] <= d) {
          pivot(P, Q, i, k);
          update(P, Q, &i, &k);
        } else if (j < nlinearity) { /* cannot pivot linearity to cobasis */
          if (zero(A[Row[i]][0])) {
            linearity[j] = 0;
          } else {
            return FALSE;
          }
        } /* end if j < nlinearity */

      } /* end of if i <= m .... */
    } /* end of for   */

    /* update linearity array to get rid of redundancies */
    i = 0;
    k = 0; /* counters for linearities         */
    while (k < nlinearity) {
      while (k < nlinearity && linearity[k] == 0)
        k++;
      if (k < nlinearity)
        linearity[i++] = linearity[k++];
    }

    nlinearity = i;
    /* lrs bug fix, 2009.6.27, nash 2015.9.16 */
    Q->nlinearity = i;

    /* column dependencies now can be recorded  */
    /* redundcol contains input column number 0..n-1 where redundancy is */
    k = 0;
    while (k < d && C[k] <= d) {
      if (C[k] <= d) /* decision variable still in cobasis */
        redundcol[nredundcol++] = C[k] - Q->hull; /* adjust for hull indices */
      k++;
    }

    /* now we know how many decision variables remain in problem */
    Q->nredundcol = nredundcol;
    Q->lastdv = d - nredundcol;
    /* 2015.9.15 bug fix : we needed first *successful* time */
    FirstTime = FALSE;
  } /* else firsttime ... we have built a dictionary from the given linearities
     */

  /* we continue from here after loading dictionary */

  /* here we save dictionary for use next time, *before* we resize */

  copy_dict(Q, P2orig, P);

  /* Remove linearities from cobasis for rest of computation */
  /* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++) { /* find cobasic index */
    k = 0;
    while (k < d && C[k] != linearity[i] + d)
      k++;
    if (k >= d) {
      /* not neccesarily an error as eg., could be repeated row/col in payoff */
    } else {
      removecobasicindex(P, Q, k);
      d = P->d;
    }
  }
  /* set index value for first slack variable */

  /* Check feasability */
  if (Q->givenstart) {
    i = Q->lastdv + 1;
    while (i <= m && !negative(A[Row[i]][0]))
      i++;
  }
  return TRUE;
} /*  end of getabasis2 */

long lrs_nashoutput(lrs_dat *Q, lrs_mp_vector output, FloatOneSumOutput *gg,
                    long player) {
  long i;
  long origin = TRUE;

  for (i = 1; i < Q->n; i++)
    if (!zero(output[i]))
      origin = FALSE;

  if (origin)
    return FALSE;

  const float den = static_cast<float>(mptoi(output[0]));
  // printf("den conversion: %ld = %f\n", mptoi(output[0]), den);

  if (player == 1) {
    // printf("p1 output dump:\n");
    // for (int j = 0; j < Q->n; ++j) {
    //   printf("%ld, ", mptoi(output[j]));
    // }
    // printf("\n");

    for (i = 1; i < Q->n; i++) {
      gg->row_strategy[i - 1] = mptoi(output[i]) / den;
    }
  } else {
    gg->value = mptoi(output[Q->n - 1]) / den;
    // printf("value set: %ld / %f = %f\n", mptoi(output[Q->n - 1]), den,
    //        gg->value);

    // printf("p2 output dump:\n");
    // for (int j = 0; j < Q->n; ++j) {
    //   printf("%ld, ", mptoi(output[j]));
    // }
    // printf("\n");

    for (i = 1; i < Q->n; i++) {
      gg->col_strategy[i - 1] = mptoi(output[i]) / den;
    }
  }
  return TRUE;
}

/*********************************************/
/* end of nash driver                        */
/*********************************************/

//==========================================================================
//   Building the problem representations (adapted from Gambit-enummixed)
//==========================================================================

//
// These two functions are based upon the program setupnash.c from the
// lrslib distribution, and the user's guide documentation.
// There are two separate functions, one for each player's problem.
// According to the user's guide, the ordering of the constraint rows
// is significant, and differs between the players; for player 1's problem
// the nonnegativity constraints come first, whereas for player 2's problem
// they appear later.  Experiments suggest this is in fact true, and
// reversing them breaks something.
//

//----------------------------------------------------------------------------------------//
void FillNonnegativityRows(lrs_dic *P, lrs_dat *Q, int firstRow, int lastRow,
                           int n) {
  const int MAXCOL = 1000; /* maximum number of columns */
  long num[MAXCOL], den[MAXCOL];
  long row, col;

  for (row = firstRow; row <= lastRow; row++) {
    num[0] = 0;
    den[0] = 1;

    for (col = 1; col < n; col++) {
      num[col] = (row - firstRow + 1 == col) ? 1 : 0;
      den[col] = 1;
    }
    lrs_set_row(P, Q, row, num, den, GE);
  }
}

//----------------------------------------------------------------------------------------//
void FillLinearityRow(lrs_dic *P, lrs_dat *Q, int m, int n) {
  const int MAXCOL = 1000; /* maximum number of columns */
  long num[MAXCOL], den[MAXCOL];
  int i;

  num[0] = -1;
  den[0] = 1;

  for (i = 1; i < n - 1; i++) {
    num[i] = 1;
    den[i] = 1;
  }

  num[n - 1] = 0;
  den[n - 1] = 1;

  lrs_set_row(P, Q, m, num, den, EQ);
}

//
// TL added this to get first row of ones. Don't know if it's needed
//----------------------------------------------------------------------------------------//
void FillFirstRow(lrs_dic *P, lrs_dat *Q, int n) {
  const int MAXCOL = 1000; /* maximum number of columns */
  long num[MAXCOL], den[MAXCOL];
  int i;

  for (i = 0; i < n; i++) {
    num[i] = 1;
    den[i] = 1;
  }
  lrs_set_row(P, Q, 0, num, den, GE);
}

void BuildRepP1Is0(lrs_dic *P, lrs_dat *Q, const FastInput *g) {
  long m = Q->m; /* number of inequalities      */
  long n = Q->n;
  FillConstraintRowsP1Is0(P, Q, g, 1);
  FillNonnegativityRows(P, Q, g->rows + 1, g->rows + g->cols, n);
  FillLinearityRow(P, Q, m, n);
  FillFirstRow(P, Q, n);
}

void BuildRepP1Is1(lrs_dic *P, lrs_dat *Q, const FastInput *g) {
  long m = Q->m; /* number of inequalities      */
  long n = Q->n;
  FillNonnegativityRows(P, Q, 1, g->rows, n);
  FillConstraintRowsP1Is1(P, Q, g, g->rows + 1); // 1 here
  FillLinearityRow(P, Q, m, n);
  FillFirstRow(P, Q, n);
}

void lrs_set_row_constraint(lrs_dic *P, lrs_dat *Q, long row, long num[],
                            long ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (GE)   - ordinary row  */
/*      = 0 (EQ)   - linearity     */
{
  lrs_mp Temp, mpone;
  lrs_mp_vector oD; /* denominator for row  */

  long i, j;

  /* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull;
  long m, d;
  lrs_alloc_mp(Temp);
  lrs_alloc_mp(mpone);
  hull = Q->hull;
  A = P->A;
  m = P->m;
  d = P->d;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector(d);
  itomp(ONE, mpone);
  itomp(ONE, oD[0]);

  i = row;
  itomp(ONE, Lcm[i]);         /* Lcm of denominators */
  itomp(ZERO, Gcd[i]);        /* Gcd of numerators */
  for (j = hull; j <= d; j++) /* hull data copied to cols 1..d */
  {
    itomp(num[j - hull], A[i][j]);
    itomp(ONE, oD[j]);
    if (!one(oD[j]))
      lcm(Lcm[i], oD[j]); /* update lcm of denominators */
    copy(Temp, A[i][j]);
    gcd(Gcd[i], Temp); /* update gcd of numerators   */
  }

  if (hull) {
    itomp(ZERO,
          A[i][0]); /*for hull, we have to append an extra column of zeroes */
    if (!one(A[i][1]) ||
        !one(oD[1])) /* all rows must have a one in column one */
      Q->polytope = FALSE;
  }
  if (!zero(A[i][hull]))    /* for H-rep, are zero in column 0     */
    Q->homogeneous = FALSE; /* for V-rep, all zero in column 1     */

  storesign(Gcd[i], POS);
  storesign(Lcm[i], POS);
  if (mp_greater(Gcd[i], mpone) || mp_greater(Lcm[i], mpone))
    for (j = 0; j <= d; j++) {
      exactdivint(A[i][j], Gcd[i], Temp); /*reduce numerators by Gcd  */
      mulint(Lcm[i], Temp, Temp);         /*remove denominators */
      exactdivint(Temp, oD[j], A[i][j]);  /*reduce by former denominator */
    }

  if (ineq == EQ) /* input is linearity */
  {
    Q->linearity[Q->nlinearity] = row;
    Q->nlinearity++;
  }

  /* 2010.4.26   Set Gcd and Lcm for the non-existant rows when nonnegative set
   */

  if (Q->nonnegative && row == m)
    for (j = 1; j <= d; j++) {
      itomp(ONE, Lcm[m + j]);
      itomp(ONE, Gcd[m + j]);
    }

  lrs_clear_mp_vector(oD, d);
  lrs_clear_mp(Temp);
  lrs_clear_mp(mpone);
} /* end of lrs_set_row_mp */

void FillConstraintRowsP1Is0(lrs_dic *P, lrs_dat *Q, const FastInput *g,
                             int firstRow) {
  const int MAXCOL = 1000; /* maximum number of columns */
  long num[MAXCOL];
  int row, s, t;

  for (row = firstRow; row < firstRow + g->rows; row++) {
    num[0] = 0;
    s = row - firstRow;
    for (t = 0; t < g->cols; t++) {
      num[t + 1] = -g->data[s * g->cols + t];
    }
    num[g->cols + 1] = 1;
    lrs_set_row_constraint(P, Q, row, num, GE);
  }
}

void FillConstraintRowsP1Is1(lrs_dic *P, lrs_dat *Q, const FastInput *g,
                             int firstRow) {
  const int MAXCOL = 1000; /* maximum number of columns */
  long num[MAXCOL];
  int row, s, t;

  for (row = firstRow; row < firstRow + g->cols; row++) {
    num[0] = 0;
    s = row - firstRow;
    for (t = 0; t < g->rows; t++) {
      num[t + 1] = g->data[t * g->cols + s] - g->den;
    }
    num[g->rows + 1] = 1;
    lrs_set_row_constraint(P, Q, row, num, GE);
  }
}
