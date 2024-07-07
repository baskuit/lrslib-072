/* lrslib.c     library code for lrs                     */

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

/*2022.1.22  lrs now has 3 modes controlled by these routines:
 lrs_run      default, V->H or H->V conversion
 redund_run   redundancy removal, redund or redund_list options
 fel_run      Fourier-Motzkin elimination, project or elimiminate options
*/

/* modified by Gary Roumanis and Skip Jordan for multithread compatability */
/* truncate needs mod to supress last pivot */
/* need to add a test for non-degenerate pivot step in reverse I guess?? */
/* Copyright: David Avis 2005,2022 avis@cs.mcgill.ca         */

#include <libgen.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <string.h>

#include "lrslib.h"

static unsigned long dict_count, dict_limit, cache_tries, cache_misses;

/* Variables and functions global to this file only */

static long lrs_checkpoint_seconds = 0;

static long lrs_global_count = 0; /* Track how many lrs_dat records are
                                     allocated */
static size_t infileLen;          /* length of cache of input file       */
static char *infile = NULL;       /* cache of input for restart          */
static char infilename[PATH_MAX];
static char outfilename[PATH_MAX];
static char tmpfilename[PATH_MAX];
static long overflow =
    0; /* =0 no overflow =1 restart overwrite =2 restart append */
static long pivoting =
    FALSE; /* =0 no overflow =1 restart overwrite =2 restart append */

static jmp_buf buf1;

static lrs_dat *lrs_global_list[MAX_LRS_GLOBALS + 1];

static lrs_dic *new_lrs_dic(long m, long d, long m_A);

static void cache_dict(lrs_dic **D_p, lrs_dat *global, long i, long j);
static long check_cache(lrs_dic **D_p, lrs_dat *global, long *i_p, long *j_p);
static void save_basis(lrs_dic *D, lrs_dat *Q);

static void lrs_dump_state();

static void pushQ(lrs_dat *global, long m, long d, long m_A);

#ifdef LRSLONG
static int tmpfd;
#endif

#ifndef TIMES
static void ptimes(void);
static double get_time(void);
#endif

char *basename(char *path);

/*******************************/
/* signals handling            */
/*******************************/
#ifndef SIGNALS
static void checkpoint(int);
static void die_gracefully(int);
static void setup_signals(void);
static void timecheck(int);
#endif

/*******************************/
/* functions  for external use */
/*******************************/

/******************************************************************/
/* lrs_run is the main reverse search part of lrs                 */
/* should be called by lrsv2_main which does setup and close also */
/******************************************************************/
long lrs_run(lrs_dic *P, lrs_dat *Q)

{

  lrs_mp_matrix Lin; /* holds input linearities if any are found             */
  long col;          /* output column index for dictionary                   */
  long startcol = 0;
  long prune = FALSE; /* if TRUE, getnextbasis will prune tree and backtrack  */

  /*********************************************************************************/
  /* Step 1: Find a starting cobasis from default of specified order */
  /*         P is created to hold  active dictionary data and may be cached */
  /*         Lin is created if necessary to hold linearity space */
  /*         Print linearity space if any, and retrieve output from first dict.
   */
  /*********************************************************************************/

  if (lrs_getfirstbasis(&P, Q, &Lin, FALSE) == 0) {
    lrs_free_dic(P, Q); /* note Q is not free here and can be reused     */
    return 1;
  }

  /* Pivot to a starting dictionary                      */
  /* There may have been column redundancy               */
  /* If so the linearity space is obtained and redundant */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */

  if (Q->homogeneous && Q->hull)
    startcol++; /* col zero not treated as redundant   */

  for (col = startcol; col < Q->nredundcol; col++) /* print linearity space */
    lrs_printoutput(Q, Lin[col]); /* Array Lin[][] holds the coeffs.     */

  if (Q->nredundcol > 0)
    lrs_clear_mp_matrix(Lin, Q->nredundcol, Q->n);

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

    // 2015.6.5   after maxcobases reached, generate subtrees that have not been
    // enumerated 2018.1.19  fix printcobasis bug when maxcobases set 2019.5.8
    // new givoutput flag to avoid printing restart cobases

    prune = lrs_checkbound(P, Q);

    if (!prune && Q->giveoutput) {
      lrs_open_outputblock(); /* keeps output together when using mplrs */

      for (col = 0; col <= P->d; col++) /* print output if any */
        if (lrs_getsolution(P, Q, Q->output, col))
          lrs_printoutput(Q, Q->output);
      pivoting = TRUE;
      lrs_close_outputblock();
    } else
      Q->giveoutput = TRUE; /* first output supressed for restart */

    /*2020.3.9  bounds on objective function check corrected */

    if ((Q->maxcobases > 0) && (Q->count[2] >= Q->maxcobases)) {
      prune = TRUE;
      if (!lrs_leaf(P, Q)) /* do not return cobases of a leaf */
        lrs_return_unexplored(P, Q);
    }

    save_basis(P, Q);

  } while (!Q->lponly && lrs_getnextbasis(&P, Q, prune)); // do ...

  if (Q->lponly)
    lrs_lpoutput(P, Q, Q->output);
  else
    lrs_printtotals(P, Q); /* print final totals, including estimates       */

  Q->m = P->m;
  lrs_free_dic(P, Q); /* note Q is not free here and can be reused     */

  return 0;
}
/*********************************************/
/* end of model test program for lrs library */
/*********************************************/

/*******************************************************/
/* redund_run is main loop for redundancy removal      */
/*******************************************************/
long redund_run(lrs_dic *P, lrs_dat *Q)

{
  lrs_mp_matrix Ain; /* holds a copy of the input matrix to output at the end */

  long ineq; /* input inequality number of current index             */
  long *redineq;

  lrs_mp_matrix Lin; /* holds input linearities if any are found             */

  long i, j, d, m;
  long nlinearity; /* number of linearities in input file                  */
  long lastdv;
  long index; /* basic index for redundancy test */
  long c1 = 0;
  long min, nin;

  /*********************************************************************************/

  /* if non-negative flag is set, non-negative constraints are not input */
  /* explicitly, and are not checked for redundancy                      */

  m = P->m_A; /* number of rows of A matrix */
  d = P->d;
  redineq = Q->redineq;
  min = Q->m;
  nin = Q->n;
  Q->Ain = lrs_alloc_mp_matrix(
      Q->m, Q->n); /* make a copy of A matrix for output later            */
  Ain = Q->Ain;

  for (i = 1; i <= m; i++) {
    for (j = 0; j <= d; j++)
      copy(Ain[i][j], P->A[i][j]);
  }

  /*********************************************************************************/
  /* Step 1: Find a starting cobasis from default of specified order */
  /*         Lin is created if necessary to hold linearity space */
  /*********************************************************************************/

  if (!lrs_getfirstbasis(&P, Q, &Lin, TRUE))
    return 1;
  /* Pivot to a starting dictionary                      */
  /* There may have been column redundancy               */
  /* If so the linearity space is obtained and redundant */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1  */

  if (Q->nredundcol > 0)
    lrs_clear_mp_matrix(Lin, Q->nredundcol, Q->n);

  /*********************************************************************************/
  /* Step 2: Test rows i where redineq[i]=1 for redundancy */
  /*********************************************************************************/

  /* note some of these may have been changed in getting initial dictionary */
  m = P->m_A;
  d = P->d;
  nlinearity = Q->nlinearity;
  lastdv = Q->lastdv;

  /* linearities are not considered for redundancy */

  for (i = 0; i < nlinearity; i++)
    redineq[Q->linearity[i]] = 2L;

  /* Q->verifyredund always false in lrs, set by mplrs to check duplicated
   * redundancy removal */
  /* Q->noredundcheck overides this to skip verification */

  if (Q->noredundcheck && Q->verifyredund)
    goto done;

  /* mplrs sets redineq[i]==-1 for guaranteed redundant inequalities */
  /* these rows must be zeroed out before testing the others         */

  if (Q->verifyredund) /* this is never run by lrs, final step of mplrs redund
                        */
  {
    for (index = lastdv + Q->redineq[0]; index <= m + d; index++) {
      ineq = Q->inequality[index - lastdv]; /* the input inequality number corr.
                                               to this index */

      if (redineq[ineq] == 1) {
        c1++;
      }
      if (redineq[ineq] == -1) {
        checkindex(P, Q,
                   -index); /* used to zero correct row of A no LP solved */
      }
    }
  }

  /* rows 0..lastdv are cost, decision variables, or linearities  */
  /* other rows need to be tested                                */

  if (Q->redineq[0] == 0) /* 2020.9.11 patch, this was set to 1 in readredund
                             but got reset somewhere */
    Q->redineq[0] = 1;

  for (index = lastdv + Q->redineq[0]; index <= m + d; index++) {
    ineq = Q->inequality[index - lastdv]; /* the input inequality number corr.
                                             to this index */
    Q->redineq[0] = ineq; /* used for restarting after arithmetic overflow    */

    if (redineq[ineq] == 1) {
      redineq[ineq] = checkindex(P, Q, index);
    }

  } /* end for index ..... */

done:

  if ((Q->mplrs && !Q->verifyredund)) {
    lrs_clear_mp_matrix(Q->Ain, min, nin);
    Q->m = P->m;
    lrs_free_dic(P, Q); /* note Q is not free here and can be reused     */
    return 0;
  }

  if (Q->fel && Q->hull)
    lrs_clear_mp_matrix(Q->Ain, min, nin);
  else
    redund_print(P, Q);

  if (Q->mplrs && !Q->noredundcheck)
    fprintf(lrs_ofp, "* %ld row(s) needed verifying\n", c1);

  if (!Q->fel)
    lrs_clear_mp_matrix(Q->Ain, min, nin);

  lrs_free_dic(P, Q);
  return 0;
}
/***********redund_run************************/

void redund_print(lrs_dic *P, lrs_dat *Q) {
  long i, j, m;
  long nlinearity; /* number of linearities in input file                  */
  long nredund;    /* number of redundant rows in input file               */
  long *redineq = Q->redineq;
  lrs_mp_matrix Ain = Q->Ain;

  m = P->m_A; /* number of rows of A matrix */
  nlinearity = Q->nlinearity;

  /* restore as mplrs loses this */
  for (i = 0; i < nlinearity; i++)
    redineq[Q->linearity[i]] = 2;

  /*
    fprintf(lrs_ofp,"\nQ->red");
    for (i = 1; i <= m; i++)
    fprintf(lrs_ofp," %ld",Q->redineq[i]);
  */

  if (!Q->hull)
    fprintf(lrs_ofp, "\nH-representation");
  else
    fprintf(lrs_ofp, "\nV-representation");

  /* linearities will be printed first in output */

  if (nlinearity > 0) {
    fprintf(lrs_ofp, "\nlinearity %ld", nlinearity);
    for (i = 1; i <= nlinearity; i++)
      fprintf(lrs_ofp, " %ld", i);
  }

  nredund = 0; /* count number of non-redundant inequalities */

  for (i = 1; i <= m; i++)
    if (redineq[i] == 0)
      nredund++;

  fprintf(lrs_ofp, "\nbegin");
  fprintf(lrs_ofp, "\n%ld %ld rational", nlinearity + nredund, Q->n);

  /* print the linearities first */

  for (i = 0; i < nlinearity; i++)
    lrs_printrow("", Q, Ain[Q->linearity[i]], Q->inputd);

  for (i = 1; i <= m; i++) {
    if (redineq[i] == 0)
      lrs_printrow("", Q, Ain[i], Q->inputd);

    /* remove redundant rows */
    /*
        if( redineq[i] != 1 && redineq[i] != -1)
         {
          row++;
          for(j=0; j<= Q->inputd;j++)
             copy(P->A[row][j],Ain[i][j]);
         }
    */
  }

  fprintf(lrs_ofp, "\nend");

  if (Q->redund)
    fprintf(lrs_ofp, "\n*Input had %ld rows and %ld columns", m, Q->n);

  redineq[0] = m - nredund - nlinearity; /* number of redundant rows */

  if (m == nredund || redineq[0] == 0) {
    if (Q->redund)
      fprintf(lrs_ofp, "\n*No redundant rows found\n");
  } else {
    j = 0;
    if (Q->redund) {
      fprintf(lrs_ofp, "\n* %ld redundant row(s) found\n", redineq[0]);
      for (i = 1; i <= m; i++)
        if (redineq[i] == 1 || redineq[i] == -1) {
          j++;
          if (j > 20) {
            j = 1;
            fprintf(lrs_ofp, "\n %ld", i);
          } else
            fprintf(lrs_ofp, " %ld", i);
        }
    }
    if (Q->noredundcheck)
      fprintf(lrs_ofp, "\n*Warning: not verified - input should be full "
                       "dimensional and duplicate free");
  }
  fprintf(lrs_ofp, "\n");
  return;
} /* end of redund_print */

/*******************/
/* lrs_printoutput */
/*  one line only   */
/*******************/
void lrs_printoutput(lrs_dat *Q, lrs_mp_vector output) {
  char *sss;
  char **ss;

  long i;
  long len = 0;

  if (Q->countonly)
    return;

  ss = (char **)malloc((1 + Q->n) * sizeof(char *));

  if (Q->hull || zero(output[0])) /*non vertex */
    for (i = 0; i < Q->n; i++) {
      ss[i] = cpmp("", output[i]);
      len = len + snprintf(NULL, 0, "%s ", ss[i]);
    }
  else
    for (i = 1; i < Q->n; i++) {
      ss[i] = cprat("", output[i], output[0]);
      len = len + snprintf(NULL, 0, "%s ", ss[i]);
    }

  sss = (char *)malloc((len + 5) * sizeof(char *));
  len = 0;

  if (Q->hull || zero(output[0])) /*non vertex */
    for (i = 0; i < Q->n; i++) {
      len = len + sprintf(sss + len, "%s ", ss[i]);
      free(ss[i]);
    }
  else { /* vertex   */
    len = sprintf(sss, " 1 ");
    for (i = 1; i < Q->n; i++) {
      len = len + sprintf(sss + len, "%s ", ss[i]);
      free(ss[i]);
    }
  }

  if (Q->mplrs)
    lrs_post_output("vertex", sss);
  else
    fprintf(lrs_ofp, "\n%s", sss);

  free(ss);
  free(sss);
}
/**************************/
/* end of lrs_printoutput */
/**************************/

/****************/
/* lrs_lpoutput */
/****************/
void lrs_lpoutput(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output) {

  if (Q->unbounded || !Q->messages)
    return;

  lrs_mp Temp1, Temp2;
  long i;

  lrs_alloc_mp(Temp1);
  lrs_alloc_mp(Temp2);

  prat("\n*Obj=", P->objnum, P->objden);
  fprintf(lrs_ofp, "    pivots=%ld ", Q->count[3]);
  fprintf(lrs_ofp, "\n");
  lrs_clear_mp(Temp1);
  lrs_clear_mp(Temp2);
}
/***********************/
/* end of lrs_lpoutput */
/***********************/
void lrs_printrow(const char *name, lrs_dat *Q, lrs_mp_vector output, long rowd)
/* print a row of A matrix in output in "original" form  */
/* rowd+1 is the dimension of output vector                */
/* if input is H-rep. output[0] contains the RHS      */
/* if input is V-rep. vertices are scaled by 1/output[1] */
{
  long i;
  fprintf(lrs_ofp, "\n%s", name);
  if (!Q->hull) /* input was inequalities, print directly */

  {
    for (i = 0; i <= rowd; i++)
      pmp("", output[i]);
    return;
  }

  /* input was vertex/ray */

  if (zero(output[1])) /*non-vertex */
  {
    for (i = 1; i <= rowd; i++)
      pmp("", output[i]);
  } else { /* vertex */
    fprintf(lrs_ofp, " 1 ");
    for (i = 2; i <= rowd; i++)
      prat("", output[i], output[1]);
  }

  return;

} /* end of lrs_printrow */

long lrs_getsolution(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output, long col)
/* check if column indexed by col in this dictionary */
/* contains output                                   */
/* col=0 for vertex 1....d for ray/facet             */
{

  long j; /* cobasic index     */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;

  if (col == ZERO) /* check for lexmin vertex */
    return lrs_getvertex(P, Q, output);

  /*  check for rays: negative in row 0 , positive if lponly */

  if (Q->lponly) {
    if (!positive(A[0][col]))
      return FALSE;
  }

  else if (!negative(A[0][col]))
    return FALSE;

  /*  and non-negative for all basic non decision variables */

  j = Q->lastdv + 1;
  while (j <= P->m && !negative(A[Row[j]][col]))
    j++;

  if (j <= P->m)
    return FALSE;
  /* 2021.5.21 check for max depth output */
  if (lexmin(P, Q, col)) {
    if (P->depth > Q->count[8])
      Q->count[8] = P->depth;
    return lrs_getray(P, Q, col, Q->n, output);
  }

  if (Q->geometric || Q->allbases || Q->lponly)
    return lrs_getray(P, Q, col, Q->n, output);

  return FALSE; /* no more output in this dictionary */

} /* end of lrs_getsolution */

void lrs_print_header(const char *name) {
  if (lrs_ofp == NULL)
    lrs_ofp = stdout;
#ifdef LRS_QUIET
  return;
#endif
  fprintf(lrs_ofp, "%s:", name);
  fprintf(lrs_ofp, TITLE);
  fprintf(lrs_ofp, VERSION);
  fprintf(lrs_ofp, "(");
  fprintf(lrs_ofp, BIT);
  fprintf(lrs_ofp, ",");
  fprintf(lrs_ofp, ARITH);
#ifdef MA
  fprintf(lrs_ofp, ",hybrid_arithmetic");
#endif
#ifdef LRSLONG
#ifndef SAFE
  fprintf(lrs_ofp, ",no_overflow_checking");
#endif
#endif
  fprintf(lrs_ofp, ")");
  if (overflow != 2) {
#ifdef GMP
    fprintf(lrs_ofp, "_gmp_v.%d.%d", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
#elif defined(FLINT)
    fprintf(lrs_ofp, "_%dbit_flint_v.%s", FLINT_BITS, FLINT_VERSION);
#endif
  }
}

long lrs_init(const char *name) /* returns TRUE if successful, else FALSE */
{

#ifndef LRS_QUIET
  if (overflow != 2)
    lrs_print_header(name);
#endif

  if (!lrs_mp_init(0, stdin, stdout)) /* initialize arithmetic */
    return FALSE;

  lrs_global_count = 0;
  lrs_checkpoint_seconds = 0;
#ifndef SIGNALS
  setup_signals();
#endif
  return TRUE;
}

void lrs_close(const char *name) {

#ifdef LRS_QUIET
  fprintf(lrs_ofp, "\n");
  if (lrs_ofp != stdout) {
    fclose(lrs_ofp);
    lrs_ofp = NULL;
  }
  return;
#endif

#ifdef LRSLONG
#ifdef SAFE
  fprintf(lrs_ofp, "\n*overflow checking on lrslong arithmetic");
#else
  fprintf(lrs_ofp,
          "\n*caution: no overflow checking on long integer arithemtic");
#endif
#endif

  fprintf(lrs_ofp, "\n*%s:", name);
  fprintf(lrs_ofp, TITLE);
  fprintf(lrs_ofp, VERSION);
  fprintf(lrs_ofp, "(");
  fprintf(lrs_ofp, BIT);
  fprintf(lrs_ofp, ",");
  fprintf(lrs_ofp, ARITH);
#ifdef MA
  fprintf(lrs_ofp, ",hybrid arithmetic");
#endif
  fprintf(lrs_ofp, ")");

#ifdef MP
  fprintf(lrs_ofp, " max decimal digits=%ld/%ld", DIG2DEC(lrs_record_digits),
          DIG2DEC(lrs_digits));
#endif

#ifndef TIMES
  ptimes();
#endif

  if (lrs_ofp != stdout) {
    fclose(lrs_ofp);
    lrs_ofp = NULL;
  }
}

/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
lrs_dat *lrs_alloc_dat(const char *name) {
  lrs_dat *Q;
  long i;

  if (lrs_global_count >= MAX_LRS_GLOBALS) {
    fprintf(stderr,
            "Fatal: Attempt to allocate more than %ld global data blocks\n",
            MAX_LRS_GLOBALS);
    return NULL;
  }

  Q = (lrs_dat *)malloc(sizeof(lrs_dat));
  if (Q == NULL)
    return Q; /* failure to allocate */

  lrs_global_list[lrs_global_count] = Q;
  Q->id = lrs_global_count;
  lrs_global_count++;
  Q->name = (char *)CALLOC((unsigned)strlen(name) + 1, sizeof(char));
  strcpy(Q->name, name);

  /* initialize variables */
  Q->mplrs = FALSE;
  Q->messages = TRUE;
// #ifdef PLRS
// Q->mplrs=TRUE;
// #endif
#ifdef LRS_QUIET
  Q->messages = FALSE;
#endif
  strcpy(Q->fname, ""); /* name of program, filled in later */
  Q->m = 0L;
  Q->n = 0L;
  Q->inputd = 0L;
  Q->deepest = 0L;
  Q->nlinearity = 0L;
  Q->nredundcol = 0L;
  Q->runs = 0L;
  Q->subtreesize = MAXD;
  Q->seed = 1234L;
  Q->totalnodes = 0L;
  for (i = 0; i < 10; i++) {
    Q->count[i] = 0L;
    Q->cest[i] = 0.0;
    if (i < 5)
      Q->startcount[i] = 0L;
  }
  Q->count[2] = 1L;      /* basis counter */
  Q->startcount[2] = 0L; /* starting basis counter */
                         /* initialize flags */
  Q->allbases = FALSE;
  Q->bound = FALSE;     /* upper/lower bound on objective function given */
  Q->countonly = FALSE; /* produce the usual output */
  Q->frequency = 0L;
  Q->dualdeg = FALSE; /* TRUE if dual degenerate starting dictionary */
  Q->geometric = FALSE;
  Q->getvolume = FALSE;
  Q->homogeneous = TRUE;
  Q->polytope = FALSE;
  Q->triangulation = FALSE;
  Q->hull = FALSE;
  Q->incidence = FALSE;
  Q->lponly = FALSE;
  Q->maxdepth = MAXD;
  Q->mindepth = -MAXD;
  Q->maxoutput = 0L;
  Q->maxcobases =
      0L; /* after maxcobases have been found unexplored subtrees reported */

  Q->redund = FALSE;
  Q->fel = FALSE;

  Q->nonnegative = FALSE;
  Q->printcobasis = FALSE;
  Q->printslack = FALSE;
  Q->truncate = FALSE; /* truncate tree when moving from opt vertex        */
  Q->extract = FALSE;
  Q->voronoi = FALSE;
  Q->maximize = FALSE;   /*flag for LP maximization                          */
  Q->minimize = FALSE;   /*flag for LP minimization                          */
  Q->givenstart = FALSE; /* TRUE if a starting cobasis is given              */
  Q->newstart = FALSE;
  Q->giveoutput = TRUE; /* set to false for first output after restart      */
  Q->verifyredund = FALSE;  /* set to true when mplrs verifies redund output  */
  Q->noredundcheck = FALSE; /* set to true when mplrs skips verifying output */
  Q->nextineq = 15; /* start redundancy testing from this row           */

  Q->facet = NULL;
  Q->redundcol = NULL;
  Q->inequality = NULL;
  Q->linearity = NULL;
  Q->vars = NULL;
  Q->startcob = NULL;
  Q->minratio = NULL;
  Q->temparray = NULL;
  Q->redineq = NULL;
  Q->Ain = NULL;
  Q->olddic = NULL;

  Q->saved_flag = 0; /* no cobasis saved initially, db */
  lrs_alloc_mp(Q->Nvolume);
  lrs_alloc_mp(Q->Dvolume);
  lrs_alloc_mp(Q->sumdet);
  lrs_alloc_mp(Q->saved_det);
  lrs_alloc_mp(Q->boundn);
  lrs_alloc_mp(Q->boundd);
  itomp(ZERO, Q->Nvolume);
  itomp(ONE, Q->Dvolume);
  itomp(ZERO, Q->sumdet);
  Q->unbounded = FALSE;

  return Q;
} /* end of allocate and initialize lrs_dat */

/* In lrs_getfirstbasis and lrs_getnextbasis we use D instead of P */
/* since the dictionary P may change, ie. &P in calling routine    */

#define D (*D_p)

long lrs_getfirstbasis(lrs_dic **D_p, lrs_dat *Q, lrs_mp_matrix *Lin,
                       long no_output)
/* gets first basis, FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is TRUE supresses output headers   */
/* 2017.12.22  could use no_output=2 to get early exit for criss-cross method */
{
  lrs_mp scale, Temp;

  long i, j, k;

  /* assign local variables to structures */

  lrs_mp_matrix A;
  long *B, *C, *Col;
  long *inequality;
  long *linearity;
  long hull = Q->hull;
  long m, d, lastdv, nlinearity, nredundcol;

  if (Q->lponly)
    no_output = TRUE;
  m = D->m;
  d = D->d;

  lastdv = Q->lastdv;

  nredundcol = 0L; /* will be set after getabasis        */
  nlinearity =
      Q->nlinearity; /* may be reset if new linearity read or in getabasis*/
  linearity = Q->linearity;

  A = D->A;
  B = D->B;
  C = D->C;
  Col = D->Col;
  inequality = Q->inequality;

  /**2020.6.17 with extract just select columns and quit */

  if (Q->extract)
    if (Q->hull || Q->nlinearity == 0) {
      extractcols(D, Q);
      return FALSE;
    }

  lrs_alloc_mp(Temp);
  lrs_alloc_mp(scale);

  if (Q->runs > 0) /* arrays for estimator */
  {
    Q->isave = (long *)CALLOC((unsigned)(m * d), sizeof(long));
    Q->jsave = (long *)CALLOC((unsigned)(m * d), sizeof(long));
  }
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
  /* for voronoi convert to h-description using the transform */
  /* a_0 .. a_d-1 -> (a_0^2 + ... a_d-1 ^2)-2a_0x_0-...-2a_d-1x_d-1 + x_d >= 0
   */
  /* note constant term is stored in column d, and column d-1 is all ones */
  /* the other coefficients are multiplied by -2 and shifted one to the right */
  if (Q->voronoi) {
    Q->hull = FALSE;
    hull = FALSE;
    for (i = 1; i <= m; i++) {
      if (zero(A[i][1])) {
        printf("\nWith voronoi option column one must be all one\n");
        return (FALSE);
      }
      copy(scale, A[i][1]); /*adjust for scaling to integers of rationals */
      itomp(ZERO, A[i][0]);
      for (j = 2; j <= d; j++) /* transform each input row */
      {
        copy(Temp, A[i][j]);
        mulint(A[i][j], Temp, Temp);
        linint(A[i][0], ONE, Temp, ONE);
        linint(A[i][j - 1], ZERO, A[i][j], -TWO);
        mulint(scale, A[i][j - 1], A[i][j - 1]);
      } /* end of for (j=1;..) */
      copy(A[i][d], scale);
      mulint(scale, A[i][d], A[i][d]);
    } /* end of for (i=1;..) */
  } /* end of if(voronoi)     */
  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      itomp(ZERO, A[0][j]);

  /* Now we pivot to standard form, and then find a primal feasible basis */
  /* Note these steps MUST be done, even if restarting, in order to get */
  /* the same index/inequality correspondance we had for the original prob. */
  /* The inequality array is used to give the insertion order */
  /* and is defaulted to the last d rows when givenstart=FALSE */

  if (Q->nonnegative) {
    /* no need for initial pivots here, labelling already done */
    Q->lastdv = d;
    Q->nredundcol = 0;
  } else {
    if (!getabasis(D, Q, inequality))
      return FALSE;
    /* bug fix 2009.12.2 */
    nlinearity =
        Q->nlinearity; /*may have been reset if some lins are redundant*/
  }

  /* 2020.2.2 */
  /* extract option asked to remove all linearities and output the reduced A
   * matrix */
  /* should be followed by redund to get minimum representation */

  if (Q->extract) {
    linextractcols(D, Q);
    return FALSE;
  }

  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = D->d;

  /********************************************************************/
  /* now we start printing the output file  unless no output requested */
  /********************************************************************/

  if (Q->count[2] == 1 && (no_output == 0)) /* don't reprint after newstart */
  {
    int len = 0;
    char *header;

    header = (char *)malloc((100 + 20 * Q->n) * sizeof(char));

    if (Q->voronoi)
      len = sprintf(header,
                    "*Voronoi Diagram: Voronoi vertices and rays are output");
    else {
      if (hull)
        len = sprintf(header, "H-representation");
      else
        len = sprintf(header, "V-representation");
    }

    /* Print linearity space                 */
    /* Don't print linearity if first column zero in hull computation */

    if (hull && Q->homogeneous)
      k = 1; /* 0 normally, 1 for homogeneous case     */
    else
      k = 0;

    if (nredundcol > k) {
      len = len + sprintf(header + len, "\nlinearity %ld ",
                          nredundcol - k); /*adjust nredundcol for homog. */
      for (i = 1; i <= nredundcol - k; i++)
        len = len + sprintf(header + len, " %ld", i);
    } /* end print of linearity space */

    len = len + sprintf(header + len, "\nbegin");
    len = len + sprintf(header + len, "\n***** %ld rational", Q->n);

    if (Q->mplrs)
      lrs_post_output("header", header);
    else if (Q->messages)
      fprintf(lrs_ofp, "\n%s", header);

    free(header);
  }

  /* end of if !no_output .......   */

  /* Reset up the inequality array to remember which index is which input
   * inequality */
  /* inequality[B[i]-lastdv] is row number of the inequality with index B[i] */
  /* inequality[C[i]-lastdv] is row number of the inequality with index C[i] */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0) /* some cobasic indices will be removed */
  {
    for (i = 0; i < nlinearity; i++) /* remove input linearity indices */
      inequality[linearity[i]] = 0;
    k = 1; /* counter for linearities         */
    for (i = 1; i <= m - nlinearity; i++) {
      while (k <= m && inequality[k] == 0)
        k++; /* skip zeroes in corr. to linearity */
      inequality[i] = inequality[k++];
    }
  }
  /* end if linearity */

  if (nredundcol > 0) {
    const unsigned int Qn = Q->n;
    *Lin = lrs_alloc_mp_matrix(nredundcol, Qn);

    for (i = 0; i < nredundcol; i++) {

      if (!(Q->homogeneous && Q->hull &&
            i == 0)) /* skip redund col 1 for homog. hull */
      {

        lrs_getray(D, Q, Col[0], D->C[0] + i - hull,
                   (*Lin)[i]); /* adjust index for deletions */
      }

      if (!removecobasicindex(D, Q, 0L)) {
        lrs_clear_mp_matrix(*Lin, nredundcol, Qn);
        return FALSE;
      }
    }

  } /* end if nredundcol > 0 */

  /*2017.12.22   If you want to do criss-cross now is the time ! */

  /* Do dual pivots to get primal feasibility */
  if (!primalfeasible(D, Q)) {
    if (!Q->mplrs)
      fprintf(lrs_ofp, "\nend");
    lrs_warning(Q, "finalwarn", "\nNo feasible solution\n");
    return FALSE;
  }

  /* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize) {
    Q->unbounded = !lrs_solvelp(D, Q, Q->maximize);
    if (Q->lponly) {

      lrs_clear_mp(Temp);
      lrs_clear_mp(scale);
      return TRUE;
    }

    else /* check to see if objective is dual degenerate */
    {
      j = 1;
      while (j <= d && !zero(A[0][j]))
        j++;
      if (j <= d)
        Q->dualdeg = TRUE;
    }
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
    j = inequality[B[i] - lastdv];
    inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
    inequality[C[0] - lastdv] = j;
    C[0] = B[i];
    B[i] = i;
    reorder1(C, Col, ZERO, d);
  }

  /* Check to see if necessary to resize */
  /* bug fix 2018.6.7 new value of d required below */
  if (Q->inputd > D->d)
    *D_p = resize(D, Q);

  lrs_clear_mp(Temp);
  lrs_clear_mp(scale);
  return TRUE;
}
/********* end of lrs_getfirstbasis  ***************/

/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/

long lrs_getnextbasis(lrs_dic **D_p, lrs_dat *Q, long backtrack)
/* gets next reverse search tree basis, FALSE if none  */
/* switches to estimator if maxdepth set               */
/* backtrack TRUE means backtrack from here            */

{
  /* assign local variables to structures */
  long i = 0L, j = 0L;
  long m = D->m;
  long d = D->d;
  long saveflag;
  long cob_est =
      0; /* estimated number of cobases in subtree from current node */

  if (backtrack && D->depth == 0)
    return FALSE; /* cannot backtrack from root      */

  if (Q->maxoutput > 0 && Q->count[0] + Q->count[1] - Q->hull >= Q->maxoutput)
    return FALSE; /* output limit reached            */

  while ((j < d) || (D->B[m] != m)) /*main while loop for getnextbasis */
  {
    if (D->depth >= Q->maxdepth) {
      if (Q->runs > 0 && !backtrack) /*get an estimate of remaining tree */
      {

        // 2015.2.9 do iterative estimation backtracking when estimate is small

        saveflag = Q->printcobasis;
        Q->printcobasis = FALSE;
        cob_est = lrs_estimate(D, Q);
        Q->printcobasis = saveflag;
        if (cob_est <= Q->subtreesize) /* stop iterative estimation */
        {
          backtrack = TRUE;
        }

      } else // either not estimating or we are backtracking

        // 2018.1.19              if (!backtrack && !Q->printcobasis)
        if (!backtrack)
          if (!lrs_leaf(D, Q)) /* 2015.6.5 cobasis returned if not a leaf */
            lrs_return_unexplored(D, Q);

      backtrack = TRUE;

      if (Q->maxdepth == 0 &&
          cob_est <= Q->subtreesize) /* root estimate only */
        return FALSE;                /* no nextbasis  */
    } // if (D->depth >= Q->maxdepth)

    /*      if ( Q->truncate && negative(D->A[0][0]))*/ /* truncate when moving
                                                           from opt. vertex */
    /*          backtrack = TRUE;    2011.7.14 */

    if (backtrack) /* go back to prev. dictionary, restore i,j */
    {
      backtrack = FALSE;

      if (check_cache(D_p, Q, &i, &j)) {
      } else {
        D->depth--;
        selectpivot(D, Q, &i, &j);
        pivot(D, Q, i, j);
        update(D, Q, &i, &j); /*Update B,C,i,j */
      }

      j++; /* go to next column */
    } /* end of if backtrack  */

    if (D->depth < Q->mindepth)
      break;

    /* try to go down tree */

    /* 2011.7.14 patch */
    while ((j < d) &&
           (!reverse(D, Q, &i, j) || (Q->truncate && Q->minratio[D->m] == 1)))
      j++;
    if (j == d)
      backtrack = TRUE;
    else
    /*reverse pivot found */
    {
      cache_dict(D_p, Q, i, j);
      /* Note that the next two lines must come _after_ the
         call to cache_dict */

      D->depth++;
      if (D->depth > Q->deepest)
        Q->deepest++;

      pivot(D, Q, i, j);
      update(D, Q, &i, &j); /*Update B,C,i,j */

      D->lexflag = lexmin(D, Q, ZERO); /* see if lexmin basis */
      Q->count[2]++;
      Q->totalnodes++;

      return TRUE; /*return new dictionary */
    }

  } /* end of main while loop for getnextbasis */
  return FALSE; /* done, no more bases */
} /*end of lrs_getnextbasis */

/*************************************/
/* print out one line of output file */
/*************************************/
long lrs_getvertex(lrs_dic *P, lrs_dat *Q, lrs_mp_vector output)
/*Print out current vertex if it is lexmin and return it in output */
/* return FALSE if no output generated  */
{
  lrs_mp_matrix A = P->A;
  long i;
  long ind;  /* output index                                  */
  long ired; /* counts number of redundant columns            */
             /* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long *B = P->B;
  long *Row = P->Row;

  long lastdv = Q->lastdv;

  long hull;
  long lexflag;

  hull = Q->hull;
  lexflag = P->lexflag;
  if (lexflag || Q->allbases) {
    ++(Q->count[1]);
    /* 2021.5.21 check for max depth output */
    if (P->depth > Q->count[8])
      Q->count[8] = P->depth;
  }

  if (Q->getvolume) {
    linint(Q->sumdet, 1, P->det, 1);
    updatevolume(P, Q);
  }
  if (Q->triangulation) /* this will print out a triangulation */
    lrs_printcobasis(P, Q, ZERO);

  /*print cobasis if printcobasis=TRUE and count[2] a multiple of frequency */
  /* or for lexmin basis, except origin for hull computation - ugly!        */

  if (Q->printcobasis)
    if ((lexflag && !hull) ||
        ((Q->frequency > 0) &&
         (count[2] == (count[2] / Q->frequency) * Q->frequency)))
      lrs_printcobasis(P, Q, ZERO);

  if (hull)
    return FALSE; /* skip printing the origin */

  if (!lexflag && !Q->allbases &&
      !Q->lponly) /* not lexmin, and not printing forced */
    return FALSE;

  /* copy column 0 to output */

  i = 1;
  ired = 0;
  copy(output[0], P->det);

  for (ind = 1; ind < Q->n; ind++) /* extract solution */

    if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
    /* column was deleted as redundant */
    {
      itomp(ZERO, output[ind]);
      ired++;
    } else
    /* column not deleted as redundant */
    {
      getnextoutput(P, Q, i, ZERO, output[ind]);
      i++;
    }

  reducearray(output, Q->n);
  if (lexflag && one(output[0]))
    ++Q->count[4]; /* integer vertex */

  /* uncomment to print nonzero basic variables

   printf("\n nonzero basis: vars");
    for(i=1;i<=lastdv; i++)
     {
      if ( !zero(A[Row[i]][0]) )
           printf(" %ld ",B[i]);
     }
  */

  /* printslack inequality indices  */

  if (Q->printslack) {
    fprintf(lrs_ofp, "\nslack ineq:");
    for (i = lastdv + 1; i <= P->m; i++) {
      if (!zero(A[Row[i]][0]))
        fprintf(lrs_ofp, " %ld ", Q->inequality[B[i] - lastdv]);
    }
  }

  return TRUE;
} /* end of lrs_getvertex */

long lrs_getray(lrs_dic *P, lrs_dat *Q, long col, long redcol,
                lrs_mp_vector output)
/*Print out solution in col and return it in output   */
/*redcol =n for ray/facet 0..n-1 for linearity column */
/*hull=1 implies facets will be recovered             */
/* return FALSE if no output generated in column col  */
{
  long i;
  long ind;  /* output index                                  */
  long ired; /* counts number of redundant columns            */
             /* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull = Q->hull;
  long n = Q->n;

  long *B = P->B;
  long *Row = P->Row;
  long lastdv = Q->lastdv;

  if (redcol == n) {
    ++count[0];
    if (Q->printcobasis)
      lrs_printcobasis(P, Q, col);
  }

  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++) /* print solution */
  {
    if (ind == 0 && !hull) /* must have a ray, set first column to zero */
      itomp(ZERO, output[0]);

    else if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
    /* column was deleted as redundant */
    {
      if (redcol == ind) /* true for linearity on this cobasic index */
        /* we print reduced determinant instead of zero */
        copy(output[ind], P->det);
      else
        itomp(ZERO, output[ind]);
      ired++;
    } else
    /* column not deleted as redundant */
    {
      getnextoutput(P, Q, i, col, output[ind]);
      i++;
    }
  }
  reducearray(output, n);
  /* printslack for rays: 2006.10.10 */
  /* printslack inequality indices  */

  if (Q->printslack) {
    fprintf(lrs_ofp, "\nslack ineq:");
    for (i = lastdv + 1; i <= P->m; i++) {
      if (!zero(P->A[Row[i]][col]))
        fprintf(lrs_ofp, " %ld ", Q->inequality[B[i] - lastdv]);
    }
  }

  return TRUE;
} /* end of lrs_getray */

void getnextoutput(lrs_dic *P, lrs_dat *Q, long i, long col, lrs_mp out)
/* get A[B[i][col] and copy to out */
{
  long row;
  long m = P->m;
  long d = P->d;
  long lastdv = Q->lastdv;
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long j;

  if (i == d && Q->voronoi)
    return; /* skip last column if voronoi set */

  row = Row[i];

  if (Q->nonnegative) /* if m+i basic get correct value from dictionary */
  /* the slack for the inequality m-d+i contains decision    */
  /* variable x_i. We first see if this is in the basis      */
  /* otherwise the value of x_i is zero, except for a ray    */
  /* when it is one (det/det) for the actual column it is in */
  {
    for (j = lastdv + 1; j <= m; j++) {
      if (Q->inequality[B[j] - lastdv] == m - d + i) {
        copy(out, A[Row[j]][col]);
        return;
      }
    }
    /* did not find inequality m-d+i in basis */
    if (i == col)
      copy(out, P->det);
    else
      itomp(ZERO, out);

  } else
    copy(out, A[row][col]);

} /* end of getnextoutput */

void lrs_printcobasis(lrs_dic *P, lrs_dat *Q, long col)
/* col is output column being printed */
{
  char *ss, *sdet, *sin_det, *sz;
  long i;
  long rflag; /* used to find inequality number for ray column */
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  lrs_mp Nvol, Dvol; /* hold rescaled det of current basis */
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long *Row = P->Row;
  long *inequality = Q->inequality;
  long *temparray = Q->temparray;
  long *count = Q->count;
  long hull = Q->hull;
  long d = P->d;
  long lastdv = Q->lastdv;
  long m = P->m;
  long firstime = TRUE;
  long nincidence; /* count number of tight inequalities */
  long len = 0;

  lrs_alloc_mp(Nvol);
  lrs_alloc_mp(Dvol);

  /* convert lrs_mp to char, compute length of string ss and malloc*/

  sdet = cpmp(" det=", P->det);

  rescaledet(P, Q, Nvol, Dvol); /* scales determinant in case input rational */
  sin_det = cprat("in_det=", Nvol, Dvol);

  sz = cprat("z=", P->objnum, P->objden);

  len = snprintf(NULL, 0, "%s%s%s", sdet, sin_det, sz);

  ss = (char *)malloc(len + (d + m) * 20);

  /* build the printcobasis string */

  len = 0;
  if (hull)
    len = len + sprintf(ss + len, "F#%ld B#%ld h=%ld vertices/rays ", count[0],
                        count[2], P->depth);
  else if (Q->voronoi)
    len = len + sprintf(ss + len, "V#%ld R#%ld B#%ld h=%ld data points ",
                        count[1], count[0], count[2], P->depth);
  else
    len = len + sprintf(ss + len, "V#%ld R#%ld B#%ld h=%ld facets ", count[1],
                        count[0], count[2], P->depth);

  rflag = (-1);
  for (i = 0; i < d; i++) {
    temparray[i] = inequality[C[i] - lastdv];
    if (Col[i] == col)
      rflag = temparray[i]; /* look for ray index */
  }
  for (i = 0; i < d; i++)
    reorder(temparray, d);
  for (i = 0; i < d; i++) {
    len = len + sprintf(ss + len, " %ld", temparray[i]);

    if (!(col == ZERO) &&
        (rflag == temparray[i])) /* missing cobasis element for ray */
      len = len + sprintf(ss + len, "*");
  }

  /* get and print incidence information */
  if (col == 0)
    nincidence = d;
  else
    nincidence = d - 1;

  for (i = lastdv + 1; i <= m; i++)
    if (zero(A[Row[i]][0]))
      if ((col == ZERO) || zero(A[Row[i]][col])) {
        nincidence++;
        if (Q->incidence) {
          if (firstime) {
            len = len + sprintf(ss + len, " :");
            firstime = FALSE;
          }
          len = len + sprintf(ss + len, " %ld", inequality[B[i] - lastdv]);
        }
      }

  len = len + sprintf(ss + len, " I#%ld", nincidence);

  sprintf(ss + len, "%s %s %s ", sdet, sin_det, sz);

  if (Q->mplrs)
    lrs_post_output("cobasis", ss);
  else
    fprintf(lrs_ofp, "\n%s", ss);

  free(ss);
  free(sdet);
  free(sin_det);
  free(sz);
  lrs_clear_mp(Nvol);
  lrs_clear_mp(Dvol);

} /* end of lrs_printcobasis */

/*********************/
/* print final totals */
/*********************/
void lrs_printtotals(lrs_dic *P, lrs_dat *Q) {
  static int first_time = 1;
  /* print warnings */

  if (first_time) {
    first_time = 0;

    if (!Q->mplrs)
      fprintf(lrs_ofp, "\nend");

    if (Q->dualdeg) {
      lrs_warning(Q, "finalwarn",
                  "*Warning: Starting dictionary is dual degenerate");
      lrs_warning(Q, "finalwarn",
                  "*Complete enumeration may not have been produced");
      if (Q->maximize)
        lrs_warning(Q, "finalwarn",
                    "*Recommendation: Add dualperturb option before maximize "
                    "in input file\n");
      else
        lrs_warning(Q, "finalwarn",
                    "*Recommendation: Add dualperturb option before minimize "
                    "in input file\n");
    }

    if (Q->unbounded) {
      lrs_warning(Q, "finalwarn",
                  "*Warning: Starting dictionary contains rays");
      lrs_warning(Q, "finalwarn",
                  "*Complete enumeration may not have been produced");
      if (Q->maximize)
        lrs_warning(Q, "finalwarn",
                    "*Recommendation: Change or remove maximize option or add "
                    "bounds\n");
      else
        lrs_warning(Q, "finalwarn",
                    "*Recommendation: Change or remove minimize option or add "
                    "bounds\n");
    }

    if (Q->truncate)
      lrs_warning(Q, "finalwarn", "*Tree truncated at each new vertex");
  }

  if (!Q->hull) {
    if (Q->allbases)
      lrs_warning(Q, "finalwarn",
                  "*Note! Duplicate vertices/rays may be present");
    else if (Q->count[0] > 1 && !Q->homogeneous)
      lrs_warning(Q, "finalwarn", "*Note! Duplicate rays may be present");
  }

  if (Q->mplrs) {
    char *vol;
    if (Q->hull && Q->getvolume) {
      rescalevolume(P, Q, Q->Nvolume, Q->Dvolume);
      vol = cprat("", Q->Nvolume, Q->Dvolume);
      lrs_post_output("volume", vol);
      free(vol);
    }
    return;
  }

  if (!Q->messages)
    return;

  long i;
  double x;
  /* local assignments */
  double *cest = Q->cest;
  long *count = Q->count;
  long *inequality = Q->inequality;
  long *linearity = Q->linearity;
  long *temparray = Q->temparray;

  long *C = P->C;

  long hull = Q->hull;
  long homogeneous = Q->homogeneous;
  long nlinearity = Q->nlinearity;
  long nredundcol = Q->nredundcol;
  long d, lastdv;
  d = P->d;
  lastdv = Q->lastdv;
  if (Q->hull) /* count[1] stores the number of linearities */
    Q->count[1] = Q->nredundcol - Q->homogeneous;

  /* warnings for lrs only */
  if (Q->maxdepth < MAXD)
    fprintf(lrs_ofp, "\n*Tree truncated at depth %lld", Q->maxdepth);
  if (Q->maxcobases > 0L)
    fprintf(lrs_ofp, "\n*maxcobases = %ld", Q->maxcobases - Q->startcount[2]);
  if (Q->maxoutput > 0L)
    fprintf(lrs_ofp, "\n*Maximum number of output lines = %ld", Q->maxoutput);

  /* next block with volume rescaling must come before estimates are printed */

  if (Q->getvolume && Q->runs == 0) {

    rescalevolume(P, Q, Q->Nvolume, Q->Dvolume);

    if (Q->polytope)
      prat("\n*Volume=", Q->Nvolume, Q->Dvolume);
    else
      prat("\n*Pseudovolume=", Q->Nvolume, Q->Dvolume);
  }

  if (hull) /* output things that are specific to hull computation */
  {
    fprintf(lrs_ofp, "\n*Totals: facets=%ld bases=%ld", count[0], count[2]);

    if (count[1] > 0) {
      fprintf(lrs_ofp, " linearities=%ld", count[1]);
      fprintf(lrs_ofp, " facets+linearities=%ld", count[1] + count[0]);
    }
    printf(" max_facet_depth=%ld", count[8]);
    if (lrs_ofp != stdout) {
      printf("\n*Totals: facets=%ld bases=%ld", count[0], count[2]);

      if (count[1] > 0) {
        printf(" linearities=%ld", count[1]);
        printf(" facets+linearities=%ld", count[1] + count[0]);
      }
      printf(" max_facet_depth=%ld", count[8]);
    }

    if (Q->runs > 0) {
      fprintf(lrs_ofp, "\n*Estimates: facets=%.0f bases=%.0f",
              count[0] + cest[0], count[2] + cest[2]);
      if (Q->getvolume) {
        rescalevolume(P, Q, Q->Nvolume, Q->Dvolume);
        rattodouble(Q->Nvolume, Q->Dvolume, &x);
        for (i = 2; i < d; i++)
          cest[3] = cest[3] / i; /*adjust for dimension */
        if (cest[3] == 0)
          prat(" volume=", Q->Nvolume, Q->Dvolume);
        else
          fprintf(lrs_ofp, " volume=%g", cest[3] + x);
      }

      fprintf(lrs_ofp, "\n*Total number of tree nodes evaluated: %ld",
              Q->totalnodes);
#ifndef TIMES
      fprintf(lrs_ofp, "\n*Estimated total running time=%.1f secs ",
              (count[2] + cest[2]) / Q->totalnodes * get_time());
#endif
    }

  } else /* output things specific to vertex/ray computation */
  {
    fprintf(lrs_ofp, "\n*Totals: vertices=%ld rays=%ld bases=%ld", count[1],
            count[0], count[2]);

    fprintf(lrs_ofp, " integer_vertices=%ld  max_vertex_depth=%ld", count[4],
            count[8]);

    if (nredundcol > 0)
      fprintf(lrs_ofp, " linearities=%ld", nredundcol);
    if (count[0] + nredundcol > 0) {
      fprintf(lrs_ofp, " vertices+rays");
      if (nredundcol > 0)
        fprintf(lrs_ofp, "+linearities");
      fprintf(lrs_ofp, "=%ld", nredundcol + count[0] + count[1]);
    }

    if (lrs_ofp != stdout) {
      printf("\n*Totals: vertices=%ld rays=%ld bases=%ld", count[1], count[0],
             count[2]);

      printf(" integer_vertices=%ld max_vertex_depth=%ld", count[4], count[8]);

      if (nredundcol > 0)
        printf(" linearities=%ld", nredundcol);
      if (count[0] + nredundcol > 0) {
        printf(" vertices+rays");
        if (nredundcol > 0)
          printf("+linearities");
        printf("=%ld", nredundcol + count[0] + count[1]);
      }
    } /* end lrs_ofp != stdout */

    if (Q->runs > 0) {
      fprintf(lrs_ofp, "\n*Estimates: vertices=%.0f rays=%.0f",
              count[1] + cest[1], count[0] + cest[0]);
      fprintf(lrs_ofp, " bases=%.0f integer_vertices=%.0f ", count[2] + cest[2],
              count[4] + cest[4]);

      if (Q->getvolume) {
        rattodouble(Q->Nvolume, Q->Dvolume, &x);
        for (i = 2; i <= d - homogeneous; i++)
          cest[3] = cest[3] / i; /*adjust for dimension */
        fprintf(lrs_ofp, " pseudovolume=%g", cest[3] + x);
      }
      fprintf(lrs_ofp, "\n*Total number of tree nodes evaluated: %ld",
              Q->totalnodes);
#ifndef TIMES
      fprintf(lrs_ofp, "\n*Estimated total running time=%.1f secs ",
              (count[2] + cest[2]) / Q->totalnodes * get_time());
#endif
    }

  } /* end of output for vertices/rays */

  fprintf(
      lrs_ofp,
      "\n*Dictionary Cache: max size= %ld misses= %ld/%ld   Tree Depth= %ld",
      dict_count, cache_misses, cache_tries, Q->deepest);
  if (lrs_ofp != stdout)
    printf(
        "\n*Dictionary Cache: max size= %ld misses= %ld/%ld   Tree Depth= %ld",
        dict_count, cache_misses, cache_tries, Q->deepest);

  return;

} /* end of lrs_printtotals */
/************************/
/*  Estimation function */
/************************/
long lrs_estimate(lrs_dic *P, lrs_dat *Q)
/*returns estimate of subtree size (no. cobases) from current node    */
/*current node is not counted.                   */
/*cest[0]rays [1]vertices [2]bases [3]volume     */
/*    [4] integer vertices                       */
{

  lrs_mp_vector
      output;        /* holds one line of output; ray,vertex,facet,linearity */
  lrs_mp Nvol, Dvol; /* hold volume of current basis */
  long estdepth = 0; /* depth of basis/vertex in subtree for estimate */
  long i = 0, j = 0, k, nchild, runcount, col;
  double prod = 0.0;
  double cave[] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double nvertices, nbases, nrays, nvol, nivertices;
  long rays = 0;
  double newvol = 0.0;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *isave = Q->isave;
  long *jsave = Q->jsave;
  double *cest = Q->cest;
  long d = P->d;
  lrs_alloc_mp(Nvol);
  lrs_alloc_mp(Dvol);
  /* Main Loop of Estimator */

  output = lrs_alloc_mp_vector(
      Q->n); /* output holds one line of output from dictionary     */

  for (runcount = 1; runcount <= Q->runs;
       runcount++) { /* runcount counts number of random probes */
    j = 0;
    nchild = 1;
    prod = 1;
    nvertices = 0.0;
    nbases = 0.0;
    nrays = 0.0;
    nvol = 0.0;
    nivertices = 0.0;

    while (nchild != 0) /* while not finished yet */
    {

      nchild = 0;
      while (j < d) {
        if (reverse(P, Q, &i, j)) {
          isave[nchild] = i;
          jsave[nchild] = j;
          nchild++;
        }
        j++;
      }

      if (estdepth == 0 && nchild == 0) {
        cest[0] = cest[0] + rays; /* may be some rays here */
        lrs_clear_mp(Nvol);
        lrs_clear_mp(Dvol);
        lrs_clear_mp_vector(output, Q->n);
        return (0L); /*subtree is a leaf */
      }

      prod = prod * nchild;
      nbases = nbases + prod;

      if (nchild > 0) /*reverse pivot found choose random child */
      {
        k = myrandom(Q->seed, nchild);
        Q->seed = myrandom(Q->seed, 977L);
        i = isave[k];
        j = jsave[k];
        estdepth++;
        Q->totalnodes++; /* calculate total number of nodes evaluated */
        pivot(P, Q, i, j);
        update(P, Q, &i, &j);   /*Update B,C,i,j */
        if (lexmin(P, Q, ZERO)) /* see if lexmin basis for vertex */
        {
          nvertices = nvertices + prod;
          /* integer vertex estimate */
          if (lrs_getvertex(P, Q, output)) {
            --Q->count[1];
            if (one(output[0])) {
              --Q->count[4];
              nivertices = nivertices + prod;
            }
          }
        }

        rays = 0;
        for (col = 1; col <= d; col++)
          if (negative(A[0][col]) && (lrs_ratio(P, Q, col) == 0) &&
              lexmin(P, Q, col))
            rays++;
        nrays = nrays + prod * rays; /* update ray info */

        if (Q->getvolume) {
          rescaledet(P, Q, Nvol,
                     Dvol); /* scales determinant in case input rational */
          rattodouble(Nvol, Dvol, &newvol);
          nvol = nvol + newvol * prod; /* adjusts volume for degree */
        }
        j = 0;
      }
    }
    cave[0] = cave[0] + nrays;
    cave[1] = cave[1] + nvertices;
    cave[2] = cave[2] + nbases;
    cave[3] = cave[3] + nvol;
    cave[4] = cave[4] + nivertices;

    /*  backtrack to root and do it again */

    while (estdepth > 0) {
      estdepth = estdepth - 1;
      selectpivot(P, Q, &i, &j);
      pivot(P, Q, i, j);
      update(P, Q, &i, &j); /*Update B,C,i,j */
      /*fprintf(lrs_ofp,"\n0  +++"); */
      j++;
    }

  } /* end of for loop on runcount */

  j = (long)cave[2] / Q->runs;

  // 2015.2.9   Do not update totals if we do iterative estimating and subtree
  // is too big
  if (Q->subtreesize == 0 || j <= Q->subtreesize)
    for (i = 0; i < 5; i++)
      cest[i] = cave[i] / Q->runs + cest[i];

  lrs_clear_mp(Nvol);
  lrs_clear_mp(Dvol);
  lrs_clear_mp_vector(output, Q->n);
  return (j);
} /* end of lrs_estimate  */

/*********************************/
/* Internal functions            */
/*********************************/
/* Basic Dictionary functions    */
/******************************* */

long reverse(lrs_dic *P, lrs_dat *Q, long *r, long s)
/*  find reverse indices  */
/* TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  long i, j, enter, row, col;

  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long d = P->d;

  enter = C[s];
  col = Col[s];
  if (!negative(A[0][col])) {
    Q->minratio[P->m] = 0; /* 2011.7.14 */
    return (FALSE);
  }

  *r = lrs_ratio(P, Q, col);
  if (*r == 0) /* we have a ray */
  {
    Q->minratio[P->m] = 0; /* 2011.7.14 */
    return (FALSE);
  }

  row = Row[*r];

  /* check cost row after "pivot" for smaller leaving index    */
  /* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
  /* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++)
    if (i != s) {
      j = Col[i];
      if (positive(A[0][j]) ||
          negative(A[row][j])) /*or else sign test fails trivially */
        if ((!negative(A[0][j]) && !positive(A[row][j])) ||
            comprod(A[0][j], A[row][col], A[0][col], A[row][j]) ==
                -1) {            /*+ve cost found */
          Q->minratio[P->m] = 0; /* 2011.7.14 */

          return (FALSE);
        }
    }
  return (TRUE);
} /* end of reverse */

long selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using lexicographic rule   */
/* returns TRUE if pivot found else FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  long j, col;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;

  /*find positive cost coef */
  while ((j < d) && (!positive(A[0][Col[j]])))
    j++;

  if (j < d) /* pivot column found! */
  {
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0)
      return (TRUE); /* unbounded if *r=0 */
  }
  return (FALSE);
} /* end of selectpivot        */
/******************************************************* */

void pivot(lrs_dic *P, lrs_dat *Q, long bas, long cob)
/* Qpivot routine for array A              */
/* indices bas, cob are for Basis B and CoBasis C    */
/* corresponding to row Row[bas] and column       */
/* Col[cob]   respectively                       */
{
  long r, s;
  long i, j;

  /* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long d, m_A;

#ifndef LRSLONG
  lrs_mp Ns, Nt;
  lrs_alloc_mp(Ns);
  lrs_alloc_mp(Nt);
#endif
  lrs_mp Ars;
  lrs_alloc_mp(Ars);

  d = P->d;
  m_A = P->m_A;
  Q->count[3]++; /* count the pivot */

  r = Row[bas];
  s = Col[cob];

  /* Ars=A[r][s]    */
  copy(Ars, A[r][s]);
  storesign(P->det, sign(Ars)); /*adjust determinant to new sign */

  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
        if (j != s) {
          /*          A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */

#ifdef LRSLONG
          qpiv(A[i][j], Ars, A[i][s], A[r][j], P->det);
#else
          mulint(A[i][j], Ars, Nt);
          mulint(A[i][s], A[r][j], Ns);
          decint(Nt, Ns);
          exactdivint(Nt, P->det, A[i][j]);
#endif
        } /* end if j ....  */

  if (sign(Ars) == POS) {
    for (j = 0; j <= d; j++) /* no need to change sign if Ars neg */
      /*   A[r][j]=-A[r][j];              */
      if (!zero(A[r][j]))
        changesign(A[r][j]);
  } /* watch out for above "if" when removing this "}" ! */
  else
    for (i = 0; i <= m_A; i++)
      if (!zero(A[i][s]))
        changesign(A[i][s]);

  /*  A[r][s]=P->det;                  */

  copy(A[r][s], P->det); /* restore old determinant */
  copy(P->det, Ars);
  storesign(P->det, POS); /* always keep positive determinant */

  /* set the new rescaled objective function value */

  mulint(P->det, Q->Lcm[0], P->objden);
  mulint(Q->Gcd[0], A[0][0], P->objnum);

  if (!Q->maximize)
    changesign(P->objnum);
  if (zero(P->objnum))
    storesign(P->objnum, POS);
  else
    reduce(P->objnum, P->objden);

#ifndef LRSLONG
  lrs_clear_mp(Ns);
  lrs_clear_mp(Nt);
#endif
  lrs_clear_mp(Ars);
} /* end of pivot */

long primalfeasible(lrs_dic *P, lrs_dat *Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule
 */
{
  long primalinfeasible = TRUE;
  long i, j;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  /*temporary: try to get new start after linearity */

  while (primalinfeasible) {
    i = lastdv + 1;
    while (i <= m && !negative(A[Row[i]][0]))
      i++;
    if (i <= m) {
      j = 0; /*find a positive entry for in row */
      while (j < d && !positive(A[Row[i]][Col[j]]))
        j++;
      if (j >= d)
        return (FALSE); /* no positive entry */
      pivot(P, Q, i, j);
      update(P, Q, &i, &j);
    } else
      primalinfeasible = FALSE;
  } /* end of while primalinfeasibile */
  return (TRUE);
} /* end of primalfeasible */

long lrs_solvelp(lrs_dic *P, lrs_dat *Q, long maximize)
/* Solve primal feasible lp by Dantzig`s rule and lexicographic ratio test */
/* return TRUE if bounded, FALSE if unbounded                              */
{
  long i, j, k = 0L;
  long notdone = TRUE;
  /* assign local variables to structures */
  long d = P->d;

  /* lponly=1 Dantzig,  =2 random,  =3 hybrid, =4 Bland */

  if (Q->lponly <= 1) /* Dantzig's rule */
    while (dan_selectpivot(P, Q, &i, &j)) {
      pivot(P, Q, i, j);
      update(P, Q, &i, &j); /*Update B,C,i,j */
    }

  if (Q->lponly == 2) /* random edge rule */
    while (ran_selectpivot(P, Q, &i, &j)) {
      pivot(P, Q, i, j);
      update(P, Q, &i, &j); /*Update B,C,i,j */
    }

  if (Q->lponly == 3) /* alternate Dantzig/randome rules */
    while (notdone) {
      if (k % 2) /* odd for dantzig even for random */
        notdone = dan_selectpivot(P, Q, &i, &j);
      else
        notdone = ran_selectpivot(P, Q, &i, &j);

      if (notdone) {
        pivot(P, Q, i, j);
        update(P, Q, &i, &j); /*Update B,C,i,j */
      }
      k++;
    }

  if (Q->lponly == 4) /* Bland's rule - used for vertex enumeration */
    while (selectpivot(P, Q, &i, &j)) {
      pivot(P, Q, i, j);
      update(P, Q, &i, &j); /*Update B,C,i,j */
    }

  if (j < d && i == 0) /* selectpivot gives information on unbounded solution */
  {
    if (Q->lponly && Q->messages)
      fprintf(lrs_ofp, "\n*Unbounded solution");
    return FALSE;
  }
  return TRUE;
} /* end of lrs_solvelp  */

long getabasis(lrs_dic *P, lrs_dat *Q, long order[])

/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */

{
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
  char mess[100];
  m = P->m;
  d = P->d;
  nlinearity = Q->nlinearity;

  for (j = 0l; j < m; j++) {
    i = 0l;
    while (i <= m && B[i] != d + order[j])
      i++;                       /* find leaving basis index i */
    if (j < nlinearity && i > m) /* cannot pivot linearity to cobasis */
    {
      if (Q->messages)
        fprintf(lrs_ofp, "\nCannot find linearity in the basis");
      return FALSE;
    }
    if (i <= m) { /* try to do a pivot */
      k = 0l;
      while (C[k] <= d && zero(A[Row[i]][Col[k]])) {
        k++;
      }
      if (C[k] <= d) {

        pivot(P, Q, i, k);
        update(P, Q, &i, &k);
      } else if (j < nlinearity) { /* cannot pivot linearity to cobasis */
        if (zero(A[Row[i]][0])) {
          if (Q->messages && overflow != 2) {
            sprintf(mess,
                    "*Input linearity in row %ld is redundant--converted to "
                    "inequality",
                    order[j]);
            lrs_warning(Q, "warning", mess);
          }
          linearity[j] = 0l;
          Q->redineq[j] = 1; /* check for redundancy if running redund */
        } else {
          lrs_warning(Q, "warning", "*No feasible solution");
          return FALSE;
        }
      }
    }
  }

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
  /* bug fix, 2009.6.27 */ Q->nlinearity = i;

  /* column dependencies now can be recorded  */
  /* redundcol contains input column number 0..n-1 where redundancy is */
  k = 0;
  while (k < d && C[k] <= d) {
    if (C[k] <= d) { /* decision variable still in cobasis */
      redundcol[nredundcol++] = C[k] - Q->hull; /* adjust for hull indices */
    }
    k++;
  }

  /* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;

  /* Remove linearities from cobasis for rest of computation */
  /* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++) { /* find cobasic index */
    k = 0;
    while (k < d && C[k] != linearity[i] + d)
      k++;
    if (k >= d) {
      lrs_warning(Q, "warning", "\nError removing linearity");
      return FALSE;
    }
    if (!removecobasicindex(P, Q, k))
      return FALSE;
    d = P->d;
  }

  /* Check feasability */
  if (Q->givenstart) {
    i = Q->lastdv + 1;
    while (i <= m && !negative(A[Row[i]][0]))
      i++;
    if (i <= m)
      fprintf(lrs_ofp, "\n*Infeasible startingcobasis - will be modified");
  }
  return TRUE;
} /*  end of getabasis */

long removecobasicindex(lrs_dic *P, lrs_dat *Q, long k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  long i, j, cindex, deloc;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  cindex = C[k];  /* cobasic index to remove              */
  deloc = Col[k]; /* matrix column location to remove     */

  for (i = 1; i <= m; i++) /* reduce basic indices by 1 after index */
    if (B[i] > cindex)
      B[i]--;

  for (j = k; j < d; j++) /* move down other cobasic variables    */
  {
    C[j] = C[j + 1] - 1; /* cobasic index reduced by 1           */
    Col[j] = Col[j + 1];
  }

  if (deloc != d) {
    /* copy col d to deloc */
    for (i = 0; i <= m; i++)
      copy(A[i][deloc], A[i][d]);

    /* reassign location for moved column */
    j = 0;
    while (Col[j] != d)
      j++;

    Col[j] = deloc;
  }

  P->d--;
  return TRUE;
} /* end of removecobasicindex */

lrs_dic *resize(lrs_dic *P, lrs_dat *Q)
/* resize the dictionary after some columns are deleted, ie. inputd>d */
/* a new lrs_dic record is created with reduced size, and items copied over */
{
  lrs_dic *P1; /* to hold new dictionary in case of resizing */

  long i, j;
  long m, d, m_A;

  m = P->m;
  d = P->d;
  m_A = P->m_A;

  /* get new dictionary record */

  P1 = new_lrs_dic(m, d, m_A);

  /* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P1->d_orig = d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  copy(P1->det, P->det);
  copy(P1->objnum, P->objnum);
  copy(P1->objden, P->objden);

  for (i = 0; i <= m; i++) {
    P1->B[i] = P->B[i];
    P1->Row[i] = P->Row[i];
  }
  for (i = 0; i <= m_A; i++) {
    for (j = 0; j <= d; j++)
      copy(P1->A[i][j], P->A[i][j]);
  }

  for (j = 0; j <= d; j++) {
    P1->Col[j] = P->Col[j];
    P1->C[j] = P->C[j];
  }

  lrs_free_dic(P, Q);

  /* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;

  return P1;
}
/********* resize                    ***************/

long restartpivots(lrs_dic *P, lrs_dat *Q)
/* facet contains a list of the inequalities in the cobasis for the restart */
/* inequality contains the relabelled inequalities after initialization     */
{
  long i, j, k;
  long *Cobasic; /* when restarting, Cobasic[j]=1 if j is in cobasis */
                 /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *inequality = Q->inequality;
  long *facet = Q->facet;
  long nlinearity = Q->nlinearity;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  Cobasic = (long *)CALLOC((unsigned)m + d + 2, sizeof(long));

  /* set Cobasic flags */
  for (i = 0; i < m + d + 1; i++)
    Cobasic[i] = 0;
  for (i = 0; i < d; i++) /* find index corresponding to facet[i] */
  {
    j = 1;
    while (facet[i + nlinearity] != inequality[j])
      j++;
    Cobasic[j + lastdv] = 1;
  }

  /* Note that the order of doing the pivots is important, as */
  /* the B and C vectors are reordered after each pivot       */

  /* Suggested new code from db starts */
  i = m;
  while (i > d) {
    while (Cobasic[B[i]]) {
      k = d - 1;
      while ((k >= 0) && (zero(A[Row[i]][Col[k]]) || Cobasic[C[k]])) {
        k--;
      }
      if (k >= 0) {
        /*db asks: should i really be modified here? (see old code) */
        /*da replies: modifying i only makes is larger, and so      */
        /*the second while loop will put it back where it was       */
        /*faster (and safer) as done below                          */
        long ii = i;
        pivot(P, Q, ii, k);
        update(P, Q, &ii, &k);
      } else {
        lrs_warning(Q, "warning",
                    "\nInvalid Co-basis - does not have correct rank");
        free(Cobasic);
        return FALSE;
      }
    }
    i--;
  }
  /* Suggested new code from db ends */

  /* check restarting from a primal feasible dictionary               */
  for (i = lastdv + 1; i <= m; i++)
    if (negative(A[Row[i]][0])) {
      lrs_warning(Q, "warning",
                  "\nTrying to restart from infeasible dictionary");
      free(Cobasic);
      return FALSE;
    }
  free(Cobasic);
  return TRUE;

} /* end of restartpivots */

long lrs_ratio(lrs_dic *P, lrs_dat *Q, long col) /*find lex min. ratio */
/* find min index ratio -aig/ais, ais<0 */
/* if multiple, checks successive basis columns */
/* recoded Dec 1997                     */
{
  long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  long firstime; /*For ratio test, true on first pass,else false */
  lrs_mp Nmin, Dmin;
  long degencount, ndegencount;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long *Col = P->Col;
  long *minratio = Q->minratio;
  long m, d, lastdv;

  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  nstart = 0;
  ndegencount = 0;
  degencount = 0;
  minratio[P->m] = 1; /*2011.7.14 non-degenerate pivot flag */

  for (j = lastdv + 1; j <= m; j++) {
    /* search rows with negative coefficient in dictionary */
    /*  minratio contains indices of min ratio cols        */
    if (negative(A[Row[j]][col])) {
      minratio[degencount++] = j;
      if (zero(A[Row[j]][0]))
        minratio[P->m] = 0; /*2011.7.14 degenerate pivot flag */
    }
  } /* end of for loop */
  if (degencount == 0)
    return (degencount); /* non-negative pivot column */

  lrs_alloc_mp(Nmin);
  lrs_alloc_mp(Dmin);
  ratiocol = 0;   /* column being checked, initially rhs */
  start = 0;      /* starting location in minratio array */
  bindex = d + 1; /* index of next basic variable to consider */
  cindex = 0;     /* index of next cobasic variable to consider */
  basicindex =
      d; /* index of basis inverse for current ratio test, except d=rhs test */
  while (degencount > 1) /*keep going until unique min ratio found */
  {
    if (B[bindex] == basicindex) /* identity col in basis inverse */
    {
      if (minratio[start] == bindex)
      /* remove this index, all others stay */
      {
        start++;
        degencount--;
      }
      bindex++;
    } else
    /* perform ratio test on rhs or column of basis inverse */
    {
      firstime = TRUE;
      /*get next ratio column and increment cindex */
      if (basicindex != d)
        ratiocol = Col[cindex++];
      for (j = start; j < start + degencount; j++) {
        i = Row[minratio[j]]; /* i is the row location of the next basic
                                 variable */
        comp = 1;             /* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
        if (firstime)
          firstime = FALSE; /*force new min ratio on first time */
        else {
          if (positive(Nmin) || negative(A[i][ratiocol])) {
            if (negative(Nmin) || positive(A[i][ratiocol]))
              comp = comprod(Nmin, A[i][col], A[i][ratiocol], Dmin);
            else
              comp = -1;
          }

          else if (zero(Nmin) && zero(A[i][ratiocol]))
            comp = 0;

          if (ratiocol == ZERO)
            comp = -comp; /* all signs reversed for rhs */
        }
        if (comp == 1) { /*new minimum ratio */
          nstart = j;
          copy(Nmin, A[i][ratiocol]);
          copy(Dmin, A[i][col]);
          ndegencount = 1;
        } else if (comp == 0) /* repeated minimum */
          minratio[nstart + ndegencount++] = minratio[j];

      } /* end of  for (j=start.... */
      degencount = ndegencount;
      start = nstart;
    } /* end of else perform ratio test statement */
    basicindex++; /* increment column of basis inverse to check next */
  } /*end of while loop */
  lrs_clear_mp(Nmin);
  lrs_clear_mp(Dmin);
  return (minratio[start]);
} /* end of ratio */

long lexmin(lrs_dic *P, lrs_dat *Q, long col)
/*test if basis is lex-min for vertex or ray, if so TRUE */
/* FALSE if a_r,g=0, a_rs !=0, r > s          */
{
  /*do lexmin test for vertex if col=0, otherwise for ray */
  long r, s, i, j;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++) {
    r = Row[i];
    if (zero(A[r][col])) /* necessary for lexmin to fail */
      for (j = 0; j < d; j++) {
        s = Col[j];
        if (B[i] > C[j]) /* possible pivot to reduce basis */
        {
          if (zero(A[r][0])) /* no need for ratio test, any pivot feasible */
          {
            if (!zero(A[r][s]))
              return (FALSE);
          } else if (negative(A[r][s]) && ismin(P, Q, r, s)) {
            return (FALSE);
          }
        } /* end of if B[i] ... */
      }
  }
  return (TRUE);
} /* end of lexmin */

long ismin(lrs_dic *P, lrs_dat *Q, long r, long s)
/*test if A[r][s] is a min ratio for col s */
{
  long i;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if ((i != r) && negative(A[i][s]) &&
        comprod(A[i][0], A[r][s], A[i][s], A[r][0])) {
      return (FALSE);
    }

  return (TRUE);
}

void update(lrs_dic *P, lrs_dat *Q, long *i, long *j)
/*update the B,C arrays after a pivot */
/*   involving B[bas] and C[cob]           */
{

  long leave, enter;
  /* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  leave = B[*i];
  enter = C[*j];
  B[*i] = enter;
  reorder1(B, Row, *i, m + 1);
  C[*j] = leave;
  reorder1(C, Col, *j, d);
  /* restore i and j to new positions in basis */
  for (*i = 1; B[*i] != enter; (*i)++)
    ; /*Find basis index */
  for (*j = 0; C[*j] != leave; (*j)++)
    ; /*Find co-basis index */
} /* end of update */

long lrs_degenerate(lrs_dic *P, lrs_dat *Q)
/* TRUE if the current dictionary is primal degenerate */
/* not thoroughly tested   2000/02/15                  */
{
  long i;
  long *Row;

  lrs_mp_matrix A = P->A;
  long d = P->d;
  long m = P->m;

  Row = P->Row;

  for (i = d + 1; i <= m; i++)
    if (zero(A[Row[i]][0]))
      return TRUE;

  return FALSE;
}

/*********************************************************/
/*                 Miscellaneous                         */
/******************************************************* */

void reorder(long a[], long range)
/*reorder array in increasing order with one misplaced element */
{
  long i, temp;
  for (i = 0; i < range - 1; i++)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }
  for (i = range - 2; i >= 0; i--)
    if (a[i] > a[i + 1]) {
      temp = a[i];
      a[i] = a[i + 1];
      a[i + 1] = temp;
    }

} /* end of reorder */

void reorder1(long a[], long b[], long newone, long range)
/*reorder array a in increasing order with one misplaced element at index newone
 */
/*elements of array b are updated to stay aligned with a */
{
  long temp;
  while (newone > 0 && a[newone] < a[newone - 1]) {
    temp = a[newone];
    a[newone] = a[newone - 1];
    a[newone - 1] = temp;
    temp = b[newone];
    b[newone] = b[newone - 1];
    b[--newone] = temp;
  }
  while (newone < range - 1 && a[newone] > a[newone + 1]) {
    temp = a[newone];
    a[newone] = a[newone + 1];
    a[newone + 1] = temp;
    temp = b[newone];
    b[newone] = b[newone + 1];
    b[++newone] = temp;
  }
} /* end of reorder1 */

void rescaledet(lrs_dic *P, lrs_dat *Q, lrs_mp Vnum, lrs_mp Vden)
/* rescale determinant to get its volume */
/* Vnum/Vden is volume of current basis  */
{
  lrs_mp gcdprod; /* to hold scale factors */
  long i;
  /* assign local variables to structures */
  long *C = P->C;
  long *B = P->B;
  long m, d, lastdv;

  lrs_alloc_mp(gcdprod);
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  itomp(ONE, gcdprod);
  itomp(ONE, Vden);

  for (i = 0; i < d; i++)
    if (B[i] <= m) {
      mulint(Q->Gcd[Q->inequality[C[i] - lastdv]], gcdprod, gcdprod);
      mulint(Q->Lcm[Q->inequality[C[i] - lastdv]], Vden, Vden);
    }
  mulint(P->det, gcdprod, Vnum);
  //  reduce (Vnum, Vden);
  lrs_clear_mp(gcdprod);
} /* end rescaledet */

void rescalevolume(lrs_dic *P, lrs_dat *Q, lrs_mp Vnum, lrs_mp Vden)
/* adjust volume for dimension */
{
  lrs_mp temp, dfactorial;
  /* assign local variables to structures */
  long lastdv = Q->lastdv;

  lrs_alloc_mp(temp);
  lrs_alloc_mp(dfactorial);

  /*reduce Vnum by d factorial  */
  getfactorial(dfactorial, lastdv);
  mulint(dfactorial, Vden, Vden);
  if (Q->hull && !Q->homogeneous) { /* For hull option multiply by d to correct
                                       for lifting */
    itomp(lastdv, temp);
    mulint(temp, Vnum, Vnum);
  }

  reduce(Vnum, Vden);
  lrs_clear_mp(temp);
  lrs_clear_mp(dfactorial);
}

void updatevolume(lrs_dic *P,
                  lrs_dat *Q) /* rescale determinant and update the volume */
{
  lrs_mp tN, tD, Vnum, Vden;
  lrs_alloc_mp(tN);
  lrs_alloc_mp(tD);
  lrs_alloc_mp(Vnum);
  lrs_alloc_mp(Vden);
  rescaledet(P, Q, Vnum, Vden);
  copy(tN, Q->Nvolume);
  copy(tD, Q->Dvolume);
  linrat(tN, tD, ONE, Vnum, Vden, ONE, Q->Nvolume, Q->Dvolume);
  lrs_clear_mp(tN);
  lrs_clear_mp(tD);
  lrs_clear_mp(Vnum);
  lrs_clear_mp(Vden);

} /* end of updatevolume */

/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/

long checkredund(lrs_dic *P, lrs_dat *Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* return 0=nonredundant -1=strict redundant(interior) 1=non-strict redundant*/
{
  lrs_mp Ns, Nt;
  long i, j;
  long r, s;

  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row, *Col;
  long d = P->d;

  lrs_alloc_mp(Ns);
  lrs_alloc_mp(Nt);
  Row = P->Row;
  Col = P->Col;
  while (selectpivot(P, Q, &i, &j)) {
    Q->count[2]++;

    /* sign of new value of A[0][0]            */
    /* is      A[0][s]*A[r][0]-A[0][0]*A[r][s] */

    r = Row[i];
    s = Col[j];

    mulint(A[0][s], A[r][0], Ns);
    mulint(A[0][0], A[r][s], Nt);

    if (mp_greater(Ns, Nt)) {
      lrs_clear_mp(Ns);
      lrs_clear_mp(Nt);
      return 0; /* non-redundant */
    }

    pivot(P, Q, i, j);
    update(P, Q, &i, &j); /*Update B,C,i,j */
  }
  lrs_clear_mp(Ns);
  lrs_clear_mp(Nt);

  if (j < d && i == 0) /* unbounded is also non-redundant */
    return 0;

  /* 2020.6.8 check for strict redundancy and return -1 if so */

  if (negative(P->A[0][0]))
    return -1;
  else
    return 1;

} /* end of checkredund  */

long checkcobasic(lrs_dic *P, lrs_dat *Q, long index)
/* TRUE if index is cobasic and nonredundant                         */
/* FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

  /* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *B, *C, *Row, *Col;
  long d = P->d;
  long m = P->m;
  long i = 0;
  long j = 0;
  long s;

  B = P->B;
  C = P->C;
  Row = P->Row;
  Col = P->Col;

  while ((j < d) && C[j] != index)
    j++;

  if (j == d)
    return FALSE; /* not cobasic index */

  /* index is cobasic */

  s = Col[j];
  i = Q->lastdv + 1;

  while ((i <= m) && (zero(A[Row[i]][s]) || !zero(A[Row[i]][0])))
    i++;

  if (i > m) {
    return TRUE;
  }

  pivot(P, Q, i, j);
  update(P, Q, &i, &j); /*Update B,C,i,j */

  return FALSE; /*index is no longer cobasic */

} /* end of checkcobasic */

long checkindex(lrs_dic *P, lrs_dat *Q, long index)
/* 0 if index is non-redundant inequality    */
/*-1 if index is strict redundant inequality */
/* 1 if index is non-strict redundant ine    */
/* 2 if index is input linearity             */
/*NOTE: row is returned all zero if redundant!! */
{
  long i, j;

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *B = P->B;
  long d = P->d;
  long m = P->m;
  long zeroonly = 0;

  if (index <
      0) /* used to zero out known redundant rows in mplrs verifyredund */
  {
    zeroonly = 1;
    index = -index;
  }

  /* each slack index must be checked for redundancy */
  /* if in cobasis, it is pivoted out if degenerate */
  /* else it is non-redundant                       */

  if (checkcobasic(P, Q, index)) {
    return ZERO;
  }
  /* index is basic   */
  j = 1;
  while ((j <= m) && (B[j] != index))
    j++;

  i = Row[j];

  /* copy row i to cost row, and set it to zero */

  for (j = 0; j <= d; j++) {
    copy(A[0][j], A[i][j]);
    changesign(A[0][j]);
    itomp(ZERO, A[i][j]);
  }
  if (zeroonly)
    return 1;

  /*2020.6.6 new test for strict redundancy */

  j = checkredund(P, Q);
  if (j != 0)
    return j;

  /* non-redundant, copy back and change sign */

  for (j = 0; j <= d; j++) {
    copy(A[i][j], A[0][j]);
    changesign(A[i][j]);
  }

  return 0;

} /* end of checkindex */

/***************************************************************/
/*                                                             */
/* f I/O routines                          */
/*                                                             */
/***************************************************************/

void lprat(const char *name, long Nt, long Dt)
/*print the long precision rational Nt/Dt without reducing  */
{
  if (Nt > 0)
    fprintf(lrs_ofp, " ");
  fprintf(lrs_ofp, "%s%ld", name, Nt);
  if (Dt != 1)
    fprintf(lrs_ofp, "/%ld", Dt);
  fprintf(lrs_ofp, " ");
} /* lprat */

long lreadrat(long *Num, long *Den)
/* read a rational string and convert to long    */
/* returns true if denominator is not one        */
{
  char in[MAXINPUT], num[MAXINPUT], den[MAXINPUT];
  if (fscanf(lrs_ifp, "%s", in) == EOF)
    return (FALSE);
  atoaa(in, num, den); /*convert rational to num/dem strings */
  *Num = atol(num);
  if (den[0] == '\0') {
    *Den = 1L;
    return (FALSE);
  }
  *Den = atol(den);
  return (TRUE);
}

void lrs_getinput(lrs_dic *P, lrs_dat *Q, long *num, long *den, long m, long d)
/* code for reading data matrix in lrs/cdd format */
{
  long j, row;

  printf("\nEnter each row: b_i  a_ij j=1..%ld", d);
  for (row = 1; row <= m; row++) {
    printf("\nEnter row %ld: ", row);
    for (j = 0; j <= d; j++) {
      lreadrat(&num[j], &den[j]);
      lprat(" ", num[j], den[j]);
    }

    lrs_set_row(P, Q, row, num, den, GE);
  }

  printf("\nEnter objective row c_j j=1..%ld: ", d);
  num[0] = 0;
  den[0] = 1;
  for (j = 1; j <= d; j++) {
    lreadrat(&num[j], &den[j]);
    lprat(" ", num[j], den[j]);
  }

  lrs_set_obj(P, Q, num, den, MAXIMIZE);
}

long readlinearity(lrs_dat *Q) /* read in and check linearity list */
{
  long i, j;
  long nlinearity;
  if (fscanf(lrs_ifp, "%ld", &nlinearity) == EOF) {
    lrs_warning(Q, "warning", "\nLinearity option invalid, no indices ");
    return (FALSE);
  }
  if (nlinearity < 1) {
    lrs_warning(Q, "warning",
                "\nLinearity option invalid, indices must be positive");
    return (FALSE);
  }

  Q->linearity = (long int *)CALLOC((nlinearity + 1), sizeof(long));

  for (i = 0; i < nlinearity; i++) {
    if (fscanf(lrs_ifp, "%ld", &j) == EOF) {
      lrs_warning(Q, "warning", "\nLinearity option invalid, missing indices");
      return (FALSE);
    }
    Q->linearity[i] = j;
  }
  for (i = 1; i < nlinearity; i++) /*sort in order */
    reorder(Q->linearity, nlinearity);

  Q->nlinearity = nlinearity;
  Q->polytope = FALSE;
  return TRUE;
} /* end readlinearity */

long readredund(lrs_dat *Q) /* read in and check linearity list */
{
  long i, j, k;
  char *mess;
  int len = 0;

  if (fscanf(lrs_ifp, "%ld", &k) == EOF) {
    lrs_warning(Q, "warning", "\nredund_list option invalid: no indices ");
    return (FALSE);
  }
  if (k < 0) {
    lrs_warning(Q, "warning",
                "\nredund_list option invalid, first index must be >= 0");
    return (FALSE);
  }

  for (i = 1; i <= Q->m;
       i++) /*reset any previous redund option except =2 values */
    if (Q->redineq[i] != 2)
      Q->redineq[i] = 0;
  Q->redineq[0] = 1;

  for (i = 0; i < k; i++) {
    if (fscanf(lrs_ifp, "%ld", &j) == EOF) {
      lrs_warning(Q, "warning",
                  "\nredund_list option invalid: missing indices");
      fflush(lrs_ofp);
      return (FALSE);
    }

    if (j < 0 || j > Q->m) {
      fprintf(lrs_ofp,
              "\nredund_list option invalid: indices not between 1 and %ld",
              Q->m);
      return (FALSE);
    }
    Q->redineq[j] = 1;
  }

  if (Q->messages && overflow != 2) {
    mess = (char *)malloc(20 * Q->m * sizeof(char));
    len = sprintf(mess, "redund_list %ld ", k);
    for (i = 1; i <= Q->m; i++)
      if (Q->redineq[i] == 1)
        len = len + sprintf(mess + len, " %ld", i);
    lrs_warning(Q, "warning", mess);
    free(mess);
  }
  return TRUE;
} /* end readredund */

long readfacets(lrs_dat *Q, long facet[])
/* read and check facet list for obvious errors during start/restart */
/* this must be done after linearity option is processed!!           */
{
  long i, j;
  char str[1000000], *p, *e;

  /* assign local variables to structures */
  long m, d;
  long *linearity = Q->linearity;
  m = Q->m;
  d = Q->inputd;

  /* modified 2018.6.7 to fix bug restarting with less than full dimension input
   */
  /* number of restart indices is not known at this point */

  j = Q->nlinearity; /* note we place these after the linearity indices */

  if (fgets(str, 1000000, lrs_ifp) == NULL)
    return FALSE; /* pick up indices from the input line */
  for (p = str;; p = e) {
    facet[j] = strtol(p, &e, 10);
    if (p == e)
      break;
    if (!Q->mplrs && overflow != 2)
      fprintf(lrs_ofp, " %ld", facet[j]);

    /* 2010.4.26 nonnegative option needs larger range of indices */
    if (Q->nonnegative)
      if (facet[j] < 1 || facet[j] > m + d) {
        fprintf(lrs_ofp,
                "\n Start/Restart cobasic indices must be in range 1 .. %ld ",
                m + d);
        return FALSE;
      }

    if (!Q->nonnegative)
      if (facet[j] < 1 || facet[j] > m) {
        fprintf(lrs_ofp,
                "\n Start/Restart cobasic indices must be in range 1 .. %ld ",
                m);
        return FALSE;
      }

    for (i = 0; i < Q->nlinearity; i++)
      if (linearity[i] == facet[j]) {
        fprintf(
            lrs_ofp,
            "\n Start/Restart cobasic indices should not include linearities");
        return FALSE;
      }
    /*     bug fix 2011.8.1  reported by Steven Wu*/
    for (i = Q->nlinearity; i < j; i++)
      /*     end bug fix 2011.8.1 */

      if (facet[i] == facet[j]) {
        fprintf(lrs_ofp, "\n  Start/Restart cobasic indices must be distinct");
        return FALSE;
      }
    j++;
  }
  return TRUE;
} /* end of readfacets */

long readvars(lrs_dat *Q, char *name) {
  /* read in and check ordered list of vars for extract/project        */
  /* extract mode: *vars is an ordered list of variables to be kept    */
  /* fel mode:     *vars is an ordered list of variables to be removed */

  long i, j, len;
  long nvars, nremove;
  long k = 0;

  long *vars;
  char *mess;
  long *var; /* binary representation of vars */

  long n = Q->n;
  Q->vars = (long int *)CALLOC((n + 3), sizeof(long));
  var = (long int *)CALLOC((n + 3), sizeof(long));

  vars = Q->vars;
  for (i = 0; i <= n + 2; i++) {
    vars[i] = 0;
    var[i] = 0;
  }

  if (fscanf(lrs_ifp, "%ld", &nvars) == EOF) {
    fprintf(lrs_ofp, "\n*%s: missing indices\n", name);
    free(var);
    return FALSE;
  }

  if (nvars > n - 1) {
    nvars = n - 1;
    fprintf(lrs_ofp, "\n*%s: too many indices, first %ld taken", name, n - 1);
  }

  for (i = 0; i < nvars; i++) {
    if (fscanf(lrs_ifp, "%ld", &j) == EOF) {
      fprintf(lrs_ofp, "\n*%s: missing indices\n", name);
      free(var);
      return FALSE;
    }
    if (j > 0 && j < n) {
      if (var[j] == 1)
        fprintf(lrs_ofp, "\n*%s: duplicate index %ld skipped", name, j);
      else {
        vars[k++] = j;
        var[j] = 1;
      }
    } else {
      fprintf(lrs_ofp, "\n*%s: index %ld out of range 1 to %ld\n", name, j,
              n - 1);
      free(var);
      return FALSE;
    }
  }

  i = 0;
  while (i < n && vars[i] != 0)
    i++;
  nvars = i;

  vars[n + 1] = nvars;

  if (Q->messages && overflow != 2) {
    mess = (char *)malloc(20 * Q->n * sizeof(char));
    len = sprintf(mess, "*%s %ld  ", name, nvars);
    for (i = 0; i < nvars; i++)
      len = len + sprintf(mess + len, "%ld ", vars[i]);
    lrs_warning(Q, "warning", mess);
    free(mess);
  }

  if (strcmp(name, "project") == 0) /* convert to project vars to remove vars */
  {
    for (i = 0; i < nvars; i++)
      vars[i] = 0;
    nremove = 0;
    for (i = 1; i < n; i++)
      if (!var[i])
        vars[nremove++] = i;
    vars[n + 1] = nremove;
    vars[n] =
        1; /* used to control column selection rule =1 min cols in new matrix */
  } /* convert project vars */

  free(var);

  if (Q->fel)
    return TRUE;

  /* if nlinearity>0 fill up list with remaining decision variables */

  if (!Q->hull && Q->nlinearity > 0)
    for (i = 1; i < n; i++) {
      j = 0;
      while (j < nvars && vars[j] != i)
        j++;
      if (j == nvars)
        vars[nvars++] = i;
    }
  return TRUE;
} /* readvars */

long extractcols(lrs_dic *P, lrs_dat *Q) {
  /* 2020.6.17*/
  /* extract option just pulls out the columns - extract mode */
  /* project option also removes redundant rows -  fel mode   */

  long i, j, m, n;
  long ncols;
  long rows;
  lrs_mp_matrix A;
  long *Col, *Row, *remain, *output, *redineq;

  lrs_dic *P1;

  Col = P->Col;
  Row = P->Row;
  remain = Q->vars;
  output = Q->temparray;
  m = P->m;
  n = Q->n;
  if (Q->fel)
    ncols = n - remain[n + 1] - 1;
  else
    ncols = remain[n + 1];

  for (j = 0; j < n; j++)
    output[j] = 0;

  for (j = 0; j < n; j++)
    output[remain[j]] = 1;

  if (Q->fel) /* complement for fel mode - don't ask! */
    for (j = 1; j < n; j++)
      output[j] = 1 - output[j];

  if (Q->fel) /* fel mode remove redundancy */
  {
    for (i = 1; i <= m; i++) /* zero out removed cols */
      for (j = 0; j < n; j++)
        if (!output[j])
          itomp(ZERO, P->A[Row[i]][Col[j]]);

    P1 = lrs_getdic(Q);
    Q->Qhead = P;
    Q->Qtail = P;

    copy_dict(Q, P1, P);
    Q->Qhead = P1;
    Q->Qtail = P1;
    Q->olddic = P; /*in case of overflow */
    A = P1->A;

    redund_run(P1, Q);
    redineq = Q->redineq;
    rows = 0;
    for (i = 1; i <= P->m_A; i++)
      if (redineq[i] == 0 || redineq[i] == 2)
        rows++;

    Q->Qhead = P;
    Q->Qtail = P;

  } /* end   if(Q->fel)...             */
  else /* initialization for extract mode */
  {
    redineq = Q->redineq;
    rows = m;
    for (i = 1; i <= m; i++)
      redineq[i] = 0;
  }

  A = P->A;
  m = Q->m;

  if (Q->hull)
    fprintf(lrs_ofp, "\nV-representation");
  else
    fprintf(lrs_ofp, "\nH-representation");

  if (Q->nlinearity > 0) {
    fprintf(lrs_ofp, "\nlinearity %ld", Q->nlinearity);
    for (i = 0; i < Q->nlinearity; i++)
      fprintf(lrs_ofp, " %ld", Q->linearity[i]);
  }

  fprintf(lrs_ofp, "\nbegin\n%ld %ld rational", rows, ncols + 1);
  for (i = 1; i <= m; i++) {
    if (redineq[i] != 1) {
      reducearray(A[Row[i]], n + Q->hull); /*we already decremented n */
      fprintf(lrs_ofp, "\n");
      if (Q->hull) {
        for (j = 0; j < n; j++)
          if (output[j]) {
            if (zero(A[Row[i]][Col[0]]))
              pmp("", A[Row[i]][Col[j]]);
            else
              prat("", A[Row[i]][Col[j]], A[Row[i]][Col[0]]);
          }
      } else /* no lifting */
      {
        pmp("", A[Row[i]][0]);
        for (j = 1; j < n; j++)
          if (output[j])
            pmp("", A[Row[i]][Col[j - 1]]);
      }
    }
  }
  fprintf(lrs_ofp, "\nend");

  fprintf(lrs_ofp, "\n*columns retained:");
  for (j = 0; j < n; j++)
    if (output[j])
      fprintf(lrs_ofp, " %ld", j);

  fprintf(lrs_ofp, "\n");

  return 0;
} /* extractcols */

long linextractcols(lrs_dic *P, lrs_dat *Q)
/* 2020.2.2 */
/* extract option to output the reduced A matrix after linearities are removed
 */
/* should be followed by redund to get minimum representation */
{
  long d, i, j, k, m, n;
  long ii, jj;
  long nlinearity = Q->nlinearity;
  lrs_mp_matrix A;
  long *B, *C, *Col, *Row, *remain;

  A = P->A;
  B = P->B;
  C = P->C;
  Col = P->Col;
  Row = P->Row;
  remain = Q->vars;

  m = P->m;
  n = Q->n;
  d = Q->inputd;

  fprintf(lrs_ofp, "\n*extract col order: ");

  for (j = 0; j < n - 1; j++)
    fprintf(lrs_ofp, "%ld ", remain[j]);

  for (k = 0; k < n - 1; k++) /* go through input order for vars to remain */
  {
    i = 1;
    while (i <= m) {
      if (B[i] == remain[k]) {
        j = 0;
        while (j + nlinearity < d && (C[j] <= d || zero(A[Row[i]][Col[j]])))
          j++;
        if (j + nlinearity < d) {
          ii = i;
          jj = j;
          pivot(P, Q, ii, jj);
          update(P, Q, &ii, &jj);
          i = 0;
        }
      } /* if B[i]...    */
      i++;
    } /* while i <= m  */
  } /* for k=0       */

  if (Q->hull)
    fprintf(lrs_ofp, "\n*columns retained:");
  else
    fprintf(lrs_ofp, "\n*columns retained: 0");
  for (j = 0; j < d - nlinearity; j++)
    fprintf(lrs_ofp, " %ld", C[j] - Q->hull);

  if (Q->hull)
    fprintf(lrs_ofp, "\nV-representation\nbegin");
  else
    fprintf(lrs_ofp, "\nH-representation\nbegin");
  fprintf(lrs_ofp, "\n%ld %ld rational", m - nlinearity,
          d - nlinearity + 1 - Q->hull);

  for (i = nlinearity + 1; i <= m; i++) {
    reducearray(A[Row[i]], n - nlinearity);
    fprintf(lrs_ofp, "\n");
    if (!Q->hull)
      pmp("", A[Row[i]][0]);
    for (j = 0; j < d - nlinearity; j++)
      pmp("", A[Row[i]][Col[j]]);
  }
  fprintf(lrs_ofp, "\nend");
  if (Q->hull)
    fprintf(lrs_ofp, "\n*columns retained:");
  else
    fprintf(lrs_ofp, "\n*columns retained: 0");
  for (j = 0; j < d - nlinearity; j++)
    fprintf(lrs_ofp, " %ld", C[j] - Q->hull);
  fprintf(lrs_ofp, "\n");

  return 0;
} /* linextractcols  */

void printA(lrs_dic *P, lrs_dat *Q) /* print the integer m by n array A
                                       with B,C,Row,Col vectors         */
{
  long i, j;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  fprintf(lrs_ofp, "\n Basis    ");
  for (i = 0; i <= m; i++)
    fprintf(lrs_ofp, "%ld ", B[i]);
  fprintf(lrs_ofp, " Row ");
  for (i = 0; i <= m; i++)
    fprintf(lrs_ofp, "%ld ", Row[i]);
  fprintf(lrs_ofp, "\n Co-Basis ");
  for (i = 0; i <= d; i++)
    fprintf(lrs_ofp, "%ld ", C[i]);
  fprintf(lrs_ofp, " Column ");
  for (i = 0; i <= d; i++)
    fprintf(lrs_ofp, "%ld ", Col[i]);
  pmp(" det=", P->det);
  fprintf(lrs_ofp, "\n");
  i = 0;
  while (i <= m) {
    for (j = 0; j <= d; j++)
      pimat(P, i, j, A[Row[i]][Col[j]], "A");
    fprintf(lrs_ofp, "\n");
    if (i == 0 && Q->nonnegative) /* skip basic rows - don't exist! */
      i = d;
    i++;
    fflush(stdout);
  }
  fflush(stdout);
}

void pimat(lrs_dic *P, long r, long s, lrs_mp Nt, const char *name)
/*print the long precision integer in row r col s of matrix A */
{
  long *B = P->B;
  long *C = P->C;
  if (s == 0)
    fprintf(lrs_ofp, "%s[%ld][%ld]=", name, B[r], C[s]);
  else
    fprintf(lrs_ofp, "[%ld]=", C[s]);
  pmp("", Nt);
}

/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

static void cache_dict(lrs_dic **D_p, lrs_dat *global, long i, long j) {

  if (dict_limit > 1) {
    /* save row, column indicies */
    (*D_p)->i = i;
    (*D_p)->j = j;

    /* Make a new, blank spot at the end of the queue to copy into */

    pushQ(global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A);

    /*2019.6.7 This ought not to happen but it does */

    if (global->Qtail == *D_p)
      return;

    copy_dict(global, global->Qtail, *D_p); /* Copy current dictionary */
  }
  *D_p = global->Qtail;
}

void copy_dict(lrs_dat *global, lrs_dic *dest, lrs_dic *src) {
  long m = src->m;
  long m_A = src->m_A; /* number of rows in A */
  long d = src->d;
  long r, s;

  if (dest == src) {
    if (global->mplrs)
      lrs_post_output("warning", "*copy_dict has dest=src -ignoring copy");
    else
      fprintf(stderr, "*copy_dict has dest=src -ignoring copy");
    return;
  }

  for (r = 0; r <= m_A; r++)
    for (s = 0; s <= d; s++)
      copy(dest->A[r][s], src->A[r][s]);

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->d_orig = src->d_orig;
  dest->m_A = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  copy(dest->det, src->det);
  copy(dest->objnum, src->objnum);
  copy(dest->objden, src->objden);

  memcpy(dest->B, src->B, (m + 1) * sizeof(long));
  memcpy(dest->C, src->C, (d + 1) * sizeof(long));
  memcpy(dest->Row, src->Row, (m + 1) * sizeof(long));
  memcpy(dest->Col, src->Col, (d + 1) * sizeof(long));
}

/*
 * pushQ(lrs_dat *globals,m,d):
 * this routine ensures that Qtail points to a record that
 * may be copied into.
 *
 * It may create a new record, or it may just move the head pointer
 * forward so that know that the old record has been overwritten.
 */
#if 0
#define TRACE(s)                                                               \
  fprintf(stderr, "\n%s %p %p\n", s, global->Qhead, global->Qtail);
#else
#define TRACE(s)
#endif

static void pushQ(lrs_dat *global, long m, long d, long m_A) {

  if ((global->Qtail->next) == global->Qhead) {
    /* the Queue is full */
    if (dict_count < dict_limit) {
      /* but we are allowed to create more */
      lrs_dic *p;

      p = new_lrs_dic(m, d, m_A);

      if (p) {

        /* we successfully created another record */

        p->next = global->Qtail->next;
        (global->Qtail->next)->prev = p;
        (global->Qtail->next) = p;
        p->prev = global->Qtail;

        dict_count++;
        global->Qtail = p;

        TRACE("Added new record to Q");

      } else {
        /* virtual memory exhausted. bummer */
        global->Qhead = global->Qhead->next;
        global->Qtail = global->Qtail->next;

        TRACE("VM exhausted");
      }
    } else {
      /*
       * user defined limit reached. start overwriting the
       * beginning of Q
       */
      global->Qhead = global->Qhead->next;
      global->Qtail = global->Qtail->next;
      TRACE("User  limit");
    }
  }

  else {
    global->Qtail = global->Qtail->next;
    TRACE("Reusing");
  }
}

lrs_dic *lrs_getdic(lrs_dat *Q)
/* create another dictionary for Q without copying any values */
/* derived from lrs_alloc_dic,  used by nash.c                */
{
  lrs_dic *p;

  long m;

  m = Q->m;

  /* nonnegative flag set means that problem is d rows "bigger"     */
  /* since nonnegative constraints are not kept explicitly          */

  if (Q->nonnegative)
    m = m + Q->inputd;

  p = new_lrs_dic(m, Q->inputd, Q->m);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  return p;
}

#define NULLRETURN(e)                                                          \
  if (!(e))                                                                    \
    return NULL;

static lrs_dic *new_lrs_dic(long m, long d, long m_A) {
  lrs_dic *p;

  NULLRETURN(p = (lrs_dic *)malloc(sizeof(lrs_dic)));

  NULLRETURN(p->B = (long int *)calloc((m + 1), sizeof(long)));
  NULLRETURN(p->Row = (long int *)calloc((m + 1), sizeof(long)));

  NULLRETURN(p->C = (long int *)calloc((d + 1), sizeof(long)));
  NULLRETURN(p->Col = (long int *)calloc((d + 1), sizeof(long)));

#if defined(GMP) || defined(FLINT)
  lrs_alloc_mp(p->det);
  lrs_alloc_mp(p->objnum);
  lrs_alloc_mp(p->objden);
#endif

  p->d_orig = d;
  p->A = lrs_alloc_mp_matrix(m_A, d);

  return p;
}

void lrs_free_dic(lrs_dic *P, lrs_dat *Q) {
  /* do the same steps as for allocation, but backwards */
  /* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic *P1;

  if (Q == NULL) {
    if (Q->mplrs)
      lrs_post_output("warning",
                      "*lrs_free_dic trying to free null Q : skipped");
    else
      fprintf(stderr, "*lrs_free_dic trying to free null Q : skipped");
    return;
  }

  if (P == NULL) {
    if (Q->mplrs)
      lrs_post_output("warning",
                      "*lrs_free_dic trying to free null P : skipped");
    else
      fprintf(stderr, "*lrs_free_dic trying to free null P : skipped");
    return;
  }
  /* repeat until cache is empty */

  do {
    /* I moved these here because I'm not certain the cached dictionaries
       need to be the same size. Well, it doesn't cost anything to be safe. db
     */

    long d = P->d_orig;
    long m_A = P->m_A;

    lrs_clear_mp_matrix(P->A, m_A, d);

    /* "it is a ghastly error to free something not assigned my malloc" KR167 */
    /* so don't try: free (P->det);                                           */

    lrs_clear_mp(P->det);
    lrs_clear_mp(P->objnum);
    lrs_clear_mp(P->objden);

    free(P->Row);
    free(P->Col);
    free(P->C);
    free(P->B);

    /* go to next record in cache if any */
    P1 = P->next;
    free(P);
    P = P1;

  } while (Q->Qhead != P);

  Q->Qhead = NULL;
  Q->Qtail = NULL;
}

void lrs_free_dic2(lrs_dic *P, lrs_dat *Q) {
  /* do the same steps as for allocation, but backwards */
  /* same as lrs_free_dic except no cache for P */
  /* I moved these here because I'm not certain the cached dictionaries
     need to be the same size. Well, it doesn't cost anything to be safe. db */

  long d = P->d_orig;
  long m_A = P->m_A;

  lrs_clear_mp_matrix(P->A, m_A, d);

  /* "it is a ghastly error to free something not assigned my malloc" KR167 */
  /* so don't try: free (P->det);                                           */

  lrs_clear_mp(P->det);
  lrs_clear_mp(P->objnum);
  lrs_clear_mp(P->objden);

  free(P->Row);
  free(P->Col);
  free(P->C);
  free(P->B);

  free(P);
}

void lrs_free_dat(lrs_dat *Q) {

  int i = 0;

  if (Q == NULL) {
    if (Q->mplrs)
      lrs_post_output("warning",
                      "*lrs_free_dat trying to free null Q : skipped");
    else
      fprintf(stderr, "*lrs_free_dat trying to free null Q : skipped");
    return;
  }

  /* most of these items were allocated in lrs_alloc_dic */

  lrs_clear_mp_vector(Q->Gcd, Q->m);
  lrs_clear_mp_vector(Q->Lcm, Q->m);
  lrs_clear_mp_vector(Q->output, Q->n);

  lrs_clear_mp(Q->sumdet);
  lrs_clear_mp(Q->Nvolume);
  lrs_clear_mp(Q->Dvolume);
  lrs_clear_mp(Q->saved_det);
  lrs_clear_mp(Q->boundd);
  lrs_clear_mp(Q->boundn);

  free(Q->facet);
  free(Q->redundcol);
  free(Q->inequality);
  free(Q->linearity);
  free(Q->vars);
  free(Q->startcob);
  free(Q->minratio);
  free(Q->redineq);
  free(Q->temparray);

  free(Q->name);
  free(Q->saved_C);

  /*2020.8.1 DA: lrs_global_list is not a stack but a list, so have to delete Q
   */

  while (i < lrs_global_count && lrs_global_list[i] != Q)
    i++;

  if (i == lrs_global_count)
    lrs_warning(Q, "warning", "lrs_free_dat(Q) not in global list - skipped");
  else
    while (i < lrs_global_count) {
      lrs_global_list[i] = lrs_global_list[i + 1];
      i++;
    }

  lrs_global_count--;
  free(Q);
}

static long check_cache(lrs_dic **D_p, lrs_dat *global, long *i_p, long *j_p) {
  /* assign local variables to structures */

  cache_tries++;

  if (global->Qtail == global->Qhead) {
    TRACE("cache miss");
    /* Q has only one element */
    cache_misses++;
    return 0;

  } else {
    global->Qtail = global->Qtail->prev;

    *D_p = global->Qtail;

    *i_p = global->Qtail->i;
    *j_p = global->Qtail->j;

    TRACE("restoring dict");
    return 1;
  }
}

lrs_dic *lrs_alloc_dic(lrs_dat *Q)
/* allocate and initialize lrs_dic */
{

  lrs_dic *p;
  long i, j;
  long m, d, m_A;

  if (Q->hull)        /* d=col dimension of A */
    Q->inputd = Q->n; /* extra column for hull */
  else
    Q->inputd = Q->n - 1;

  m = Q->m;
  d = Q->inputd;
  m_A = m; /* number of rows in A */

  /* nonnegative flag set means that problem is d rows "bigger"     */
  /* since nonnegative constraints are not kept explicitly          */

  if (Q->nonnegative)
    m = m + d;

  p = new_lrs_dic(m, d, m_A);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  dict_count = 1;
  dict_limit = 50;
  cache_tries = 0;
  cache_misses = 0;

  /* Initializations */

  p->d = p->d_orig = d;
  p->m = m;
  p->m_A = m_A;
  p->depth = 0L;
  p->lexflag = TRUE;
  itomp(ONE, p->det);
  itomp(ZERO, p->objnum);
  itomp(ONE, p->objden);

  /*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
  /*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      itomp(ZERO, p->A[i][j]);

  if (Q->nlinearity == ZERO) /* linearity may already be allocated */
    Q->linearity = (long int *)CALLOC((m + d + 1), sizeof(long));

  Q->inequality = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->facet = (long int *)CALLOC((unsigned)m + d + 1, sizeof(long));
  Q->redundcol = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->minratio = (long int *)CALLOC((m + d + 1), sizeof(long));
  /*  2011.7.14  minratio[m]=0 for degen =1 for nondegen pivot*/
  Q->redineq = (long int *)CALLOC((m + d + 1), sizeof(long));
  Q->temparray = (long int *)CALLOC((unsigned)m + d + 1, sizeof(long));

  Q->inequality[0] = 2L;
  Q->Gcd = lrs_alloc_mp_vector(m);
  Q->Lcm = lrs_alloc_mp_vector(m);
  Q->output = lrs_alloc_mp_vector(Q->n);
  Q->saved_C = (long int *)CALLOC(d + 1, sizeof(long));

  Q->lastdv = d; /* last decision variable may be decreased */
                 /* if there are redundant columns          */

  for (i = 0; i < m + d + 1; i++) {
    Q->redineq[i] = 1;
    Q->inequality[i] = 0;
  }

  /*initialize basis and co-basis indices, and row col locations */
  /*if nonnegative, we label differently to avoid initial pivots */
  /* set basic indices and rows */
  if (Q->nonnegative)
    for (i = 0; i <= m; i++) {
      p->B[i] = i;
      if (i <= d)
        p->Row[i] = 0; /* no row for decision variables */
      else
        p->Row[i] = i - d;
    }
  else
    for (i = 0; i <= m; i++) {
      if (i == 0)
        p->B[0] = 0;
      else
        p->B[i] = d + i;
      p->Row[i] = i;
    }

  for (j = 0; j < d; j++) {
    if (Q->nonnegative)
      p->C[j] = m + j + 1;
    else
      p->C[j] = j + 1;
    p->Col[j] = j + 1;
  }
  p->C[d] = m + d + 1;
  p->Col[d] = 0;
  return p;
} /* end of lrs_alloc_dic */

/*
   this routine makes a copy of the information needed to restart,
   so that we can guarantee that if a signal is received, we
   can guarantee that nobody is messing with it.
   This as opposed to adding all kinds of critical regions in
   the main line code.

   It is also used to make sure that in case of overflow, we
   have a valid cobasis to restart from.
 */
static void save_basis(lrs_dic *P, lrs_dat *Q) {
  int i;
  /* assign local variables to structures */
  long *C = P->C;
  long d;

#ifndef SIGNALS
  sigset_t oset, blockset;
  sigemptyset(&blockset);
  sigaddset(&blockset, SIGTERM);
  sigaddset(&blockset, SIGHUP);
  sigaddset(&blockset, SIGUSR1);

  errcheck("sigprocmask", sigprocmask(SIG_BLOCK, &blockset, &oset));
#endif
  d = P->d;

  Q->saved_flag = 1;

  for (i = 0; i < 5; i++)
    Q->saved_count[i] = Q->count[i];

  for (i = 0; i < d + 1; i++)
    Q->saved_C[i] = C[i];

  copy(Q->saved_det, P->det);

  Q->saved_d = P->d;
  Q->saved_depth = P->depth;

#ifndef SIGNALS
  errcheck("sigprocmask", sigprocmask(SIG_SETMASK, &oset, 0));
#endif
}

/* digits overflow is a call from lrs_mp package */

void digits_overflow() {
  fprintf(lrs_ofp, "\noverflow at digits=%ld", DIG2DEC(lrs_digits));
  fprintf(lrs_ofp, "\nrerun with option: digits n, where n > %ld\n",
          DIG2DEC(lrs_digits));
  lrs_dump_state();

  notimpl("");
}

static void lrs_dump_state() {
  long i;

  fprintf(lrs_ofp, "\n\nlrs_lib: checkpointing:\n");

#ifdef MP
  fprintf(stderr, "lrs_lib: Current digits at %ld out of %ld\n",
          DIG2DEC(lrs_record_digits), DIG2DEC(lrs_digits));
#endif

  for (i = 0; i < lrs_global_count; i++) {
    print_basis(lrs_ofp, lrs_global_list[i]);
  }
  fprintf(lrs_ofp, "lrs_lib: checkpoint finished\n");
}

/* print out the saved copy of the basis */
void print_basis(FILE *fp, lrs_dat *global) {
  int i;
  /* assign local variables to structures */
  fprintf(fp, "lrs_lib: State #%ld: (%s)\t", global->id, global->name);

  if (global->saved_flag) {

    /* legacy output which is not actually correct for V-representations as V#
     * is not used */
    /*
          fprintf (fp, "V#%ld R#%ld B#%ld h=%ld facets ",
                   global->saved_count[1],
                   global->saved_count[0],
                   global->saved_count[2],
                   global->saved_depth);
          for (i = 0; i < global->saved_d; i++)
            fprintf (fp, "%ld ",
                     global->inequality[global->saved_C[i] - global->lastdv]);
          pmp (" det=", global->saved_det);
          fprintf (fp, "\n");
    */

    if (global->hull)
      fprintf(fp, "\nrestart %ld %ld %ld ", global->saved_count[0],
              global->saved_count[2], global->saved_depth);
    else
      fprintf(fp, "\nrestart %ld %ld %ld %ld ", global->saved_count[1],
              global->saved_count[0], global->saved_count[2],
              global->saved_depth);

    for (i = 0; i < global->saved_d; i++)
      fprintf(fp, "%ld ",
              global->inequality[global->saved_C[i] - global->lastdv]);
    if (global->saved_count[4] > 0)
      fprintf(fp, "\nintegervertices %ld", global->saved_count[4]);
    fprintf(fp, "\n");

  } else {
    fprintf(fp, "lrs_lib: Computing initial basis\n");
  }

  fflush(fp);
}

#ifndef SIGNALS

/*
   If given a signal
   USR1            print current cobasis and continue
   TERM            print current cobasis and terminate
   INT (ctrl-C) ditto
   HUP                     ditto
 */
static void setup_signals() {
  errcheck("signal", signal(SIGTERM, die_gracefully));
  errcheck("signal", signal(SIGALRM, timecheck));
  errcheck("signal", signal(SIGHUP, die_gracefully));
  errcheck("signal", signal(SIGINT, die_gracefully));
  errcheck("signal", signal(SIGUSR1, checkpoint));
}

static void timecheck(int) {
  lrs_dump_state();
  errcheck("signal", signal(SIGALRM, timecheck));
  alarm(lrs_checkpoint_seconds);
}

static void checkpoint(int) {
  lrs_dump_state();
  errcheck("signal", signal(SIGUSR1, checkpoint));
}

static void die_gracefully(int) {
  lrs_dump_state();

  exit(1);
}

#endif

#ifndef TIMES
/*
 * Not sure about the portability of this yet,
 *              - db
 */
#include <sys/resource.h>
#define double_time(t) ((double)(t.tv_sec) + (double)(t.tv_usec) / 1000000)

static void ptimes() {
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  fprintf(
      lrs_ofp,
      "\n*%0.3fu %0.3fs %ldKb %ld flts %ld swaps %ld blks-in %ld blks-out \n",
      double_time(rusage.ru_utime), double_time(rusage.ru_stime),
      rusage.ru_maxrss, rusage.ru_majflt, rusage.ru_nswap, rusage.ru_inblock,
      rusage.ru_oublock);
  if (lrs_ofp != stdout)
    printf(
        "\n*%0.3fu %0.3fs %ldKb %ld flts %ld swaps %ld blks-in %ld blks-out \n",
        double_time(rusage.ru_utime), double_time(rusage.ru_stime),
        rusage.ru_maxrss, rusage.ru_majflt, rusage.ru_nswap, rusage.ru_inblock,
        rusage.ru_oublock);
}

static double get_time() {
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  return (double_time(rusage.ru_utime));
}

#endif

/* Routines based on lp_solve */

void lrs_set_row(lrs_dic *P, lrs_dat *Q, long row, long num[], long den[],
                 long ineq)
/* convert to lrs_mp then call lrs_set_row */
{
  lrs_mp_vector Num, Den;
  long d;
  long j;

  d = P->d;

  Num = lrs_alloc_mp_vector(d + 1);
  Den = lrs_alloc_mp_vector(d + 1);

  for (j = 0; j <= d; j++) {
    itomp(num[j], Num[j]);
    itomp(den[j], Den[j]);
  }

  lrs_set_row_mp(P, Q, row, Num, Den, ineq);

  lrs_clear_mp_vector(Num, d + 1);
  lrs_clear_mp_vector(Den, d + 1);
}

void lrs_set_row_mp(lrs_dic *P, lrs_dat *Q, long row, lrs_mp_vector num,
                    lrs_mp_vector den, long ineq)
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
    copy(A[i][j], num[j - hull]);
    copy(oD[j], den[j - hull]);
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

void lrs_set_obj(lrs_dic *P, lrs_dat *Q, long num[], long den[], long max) {
  long i;

  if (max == MAXIMIZE)
    Q->maximize = TRUE;
  else {
    Q->minimize = TRUE;
    for (i = 0; i <= P->d; i++)
      num[i] = -num[i];
  }

  lrs_set_row(P, Q, 0L, num, den, GE);
}

void lrs_set_obj_mp(lrs_dic *P, lrs_dat *Q, lrs_mp_vector num,
                    lrs_mp_vector den, long max) {
  long i;

  if (max == MAXIMIZE)
    Q->maximize = TRUE;
  else {
    Q->minimize = TRUE;
    for (i = 0; i <= P->d; i++)
      changesign(num[i]);
  }

  lrs_set_row_mp(P, Q, 0L, num, den, GE);
}

long lrs_solve_lp(lrs_dic *P, lrs_dat *Q)
/* user callable function to solve lp only */
{
  lrs_mp_matrix Lin; /* holds input linearities if any are found             */
  long col;

  Q->lponly = TRUE;

  if (!lrs_getfirstbasis(&P, Q, &Lin, FALSE))
    return FALSE;

  /* There may have been column redundancy                */
  /* If so the linearity space is obtained and redundant  */
  /* columns are removed. User can access linearity space */
  /* from lrs_mp_matrix Lin dimensions nredundcol x d+1   */

  for (col = 0; col < Q->nredundcol; col++) /* print linearity space */
    lrs_printoutput(Q, Lin[col]); /* Array Lin[][] holds the coeffs.     */

  return TRUE;
} /* end of lrs_solve_lp */

long dan_selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns TRUE if pivot found else FALSE                        */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long j, k, col;
  lrs_mp coeff;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  /*  printf("\n*dantzig"); */
  lrs_alloc_mp(coeff);
  *r = 0;
  *s = d;
  j = 0;
  k = 0;

  itomp(0, coeff);
  /*find positive cost coef */
  while (k < d) {
    if (mp_greater(A[0][Col[k]], coeff)) {
      j = k;
      copy(coeff, A[0][Col[j]]);
    }
    k++;
  }

  if (positive(coeff)) /* pivot column found! */
  {
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0) {
      lrs_clear_mp(coeff);
      return (TRUE); /* pivot found */
    }
  }
  lrs_clear_mp(coeff);
  return (FALSE);
} /* end of dan_selectpivot        */

long ran_selectpivot(lrs_dic *P, lrs_dat *Q, long *r, long *s)
/* select pivot indices using random edge rule                   */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long i, j, k, col, t;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;
  long *perm;

  perm = (long *)calloc((d + 1), sizeof(long));
  *r = 0;
  *s = d;
  k = 0;
  /*  printf("\n*random edge"); */

  /* generate random permutation of 0..d-1 */
  for (i = 0; i < d; i++)
    perm[i] = i;

  for (i = 0; i < d; i++) {
    j = rand() % (d - i) + i;
    t = perm[j];
    perm[j] = perm[i];
    perm[i] = t; // Swap i and j
  }

  /*find first positive cost coef according to perm */
  while (k < d && !positive(A[0][Col[perm[k]]]))
    k++;

  if (k < d) /* pivot column found! */
  {
    j = perm[k];
    *s = j;
    col = Col[j];

    /*find min index ratio */
    *r = lrs_ratio(P, Q, col);
    if (*r != 0) {
      free(perm);
      return (TRUE); /* pivot found */
    }
  }
  free(perm);
  return (FALSE);
} /* end of ran_selectpivot        */

long phaseone(lrs_dic *P, lrs_dat *Q)
/* Do a dual pivot to get primal feasibility (pivot in X_0)*/
/* Bohdan Kaluzny's handiwork                                    */
{
  long i, j, k;
  /* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d;
  lrs_mp b_vector;
  lrs_alloc_mp(b_vector);
  m = P->m;
  d = P->d;
  i = 0;
  k = d + 1;

  itomp(0, b_vector);

  fprintf(lrs_ofp, "\nLP: Phase One: Dual pivot on artificial variable");

  /*find most negative b vector */
  while (k <= m) {
    if (mp_greater(b_vector, A[Row[k]][0])) {
      i = k;
      copy(b_vector, A[Row[i]][0]);
    }
    k++;
  }

  if (negative(b_vector)) /* pivot row found! */
  {
    j = 0; /*find a positive entry for in row */
    while (j < d && !positive(A[Row[i]][Col[j]]))
      j++;
    if (j >= d) {
      lrs_clear_mp(b_vector);
      return (FALSE); /* no positive entry */
    }
    pivot(P, Q, i, j);
    update(P, Q, &i, &j);
  }
  lrs_clear_mp(b_vector);
  return (TRUE);
}

long lrs_set_digits(long dec_digits) {
  /* convert user specified decimal digits to mp digits */

  if (dec_digits > 0)
    lrs_digits = DEC2DIG(dec_digits);
  if (lrs_digits > MAX_DIGITS) {
    fprintf(lrs_ofp,
            "\nDigits must be at most %ld\nChange MAX_DIGITS and recompile",
            DIG2DEC(MAX_DIGITS));
    fflush(stdout);
    return (FALSE);
  }
  return (TRUE);
}

long lrs_checkbound(lrs_dic *P, lrs_dat *Q) {
  /* check bound on objective and return TRUE if exceeded */

  if (!Q->bound)
    return FALSE;

  if (Q->maximize && comprod(Q->boundn, P->objden, P->objnum, Q->boundd) == 1) {
    return TRUE;
  }
  if (Q->minimize &&
      comprod(Q->boundn, P->objden, P->objnum, Q->boundd) == -1) {
    return TRUE;
  }
  return FALSE;
}

long lrs_leaf(lrs_dic *P, lrs_dat *Q) {
  /* check if current dictionary is a leaf of reverse search tree */
  long col = 0;
  long tmp = 0;

  while (col < P->d && !reverse(P, Q, &tmp, col))
    col++;
  if (col < P->d)
    return 0; /* dictionary is not a leaf */
  else
    return 1;
}

/* prevent output flushes in mplrs */
void lrs_open_outputblock(void) {}

/* re-enable output flushes in mplrs */
void lrs_close_outputblock(void) {}

void lrs_post_output(const char *type, const char *data) {}

void lrs_return_unexplored(
    lrs_dic *P, lrs_dat *Q) /* send cobasis data for unexplored nodes */

{}

#ifdef MP
void lrs_overflow(int parm) { lrs_exit(parm); }
#endif

#ifdef LRSLONG

/* replace by user overflow routine if not using lrsv2_main() */
void lrs_overflow(int parm) { lrsv2_overflow(parm); }

void lrsv2_overflow(int parm) {
  lrs_dat *Q;
  lrs_dic *P;
  char *restart;
  char *part;

  int i;
  int try_restart = FALSE;

  if (lrs_global_list[0] == NULL) {
    fprintf(stderr, "*lrs_overflow has null Q ");
    lrs_exit(parm);
  }

  /* db's cunningly hidden locations */
  Q = lrs_global_list[lrs_global_count - 1];
  P = Q->Qhead;

  /* non mplrs overflow handling             */
  /* lrs, redund,fel restarted at the moment */

#ifdef MA
  if (strcmp(Q->fname, "lrs") == 0 || strcmp(Q->fname, "lrsmp") == 0 ||
      Q->redund || Q->fel)
    try_restart = TRUE;
#endif

  if (lrs_ifp != NULL)
    fclose(lrs_ifp);

  if (!try_restart) /* hard exit */
  {
    if (strcmp(BIT, "64bit") == 0) {
      fprintf(
          stderr,
          "\n*64bit integer overflow: try running 128bit or gmp versions\n");
      if (lrs_ofp != stdout)
        fprintf(
            lrs_ofp,
            "\n*64bit integer overflow: try running 128bit or gmp versions\n");
    } else {
      fprintf(stderr, "\n*128bit integer overflow: try running gmp version\n");
      if (lrs_ofp != stdout)
        fprintf(lrs_ofp,
                "\n*128bit integer overflow: try running gmp version\n");
    }
    lrs_exit(parm);
  }

  /* try to restart */
  if (overflow == 0) /*  first overflow */
  {
    if (*tmpfilename != '\0') /* we made a temporary file for stdin  */
      if (remove(tmpfilename) != 0)
        fprintf(lrs_ofp, "\nCould not delete temporary file");
    strncpy(tmpfilename, "/tmp/lrs_restartXXXXXX", PATH_MAX);
    /* XXX in principle this file descriptor should be used instead of the name
     */
    tmpfd = mkstemp(tmpfilename);
  } else
    strcpy(tmpfilename, infilename);

  if (!pivoting || Q->redund || Q->getvolume || Q->fel ||
      Q->extract) /* we make restart from original input   */
  {
    overflow = 1L;
    lrs_cache_to_file(tmpfilename, " ");
  } else {
    restart = (char *)malloc(Q->saved_d * 20 + 100);
    part = (char *)malloc(Q->saved_d * 20 + 100);
    overflow = 2L;
    if (Q->hull)
      sprintf(restart, " %ld %ld %ld ", Q->saved_count[2], Q->saved_count[0],
              Q->saved_depth);
    else
      sprintf(restart, " %ld %ld %ld %ld ", Q->saved_count[1],
              Q->saved_count[0], Q->saved_count[2], Q->saved_depth);

    for (i = 0; i < Q->saved_d; i++) {
      sprintf(part, "%ld ", Q->inequality[Q->saved_C[i] - Q->lastdv]);
      strcat(restart, part);
    }
    sprintf(part, "\nintegervertices %ld", Q->saved_count[4]);
    strcat(restart, part);

    lrs_cache_to_file(tmpfilename, restart);
    free(restart);
    free(part);
  }

  if (Q->fel || Q->redund)
    if (Q->Ain != NULL)
      lrs_clear_mp_matrix(Q->Ain, Q->m, Q->n);

  Q->m = P->m;

  lrs_free_dic(P, Q); /* note Q is not freed here and is needed again  */

  if (Q->fel && !Q->hull)
    lrs_free_dat(Q); /* in this case we free Q as it was alloc'ed in fel_run */

  if (lrs_ofp != NULL && lrs_ofp != stdout) {
    fclose(lrs_ofp);
    lrs_ofp = NULL;
  }
  close(tmpfd);

  longjmp(buf1, 1); /* return to lrsv2_main */

  lrs_exit(parm); /* should not happen */
}
#endif

void lrs_exit(int i) {
  fflush(stdout);
  exit(i);
}

void lrs_free_all_memory(lrs_dic *P, lrs_dat *Q) {

  if (Q->runs > 0) {
    free(Q->isave);
    free(Q->jsave);
  }
  if (P != NULL) /* may not have allocated P yet */
  {
    long savem = P->m;  /* need this to clear Q*/
    lrs_free_dic(P, Q); /* deallocate lrs_dic */
    Q->m = savem;
  }
  lrs_free_dat(Q); /* deallocate lrs_dat */
#ifdef LRSLONG
  free(infile); /* we cached input file for possible restart */
#endif
  return;
}

long lrs_stdin_to_file(char *filename) {
  FILE *fptr1, *fptr2;
  char c;

  fptr1 = stdin;
  fptr2 = fopen(filename, "w");
  if (fptr2 == NULL) {
    printf("Cannot open file %s \n", filename);
    exit(0);
  }

  c = fgetc(fptr1);
  while (c != EOF) {
    fputc(c, fptr2);
    c = fgetc(fptr1);
  }

  fclose(fptr2);
  fptr2 = NULL;

  return 0;
}

long lrs_file_to_cache(FILE *ifp) {
  long ret;

  if (ifp != NULL)
    if (fseek(ifp, 0L, SEEK_END) == 0) {
      ret = ftell(ifp);
      if (ret == -1) {
        fputs("*Error reading file", stderr);
        return 1;
      }
      infileLen = ret;
      infile = (char *)malloc(sizeof(char) * (infileLen + 1));

      if (fseek(ifp, 0L, SEEK_SET) != 0) {
        fputs("*Error resetting input file", stderr);
        return 1;
      }
      infileLen = fread(infile, sizeof(char), infileLen, ifp);
      if (ferror(ifp) != 0) {
        fputs("*Error reading input file", stderr);
        return 1;
      } else
        infile[infileLen++] = '\0'; /* Just to be safe. */
    }
  rewind(ifp);
  return 0;
}

long lrs_cache_to_file(char *name, const char *restart) {
  FILE *ofp = fopen(name, "wb");

  if (ofp == NULL) {
    printf("*Error opening output file %s", name);
    return 1;
  }
  fwrite(infile, sizeof(char), infileLen, ofp);

  if (lrs_global_list[0]->count[2] > 1L && overflow == 2)
    fprintf(ofp, "\nrestart %s", restart);

  fclose(ofp);
  return 0;
}

void lrs_warning(lrs_dat *Q, char *type, char *ss) {
  if (Q->messages) {
    if (Q->mplrs)
      lrs_post_output(type, ss);
    else {
      fprintf(lrs_ofp, "\n%s", ss);
      if (lrs_ofp != stdout)
        fprintf(stderr, "\n%s", ss);
    }
  }
}

/*********************************************************************************/
/*2022.1.18  For V-rep with lponly just find optimizing input rows */
/*********************************************************************************/
long lrs_check_inequality(lrs_dic *P, lrs_dat *Q) {
  lrs_mp_matrix A = P->A;
  lrs_mp tmp, opt, total;
  long m, d, i, j, count;

  lrs_alloc_mp(tmp);
  lrs_alloc_mp(total);
  lrs_alloc_mp(opt);

  fprintf(lrs_ofp, "\n");

  m = P->m;
  d = P->d;
  itomp(0, opt);

  if (Q->nonnegative) /* skip basic rows - don't exist! */
    m = m - d;

  for (i = 1; i <= m; i++) {
    itomp(0, total);
    for (j = 1; j <= d; j++) {
      mulint(A[0][j], A[i][j], tmp);
      linint(total, 1, tmp, 1);
    }
    if (i == 1)
      copy(opt, total);
    else if (mp_greater(total, opt))
      copy(opt, total);
  }
  fprintf(lrs_ofp, "\n*optimum rows:");
  count = 0;
  for (i = 1; i <= m; i++) /* once more to print optima */
  {
    itomp(0, total);
    for (j = 1; j <= d; j++) {
      mulint(A[0][j], A[i][j], tmp);
      linint(total, 1, tmp, 1);
    }
    if (!mp_greater(opt, total)) {
      count++;
      fprintf(lrs_ofp, " %ld", i);
    }
  }

  if (Q->minimize) /* lrs only maximizes */
  {
    changesign(opt);
    prat("\n*min value:", opt, P->det);
  } else
    pmp("\n*max value:", opt);
  fprintf(lrs_ofp, " obtained by %ld rows", count);
  fprintf(lrs_ofp, "\n");
  lrs_clear_mp(tmp);
  lrs_clear_mp(opt);

  return 0;
}
