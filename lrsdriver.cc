/* This file contains functions and variables that should not be duplicated per
 * arithmetic */

#include "lrsdriver.h"
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Globals; these need to be here, rather than lrsdriver.h, so they are
   not multiply defined. */

FILE *lrs_cfp; /* output file for checkpoint information       */
FILE *lrs_ifp; /* input file pointer       */
FILE *lrs_ofp; /* output file pointer      */