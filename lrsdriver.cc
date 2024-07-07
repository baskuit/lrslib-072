/* This file contains functions and variables that should not be duplicated per arithmetic */

#include <stdio.h>
#include <string.h>
#include <setjmp.h>
#include <stdlib.h>
#include <limits.h>
#include "lrsdriver.h"

/* Globals; these need to be here, rather than lrsdriver.h, so they are
   not multiply defined. */

FILE *lrs_cfp;			/* output file for checkpoint information       */
FILE *lrs_ifp;			/* input file pointer       */
FILE *lrs_ofp;			/* output file pointer      */