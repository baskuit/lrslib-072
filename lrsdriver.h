/* This file contains functions and variables that should not be duplicated per arithmetic */

#ifndef LRS_DRIVER_H
#define LRS_DRIVER_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lrsrestart.h"

struct lrs_dic_struct;
typedef struct lrs_dic_struct lrs_dic;

struct lrs_dat;
typedef struct lrs_dat lrs_dat;

extern FILE *lrs_cfp;			/* output file for checkpoint information       */
extern FILE *lrs_ifp;			/* input file pointer       */
extern FILE *lrs_ofp;			/* output file pointer      */

#endif
