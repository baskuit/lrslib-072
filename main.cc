/*********************************************************/
/* lrsnash is driver for computing all nash equilibria   */
/* for two person games given by mxn payoff matrices A,B */
/* Usage: lrsnash [options] game [output file]           */
/*                                                       */
/* game is a file containing:                            */
/* m n                                                   */
/* A          (row by row)                               */
/* B          (row by row)                               */
/*                                                       */
/* Derived from nash.c in lrslib-060 by                  */
/* Terje Lensberg, October 26, 2015:                     */
/* Simplified API via lrs_solve_nash(game *g)            */
/*********************************************************/
/*
Compile:
gcc -O3 -o lrsnash lrsnash.c lrsnashlib.c lrslib.c lrsgmp.c -lgmp -DGMP
*/

char Usage[] =
    "usage (standard): %s [options...] <input_file...>\n"
    "usage (legacy):   %s <input_file1> <input_file2> [<output_file>]\n"
    "type %s -h for more information\n\n";

char Helptext[] =
    "\nusage (standard): %s [options...] <input_file...>\n"
    "  Input file structure: Input files to setupnash\n"
    "  Input files can be specified separately, or by using wildcards, as in "
    "'game*'\n"
    "  Options:\n"
    "    -p, --printgame       Prints the payoff matrix for the game\n"
    "    -s, --standard        Promise that input files have standard "
    "structure\n"
    "    -o, --outfile <file>  Send output to <file>\n"
    "    -h, --help            Prints this text\n"
    "     Short options can be grouped, as in '-ps' and '-do out.txt'\n"
    "usage (legacy): %s <input_file1> <input_file2> [<output_file>]\n"
    "  Input file structure: Output files from setupnash\n"
    "  Passing options with legacy input files produces an error\n"
    "  (options must be specified in the input files)\n\n";

char LegacyMsg[] =
    "\nProcessing legacy input files. Alternatively, you may skip\n"
    "setupnash and pass its input file to this program.\n";

#include "lrsdriver.h"
#include "lrslib.h"
#include "lib.h"
#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//========================================================================
// Games
//========================================================================

//----------------------------------------------------------------------------------------//
// Reading games
//----------------------------------------------------------------------------------------//

char *Outfile = NULL;

//----------------------------------------------------------------------------------------//
int openIO(void) {
  if (!lrs_init("*lrsnash:"))
    return FALSE;
  fprintf(stderr, "\n");
  if (Outfile != NULL) {
    if ((lrs_ofp = fopen(Outfile, "w")) == NULL) {
      fprintf(stderr, "\nBad output file name\n");
      return FALSE;
    }
  }
  return TRUE;
}

void closeIO(void) {
  if (lrs_ofp != stdout)
    fprintf(stdout, "\n");
  lrs_close("lrsnash:");
}

#define RATWARN(name, filename)                                                \
  fprintf(stderr,                                                              \
          "\nWarning: String '%s' is not a rational number in file %s.\n",     \
          name, filename);
#define RECWARN(filename)                                                      \
  fprintf(stderr, "\nWarning: Excess data in file %s.\n", filename);
#define ERREXIT                                                                \
  if (lrs_ofp != NULL)                                                         \
    closeIO();                                                                 \
  exit(1);
#define FILEERROR(name)                                                        \
  {                                                                            \
    fprintf(stderr, "\nError: Cannot find input file '%s'. \
  Execution halted\n",                                                         \
            name);                                                             \
    ERREXIT                                                                    \
  }
#define READERROR(name)                                                        \
  {                                                                            \
    fprintf(stderr, "\nError: Premature end of input file '%s'. \
  Execution halted\n",                                                         \
            name);                                                             \
    ERREXIT                                                                    \
  }
#define SIZEERROR(name)                                                        \
  {                                                                            \
    fprintf(                                                                   \
        stderr,                                                                \
        "\nError: Number of strategies exceeds maximum (%d) in input file '%s'. \
  Execution halted\n",                                                         \
        MAXSTRAT, name);                                                       \
    ERREXIT                                                                    \
  }

//----------------------------------------------------------------------------------------//
// Simple function to convert string to (num, den)
int tl_readrat(long *num, long *den, char *str) {
  char *div = strchr(str, '/');
  if (div == NULL) {
    *num = atol(str);
    *den = 1;
  } else if (div == str || *(div + 1) == 0) { //  str = '/x' or str = 'x/'
    return FALSE;
  } else {
    *div = 0; // Note: 'str' is modified here
    *num = atol(str);
    *den = atol(div + 1);
  }
  return TRUE;
}

//----------------------------------------------------------------------------------------//
int readGame(game *g, const char *filename) {
  FILE *IN;
  long pos, s, t, nr, nc;
  char in[MAXINPUT];
  strcpy(((gInfo *)g->aux)->name, filename);
  if ((IN = fopen(filename, "r")) == NULL)
    FILEERROR(filename);
  if (fscanf(IN, "%ld %ld", &nr, &nc) < 2)
    READERROR(filename);
  if (nr > MAXSTRAT || nc > MAXSTRAT)
    SIZEERROR(filename);
  g->nstrats[ROW] = nr;
  g->nstrats[COL] = nc;
  initFwidth(g);
  // Read payoffs
  for (pos = 0; pos < 2; pos++) {
    for (s = 0; s < nr; s++) {
      for (t = 0; t < nc; t++) {
        if (fscanf(IN, "%s", in) < 1)
          READERROR(filename);
        updateFwidth(g, t, pos, in);
        if (!tl_readrat(&g->payoff[s][t][pos].num, &g->payoff[s][t][pos].den,
                        in))
          RATWARN(in, filename);
      }
    }
  }
  if (fscanf(IN, "%s", in) == 1) // Too many payoff entries
    RECWARN(filename);
  fclose(IN);
  return TRUE;
}

//----------------------------------------------------------------------------------------//

//========================================================================
// Command line processing
//========================================================================

// Flags to be set from command line options
static long Print_game_flag;
static long Standard_input_flag;

void printUsage(const char *progname) {
  fprintf(stderr, Usage, progname, progname, progname);
}

void printInfo(const char *progname) {
  fprintf(stderr, Helptext, progname, progname);
}

//----------------------------------------------------------------------------------------//
// Collects flags and reads list of games
int getArgs(int argc, char **argv) {
  int c, error = FALSE;
  const char shortOptions[] = ":vdpsho:";

  if (argc < 2) {
    printUsage(argv[0]);
    return FALSE;
  }

  while (1) {
    static struct option long_options[] = {
        {"printgame", no_argument, 0, 'p'},
        {"standard", no_argument, 0, 's'},
        {"outfile", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'}
        //      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, shortOptions, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case '?':
      fprintf(stderr, "\nError: Unknown option '-%c'.\n", optopt);
      error = TRUE;
      break;

    case ':':
      fprintf(stderr, "\nError: Missing argument to option '-%c'.\n", optopt);
      error = TRUE;
      break;

    case 'd':
      Debug_flag = TRUE;
      break;

    case 'p':
      Print_game_flag = TRUE;
      break;

    case 's':
      Standard_input_flag = TRUE;
      break;

    case 'h':
      printInfo(argv[0]);
      return FALSE;
      break;

    case 'o':
      Outfile = optarg;
      break;

    default:
      abort();
    }
  }
  if (error) {
    fprintf(stderr, "Execution halted\n");
    return FALSE;
  }
  return TRUE;
}

//========================================================================
// Main()
//========================================================================
int main(int argc, char **argv) {
  game Game; // Storage for one game
  game *g = &Game;
  gInfo GI; // Storage for auxiliary information about the game
  g->aux = &GI;

  if (!getArgs(argc,
               argv)) // Read options and input file names. When we get here:
    return 1;         // optind is a global integer supplied by getopt, and
                      // argv[optind] is the first non-option argument in argv

  if (!openIO())
    return 1;
  while (optind < argc) { // Handle standard input file[s]
    if (readGame(g, argv[optind++])) {
      if (Print_game_flag)
        printGame(g);
      for (int i = 0; i < (1 << 0); ++i) {
        SolveInput input;
        input.rows = g->nstrats[0];
        input.cols = g->nstrats[1];
        input.den = 64;
        input.data = (int *)calloc(input.rows * input.cols, sizeof(int));
        for (int ii = 0; ii < input.rows; ++ii) {
          for (int jj = 0; jj < input.cols; ++jj) {
            input.data[ii * input.cols + jj] = g->payoff[ii][jj][0].num;
          }
        }
        lrs_solve_nash(&input);
        free(input.data);
      }
    }
  }
  closeIO();

  return 0;
}
