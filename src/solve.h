/*******************************************************/
/* lrsnashlib is a library of routines for computing   */
/* computing all nash equilibria for two person games  */
/* given by mxn payoff matrices A,B                    */
/*                                                     */
/*                                                     */
/* Main user callable function is                      */
/*         lrs_solve_nash(game *g)                     */
/*                                                     */
/* Sample driver: lrsnash.c                            */
/* Derived from nash.c in lrslib-060                   */
/* by Terje Lensberg, October 26, 2015:                */
/*******************************************************/

/* minor mod 2018.2.12 to set Q->fname="nash"          */

/*************/
/* Games     */
/*************/

#define MAXSTRAT 200
#define ROW 0
#define COL 1

typedef struct
{
	long num;
	long den;
} ratnum;

typedef struct
{
	long nstrats[2];
#ifdef SURSKIT
#ifdef GMP
	mpq_t *row_payoff; // must be dynamically allocated
	mpq_t *col_payoff;
#else
	ratnum *row_payoff;
	ratnum *col_payoff;
#endif
#else
	ratnum payoff[MAXSTRAT][MAXSTRAT][2];
#endif
	// For auxiliary information
	void *aux;
} game;

typedef struct
{
	char name[100];
	int fwidth[MAXSTRAT][2]; // Column field widths (for output)
} gInfo;

int lrs_solve_nash(game *g);
int solve(game *g, lrs_mp_vector row_data, lrs_mp_vector col_data);
int solve_pointer(game *g, const mpq_t** row_payoff_data, const mpq_t** col_payoff_data, lrs_mp_vector row_data, lrs_mp_vector col_data);
int solve_pointer_constant_sum(game *g, const mpq_t** row_payoff_data, lrs_mp_vector row_data, lrs_mp_vector col_data, int payoff_sum_num, int payoff_sum_den);
int solve_float (game *g, const long* payoff_data, long den, lrs_mp_vector row_data, lrs_mp_vector col_data);

long nash2_main(lrs_dic *P1, lrs_dat *Q1, lrs_dic *P2orig,
				lrs_dat *Q2, long *numequilib, lrs_mp_vector output, long linindex[]);
long nash2_main_(lrs_dic *P1, lrs_dat *Q1, lrs_dic *P2orig,
				lrs_dat *Q2, long *numequilib, lrs_mp_vector output, long linindex[], lrs_mp_vector data);
/* lrs driver, argv[2]= 2nd input file for nash equilibria */

long lrs_getfirstbasis2(lrs_dic **D_p, lrs_dat *Q, lrs_dic *P2orig, lrs_mp_matrix *Lin, long no_output,
						long linindex[]);

long getabasis2(lrs_dic *P, lrs_dat *Q, lrs_dic *P2orig, long order[], long linindex[]);

long lrs_nashoutput(lrs_dat *Q, lrs_mp_vector output, long player);
long lrs_nashoutput_(lrs_dat *Q, lrs_mp_vector output, long player, lrs_mp_vector data);
/* returns TRUE and prints output if not the origin */

int lrs_solve_nash_legacy(int argc, char *argv[]);

void BuildRep(lrs_dic *P, lrs_dat *Q, const game *g, int p1, int p2);
void BuildRepPointer(lrs_dic *P, lrs_dat *Q, const game *g, const mpq_t **row_payoff_data, const mpq_t **col_payoff_data, int p1, int p2);
void BuildRepPointerConstantSum(lrs_dic *P, lrs_dat *Q, const game *g, const mpq_t **row_payoff_data, int p1, int p2, int payoff_sum_num, int payoff_sum_den);
void BuildRepFloat(lrs_dic *P, lrs_dat *Q, const game *g, const long *payoff_data, long den, int p1, int p2);
void FillFirstRow(lrs_dic *P, lrs_dat *Q, int n);
void FillLinearityRow(lrs_dic *P, lrs_dat *Q, int m, int n);
void FillConstraintRows(lrs_dic *P, lrs_dat *Q, const game *g, int p1, int p2, int firstRow);
void FillNonnegativityRows(lrs_dic *P, lrs_dat *Q, int firstRow, int lastRow, int n);
void printGame(game *g);
void setFwidth(game *g, int len);
void resetNashSolver(); /* Call this function for every new game to be solved */
void initFwidth(game *g);
void updateFwidth(game *g, int col, int pos, char *str);

extern __thread long Debug_flag;
extern __thread long Verbose_flag;
