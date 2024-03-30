#include <gmp.h>

#ifndef LIBLRSNASHGMP
#define LIBLRSNASHGMP

#ifdef __cplusplus
extern "C"
{
#endif

    extern void solve_gmp(
        int rows, int cols,
        const mpq_t *row_payoff_data, const mpq_t *col_payoff_data,
        mpz_t *row_solution_data, mpz_t *col_solution_data);

    extern void solve_gmp_pointer(
        int rows, int cols,
        const mpq_t **row_payoff_data, const mpq_t **col_payoff_data,
        mpz_t *row_solution_data, mpz_t *col_solution_data);

    extern void solve_gmp_pointer_constant_sum(
        int rows, int cols,
        const mpq_t **row_payoff_data,
        mpz_t *row_solution_data, mpz_t *col_solution_data,
        int payoff_sum_num, int payoff_sum_den);

    extern void solve_gmp_float(
        int rows, int cols,
        const long *payoff_data, long den,
        mpz_t *row_solution_data, mpz_t *col_solution_data);

    extern mpz_t *alloc(size_t n);
    extern void dealloc(mpz_t *x, size_t n);

#ifdef __cplusplus
}
#endif

#endif // LIBLRSNASHGMP