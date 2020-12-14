#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* include inline debug statements? */
#define DEBUG 0

/* at what length does quicksort switch to insertion sort? */
#define QUICKSORT_INSERTION_CUTOFF       20
#define QUICKSORT3WAY_INSERTION_CUTOFF   20


#ifdef LONG_VECTOR_SUPPORT
#  define MY_SIZE_T R_xlen_t
#  define MY_LENGTH xlength
#else
#  define MY_SIZE_T int
#  define MY_LENGTH length
#endif


/* FUNCTION PROTOTYPE DECLARATION *********************************/

static void
fr_quicksort_double_i_(const double * a, MY_SIZE_T indx[], const MY_SIZE_T n);

SEXP fastrank_num_avg_(SEXP s_x);




/* FUNCTION REGISTRATION *********************************/

static R_CallMethodDef callMethods[] = {
    {"fastrank_num_avg_", (DL_FUNC) &fastrank_num_avg_, 1},
    {NULL,                NULL,                         0}
};

void R_init_fastrank(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}



/* SORTING *********************************
 *
 * See also BENCHMARKING.md, which I'll create after I've completed all the
 * benchmarking I intend to do.
 *
 * Quicksort from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort
 * and especially Wikipedia modified two ways:
 *
 * 1. Return a vector of indices and not modify the array of values a[],
 * which requires that the indx[] vector is pre-allocated and filled 0..n-1
 *
 * 2. Use insertion sort for vectors length <= QUICKSORT_INSERTION_CUTOFF
 *
 * Note the papers and resources (esp Sedgwick) I have collected, the LESSER
 * loops should probably be changed to <= to avoid pathological behaviour.
 * Also if a[indx[i]] == a[indx[j]] then they should only be swapped if
 * j < i.  That might be enough to make it stable...
 *
 * Also note the Dual Pivot Quicksort papers, that might be the truly
 * best way to go.
 */


#define SWAP(__T, __A, __B) { __T t = __A; __A = __B; __B = t; }

#undef __TYPE
#undef __LESSER
#define FR_quicksort_body(__TYPE, __LESSER) \
    MY_SIZE_T i; /* used as param outside of body */ \
    { \
    __TYPE pvt; \
    MY_SIZE_T j, it; \
    if (n <= QUICKSORT_INSERTION_CUTOFF) { \
        for (i = 1; i < n; ++i) { \
            it = indx[i]; \
            for (j = i; j > 0 && __LESSER(a[it], a[indx[j - 1]]); --j) { \
                indx[j] = indx[j - 1]; \
            } \
            indx[j] = it; \
        } \
        return; \
    } \
    pvt = a[indx[n / 2]]; \
    for (i = 0, j = n - 1; ; i++, j--) { \
        while (__LESSER(a[indx[i]], pvt)) i++; \
        while (__LESSER(pvt, a[indx[j]])) j--; \
        if (i >= j) break; \
        SWAP(MY_SIZE_T, indx[i], indx[j]); \
    } \
    }


#undef LESSER
#undef EQUAL
#define LESSER(__A, __B) (__A < __B)
#define EQUAL(__A, __B) (__A == __B)
static void
fr_quicksort_double_i_ (const double * a,
                        MY_SIZE_T indx[],
                        const MY_SIZE_T n) {

    FR_quicksort_body(double, LESSER)

    fr_quicksort_double_i_(a, indx,     i    );
    fr_quicksort_double_i_(a, indx + i, n - i);
}





/* Rank a numeric vector giving ties their average rank */
SEXP fastrank_num_avg_(SEXP s_x) {
    MY_SIZE_T n = MY_LENGTH(s_x);
    double *x = REAL(s_x);
    if (DEBUG) {
        Rprintf("    x:  ");
        for (int i = 0; i < n; ++i) Rprintf("%.3f ", x[i]);
        Rprintf("\n");
    }

    /* double because "average" */
    SEXP s_ranks = PROTECT(allocVector(REALSXP, n));
    double *ranks = REAL(s_ranks);
    MY_SIZE_T *indx = (MY_SIZE_T *) R_alloc(n, sizeof(MY_SIZE_T));
    /* pre-fill indx with index from 0..n-1 */
    for (MY_SIZE_T i = 0; i < n; ++i) indx[i] = i;
    fr_quicksort_double_i_(x, indx, n);
    if (DEBUG) {
        Rprintf(" indx:   ");
        for (int i = 0; i < n; ++i) Rprintf("%d    ", indx[i]);
        Rprintf("\n");
    }

    MY_SIZE_T ib = 0;
    double b = x[indx[0]];
    MY_SIZE_T i;
    for (i = 1; i < n; ++i) {
        if (x[indx[i]] != b) { /* consecutive numbers differ */
            if (ib < i - 1) {  /* average of sum of ranks */
                double rnk = (i - 1 + ib + 2) / 2.0;
                for (MY_SIZE_T j = ib; j <= i - 1; ++j)
                    ranks[indx[j]] = rnk;
            } else {
                ranks[indx[ib]] = (double)(ib + 1);
            }
            b = x[indx[i]];
            ib = i;
        }
    }
    /* now check leftovers */
    if (ib == i - 1)  /* last two were unique */
        ranks[indx[ib]] = (double)i;
    else {  /* ended with ties */
        double rnk = (i - 1 + ib + 2) / 2.0;
        for (MY_SIZE_T j = ib; j <= i - 1; ++j)
            ranks[indx[j]] = rnk;
    }
    if (DEBUG) {
        Rprintf("ranks:  ");
        for (int i = 0; i < n; ++i) Rprintf("%.1f   ", ranks[i]);
        Rprintf("\n");
    }

    UNPROTECT(1);
    return s_ranks;
}
