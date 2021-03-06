/*
 * lbfgsb_wrapper.c -
 *
 *	C wrapper for L-BFGS-B FORTRAN routines.
 *
 */


#include <string.h>
/*#include <stdio.h>*/

/* The following FORTRAN-C data types equivalences are assumed:
 *	  DOUBLE PRECISION = double
 *	           INTEGER = integer_t (define macro below)
 *	           LOGICAL = logical_t (define macro below)
 *
 * Define the following macros to match data types from your FORTRAN
 * compiler (note: Yorick assumes that integer_t = long):
 */

#define integer_t long
#define logical_t int
#define ftnlen_t  integer_t /* length of strings */
#define real_t    double    /* FORTRAN "DOUBLE PRECISION" type */

#define LBFGSB_TASK_LENGTH 60
#define TRUE  1
#define FALSE 0

extern integer_t lbfgsb_wrapper(integer_t n, integer_t m, real_t x[],
				real_t *f, real_t g[], const real_t l[],
				const real_t u[], const integer_t bnd[],
				real_t factr, real_t pgtol,
				char task[], char csave[], integer_t isave[],
				real_t dsave[], integer_t iprint);
/*
 * lbfgs_wrapper -
 *
 *   This subroutine partitions the working arrays wa and iwa, and
 *      then uses the limited memory BFGS method to solve the bound
 *      constrained optimization problem by calling mainlb.
 *      (The direct method will be used in the subspace minimization.)
 *   
 *   N is an integer variable.
 *      On entry N is the dimension of the problem.
 *      On exit N is unchanged.
 *   
 *   M is an integer variable.
 *      On entry M is the maximum number of variable metric corrections
 *        used to define the limited memory matrix.
 *      On exit M is unchanged.
 *   
 *   X is a double precision array of dimension N.
 *      On entry X is an approximation to the solution.
 *      On exit X is the current approximation.
 *   
 *   F is a double precision variable.
 *      on first entry F is unspecified.
 *      On final exit F is the value of the function at X.
 *   
 *   G is a double precision array of dimension N.
 *      On first entry G is unspecified.
 *      On final exit G is the value of the gradient at X.
 *   
 *   L is a double precision array of dimension N.
 *      On entry L is the lower bound on X.
 *      On exit L is unchanged.
 *   
 *   U is a double precision array of dimension N.
 *      On entry U is the upper bound on X.
 *      On exit U is unchanged.
 *   
 *   BND is an integer array of dimension N.
 *      On entry BND represents the type of bounds imposed on the
 *        variables, and must be specified as follows:
 *        BND(i)=0 if X(i) is unbounded,
 *               1 if X(i) has only a lower bound,
 *               2 if X(i) has both lower and upper bounds, and
 *               3 if X(i) has only an upper bound.
 *      On exit BND is unchanged.
 *
 *   FACTR is a double precision variable.
 *      On entry FACTR >= 0 is specified by the user.  The iteration
 *        will stop when
 *   
 *        (F^k - F^{k+1})/max{|F^k|,|F^{k+1}|,1} <= FACTR*EPSMCH
 *   
 *        where EPSMCH is the machine precision, which is automatically
 *        generated by the code. Typical values for FACTR: 1.d+12 for
 *        low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
 *        high accuracy.
 *      On exit FACTR is unchanged.
 *   
 *   PGTOL is a double precision variable.
 *      On entry PGTOL >= 0 is specified by the user.  The iteration
 *        will stop when
 *   
 *                max{|PG(i)|; i = 1, ..., n} <= PGTOL
 *   
 *        where PG(i) is the i-th component of the projected gradient.
 *      On exit PGTOL is unchanged.
 *   
 *   TASK is a working array of characters of length 60 indicating
 *      the current job when entering and quitting this subroutine.
 *      TASK must be initialized to "START" upon first call, then TASK
 *      must not be modified between calls.
 *   
 *   CSAVE is a working array of characters of length 60.  CSAVE must not
 *      be modified between calls. 
 *   
 *   ISAVE is an integer working array of length 3*N + 27.
 *      ISAVE must not be modified between calls.  On exit with
 *      TASK = "NEW_X", the following information is available:
 *        ISAVE[0] = the total number of intervals explored in the
 *                   search of Cauchy points;
 *        ISAVE[4] = the total number of skipped BFGS updates before
 *                   the current iteration;
 *        ISAVE[8] = the number of current iteration;
 *        ISAVE[9] = the total number of BFGS updates prior the current
 *                   iteration;
 *        ISAVE[11] = the number of intervals explored in the search of
 *                    Cauchy point in the current iteration;
 *        ISAVE[12] = the total number of function and gradient
 *                    evaluations;
 *        ISAVE[14] = the number of function value or gradient
 *                    evaluations in the current iteration;
 *        if ISAVE[15] = 0, then the subspace argmin is within the box;
 *        if ISAVE[15] = 1, then the subspace argmin is beyond the box;
 *        ISAVE[16] = number of free variables in the current iteration;
 *        ISAVE[17] = number of active constraints in the current iteration;
 *        N + 1 - ISAVE[18] = number of variables leaving the set of active
 *                            constraints in the current iteration;
 *        ISAVE[19] = number of variables entering the set of active
 *                    constraints in the current iteration.
 *        if ISAVE[23] = TRUE, then the initial X has been replaced by
 *                             its projection in the feasible set;
 *        if ISAVE[24] = TRUE, then the problem is constrained;
 *        if ISAVE[25] = TRUE, then each variable has upper and lower bounds;
 *   
 *   DSAVE is a double precision working array of length
 *           (2*M + 4)*N + (11*M + 8)*M + 29
 *      DSAVE must not be modified between calls.  On exit with
 *      TASK = "NEW_X", the following information is available:
 *	  DSAVE[0] = current 'THETA' in the BFGS matrix;
 *	  DSAVE[1] = F(X) in the previous iteration;
 *	  DSAVE[2] = FACTR*EPSMCH;
 *	  DSAVE[3] = 2-norm of the line search direction vector;
 *	  DSAVE[4] = the machine precision epsmch generated by the code;
 *	  DSAVE[6] = the accumulated time spent on searching for Cauchy points;
 *	  DSAVE[7] = the accumulated time spent on subspace minimization;
 *	  DSAVE[8] = the accumulated time spent on line search;
 *	  DSAVE[10] = the slope of the line search function at the current
 *        	      point of line search;
 *	  DSAVE[11] = the maximum relative step length imposed in line search;
 *	  DSAVE[12] = the infinity norm of the projected gradient;
 *	  DSAVE[13] = the relative step length in the line search;
 *	  DSAVE[14] = the slope of the line search function at the starting
 *        	      point of the line search;
 *	  DSAVE[15] = the square of the 2-norm of the line search direction
 *	  	      vector.
 *
 *   IPRINT is an integer variable that must be set by the user.
 *      It controls the frequency and type of output generated:
 *       IPRINT<0    no output is generated;
 *       IPRINT=0    print only one line at the last iteration;
 *       0<IPRINT<99 print also F and |proj G| every IPRINT iterations;
 *       IPRINT=99   print details of every iteration except N-vectors;
 *       IPRINT=100  print also the changes of active set and final X;
 *       IPRINT>100  print details of every iteration including X and G;
 *      When IPRINT > 0, the file iterate.dat will be created to
 *                       summarize the iteration.
 *   
 *
 *   References:
 *   
 *      [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
 *      memory algorithm for bound constrained optimization'',
 *      SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
 *   
 *      [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
 *      limited memory FORTRAN code for solving bound constrained
 *      optimization problems'', Tech. Report, NAM-11, EECS Department,
 *      Northwestern University, 1994.
 *   
 *      (Postscript files of these papers are available via anonymous
 *       ftp to ece.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
 *   
 *                          *  *  *
 *   
 *   NEOS, November 1994. (Latest revision April 1997.)
 *   Optimization Technology Center.
 *   Argonne National Laboratory and Northwestern University.
 *   Written by
 *                       Ciyou Zhu
 *   in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
 *   
 *                          *  *  *
 *
 */


#if 0
static char *str_copy(char *dst, const char *src, unsigned int len);
static char *str_copy(char *dst, const char *src, unsigned int len)
{
  unsigned int i;
  for (i=0 ; i<len && src[i] ; ++i) dst[i] = src[i];
  while (i<len) dst[i++] = ' ';
  dst[len] = 0;
  return dst;
}
#endif

static char *str_trim(char *s, unsigned int len);
static char *str_trim(char *s, unsigned int len)
{
  unsigned int i, j=len;
  for (i=0 ; i<len ; ++i) {
    if (s[i] == ' ') {
      if (i < j) j = i;
    } else {
      j = len;
    }
  }
  s[j] = 0;
  return s;
}

static char *str_trim(char *s, unsigned int len);
static char *str_pad(char *s, unsigned int len)
{
  unsigned int i;
  for (i=0 ; i<len && s[i]; ++i)
    ;
  while (i<len) s[i++] = ' ';
  s[len] = 0;
  return s;
}

integer_t lbfgsb_wrapper(integer_t n, integer_t m, real_t x[], real_t *f,
			 real_t g[], const real_t l[], const real_t u[],
			 const integer_t bnd[], real_t factr, real_t pgtol,
			 char task[], char csave[], integer_t isave[],
			 real_t dsave[], integer_t iprint)
{
  extern void mainlb_(integer_t *n, integer_t *m, real_t *x,
		      const real_t *l, const real_t *u, const integer_t *bnd,
		      real_t *f, real_t *g, real_t *factr, real_t *pgtol,
		      real_t *ws, real_t *wy, real_t *sy, real_t *ss,
		      real_t *wt, real_t *wn, real_t *snd,
		      real_t *z, real_t *r, real_t *d, real_t *t, real_t *wa,
		      integer_t *index, integer_t *iwhere, integer_t *indx2,
		      char *task, integer_t *iprint, char *csave,
		      logical_t *lsave, integer_t *isave, real_t *dsave,
		      ftnlen_t task_len, ftnlen_t csave_len);

 /* Notes about storage in working arrays:
  *
  *   - ISAVE[0..22]             -> ISAVE  in MAINLB (23 elements)
  *   - ISAVE[23..26]            -> LSAVE  in MAINLB  (4 elements)
  *   - ISAVE[26 + (1..N)]       -> INDEX  in MAINLB  (N elements)
  *   - ISAVE[26 + N + (1..N)]   -> IWHERE in MAINLB  (N elements)
  *   - ISAVE[26 + 2*N + (1..N)] -> INDX2  in MAINLB  (N elements)
  *
  *   - DSAVE[0..28]             -> DSAVE  in MAINLB (29 elements)
  *     index 27 of ISAVE;
  *   - arrays WS, ... get stored at index 29 of array DSAVE
  */
  integer_t l1   = m*n;
  integer_t l2   = m*m;
  integer_t l3   = 4*m*m;
  integer_t lws  = 29; /* see note */
  integer_t lwy  = lws + l1;
  integer_t lsy  = lwy + l1;
  integer_t lss  = lsy + l2;
  integer_t lwt  = lss + l2;
  integer_t lwn  = lwt + l2;
  integer_t lsnd = lwn + l3;
  integer_t lz   = lsnd + l3;
  integer_t lr   = lz + n;
  integer_t ld   = lr + n;
  integer_t lt   = ld + n;
  integer_t lwa  = lt + n;
  int c;
  logical_t lsave[4];

  /* copy LOGICAL values */
  lsave[0] = (isave[23] != 0);
  lsave[1] = (isave[24] != 0);
  lsave[2] = (isave[25] != 0);
  lsave[3] = (isave[26] != 0);

  /* call FORTRAN MAINLB subroutine with character arrays fixed */
  mainlb_(&n, &m, x, l, u, bnd, f, g, &factr, &pgtol,
	  &dsave[lws], &dsave[lwy],  &dsave[lsy], &dsave[lss], &dsave[lwt],
	  &dsave[lwn], &dsave[lsnd], &dsave[lz],  &dsave[lr],  &dsave[ld],
	  &dsave[lt], &dsave[lwa],
	  &isave[27], &isave[27 + n], &isave[27 + 2*n],
	  str_pad(task, LBFGSB_TASK_LENGTH), &iprint,
	  str_pad(csave, LBFGSB_TASK_LENGTH), lsave/*&isave[23]*/,
	  &isave[0], &dsave[0],
	  (ftnlen_t)LBFGSB_TASK_LENGTH, (ftnlen_t)LBFGSB_TASK_LENGTH);

  /* copy back LOGICAL values */
  isave[23] = lsave[0];
  isave[24] = lsave[1];
  isave[25] = lsave[2];
  isave[26] = lsave[3];
  
  /* fix character arrays and compute return value (JOB) */
  str_trim(task,  LBFGSB_TASK_LENGTH);
  str_trim(csave, LBFGSB_TASK_LENGTH);
  c = task[0];
  if (c == 'F' && task[1] == 'G')                     return 1;
  if (c == 'N' && ! strncmp(task, "NEW_X", 5))        return 2;
  if (c == 'C' && ! strncmp(task, "CONVERGENCE", 11)) return 3;
  if (c == 'W' && ! strncmp(task, "WARNING", 7))      return 4;
  if (c == 'S' && ! strncmp(task, "START", 5))        return 0;
  /*if (c == 'E' && ! strncmp(task, "ERROR", 5))*/    return 5;
}
