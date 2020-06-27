/* This code is designed to solve min 0.5*<X-G, X-G> s.t. X_ii =1, i=1,2,...,n		*/
/* and X>=tau*I (symmetric and positive semi-definite)								*/
/* based on the algorithm  in "A Quadratically Convergent Newton Method for			*/
/* Computing the Nearest Correlation Matrix											*/
/* By Houduo Qi and Defeng Sun														*/
/* SIAM J. Matrix Anal. Appl. 28 (2006) 360--385.									*/

/* This particular C implementation is a result of the summer research project		*/
/* of Pawel Zaczkowski, Imperial College London,									*/
/* with Professor Defeng Sun, National University of Singapore						*/
 
/* Last modified date:  August 11, 2010												*/
/* The  input argument is the given symmetric G										*/
/* The outputs are the optimal primal and dual solutions							*/
/* Diagonal Preconditioner is added													*/

/* Please send your comments and suggestions to										*/
/* pawelzaczkowski@cantab.net or matsundf@nus.edu.sg								*/

/* Warning: Accuracy may not be guaranteed!!!!!										*/

#include <Accelerate/Accelerate.h> 
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

#define PERTURBATION 1.0e-9

/* declare a structure for keeping dimensions and entries of matrices together		*/
struct matrix
{
	double *entries;
	int rows;
	int columns;
};

/* function declarations */

void Correlation_Newton(struct matrix* G, struct matrix* X, double* y);
/* PURPOSE: calculating the nearest correlation matrix of G							*/
/* INPUT:	struct matrix* G														*/
/* OUTPUT:	struct matrix* X														*/
/* OUTPUT:	double* y																*/
/* based on the algorithm  in "A Quadratically Convergent Newton Method for			*/
/* Computing the Nearest Correlation Matrix											*/
/* By Houduo Qi and Defeng Sun														*/
/* SIAM J. Matrix Anal. Appl. 28 (2006) 360--385.									*/

void pre_cg(double* b,double tol,int maxit,double* c,struct matrix* Omega12,struct matrix* P, double* p, int* flag, double* relres, int* iterk);
/* PURPOSE:	PCG method																*/ 
/* INPUT:   double* b																*/ 
/* INPUT:   double tol																*/
/* INPUT:   int maxit																*/
/* INPUT:   double* c																*/
/* INPUT:   struct matrix* Omega12													*/
/* INPUT:   struct matrix* P														*/
/* OUTPUT:  double* p																*/
/* OUTPUT:  int* flag																*/
/* OUTPUT:  double* relres															*/
/* OUTPUT:  int* iterk																*/
/* This is exactly the algorithm by  Hestenes and Stiefel (1952)					*/
/* An iterative method to solve A(x) =b												*/
/* The symmetric positive definite matrix M is a preconditioner for A.				*/
/* See Pages 527 and 534 of Golub and va Loan (1996)								*/

void precond_matrix (struct matrix* Omega12, struct matrix * P, double* c);
/* PURPOSE: generating the diagonal preconditioner									*/
/* INPUT:	struct matrix* Omega12													*/
/* INPUT:	struct matrix* P														*/
/* OUTPUT:	double* c																*/

void gradient(double *y, double* lambda, struct matrix* P, double* b0, double* f, double* Fy);
/* PURPOSE: generating F(y)															*/
/* INPUT:	double* y																*/
/* INPUT:	double* lambda															*/
/* INPUT:	struct matrix* P														*/
/* INPUT:	double* b0																*/
/* OUTPUT:	double* f																*/
/* OUTPUT:	double* Fy																*/

void omega_mat(double* lambda, int n, struct matrix * Omega12);
/* PURPOSE: generating the essential part of the first-order difference d			*/
/* INPUT:	double* lambda															*/
/* INPUT:	int n																	*/
/* OUTPUT:	struct matrix* Omega12													*/

void Jacobian_matrix (double* x, struct matrix* Omega12, struct matrix* P, struct matrix* Ax);
/* PURPOSE: generating Jacobian matrix												*/
/* INPUT:	double* x																*/
/* INPUT:	struct matrix* Omega12													*/
/* INPUT:	struct matrix* P														*/
/* OUTPUT:	struct matrix* Ax														*/

void printMatrix(struct matrix* M);
/* PURPOSE: prints a given matrix on screen											*/
/* INPUT:   struct matrix* M														*/ 
 
double max(double a, double b);
double norm(double *a, int n);