#include "main.h"

void pre_cg(double* b,double tol,int maxit,double* c,struct matrix* Omega12,struct matrix* P, double* p, int* flag, double* relres, int* iterk)
{
	int n = P->rows;
	int i;
	double* r = (double*) malloc (sizeof(double) *n);

	for(i=0; i<n; i++) {
		r[i] = b[i];
	}

	double n2b = norm(b, n);

	double tolb = tol * n2b;

	for(i=0; i<n; i++) p[i] = 0;

	*flag = 1;
	*iterk = 0;
	*relres = 1000;

	double* z = (double*) malloc(sizeof(double) * n);
	for(i=0; i<n; i++) {
		z[i] = r[i] / c[i];
	}

	double rz1 = 0;
	for(i=0; i<n; i++) {
		rz1 += r[i]*z[i];
	}

	double rz2 = 1;

	double *d = (double*) malloc(sizeof(double) * n);
	for(i=0; i<n; i++) {
		d[i] = z[i];
	}

	int k;

	struct matrix w;
	w.rows = n;
	w.columns = 1;
	w.entries = (double*) malloc(sizeof(double)*n);

	for(k=1; k<=maxit; k++)
	{
		if(k>1) {
			for(i=0; i<n; i++) {
				d[i] = z[i] + d[i]*rz1/rz2;
			}
		}

		Jacobian_matrix(d, Omega12, P, &w);

		double denom = 0;
		for (i=0; i<n; i++) {
			denom += d[i] * w.entries[i];
		}

		*iterk = k;

		double normr = norm(r, n);

		*relres = normr / n2b;

		if(denom<=0)
		{
			double normd = norm(d, n);
			for(i=0; i<n; i++) {
				p[i] = d[i] / normd;
			}
			break;
		}
		else
		{
			double alpha = rz1/denom;
			for(i=0; i<n; i++)
			{
				p[i] = p[i] + alpha * d[i];
				r[i] = r[i] - alpha * w.entries[i];
			}
		}

		for(i=0; i<n; i++) {
			z[i] = r[i] / c[i];
		}

		normr = norm(r, n);

		if(normr <= tolb)
		{
			*iterk = k;
			*relres = normr / n2b;
			*flag = 0;
			break;
		}
		rz2 = rz1;

		rz1 = 0;
		for(i=0; i<n; i++) {
			rz1 += r[i] * z[i];
		}
	}
	free(w.entries);
	free(r);
	free(z);
	free(d);
}

void Correlation_Newton(struct matrix* G, struct matrix* X, double* y)
{
	int i,j;
	int n = G->rows;

	double tau = 0;

	double* b;
	b = (double*) malloc(sizeof(double) * n);
	for (i=0; i<n; i++) b[i] = 1;

	double *b0;
	b0 = (double*) malloc (sizeof(double) * n);
	for (i=0; i<n; i++) b0[i] = b[i];

	//make G symmetric
	// Access upper triangle
	for (i=0; i<n; i++) {
		for(j=i+1; j<n; j++) {
			G->entries[j*n+i] = (G->entries[j*n+i]+G->entries[j+i*n])/2; // set upper tri to average of upper/lower
			G->entries[j+i*n] = G->entries[j*n+i]; // set lower tri to upper tri
		}
	}

	// Access diagonals
	for(i=0; i<n; i++) {
		G->entries[i*n+i] = G->entries[i*n+i] - tau;
	}

	for(i=0; i<n; i++) b0[i] = b0[i] - tau;

	double Res_b[300];
	for(i=0; i<300; i++) Res_b[i] = 0;

	//initial point
	for(i=0; i<n; i++) y[i] = 0;

	double* Fy;
	Fy = (double *) malloc (sizeof(double) * n);
	for(i=0; i<n; i++) Fy[i] = 0;

	int k=0;
	int f_eval = 0;

	int Iter_Whole = 200;
	int iter_Inner = 20;
	int maxit = 200;
	int iterk = 0;
	double tol = 1.0e-2;

	double error_tol = 1.0e-6;
	double sigma_1 = 1.0e-4;

	double * x0;
	x0 = (double *) malloc(sizeof(double) *n);
	for (i=0; i<n; i++) x0[i] = y[i];

	double *c = (double *) malloc (sizeof(double) *n);
	for (i=0; i<n; i++) c[i] =1 ;

	double *d = (double*) malloc(sizeof(double) *n );
	for(i=0; i<n; i++) d[i] = 0;

	double val_G = 0;
	for (i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			val_G += G->entries[j*n+i] * G->entries[j*n+i] / 2.;
		}
	}

	X->entries = (double*) malloc(sizeof(double) *n *n);
	// X = G + diag(y)
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			if(i == j) {
				X->entries[n*j+i] = G->entries[n*j+i] + y[i]; //
			} else {
				X->entries[n*j+i] = G->entries[n*j+i];
			}
		}
	}

	double *lambda = (double*) malloc(sizeof(double) * n);
	int info;
	char jobz = 'V';
	char uplo = 'L';
	int lwork = 1+6*n+2*n*n;
	int liwork = 3+5*n;

	//copy G as evec decomposition kills it..
	struct matrix P;
	P.rows = n;
	P.columns = n;
	P.entries = (double*) malloc (sizeof(double) * n * n);

	for (j=0; j<n; j++) {
		for (i=0; i<n; i++) {
			P.entries[j*n+i] = X->entries[j*n+i];
		}
	}

	double* a = (double*) malloc(sizeof(double) * lwork);
	int* aa = (int*) malloc(sizeof(int) * liwork);

	//G contains the evecs now, lambda contains evals; if info = 0, evals are in ascending order
	dsyevd_(&jobz, &uplo, &G->rows, P.entries, &G->rows, lambda, a, &lwork, aa, &liwork, &info);

	/* we want the evals in descending order; dsyev returns them in ascending order at best.
	 If that's the case, we simply revert the order.
	 In case they are scattered, we would need to sort it together with matrix P attached to it. This is a bit messy, tbd.. */
	if (info == 0)
	{
		for(j=0; j<n/2; j++)
		{
			if(lambda[j] < lambda[n-1-j])
			{
				double tmp = lambda[j];
				lambda[j] = lambda[n-1-j];
				lambda[n-1-j] = tmp;
				for(i=0; i<n; i++)
				{
					tmp = P.entries[n*j + i];
					P.entries[n*j+i] = P.entries[n*(n-1-j)+i];
					P.entries[n*(n-1-j)+i] = tmp;
				}
			}
		}
	}
	else
	{
		printf("Evals not in ascending order!!!\n");
		return;
	}

	double f0;
	gradient(y, lambda, &P, b0, &f0, Fy);

	double f = f0;

	f_eval = f_eval + 1;
	for(i=0; i<n; i++) {
		b[i] = b0[i] - Fy[i];
	}
	double norm_b = norm(b, n);

	double Initial_f = val_G - f0;

	struct matrix Omega12;
	omega_mat(lambda, n, &Omega12);
	for (i=0; i<n; i++) x0[i] = y[i];

	while(norm_b>error_tol && k< Iter_Whole)
	{
		precond_matrix(&Omega12, &P, c);

		int flag;
		double relres;

		pre_cg(b, tol, maxit, c, &Omega12, &P, d, &flag, &relres, &iterk);


		double slope = 0;
		for(i=0; i<n; i++) {
			slope += (Fy[i] - b0[i])*d[i];
		}

		for (i=0; i<n; i++) {
			y[i] = x0[i] + d[i];
		}

		// X = G + diag(y)
		for(i=0; i<n; i++) {
			for(j=0; j<n; j++) {
				if (i == j) {
					X->entries[n*j+i] = G->entries[n*j+i] + y[i];
				} else {
					X->entries[n*j+i] = G->entries[n*j+i];
				}
			}
		}

		// Copy X to P
		for (j=0; j<n; j++) {
			for (i=0; i<n; i++) {
				P.entries[j*n+i] = X->entries[j*n+i];
			}
		}

		//G contains the evecs now, lambda contains evals; if info = 0, evals are in ascending order
		dsyevd_(&jobz, &uplo, &G->rows, P.entries, &G->rows, lambda, a, &lwork, aa, &liwork, &info);
		//dsyev_(&jobz, &uplo, &G->rows, P.entries, &G->rows, lambda, a, &lwork, &info);

		if (info == 0)
		{
			for(j=0; j<n/2+1; j++)
			{
				if(lambda[j] < lambda[n-1-j])
				{
					double tmp = lambda[j];
					lambda[j] = lambda[n-1-j];
					lambda[n-1-j] = tmp;
					for(i=0; i<n; i++)
					{
						tmp = P.entries[n*j + i];
						P.entries[n*j+i] = P.entries[n*(n-1-j)+i];
						P.entries[n*(n-1-j)+i] = tmp;
					}
				}
			}
		}
		else
		{
			printf("Evals not in ascending order!!!\n");
			return;
		}

		gradient(y, lambda, &P, b0, &f, Fy);

		int k_inner = 0;

		while(k_inner<=iter_Inner && f>f0 + sigma_1*pow(0.5, k_inner) * slope	+ 10e-6)
		{
			k_inner++;
			for (i=0; i<n; i++) {
				y[i] = x0[i] + pow(.5, k_inner) * d[i];
			}

			for(i=0; i<n; i++) {
				for(j=0; j<n; j++) {
					if (i==j) {
						X->entries[n*j+i] = G->entries[n*j+i] + y[i];
					} else {
						X->entries[n*j+i] = G->entries[n*j+i];
					}
				}
			}

			for (j=0; j<n; j++) {
				for (i=0; i<n; i++) {
					P.entries[j*n+i] = X->entries[j*n+i];
				}
			}

			//G contains the evecs now, D contains evals; if info = 0, evals are in ascending order
			dsyevd_(&jobz, &uplo, &G->rows, P.entries, &G->rows, lambda, a, &lwork, aa, &liwork, &info);
			//dsyev_(&jobz, &uplo, &G->rows, P.entries, &G->rows, lambda, a, &lwork, &info);

			eig_time = eig_time + (double) clock()/(double) CLOCKS_PER_SEC - eig_time0;

			if (info == 0)
			{
				for(j=0; j<n/2+1; j++)
				{
					if(lambda[j] < lambda[n-1-j])
					{
						double tmp = lambda[j];
						lambda[j] = lambda[n-1-j];
						lambda[n-1-j] = tmp;
						for(i=0; i<n; i++)
						{
							tmp = P.entries[n*j + i];
							P.entries[n*j+i] = P.entries[n*(n-1-j)+i];
							P.entries[n*(n-1-j)+i] = tmp;
						}
					}
				}
			}
			else
			{
				printf("Evals not in ascending order!!!\n");
				return;
			}
			gradient(y, lambda, &P, b0, &f, Fy);
		}

		f_eval = f_eval + k_inner + 1;

		for(i=0; i<n; i++) {
			x0[i] = y[i];
		}
		f0 = f;

		k++;

		for(i=0; i<n; i++) {
			b[i] = b0[i] - Fy[i];
		}

		norm_b = norm(b,n);

		Res_b[k] = norm_b;

		if(Omega12.rows != 0) {
			free(Omega12.entries);
		}
		omega_mat(lambda, n, &Omega12);
	}

	int r = 0;
	while(lambda[r]>0 && r<n) r++;
	// TODO
	if (r == 0) {
		for(i=0; i<n*n; i++) {
			X->entries[i] = 0; // case 0
		}
	} else if (r == n) {
		// do nothing
	} else {
		if((double) r<=(double) n /2.)
		{
			if(r > 1) // case 2 (r > 1)
			{
				for (j=0; j<r; j++)	{
					for(i=0; i<n; i++) {
						P.entries[i+j*n] *= sqrt(lambda[j]);
					}
				}
				cblas_dgemm(
					CblasColMajor, CblasNoTrans, CblasTrans,
					// A[MxK], B[KxN], C[MxN]
					n, n, r,      // M, N, K
					1,            // alpha
					P.entries, n, // A[n,r] "first r columns"
					P.entries, n, // B'[r,n]
					0,            // beta
					X->entries, n // Matrix C[n,n]
				);
			}
			else // case 1 (r == 1)
			{
				cblas_dgemm(
					CblasColMajor, CblasNoTrans, CblasTrans,
					// A[MxK], B[KxN], C[MxN]
					n, n, r,             // M, N, K
					lambda[0]*lambda[0], // alpha
					P.entries, n,        // Matrix A[n,r] "first column"
					P.entries, n,        // Matrix B[r,n]
					0,                   // beta
					X->entries, n        // Matrix C[n,n]
				);
			}
		}
		else // case 3 (r > n/2)
		{
			for (j=0; j<n-r; j++) {
				for(i=0; i<n; i++) {
					P.entries[i+(j+r)*n] *= sqrt(-lambda[r+j]);
				}
			}
			cblas_dgemm(
				// C = A*B + C
				CblasColMajor, CblasNoTrans, CblasTrans,
				// A[MxK], B[KxN], C[MxN]
				n, n, n-r,        // M, N, K
				1,                // alpha
				P.entries+n*r, n, // Matrix A[n,n-r] "last n-r columns"
				P.entries+n*r, n, // Matrix B[n-r,n]
				1,                // beta
				X->entries, n     // Matrix C[n,n]
			);
		}
	}

	free(b);
	free(b0);
	free(Fy);
	free(c);
	free(d);
	free(lambda);
	free(P.entries);
	free(a);
	free(aa);

	double Final_f = val_G - f;
	double val_obj = 0;
	for(i=0; i<n*n; i++) {
		val_obj += pow(X->entries[i] - G->entries[i], 2)/2;
	}

	for(i=0; i<n; i++)
		X->entries[i+n*i] += tau;
}

void precond_matrix (struct matrix* Omega12, struct matrix * P, double* c)
{
	int r = Omega12->rows;
	int s = Omega12->columns;
	int n = P->rows;
	int i,j;

	for(i=0; i<n; i++) c[i] = 1.;

	if(r==n) return;

	//H=P.*P;
	struct matrix H;
	H.rows = n;
	H.columns = n;
	H.entries = (double*) malloc(sizeof(double) *n *n);
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			H.entries[j*n + i ] = P->entries[j*n+i]*P->entries[j*n+i];
		}
	}

	if(r>0)
	{
		if((double) r < (double) n /2.)
		{
			struct matrix H12;
			H12.rows = n;
			H12.columns = s;
			H12.entries = (double*) malloc (sizeof(double) * n *s);

			//H12 = H(1:r,:)'*Omega12;
			cblas_dgemm(
				CblasColMajor, CblasNoTrans, CblasNoTrans,
				// A[MxK], B[KxN], C[MxN]
				n, s, r,				// M,N,K
				1,                      //
				H.entries, n,			// A[n,r]
				Omega12->entries, r,	// B[r,s]
				0,						//
				H12.entries, n 			// C[n,s]
			);

			for (i=0; i<n; i++)
			{
				// c(i) = sum(H(1:r,i))*(d'*H(1:r,i));
				c[i] = 0;
				for(j=0; j<r; j++) {
					c[i] += H.entries[i+j*n];
				}
				c[i] *= c[i];

				//c(i) = c(i) +2.0*(H12(i,:)*H(r+1:n,i));
				for(j=0; j<s; j++) {
					c[i] += 2.0 * H12.entries[i+j*n]*H.entries[i+(r+j)*n];
				}

				if(c[i] < 1.0e-8) c[i] = 1.0e-8;
			}
			free(H12.entries);
		}
		else
		{
			struct matrix H12;
			H12.rows = r;
			H12.columns = n;
			H12.entries = (double*) malloc (sizeof(double) * r *n);

			struct matrix Omega;
			Omega.rows = Omega12->rows;
			Omega.columns = Omega12->columns;
			Omega.entries = (double*) malloc(sizeof(double) * Omega.rows * Omega.columns);

			for(i=0; i<r; i++) {
				for(j=0; j<s; j++) {
					Omega.entries[i+j*r] = 1. - Omega12->entries[i+j*r];
				}
			}

			//H12 = Omega12*H(r+1:n,:);
			cblas_dgemm(
				CblasColMajor, CblasNoTrans, CblasTrans,
				// A[MxK], B[KxN], C[MxN]
				r, n, s,			// M,N,K
				1,					//
				Omega.entries, r,	// A[r,s]
				H.entries+r*n, n,	// B[n,s]' -> B[s,n]
				0,					//
				H12.entries, r 		// C[r,n]
			);

			for (i=0; i<n; i++)
			{
				// c(i) = sum(H(r+1:n,i))*(d'*H(r+1:n,i));
				double temp = 0;
				c[i]=0;
				for(j=0; j<s; j++) {
					c[i] += H.entries[i+(j+r)*n];
				}
				c[i] *= c[i];

				// c(i) = c(i) + 2.0*(H(1:r,i)'*H12(:,i));
				for(j=0; j<r; j++) {
					c[i] += 2.0 * H.entries[i+n*j]*H12.entries[i*r+j];
				}

				double alpha = 0;
				for (j=0; j<n; j++) {
					alpha += H.entries[i+j*n];
				}

				temp = 0;
				for(j=0; j<n; j++) {
					temp += H.entries[i+j*n];
				}

				c[i] = -c[i] + alpha*temp;

				if(c[i] < 1.0e-8) c[i] = 1.0e-8;
			}
			free (H12.entries);
			free(Omega.entries);
		}

	}
	free(H.entries);
}

void gradient(double *y, double* lambda, struct matrix* P, double* b0, double* f, double* Fy)
{
	*f = 0;
	int n = P->rows;

	struct matrix Q;
	Q.rows = P->rows;
	Q.columns = P->columns;
	Q.entries = (double*) malloc(sizeof(double) * Q.rows * Q.columns);

	int i,j;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{	// Let Q[j,i] = P[j,i] then let Q[:,i] *= sqrt(lambda_nonneg[i])
			Q.entries[j+i*n] = P->entries[j+i*n] * sqrt(max(lambda[i], 0));
		}
	}


	for (i=0; i<n; i++)
	{
		Fy[i] = 0;
		for (j=0; j<n; j++) {
			Fy[i] += Q.entries[i+j*n]*Q.entries[i+j*n];
		}
	}

	for (i=0; i<n; i++)
	{
		*f += pow(max(lambda[i], 0), 2.);
	}

	double temp=0;
	for (i=0; i<n; i++) {
		temp+= b0[i]*y[i];
	}
	*f = *f/2. - temp;
	free(Q.entries);
}

// assume lambda is already sorted in descending order
void omega_mat(double* lambda, int n, struct matrix * Omega12)
{
	int r=0;
	while(lambda[r]>0 && r<n) r++;
	int s = n - r;

	if (r==0)
	{
		Omega12->rows = 0;
		Omega12->columns = 0;
		Omega12->entries = NULL;
	}
	else if (r == n)
	{
		int i=0;
		Omega12->rows = n;
		Omega12->columns = n;
		Omega12->entries = (double*) malloc (sizeof(double)*n*n);
		for (i=0; i<n*n; i++) Omega12->entries[i] = 1;
	}
	else
	{
		int i,j;
		Omega12->rows = r;
		Omega12->columns = s;
		Omega12->entries = (double*) malloc(sizeof(double)*r*s);
		for(j=0; j<s; j++) {
			for(i=0; i<r; i++) {
				Omega12->entries[j*r+i] = lambda[i] / (lambda[i] - lambda[r+j]);
			}
		}
	}
}

void Jacobian_matrix (double* x, struct matrix* Omega12, struct matrix* P, struct matrix* Ax)
{
	int i,j;
	int n = P->rows;
	int r=Omega12->rows;
	int s=Omega12->columns;

	for (i=0; i<n; i++) Ax->entries[i] = 0;

	struct matrix H1;
	H1.rows = n;
	H1.columns = r;
	H1.entries = (double *) malloc (sizeof(double) * n*r);

	struct matrix temp;

	struct matrix Omega;
	Omega.rows = r;
	Omega.columns = s;
	Omega.entries = (double*) malloc(sizeof(double) * r * s);

	if (r>0)
	{
		if (r< ((double) n)/2.0) // r < s
		{
			temp.rows = r;
			temp.columns = s;
			temp.entries = (double *) malloc(sizeof(double) * r *s);

			for(i = 0; i < n; i++) {
				for(j = 0; j < r; j++) {
					H1.entries[j*n + i] = x[i] * P->entries[j*n + i];
				}
			}

			//Omega12 = Omega12.*(H1'*P(:,r+1:n))
			// cblas(_, _, _, n, n, _, alpha, , P, _, P, _, beta, X, n)
			// A[MxK], B[KxN], C[MxN]
			// C := aAB + bC
			cblas_dgemm(
				CblasColMajor, CblasTrans, CblasNoTrans,
				// A[MxK], B[KxN], C[MxN]
				r, s, n,           // M, N, K
				1,                 // alpha
				H1.entries, n,     // Matrix A[n, r]' -> A[r,n]
				P->entries+n*r, n, // Matrix B[n,s]
				0,                 // beta
				temp.entries, r    // Matrix C[r,s]
			);
			for (i=0; i<r*s; i++) {
				Omega.entries[i] = Omega12->entries[i] * temp.entries[i];
			}

			free(temp.entries);

			struct matrix HT;
			HT.rows = n;
			HT.columns = n;
			HT.entries = (double *) malloc(sizeof(double) * n*n);
			for (i=0; i<n*n; i++) {
				HT.entries[i] = 0; //make sure it's 0s
			}

			// P1^T * H1
			temp.rows = r;
			temp.columns =r;
			temp.entries = (double*) malloc(sizeof(double) * r * r);
			// [r,r] = [r,n]x[n,r]
			cblas_dgemm(
				CblasColMajor, CblasTrans ,CblasNoTrans,
				// A[MxK], B[KxN], C[MxN]
				r, r, n,        // M, N, K
				1,              // alpha
				P->entries, n,  // Matrix A[n,r]' -> [r,n]
				H1.entries, n,  // Matrix B[n,r]
				0,              // beta
				temp.entries, r // Matrix C[r,r]
			);

			// HT += P1 * P1^T * H1
			// [n,n] = [n,r]x[r,r]
			cblas_dgemm(
				CblasColMajor, CblasNoTrans, CblasNoTrans,
				// A[MxK], B[KxN], C[MxN]
				n, r, r,         // M, N, K
				1,               // alpha
				P->entries, n,   // Matrix A[n,r]
				temp.entries, r, // Matrix B[r,r]
				0,               // beta
				HT.entries, n    // Matrix C[n,r]
			);

			free(temp.entries);

			//HT = P1 * P1^T * H1 + P2 * Omega^T
			cblas_dgemm(
				CblasColMajor, CblasNoTrans, CblasTrans,
				// A[MxK], B[KxN], C[MxN]
				n, r, s, 			// M, N, K
				1,					// alpha
				P->entries+n*r, n,	// Matrix A[n,s]
				Omega.entries, r,	// Matrix B[r,s]' -> B[s,r]
				1,					// beta
				HT.entries, n 		// Matrix C[n,r]
			);

			//HT = P1 * P1^T * H1 + P2 * Omega^T ; P1 * Omega^T
			cblas_dgemm(
				CblasColMajor, CblasNoTrans, CblasNoTrans,
				// A[MxK], B[KxN], C[MxN]
				n, s, r,			// M, N, K
				1,					// alpha
				P->entries, n,		// Matrix A[n,r]
				Omega.entries, r,	// Matrix B[r,s]
				0,					// beta
				HT.entries+n*r, n 	// Matrix C[n,s]
			);


			for (i=0; i<n; i++)
			{
				for (j=0; j<n; j++)
					Ax->entries[i] += P->entries[j*n+i]*HT.entries[j*n+i];
				Ax->entries[i] += PERTURBATION*x[i];
			}
			free(HT.entries);
			free(H1.entries);
		}
		else
		{
			if (r == n){
				for (i=0; i<n; i++) {
					Ax->entries[i] = x[i] * (1.+PERTURBATION);
				}
			}
			else // r >= s
			{
				H1.rows = n;
				H1.columns = s;
				free(H1.entries);
				H1.entries = (double *) malloc (sizeof(double) * n *s);

				for(i=0; i<n; i++) {
					for(j=0; j<s; j++) {
						H1.entries[j*n + i] = x[i] * P->entries[(j+r)*n + i];  //H1 = diag(x) * P2
					}
				}

				//Omega12 = ones(r,s)-Omega12;
				for(i=0; i<r; i++) {
					for (j=0; j<s; j++) {
						Omega.entries[j*r+i] = 1. - Omega12->entries[j*r+i];
					}
				}


				temp.rows = r;
				temp.columns = s;
				temp.entries = (double *) malloc(sizeof(double) * r *s);

				//Omega12 = Omega12.*((P(:,1:r))'*H2);
				cblas_dgemm(
					CblasColMajor, CblasTrans, CblasNoTrans,
					// A[MxK], B[KxN], C[MxN]
					r, s, n,			// M, N, K
					1,					// alpha
					P->entries, n,		// Matrix P[n,r]' -> P[r,n]
					H1.entries, n,      // Matrix H1[n,s]
					0,					// beta
					temp.entries, r 	// Matrix M[r,s]
				);
				for (i=0; i<r*s; i++) {
					Omega.entries[i] *= temp.entries[i];
				}
				free(temp.entries);

				struct matrix HT;
				HT.rows = n;
				HT.columns = n;
				HT.entries = (double *) malloc(sizeof(double) * n*n);
				for (i=0; i<n*n; i++) HT.entries[i] = 0; //make sure it's 0s

				//HT += P2 * Omega^T
				cblas_dgemm(
					CblasColMajor, CblasNoTrans, CblasTrans,
					// A[MxK], B[KxN], C[MxN]
					n, r, s, 			// M, N, K
					1,					//
					P->entries+n*r, n,	// A[n,s]
					Omega.entries, r,	// B'[s,r]
					0,					//
					HT.entries, n 		// C[n,r]
				);

				cblas_dgemm(
					CblasColMajor, CblasTrans, CblasNoTrans,
					// A[MxK], B[KxN], C[MxN]
					s, s, n,			// M, N, K
					1,					//
					H1.entries, n,		// A'[s,n]
					P->entries+n*r, n,	// B[n,s]
					0,					//
					temp.entries, s 	// C[s,s]
				);

				cblas_dgemm(
					CblasColMajor, CblasNoTrans, CblasNoTrans,
					// A[MxK], B[KxN], C[MxN]
					n, s, s, 			// M, N, K
					1,					//
					P->entries+n*r, n,	// A[n,s]
					temp.entries, s,	// B[s,s]
					0,					//
					HT.entries+r*n, n 	// C[n,s]
				);


				cblas_dgemm(
					// C = AB + C
					CblasColMajor, CblasNoTrans, CblasNoTrans,
					// A[MxK], B[KxN], C[MxN]
					n, s, r,			// M, N, K
					1,					//
					P->entries, n,		// A[n,r]
					Omega.entries, r,	// B[r,s]
					1,					//
					HT.entries+r*n, n 	// C[n,s]
				);

				for (i=0; i<n; i++)
				{
					Ax->entries[i] = 0;
					for (j=0; j<n; j++)
						Ax->entries[i] -= P->entries[j*n+i]*HT.entries[j*n+i];
					Ax->entries[i] += x[i] *(1. + PERTURBATION);
				}
				free(HT.entries);
				free(H1.entries);
			}
		}
	}
	free(Omega.entries);
}

int main(int argc, const char * argv[])
{
	int i;
	int N = 3000;
	struct matrix G;
	G.rows = N;
	G.columns = N;
	G.entries = (double*) malloc(sizeof(double) * N * N);
	srand((int) time(NULL));
	for (i=0; i<N*N; i++) {
		G.entries[i] = 2.* (double) rand()/ (double) RAND_MAX-1.;
	}

	struct matrix X;
	X.rows = N;
	X.columns = N;

	double y[N];

	Correlation_Newton(&G, &X, y);

	//printMatrix(&X);

	return 0;
}
