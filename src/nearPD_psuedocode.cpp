// main()
    // G[n,n]                                                                  (input correlation matrix)
    // X[n,n] output correlation matrix
    // y[n]   ??
    // nearPD(G, X, y)


// nearPD(G[n,n], X[n,n], y[n])
    // tau = 0
    // b[n] = 1
    // b0[n] = b[n]
    // Make G symmetric
    // Set diagonals of G to G - tau
        // G.diag() -= tau
    // b0 = b0 - tau
    // Res_b[300] = 0
    // y[n] = 0
    // Fy[n] = 0
    // k = 0
    // f_eval = 0
    // iter_outer = 200
    // iter_inner = 20
    // maxit = 200
    // iterk = 0
    // tol = 1.0e-2
    // err_tol = 1.0e-6
    // sigma1 = 1.0e-4
    // x0[n] = y[n]
    // c[n] = 1
    // d[n] = 0
    // val_G = sum( G[i,j]*G[i,j] / 2 )
    // X = G + diagmat(y)
    // lambda[n]
    // P[n,n]
    // lambda[n], P[n,n] = eigenvalues, eigenvectors of X
    // Sort the eigenvalues in descending order and rearrange columns of P accordingly
        // lambda[n] = reverse(lambda[n])
        // P = fliplr(P[n,n])
    // f0
    // gradient(y, lambda, P, b0, f0, Fy)                                      (modifies f0 and Fy)
    // f = f0
    // f_eval += 1
    // b[n] = b0[n] - Fy[n]
    // normb = norm(b, 2)
    // f_init = val_G - f0
    // Initialize matrix Omega12[r,s] or Omega12[n,n]                          (num of pos eigval by num of nonpos eigval)
        // omega_mat(lambda, n, Omega12)                                       (Omega12 modified)
    // x0[n] = y[n]
    // while (normb > err_tol) and (k < iter_outer)
        // precondition matrix
            // precond_matrix(Omega12[r,s], P[n,n], c[n])                      (modifies c)
        // flag
        // relres
        // pre_cg(b[n], tol, maxit, c[n], Omega12[r,s], P[n,n], d[n],
        //        flag, relres, iterk)                                         (modifies d[n], flag, relres, iterk)
        // slope = sum( (Fy[n] .- b0[n]) .* d[n] )
        // y[n] = x0[n] + d[n]
        // X[n,n] = G[n,n] + diagmat(y[n])
        // P[n,n] = X[n,n]
        // lambda[n], P[n,n] = eigenvalues, eigenvectors of X
        // Sort the eigenvalues in descending order and rearrange columns of P accordingly
            // lambda[n] = reverse(lambda[n])
            // P = fliplr(P[n,n])
        // gradient(y, lambda, P, b0, f0, Fy)                                  (modifies f0 and Fy)
        // k_inner = 1
        // while (k_inner <= iter_inner) and (f > f0 + sigma1*pow(0.5, k_inner)*slope + 1.0e-5)
            // y[n] = x0[n] .+ pow(0.5, k_inner) .* d[n]
            // X[n,n] = G[n,n] + diagmat(y[n])
            // P[n,n] = X[n,n]
            // lambda[n], P[n,n] = eigenvalues, eigenvectors of X
            // Sort the eigenvalues in descending order and rearrange columns of P accordingly
                // lambda[n] = reverse(lambda[n])
                // P = fliplr(P[n,n])
            // gradient(y, lambda, P, b0, f0, Fy)                              (modifies f0 and Fy)
            // k_inner++
        // f_eval += k_inner + 1
        // x0[n] = y[n]
        // f0 = f
        // b[n] = b0[n] .- Fy[n]
        // normb = norm(b, 2)
        // Res_b[k+1] = normb
        // omega_mat(lambda, n, Omega12)                                       (modifies Omega12)
        // k++
    // r = 0
    // lambda is assumed to be in descending order
    // Find the number of strictly positive eigenvalues
        // while(lambda[r] > 0 && r < n) r++                                   (there has to be a better way to do this)
    // s = n - r
    // if (r == 0)
        // X[n,n] = 0
    // else if (r == n)
        // do nothing
    // else if (r == 1)
        // X[n,n] = lambda[0].*lambda[0] .* ( P[n,r] * P[n,r]')                (first column of P)
    // else if (r <= s)
        // for i in [0:n), j in [0:r)
            // P[i,j] *= sqrt(lambda[j])                                       (first r columns)
        // X[n,n] = P[n,r] * P[n,r]'                                           (first r columns)
    // else (r > s)
        // for i in [0:n), j in [0:s)
            // P[i,j+r] *= sqrt( -1.0 * lambda[r+j] )                          (last s eigen-values, vectors)
        // X[n,n] += P[n,s] * P[n,s]'                                          (last s columns)
    // Final_f = val_G - f
    // val_obj = sum( pow(X[n,n] .- G[n,n], 2) / 2.0 )
    // X.diag() += tau                                                         (recall that we subtract tau from G.diag() at the beginning)
    // return X


// gradient(y[n], lambda[n], P[n,n], b0[n], f, Fy[n])                          (modifies f and Fy)
    // Q[n,n]
    // lambda_nonneg[n] = max(lambda[n], 0)
    // Let Q[i,j] = P[i,j]
    // let Q[:,j] *= sqrt(lambda_nonneg[j])
    // Fy[i] = sum( Q[i,:] .* Q[i,:] )
    // f = sum( lambda_nonneg[i].^2 ) / 2.0 - sum(b0[i] .+ y[i])


// omega_mat(lambda[n], n, Omega12[?,?])
    // r = 0
    // lambda is assumed to be in descending order
    // Find the number of strictly positive eigenvalues
        // while(lambda[r] > 0 && r < n) r++                                   (there has to be a better way to do this)
    // s = n - r                                                               (number of non-positive eigenvalues)
    // if (r == 0) then Omega12[0,0]                                           (empty matrix)
    // else if (r == n) then Omega12[n,n]                                      (all eigenvalues are positive)
        // Omega12[i,j] = 1
    // else Omega12[r, s]                                                      (if r is the number of positive eigenvalues, then n-r is the number of non-positive values)
        // Think of lambda[i] as the positive eigenvalues and lambda[r+j] as the negative ones
        // Omega12[i,j] = lambda[i] / (lambda[i] - lambda[r+j])                ( positive / (poitive - negative) )
    // return Omega12


// precond_matrix(Omega12[r,s], P[n,n], c[n])                                  (modifies c)
    // c[n] = 1
    // if (r == n) do not modify anything (return early)
    // if (r == 0) then nothing happens. Does this mean the whole algo fails?
    // H[n,n] = P[i,j] * P[i,j]                                                (TODO: Verify this assignment. The original code says: H=P'; H=H.*H;)
    // if (r < n/2.0)
        // H12[n,s] = H[n,r] * Omega12[r,s]                                    (use first r columns of H)
        // for (i = 0; i < n; i++)
            // c[i] = 0
            // for (j = 0; j < r; j++)
                // c[i] += H[i,j]                                              (rowsum of first r columns)
            // c[i] *= c[i]
            // for (j = 0; j < s; j++)
                // c[i] += 2.0 * H12[i,j]*H[i,r+j]
            // c[i] = max(c[i], 1.0e-8)
    // else
        // Omega[r,s] = 1.0 - Omega12
        // H12[r,n] = Omega[r,s] * H[n,s]'                                     (use last s columns of H)
        // for (i = 0; i < n; i++)
            // c[i] = 0
            // for (j = 0; j < s; j++)
                // c[i] += H[i, j+r]                                           (rowsum of last s columns)
            // c[i] *= c[i]
            // for (j = 0; j < r; j++)
                // c[i] += 2.0 * H[i,j] * H12[j,i]
            // alpha = sum(H[i,:])                                             (rowsum of H)
            // c[i] = -1.0 * c[i] + alpha*alpha
            // c[i] = max(c[i], 1.0e-8)


// pre_cg(b[n], tol, maxit, c[n], Omega12[r,s], P[n,n], p[n],
//        flag, relres, iterk)                                                 (modifies p, flag, relres, iterk)
    // r[n] = b[n]
    // n2b = norm(b[n], 2)
    // tolb = tol * n2b
    // p[n] = 0
    // z[n] = r[n] ./ c[n]
    // rz1 = sum(r[n] .* z[n])
    // rz2 = 1.0
    // d[n] = z[n]
    // w[n]
    // for (k = 1; k <= maxit; k++)
        // if (k > 1)
            // d[n] = z[n] + d[n] * rz1/rz2
        // Compute the Jacobian matrix of [P?] and store in w
            // Jacobian(d, Omega12, P, w)                                      (modifies w)
        // denom = sum(d[n] .* w[n])
        // iterk = k
        // normr = norm(r, 2)
        // relres = normr / n2b
        // if (denom <= 0)
            // normd = norm(d, 2)
            // p[n] = d[n] / normd
            // return
        // else
            // alpha = rz1 / denom
            // p[n] += alpha .* d[n]
            // r[n] -= alpha .* w[n]
        // z[n] = r[n] / c[n]
        // normr = norm(r, 2)
        // if (normr <= tolb)
            // iterk = k
            // relres = normr / n2b
            // flag = 0
            // return
        // else
            // rz2 = 0.0
            // rz1 = sum(r[n] .* z[n])


// Jacobian_Matrix(x[n], Omega12[r,s], P[n,n], Ax[n,1])                        (modifies Ax)
    // MAYBE CONSIDER PARTITIONING P INTO P1[n,r] and P2[n,s] where [P1 P2] = P
    // To save memory, see if you can create references to the two partitions rather than copies
        // P1[n,r] = P.head_cols(r)
        // P2[n,s] = P.tail_cols(s)
    // Ax[n,1] = 0
    // Omega[r,s]
    // if (r == 0)
    // else if (r == n)
        // Ax[n] = x[n] * (1.0 + PERTURBATION)
    // else if (r < s)
        // H1[n,r]
        // M[r,s]
        // for i in [0:n), j in [0:r)
            // H1[i,j] = x[i] * P1[i,j]
        // M[r,s] = H1[n,r]' * P2[n,s]
        // Omega[r,s] = Omega12[i,j] * M[i,j]
        // HT[n,n] = 0
        // HT[n,r] += P1[n,r] * P1[n,r]' * H1[n,r]                             (use first r columns of HT)
        // HT[n,r] += P2[n,s] * Omega[r,s]'                                    (use first r columns of HT)
        // HT[n,s] += P1[n,r] * Omega[r,s]                                     (use last s columns of HT)
        // for i in [0:n)
            // Ax[i] = rowsum(P[i,:] .* HT[i,:]) + x[i]*PERTURBATION
    // else (r >= s)
        // H1[n,s]
        // for i in [0:n), j in [0:s)
            // H1[i,j] = x[i] * P2[i,j]
        // Omega[r,s] = 1.0 - Omega12[r,s] + P1[n,r]' * H1[n,s]
        // HT[n,n] = 0
        // HT[n,r] += P2[n,s] * Omega[r,s]'                                    (use first r columns of HT)
        // HT[n,s] += P2[n,s] * H1[n,s]' * P2[n,s]                             (use last s columns of H1, HT)
        // HT[n,s] += P1[n,r] * Omega[r,s]                                     (use last s columns of HT)
        // for i in [0:n)
            // Ax[i] = rowsum(P[i,:] .* HT[i,:]) + x[i] * (1.0 + PERTURBATION)
