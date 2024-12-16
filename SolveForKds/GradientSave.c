void GradientEquil(double *kds, double *l, double *bg, double *dil, double *A, double *a,
                   double *sXc, int *n_x, double *gradient,
                   double *dmdb, double*dLdm)
{
	int i, j;
	double x[*n_x];
	double f[*n_x];
	double m[*n_x];
	double c[*n_x];
	double delxdela[*n_x];
	double delXdela;
	double delxdelkds[*n_x];
	// double dmdb[*n_x];
	double dmdx[*n_x];
	double dmdxij;
	double dxdA[*n_x];
	double dXdkds[*n_x];
	// double dLdm[*n_x];
	double X, F, M, C, a_j, A_j;
	double dLdbg;
	double dLdA;
	double dLdkds[*n_x];
	dLdbg = 0.0;
	dLdA = 0.0;
	// result[(*n_x) + 2] = 0;	

	for(j = 0.0; j < 5; j++) {
		F = 0.0;
		M = 0.0;
		C = 0.0;
		a_j = a[j];
		A_j = A[j];
		delXdela = 0;
		for(i = 0; i < *n_x; i++) {
			dLdkds[i] = 0.0;
			x[i] = l[i]*a_j/(a_j + kds[i]);
			f[i] = l[i] - x[i];
			F += f[i];
			c[i] = sXc[j*(*n_x) + i];
			C += c[i];
			delxdela[i] = f[i]/(a_j + kds[i]);
			delXdela += delxdela[i];
			delxdelkds[i] = -x[i]/(a_j + kds[i]);
		}
		dmdxij = (F - (*bg))/(float)F;
		for(i = 0; i < *n_x; i++) {
			m[i] = x[i] + (*bg)*f[i]/F;
			M += m[i];
			dmdb[j*(*n_x) + i] = f[i]/(float)F;
			dmdx[i] = (*bg)*f[i]/((float)F*(float)F);
			dxdA[i] = delxdela[i]/(1 + delXdela);
			dXdkds[i] = delxdelkds[i]/(1 + delXdela);
		}
		double sum_dLdm_dmdx;
		double sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela;
		sum_dLdm_dmdx = 0.0;
		sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela = 0.0;
		for(i = 0; i < *n_x; i++) {
			dLdm[j*(*n_x) + i] = C/M - c[i]/m[i];
			gradient[*(n_x)]     += dLdm[j*(*n_x) + i]*dmdb[j*(*n_x) + i];
			gradient[*(n_x) + 1] += dLdm[j*(*n_x) + i]*dmdx[i]*delXdela/(1 + delXdela)*dil[j];
			gradient[*(n_x) + 1] += dLdm[j*(*n_x) + i]*dmdxij*dxdA[i]*dil[j];
			gradient[i]          += dLdm[j*(*n_x) + i]*dmdxij*delxdelkds[i];
			sum_dLdm_dmdx += dLdm[j*(*n_x) + i]*dmdx[i];
			sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela += dLdm[j*(*n_x) + i]*(dmdx[i]*delXdela + dmdxij*delxdela[i]);
		}
		for(i = 0; i < *n_x; i++) {
			gradient[i] += delxdelkds[i]*sum_dLdm_dmdx;
			gradient[i] -= dXdkds[i]*sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela;
		}
	}
}
