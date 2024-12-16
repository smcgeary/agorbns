/* kinetics_simple.c__________________________________________________________*/
#include <R.h>
#include <Rmath.h>

/* Function to compute derivatives of the ODE:................................*/
void FreeAgo(double *ax, double *bx, double *A, double *l, int *n_x, double *kds,
            double *tol, double *result)		/* An estimate to the root	*/
{
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/

  a = *ax;  b = *bx;
  int i;
  fa = *A;
  for(i = 0; i < *n_x; i++) {
    fa -= a*l[i]/(a + kds[i]);
  }
  fa -= a;

  fb = *A;
  for(i = 0; i < *n_x; i++) {
    fb -= b*l[i]/(b + kds[i]);
  }
  fb -= b;
  c = a;   fc = fa;
  /*   0    a/c-----------------------b (A)   */
  /*   A   f(a)/fc)     -[A*l/(A+kds)] f(b)   */
  for(;;) {
    double prev_step = b-a;		/* Distance from the last but one to the last approximation	*/
    double tol_act;			      /* Actual tolerance		*/
    double p, q;      			      /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment	*/
    double new_step;          /* Step at this iteration       */
    if(fabs(fc) < fabs(fb)) { /* Swap data for b to be the best approximation    */
    	a = b;  b = c;  c = a;
      fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*3.000214e-13*fabs(b) + *tol/2;
    new_step = (c-b)/2;
    if(fabs(new_step) <= tol_act || fb == (double)0 ) {
      *result = b;
      return;       /* Acceptable approx. is found  */      
    }
    /* Decide if the interpolation can be tried	*/
    /* If prev_step was large enough and was in true direction, Interpolatiom may be tried */
    if(fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) {	
			register double t1, cb, t2;
    	cb = c-b;
    	if(a == c) {			/* If we have only two distinct	points linear interpolation can only be applied   */
        t1 = fb/fa;
        p = cb*t1;
        q = 1.0 - t1;
     	} else {				/* Quadric inverse interpolation*/
        q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
        p = t2*(cb*q*(q - t1) - (b - a)*(t1 - 1.0));
        q = (q - 1.0)*(t1 - 1.0)*(t2 - 1.0);
      }
    	if(p > (double)0) {		/* p was calculated with the opposite sign; make p positive and assign possible minus to q*/
        q = -q;
    	} else {
        p = -p;
      }
      if(p < (0.75*cb*q - fabs(tol_act*q)/2) && p < fabs(prev_step*q/2)) { /* If b+p/q falls in [b,c] and isn't too large it is accepted */
        new_step = p/q;
  		}/* If p/q is too large then the bissection procedure can reduce [b,c] range to more extent */
    }
    if(fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
      if(new_step > (double)0) {	/* than tolerance		*/
      	new_step = tol_act;
      } else {
      	new_step = -tol_act;
      }
    }
    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;
    fb = *A;
    for(i = 0; i < *n_x; i++) {
      fb -= b*l[i]/(b + kds[i]);
    }
    fb -= b;
    if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {  /* Adjust c for it to have a sign opposite to that of b */
      c = a;  fc = fa;
    }
  }
}



void CostEquil(double *kds, double *l, double *L, double *bg, double *dil,
               double *A, double *a, double *sXc, int *n_x, double *cost)
{
	int i, j;
	double x[*n_x];
	double f[*n_x];
	double m[*n_x];
	double c[*n_x];
	double X, F, M, C, a_j, A_j;
	for(j = 0; j < 5; j++) {
		F = 0;
		M = 0;
		C = 0;
		a_j = a[j];
		A_j = A[j];
		for(i = 0; i < *n_x; i++) {
			x[i] = l[i]*a_j/(a_j + kds[i]);
			f[i] = l[i] - x[i];
			c[i] = sXc[j*(*n_x) + i];
			m[i] = x[i] + (*bg)*f[i]/(*L);
			cost += -c[i]*log(m[i]/(A_j - a_j + (*bg)))
		}
	}
}

void GradientEquil(double *kds, double *l, double *L, double *bg, double *dil, double *A, double *a,
                   double *sXc, int *n_x, double *gradient)
{
	int i, j;
	double x[*n_x];
	double f[*n_x];
	double m[*n_x];
	double c[*n_x];
	double delxdela[*n_x];
	double delXdela;
	double delxdelkds[*n_x];
	double dmdb[*n_x];
	double dLdm[*n_x];
	double dmdx[*n_x];
	double dmdxij;
	double dxdA[*n_x];
	double dXdkds[*n_x];
	double X, F, M, C, a_j, A_j;
	double sum_dLdm_dmdx;
	double sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela;

	for(j = 0; j < 5; j++) {
		sum_dLdm_dmdx = 0.0;
		sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela = 0.0;
		F = 0;
		M = 0;
		C = 0;
		a_j = a[j];
		A_j = A[j];
		delXdela = 0;
		for(i = 0; i < *n_x; i++) {
			x[i] = l[i]*a_j/(a_j + kds[i]);
			f[i] = l[i] - x[i];
			F += f[i];
			m[i] = x[i] + (*bg)*f[i]/(F);
			c[i] = sXc[j*(*n_x) + i];
			C += c[i];
			delxdela[i] = f[i]/(a_j + kds[i]);
			delXdela += delxdela[i];
			delxdelkds[i] = -x[i]/(a_j + kds[i]);
		}
		dmdxij = (F - (*bg))/F;
		for(i = 0; i < *n_x; i++) {
			m[i] = x[i] + (*bg)*f[i]/F;
			M += m[i];
			dmdb[i] = f[i]/F;
			dmdx[i] = (*bg)*f[i]/(F*F);
			dxdA[i] = delxdela[i]/(1 + delXdela);
			dXdkds[i] = delxdelkds[i]/(1 + delXdela);
		}
		for(i = 0; i < *n_x; i++) {
			dLdm[i] = C/M - c[i]/m[i];
			gradient[i]          += log(10)*kds[i]*dLdm[i]*dmdxij*delxdelkds[i];
			gradient[*(n_x)]     += log(10)*(*bg)*dLdm[i]*dmdb[i];
			gradient[*(n_x) + 1] += log(10)*A[j]*dLdm[i]*dmdx[i]*delXdela/(1 + delXdela)*dil[j];
			gradient[*(n_x) + 1] += log(10)*A[j]*dLdm[i]*dmdxij*dxdA[i]*dil[j];
			sum_dLdm_dmdx += dLdm[i]*dmdx[i];
			sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela += dLdm[i]*(dmdx[i]*delXdela + dmdxij*delxdela[i]);
		}
		for(i = 0; i < *n_x; i++) {
			gradient[i] += log(10)*kds[i]*delxdelkds[i]*sum_dLdm_dmdx;
			gradient[i] -= log(10)*kds[i]*dXdkds[i]*sum_dLdm_dmdx_delXdela_and_dmdxij_delxdela;
		}
	}
}







/* END file mymod.c */