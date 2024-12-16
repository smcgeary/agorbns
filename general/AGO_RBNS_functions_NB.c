/* kinetics_simple.c__________________________________________________________*/
#include <R.h>
#include <Rmath.h>
#include <omp.h>


double
FreeAgo(A, l, pars, n_x)

double A;
double *l;
double *pars;
int n_x;

{
	double  xa,  xb,  xc;
	int t = 0;
	xa  = 0.0; xb = A;
	double fa, fb, fc;
	fa = xa - A; fb = xb - A;
	int i;
	for (i = 0; i < n_x; i++) {
		fa += xa*l[i]/(xa + exp(pars[i]));
		fb += xb*l[i]/(xb + exp(pars[i]));
	}
	xc = xa;   fc = fa;

	for(;;) {
		double prev_step = xb - xa;   /* Distance from the last but one to the last approximation */
		double tol_act;           /* Actual tolerance   */
		double p, q;                  /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment */
		double new_step;          /* Step at this iteration       */
		if(fabs(fc) < fabs(fb)) { /* Swap data for b to be the best approximation    */
			xa = xb;  xb = xc;  xc = xa;
			fa = fb;  fb = fc;  fc = fa;
		}
		tol_act  = (2.220446e-16)*(2*fabs(xb) + 0.5);
		new_step = (xc - xb)/2;
		if(fabs(new_step) <= tol_act || fb == (double)0 ) {
			return xb;       /* Acceptable approx. is found  */      
		}
		/* Decide if the interpolation can be tried */
		/* If prev_step was large enough and was in true direction, Interpolatiom may be tried */
		if(fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) { 
			register double t1, cb, t2;
			cb = xc - xb;
			if (xa == xc) {      /* If we have only two distinct points linear interpolation can only be applied   */
				t1 = fb/fa;
				p  = cb*t1;
				q  = 1.0 - t1;
			} else {        /* Quadric inverse interpolation*/
				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
				p = t2*(cb*q*(q - t1) - (xb - xa)*(t1 - 1.0));
				q = (q - 1.0)*(t1 - 1.0)*(t2 - 1.0);
			}
			if(p > (double)0) {   /* p was calculated with the opposite sign; make p positive and assign possible minus to q*/
				q = -q;
			} else {
				p = -p;
			}
			if(p < (0.75*cb*q - fabs(tol_act*q)/2) && p < fabs(prev_step*q/2)) { /* If b+p/q falls in [b,c] and isn't too large it is accepted */
				new_step = p/q;
			}/* If p/q is too large then the bissection procedure can reduce [b,c] range to more extent */
		}
		if(fabs(new_step) < tol_act) {  /* Adjust the step to be not less*/
			if(new_step > (double)0) {  /* than tolerance   */
				new_step = tol_act;
		 	} else {
				new_step = -tol_act;
			}
		}
		xa = xb;  fa = fb; xb += new_step;
		fb = xb - A;
		for(i = 0; i < n_x; i++) {
			fb += xb*l[i]/(xb + exp(pars[i]));
		}
		if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {  /* Adjust c for it to have a sign opposite to that of b */
			xc = xa;  fc = fa;
		}
	}
}

void
CostEquil(pars, y, dil, l, L, Y, n_i, n_j, cost)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
double *cost;

{
	double b = exp(pars[(*n_i)]);
	double kd, A_j, a_j, F, X;
	int i, j;
	for(j = 0; j < *n_j; j++) {
		A_j = exp(pars[(*n_i) + 1])*dil[j];
		a_j = FreeAgo(A_j, l, pars, *n_i);
		F   = (*L) - A_j + a_j;
		X   = A_j - a_j + b;
		for(i = 0; i < *n_i; i++) {
			kd     = exp(pars[i]);
			*cost -= y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*(a_j + b*kd/F));
		}
	}
}


void
GradEquil(pars, y, dil, l, L, Y, n_i, n_j, grad)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
double *grad;

{
	double b = exp(pars[(*n_i)]);
	double kd, A_j, a_j, F, X;
	double Ci, C1, C2, C3, D1, D2, GA_1, GA_2;
	int i, j, y_j;
	for(j = 0; j < *n_j; j++) {
		y_j = (*n_i)*j;
		A_j = exp(pars[(*n_i) + 1])*dil[j];
		a_j   = FreeAgo(A_j, l, pars, *n_i);
		F = (*L) - A_j + a_j;
		X = A_j - a_j + b;
		D1 = F - b;
		D2 = b/F;
		C1 = 0; C3 = 0; GA_1 = 0; GA_2 = 0;
		for(i = 0; i < *n_i; i++) {
			kd    = exp(pars[i]);
			Ci    = 1/(a_j + kd);
			C1   += l[i]*exp(pars[i])*Ci*Ci;
		  	C2    = l[i]*Ci*(Y[j]/X - y[y_j+i]/(l[i]*Ci*(a_j + D2*kd)))/F;
		  	C3   += C2*kd*(D2 - D1*Ci);
		  	GA_1 += C2*kd*Ci;
		  	GA_2 += C2*kd;
		  	grad[i     ] -=     C2*D1*a_j*Ci;
		  	grad[*(n_i)] +=     C2*kd;
		}
		grad[*(n_i) + 1] += A_j/(1 + C1)*(GA_1*D1 + GA_2*D2*C1);
		for(i = 0; i < *n_i; i++) {
			Ci       = 1/(a_j + exp(pars[i]));
		 	grad[i] -= l[i]*a_j*Ci*Ci*C3/(1 + C1);
		}
	}
	for(i = 0; i < *n_i; i++) {
		grad[i] = exp(pars[i])*grad[i];
	}
	grad[*(n_i) - 1] = 0.0;
	grad[*(n_i)    ] = b*grad[*(n_i)];	
}

