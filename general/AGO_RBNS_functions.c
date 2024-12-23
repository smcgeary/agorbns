 /*kinetics_simple.c__________________________________________________________*/
#include <R.h>
#include <Rmath.h>
#include <omp.h>


// double
// FreeAgo(A, l, pars, n_x)

// double A;
// double *l;
// double *pars;
// int n_x;

// {
// 	double  xa,  xb,  xc;
// 	int t = 0;
// 	A = 1;
// 	xa  = 0.0; xb = A;
// 	double fa, fb, fc;
// 	fa = xa - A; fb = xb - A;
// 	int i;
// 	for (i = 0; i < n_x; i++) {
// 		fa += xa*l[i]/(xa + exp(pars[i]));
// 		fb += xb*l[i]/(xb + exp(pars[i]));
// 	}
// 	xc = xa;   fc = fa;

// 	for(;;) {
// 		double prev_step = xb - xa;   /* Distance from the last but one to the last approximation */
// 		double tol_act;           /* Actual tolerance   */
// 		double p, q;                  /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment */
// 		double new_step;          /* Step at this iteration       */
// 		if(fabs(fc) < fabs(fb)) { /* Swap data for b to be the best approximation    */
// 			xa = xb;  xb = xc;  xc = xa;
// 			fa = fb;  fb = fc;  fc = fa;
// 		}
// 		tol_act  = (2.220446e-16)*(2*fabs(xb) + 0.5);
// 		new_step = (xc - xb)/2;
// 		if(fabs(new_step) <= tol_act || fb == (double)0 ) {
// 			return xb;       /* Acceptable approx. is found  */      
// 		}
// 		/* Decide if the interpolation can be tried */
// 		/* If prev_step was large enough and was in true direction, Interpolatiom may be tried */
// 		if(fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb)) { 
// 			register double t1, cb, t2;
// 			cb = xc - xb;
// 			if (xa == xc) {      /* If we have only two distinct points linear interpolation can only be applied   */
// 				t1 = fb/fa;
// 				p  = cb*t1;
// 				q  = 1.0 - t1;
// 			} else {        /* Quadric inverse interpolation*/
// 				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
// 				p = t2*(cb*q*(q - t1) - (xb - xa)*(t1 - 1.0));
// 				q = (q - 1.0)*(t1 - 1.0)*(t2 - 1.0);
// 			}
// 			if(p > (double)0) {   /* p was calculated with the opposite sign; make p positive and assign possible minus to q*/
// 				q = -q;
// 			} else {
// 				p = -p;
// 			}
// 			if(p < (0.75*cb*q - fabs(tol_act*q)/2) && p < fabs(prev_step*q/2)) { /* If b+p/q falls in [b,c] and isn't too large it is accepted */
// 				new_step = p/q;
// 			}/* If p/q is too large then the bissection procedure can reduce [b,c] range to more extent */
// 		}
// 		if(fabs(new_step) < tol_act) {  /* Adjust the step to be not less*/
// 			if(new_step > (double)0) {  /* than tolerance   */
// 				new_step = tol_act;
// 		 	} else {
// 				new_step = -tol_act;
// 			}
// 		}
// 		xa = xb;  fa = fb; xb += new_step;
// 		fb = xb - A;
// 		for(i = 0; i < n_x; i++) {
// 			fb += xb*l[i]/(xb + exp(pars[i]));
// 		}
// 		if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {  /* Adjust c for it to have a sign opposite to that of b */
// 			xc = xa;  fc = fa;
// 		}
// 	}
// }

double
FreeAgo(A, l, pars, n_x)

double A;
double *l;
double *pars;
int n_x;

{
	double n_steps = 0.0;
	double  xa,  xb,  xc;
	int t = 0;
	// A = 1;
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
		n_steps = n_steps + 1.0;
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
FreeAgoR(A, l, kds, n_x, a_f, n_steps, xbs)

double *A;
double *l;
double *kds;
int *n_x;
double *a_f;
int *n_steps;
double *xbs;

{
	*n_steps = 0;
	double  xa,  xb,  xc;
	int t = 0;
	xa  = 0.0; xb = *A;
	double fa, fb, fc;
	fa = xa - *A; fb = xb - *A;
	int i;
	for (i = 0; i < *n_x; i++) {
		fa += xa*l[i]/(xa + kds[i]);
		fb += xb*l[i]/(xb + kds[i]);
	}
	xc = xa;   fc = fa;

	for(;;) {
		*n_steps += 1;
		double prev_step = xb - xa;   /* Distance from the last but one to the last approximation */
		double tol_act;           /* Actual tolerance   */
		double p, q;                  /* Interpolation step is calculated in the form p/q; division operations is delayed until the last moment */
		double new_step;          /* Step at this iteration       */
		if(fabs(fc) < fabs(fb)) { /* Swap data for b to be the best approximation    */
			xa = xb;  xb = xc;  xc = xa;
			fa = fb;  fb = fc;  fc = fa;
		}
		tol_act  = (2.220446e-16)*(2*fabs(xb) + 0.5);
		xbs[*n_steps] = xb;
		new_step = (xc - xb)/2;
		if(fabs(new_step) <= tol_act || fb == (double)0 ) {
			*a_f = xb;
			return;       /* Acceptable approx. is found  */      
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
		fb = xb - *A;
		for(i = 0; i < *n_x; i++) {
			fb += xb*l[i]/(xb + kds[i]);
		}
		if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {  /* Adjust c for it to have a sign opposite to that of b */
			xc = xa;  fc = fa;
		}
	}
}


void
FreeAgoAlt(pars, l, L, n_i, a_f)

double *pars;
double *l;
double *L;
int *n_i;
double *a_f;

{
	double A_j;
	A_j = exp(pars[(*n_i) + 1]);
	*a_f = FreeAgo(A_j, l, pars, *n_i);

}


void
CostEquilOneMirna(pars, y, dil, l, L, Y, n_i, n_j, cost)

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
CostEquilDoubleSite(pars, y_1s, y_2s, dil, l_1s, l_2s, l, L, Y, n_i, n_j, cost,
                    a_j_check, x_1s_check, x_2s_check, F_j_check, b_1s_check,
                    b_2s_check, cost_1s_check, cost_2s_check, m_1s_vec,
                    m_2s_vec, X_j_check)

double *pars;
double *y_1s;
double *y_2s;
double *dil;
double *l_1s;
double *l_2s;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
double *cost;
double *a_j_check;
double *x_1s_check;
double *x_2s_check;
double *F_j_check;
double *b_1s_check;
double *b_2s_check;
double *cost_1s_check;
double *cost_2s_check;
double *m_1s_vec;
double *m_2s_vec;
double *X_j_check;

{
	double b = exp(pars[(*n_i)]);
	double kd_1, kd_2, A_j, a_j, F, X, x_1s, x_2s, b_1s, b_2s, m_norm_1s, m_norm_2s;
	int i, i_1s, i_2s, j;
	int n_2s = *n_i - 1;
	for(j = 0; j < *n_j; j++) {
		double cost_j = 0;
		A_j = exp(pars[(*n_i) + 1])*dil[j];
		a_j = FreeAgo(A_j, l, pars, *n_i);
		a_j_check[j] = a_j;
		/* This simplification doesn't work anymore, due to double sites */
		X = 0;
		/* Music making x and F */
		for(i_1s = 0; i_1s < n_2s; i_1s++) {
			kd_1     = exp(pars[i_1s]);
			x_1s     = l_1s[i_1s]*a_j/(a_j + kd_1);

			x_1s_check[(*n_i)*j + i_1s] = x_1s;
			X += x_1s;
			
			for(i_2s = 0; i_2s < n_2s; i_2s++) {
				kd_2 = exp(pars[i_2s]);
				x_2s = l_2s[n_2s*i_1s + i_2s]*(1 - kd_1*kd_2/(a_j + kd_1)/(a_j + kd_2));
				x_2s_check[n_2s*n_2s*j + n_2s*i_1s + i_2s] = x_2s;
				X += x_2s;
			}
		}
		kd_1     = exp(pars[*n_i - 1]);
		x_1s     = l_1s[*n_i - 1]*a_j/(a_j + kd_1);
		x_1s_check[(*n_i)*(j + 1) - 1] = x_1s;

		X += x_1s;
		F = (*L) - X;
		F_j_check[j] = (*L) - X;
		X_j_check[j] = X;
		/* Loop making background */ 
		for(i_1s = 0; i_1s < n_2s; i_1s++) {
			i = (*n_i)*j + i_1s;
			x_1s = x_1s_check[i];
			b_1s = (l_1s[i_1s] - x_1s)/F*b;
			b_1s_check[i] = b_1s;
			m_norm_1s = (x_1s + b_1s)/(X + b);
			m_1s_vec[i] = m_norm_1s*Y[j];
			*cost -= y_1s[i]*log((x_1s + b_1s)/(X + b));
			cost_1s_check[i] = -y_1s[i]*log((x_1s + b_1s)/(X + b));
			for(i_2s = 0; i_2s < n_2s; i_2s++) {
				i = n_2s*n_2s*j + n_2s*i_1s + i_2s;
				x_2s = x_2s_check[i];
				b_2s = (l_2s[n_2s*i_1s + i_2s] - x_2s)/F*b;
				b_2s_check[i] = b_2s;
				m_norm_2s = (x_2s + b_2s)/(X + b);
				m_2s_vec[i] = m_norm_2s*Y[j];
				// *cost -= y_2s[i]*log((x_2s + b_2s)/(X + b));
				cost_2s_check[i] = y_2s[i]*log((x_2s + b_2s)/(X + b));
			}
		}
		i = (*n_i)*(j + 1) - 1;
		x_1s = x_1s_check[i];
		b_1s = (l_1s[*n_i - 1] - x_1s)/F*b;		
		b_1s_check[i] = b_1s;
		m_norm_1s = (x_1s + b_1s)/(X + b);
		m_1s_vec[i] = m_norm_1s*Y[j];
		*cost -= y_1s[i]*log((x_1s + b_1s)/(X + b));
		cost_1s_check[i] = y_1s[i]*log((x_1s + b_1s)/(X + b));
	}
}




void
GradEquilOneMirna(pars, y, dil, l, L, Y, n_i, n_j, grad)

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


void
CostEquilMultiMirnas(pars, y, dil, l, L, Y, n_i, n_j, n_mir, norm_constant,
                     cost, cost_each, A_j_check, a_j_check, k_i_check)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
double *norm_constant;
double *cost;
double *cost_each;
double *A_j_check;
double *a_j_check;
double *k_i_check;

{
	double b, A;
	double kd, A_j, a_j, F, X;
	int i, j=0, j_mir, k;
	int n_j_mir;
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_i) + k]);
		A = exp(pars[(*n_i) + *n_mir + k]);
		n_j_mir = n_j[k];
		for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
			A_j = A*dil[j];
			a_j = FreeAgo(A_j, l, pars, *n_i);
			A_j_check[j] = A_j;
			a_j_check[j] = a_j;
			F   = (*L) - A_j + a_j;
			X   = A_j - a_j + b;
			for(i = 0; i < *n_i; i++) {
				kd     = exp(pars[i]);
				k_i_check[i] = kd;
				*cost -= y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*
				                                 (a_j + b*kd/F));
				*cost += lgamma(y[(*n_i)*j + i] + 1);
				cost_each[(*n_i)*j + i] = l[i]/(a_j + kd)/X*(a_j + b*kd/F);
				// cost_each[(*n_i)*j + i] = -1.0*y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*
				//                                  (a_j + b*kd/F));
			}
			*cost -= lgamma(Y[j] + 1);
			j += 1;
		}
	}
	*cost = *cost/(*norm_constant);
}

void
CostEquilMultiMirnasPrior(pars, y, dil, l, L, Y, n_i, n_j, n_mir, norm_constant,
                     	  lambda, lambda0, reg_mean, reg_inds, reg_vals, n_reg,
                     	  cost, cost_each, A_j_check, a_j_check, k_i_check,
                     	  residual_check, reg_ind_check)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
double *norm_constant;
double *lambda;
double *lambda0;
double *reg_mean;
int *reg_inds;
double *reg_vals;
int *n_reg;
double *cost;
double *cost_each;
double *A_j_check;
double *a_j_check;
double *k_i_check;
double *residual_check;
int *reg_ind_check;

{
	double b, A;
	double kd, A_j, a_j, F, X;
	int i, j=0, j_mir, k;
	int n_j_mir;
	int reg_ind;
	double reg_resid;
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_i) + k]);
		A = exp(pars[(*n_i) + *n_mir + k]);
		n_j_mir = n_j[k];
		for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
			A_j = A*dil[j];
			a_j = FreeAgo(A_j, l, pars, *n_i);
			A_j_check[j] = A_j;
			a_j_check[j] = a_j;
			F   = (*L) - A_j + a_j;
			X   = A_j - a_j + b;
			for(i = 0; i < *n_i; i++) {
				kd     = exp(pars[i]);
				k_i_check[i] = kd;
				*cost -= y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*
				                                 (a_j + b*kd/F));
				*cost += lgamma(y[(*n_i)*j + i] + 1);
				cost_each[(*n_i)*j + i] = l[i]/(a_j + kd)/X*(a_j + b*kd/F);
				// cost_each[(*n_i)*j + i] = -1.0*y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*
				//                                  (a_j + b*kd/F));
			}
			*cost -= lgamma(Y[j] + 1);
			j += 1;
		}
	}
	// *cost = *cost/(*norm_constant);
	for(i = 0; i < *n_reg; i++) {
		reg_ind = reg_inds[i] - 1;
		reg_ind_check[i] = reg_ind;
		reg_resid = pars[reg_ind] - reg_vals[i] - (*reg_mean);
		reg_resid = (*lambda)*reg_resid*reg_resid + (*lambda0);
		residual_check[i] = reg_resid;
		*cost += reg_resid;
	}
}


void
CostEquilMultiMirnas_8mer(pars, y, dil, l, l_adj, L, Y, n_i, n_j, n_mir, i_6m8,
                          cost, cost_each, A_j_check)

double *pars;
double *y;
double *dil;
double *l;
double *l_adj;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
int *i_6m8;
double *cost;
double *cost_each;
double *A_j_check;

{
	double b, A;
	double kd, A_j, a_j, f, x, F, X;
	double kd_6m8=exp(pars[*i_6m8]);
	int i, j=0, j_mir, k;
	int n_j_mir;
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_i) + k]);
		A = exp(pars[(*n_i) + *n_mir + k]);
		n_j_mir = n_j[k];
		for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
			A_j = A*dil[j];
			A_j_check[j] = A_j;
			a_j = FreeAgo(A_j, l_adj, pars, *n_i);
			F   = *L;
			X   = b;
			for(i = 0; i < *n_i; i++) {
				kd     = exp(pars[i]);
				x      = l[i] * (1 - kd*kd_6m8/((a_j + kd)*(a_j + kd_6m8)));
				X += x;
				F -= x;
			}
			for(i = 0; i < *n_i; i++) {
				kd     = exp(pars[i]);
				x      = l[i] * (1 - kd*kd_6m8/((a_j + kd)*(a_j + kd_6m8)));
				f      = l[i] - x;
				*cost -= y[(*n_i)*j + i]*log((x + b*f/F)/X);
				cost_each[(*n_i)*j + i] = -y[(*n_i)*j + i]*log((x + b*f/F)/X);
			}
			j += 1;
		}
	}
}




void
CostEquilMultiMirnas12mers(pars, y, dil, l_all, L, Y, n_i, n_j, n_mir, n_pos,
                           cost, cost_each)

double *pars;
double *y;
double *dil;
double *l_all;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
int *n_pos;
double *cost;
double *cost_each;
// double *model_check;
// double *A_j_check;
// int *j_check;
// double *dil_check;
// double *model_check;

{
	double b, A;
	double kd, A_j, a_j, F, X;
	int i, j=0, j_mir, k, k_pos;
	int n_j_mir;
	double l[(*n_i)], logkds[(*n_i)];
	/* Loop over mirnas: */
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_pos)*(*n_i) + k]);
		A = exp(pars[(*n_pos)*(*n_i) + *n_mir + k]);
		/* Loop over nucleotide position: */
		for(k_pos = 0; k_pos < *n_pos; k_pos++) {
			/* Define log kds and input: */
			n_j_mir = n_j[k*(*n_pos) + k_pos];
			for(i = 0; i < *n_i; i++) {
				l[i] = l_all[k_pos*(*n_i) + i];
				logkds[i] = pars[k_pos*(*n_i) + i];
			}
			/* loop over the columns */
			for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
				A_j = A*dil[j];
				// j_check[j] = n_j_mir;
				// A_j_check[j] = A_j;
				// dil_check[j] = dil[j];

				// A_j_check[j] = A_j;
				a_j = FreeAgo(A_j, l, logkds, *n_i);
				F   = (*L) - A_j + a_j;
				X   = A_j - a_j + b;
				for(i = 0; i < *n_i; i++) {
					kd     = exp(logkds[i]);
					// model_check[j*(*n_i) + i] = l[i]/(a_j + kd)/X*(a_j + b*kd/F)*Y[j];
					*cost -= y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*(a_j + b*kd/F));
					cost_each[j*(*n_i) + i] = y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*(a_j + b*kd/F));
				}
				j += 1;
			}
		}
	}
}

void
GradEquilMultiMirnas12mers(pars, y, dil, l_all, L, Y, n_i, n_j, n_mir, n_pos,
                           zero_grad, n_z, grad)

double *pars;
double *y;
double *dil;
double *l_all;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
int *n_pos;
int *zero_grad;
int *n_z;
double *grad;

{
	double A, b;
	double kd, A_j, a_j, F, X;
	double Ci, C1, C2, C3, D1, D2, GA_1, GA_2;
	int i, j=0, j_mir, i_mir, i_pos, y_j;
	int n_j_mir;
	double l[(*n_i)], logkds[(*n_i)];
	/* Loop over mirnas: */
	for(i_mir = 0; i_mir < *n_mir; i_mir++) {
		b = exp(pars[(*n_pos)*(*n_i) + i_mir]);
		A = exp(pars[(*n_pos)*(*n_i) + *n_mir + i_mir]);
		for(i_pos = 0; i_pos < *n_pos; i_pos++) {
			/* Define log kds and input: */
			n_j_mir = n_j[i_mir*(*n_pos) + i_pos];
			for(i = 0; i < *n_i; i++) {
				l[i] = l_all[i_pos*(*n_i) + i];
				logkds[i] = pars[i_pos*(*n_i) + i];
			}
			/* loop over the columns */
			for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
				y_j = (*n_i)*j;
				A_j = A*dil[j];
				a_j   = FreeAgo(A_j, l, logkds, *n_i);
				F = (*L) - A_j + a_j;
				X = A_j - a_j + b;
				D1 = F - b;
				D2 = b/F;
				C1 = 0; C3 = 0; GA_1 = 0; GA_2 = 0;
				for(i = 0; i < *n_i; i++) {
					kd    = exp(logkds[i]);
					Ci    = 1/(a_j + kd);
					C1   += l[i]*kd*Ci*Ci;
			  	C2    = l[i]*Ci*(Y[j]/X - y[y_j+i]/(l[i]*Ci*(a_j + D2*kd)))/F;
			  	C3   += C2*kd*(D2 - D1*Ci);
			  	GA_1 += C2*kd*Ci;
			  	GA_2 += C2*kd;
			  	/* Kd all */
			  	grad[i_pos*(*n_i) + i] -=     C2*D1*a_j*Ci;
			  	/* bg */
			  	grad[(*n_pos)*(*n_i) + i_mir] +=     C2*kd;
				}
				grad[(*n_pos)*(*n_i) + (*n_mir) + i_mir] += A_j/(1 + C1)*(GA_1*D1 + GA_2*D2*C1);
				for(i = 0; i < *n_i; i++) {
					Ci       = 1/(a_j + exp(logkds[i]));
					/* Kd self */
				 	grad[i_pos*(*n_i) + i] -= l[i]*a_j*Ci*Ci*C3/(1 + C1);
				}
				j += 1;
			}
			/* bg */
		}
		grad[(*n_pos)*(*n_i) + i_mir] = b*grad[(*n_pos)*(*n_i) + i_mir];
	}
	for(i = 0; i < (*n_pos)*(*n_i); i++) {
		grad[i] = exp(pars[i])*grad[i];
	}
	for(i = 0; i < *n_z; i++) {
		grad[zero_grad[i] - 1] = 0;
	}
}





void
GradEquilMultiMirna(pars, y, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z,
                    	 fixed, norm_constant, grad)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
int *zero_grad;
int *n_z;
int *fixed;
double *norm_constant;
double *grad;

{
	double A, b;
	double kd, A_j, a_j, F, X;
	double Ci, C1, C2, C3, D1, D2, GA_1, GA_2;
	int i, j=0, j_mir, k, y_j;
	int n_j_mir;
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_i) + k]);
		A = exp(pars[(*n_i) + *n_mir + k]);
		n_j_mir = n_j[k];
		for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
			y_j = (*n_i)*j;
			A_j = A*dil[j];
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
			  	grad[*(n_i) + k] +=     C2*kd;
			}
			grad[*(n_i) + *n_mir + k] += A_j/(1 + C1)*(GA_1*D1 + GA_2*D2*C1);
			for(i = 0; i < *n_i; i++) {
				Ci       = 1/(a_j + exp(pars[i]));
			 	grad[i] -= l[i]*a_j*Ci*Ci*C3/(1 + C1);
			}
			j += 1;
		}
		grad[*(n_i) + k] = b*grad[*(n_i) + k];
	}
	for(i = 0; i < *n_i; i++) {
		grad[i] = exp(pars[i])*grad[i];
	}
	for(i = 0; i < *n_z; i++) {
		grad[zero_grad[i] - 1] = 0;
	}
	// grad[*(n_i) - 1] = 0;
	if(*fixed == 1) {
		grad[*(n_i)    ] = 0;
		grad[*(n_i) + 1] = 0;
	}
	for(i = 0; i < *n_i + 2*(*n_mir); i++) {
		grad[i] = grad[i]/(*norm_constant);
	}
}



void
GradEquilMultiMirnaPrior(pars, y, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z,
                    	 fixed, norm_constant, lambda, lambda0, reg_mean,
                    	 reg_inds, reg_vals, n_reg, grad)

double *pars;
double *y;
double *dil;
double *l;
double *L;
double *Y;
int *n_i;
int *n_j;
int *n_mir;
int *zero_grad;
int *n_z;
int *fixed;
double *norm_constant;
double *lambda;
double *lambda0;
double *reg_mean;
int *reg_inds;
double *reg_vals;
int *n_reg;
double *grad;

{
	double A, b;
	double kd, A_j, a_j, F, X;
	double Ci, C1, C2, C3, D1, D2, GA_1, GA_2;
	int i, j=0, j_mir, k, y_j;
	int n_j_mir;
	int reg_ind;
	double reg_resid;
	for(k = 0; k < *n_mir; k++) {
		b = exp(pars[(*n_i) + k]);
		A = exp(pars[(*n_i) + *n_mir + k]);
		n_j_mir = n_j[k];
		for(j_mir = 0; j_mir < n_j_mir; j_mir++) {
			y_j = (*n_i)*j;
			A_j = A*dil[j];
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
			  	grad[*(n_i) + k] +=     C2*kd;
			}
			grad[*(n_i) + *n_mir + k] += A_j/(1 + C1)*(GA_1*D1 + GA_2*D2*C1);
			for(i = 0; i < *n_i; i++) {
				Ci       = 1/(a_j + exp(pars[i]));
			 	grad[i] -= l[i]*a_j*Ci*Ci*C3/(1 + C1);
			}
			j += 1;
		}
		grad[*(n_i) + k] = b*grad[*(n_i) + k];
	}
	for(i = 0; i < *n_i; i++) {
		grad[i] = exp(pars[i])*grad[i];
	}
	for(i = 0; i < *n_z; i++) {
		grad[zero_grad[i] - 1] = 0;
	}
	// grad[*(n_i) - 1] = 0;
	if(*fixed == 1) {
		grad[*(n_i)    ] = 0;
		grad[*(n_i) + 1] = 0;
	}
	// for(i = 0; i < *n_i + 2*(*n_mir); i++) {
	// 	grad[i] = grad[i]/(*norm_constant);
	// }
	// for(i = 0; i < *n_reg; i++) {
	// 	reg_ind = reg_inds[i] - 1;
	// 	reg_resid = pars[reg_ind] - reg_vals[i];
	// 	grad[reg_ind] += 2*(*lambda)*reg_resid;
	// }
	for(i = 0; i < *n_reg; i++) {
		reg_ind = reg_inds[i] - 1;
		reg_resid = pars[reg_ind] - reg_vals[i] - (*reg_mean);
		reg_resid = 2*(*lambda)*reg_resid;
		// *cost += reg_resid;
		grad[reg_ind] += reg_resid;
	}

}



// void CostEquil(double *kds, double *l, double *L, double *bg, double *dil,
// 			   double *A, double *a, double *sXc, int *n_x, int *n_col,
// 			   double *cost)
// {
//   int i, j;
//   double x[*n_x];
//   double f[*n_x];
//   double m[*n_x];
//   double c[*n_x];
//   double X, F, M, C, a_j, A_j;
//   for(j = 0; j < *n_col; j++) {
// 	F = 0;
// 	M = 0;
// 	C = 0;
// 	a_j = a[j];
// 	A_j = A[j];
// 	for(i = 0; i < *n_x; i++) {
// 	  x[i] = l[i]*a_j/(a_j + kds[i]);
// 	  f[i] = l[i] - x[i];
// 	  c[i] = sXc[j*(*n_x) + i];
// 	  m[i] = x[i] + (*bg)*f[i]/(*L - A_j + a_j);
// 	  *cost += -c[i]*log(m[i]/(A_j - a_j + (*bg)));
// 	}
//   }
// }



// void GradientEquil(double *kds, double *l, double *L, double *bg, double *dil,
// 				   double *A, double *a, double *sXc, int *n_x, int *n_col,
// 				   double *dgdp)
// {
//   int i, j;
//   double x[*n_x];
//   double f[*n_x];
//   double m[*n_x];
//   double c[*n_x];
//   double delxdela[*n_x];
//   double delXdela;
//   double delxdelkds[*n_x];
//   double dmdb[*n_x];
//   double dLdm[*n_x];
//   double dmdx[*n_x];
//   double dmdxij;
//   double dxdA[*n_x];
//   double dXdkds[*n_x];
//   double dgdKd_temp[*n_x];
//   double X, F, M, C, a_j, A_j;
//   double temp;
//   for(j = 0; j < *n_col; j++) {
// 	temp = 0.0;
// 	C = 0;
// 	a_j = a[j];
// 	A_j = A[j];
// 	F = (*L) - A_j + a_j;
// 	M = A_j - a_j + (*bg);
// 	dmdxij = (F - (*bg))/F;
// 	delXdela = 0;
// 	for(i = 0; i < *n_x; i++) {
// 	  x[i] = l[i]*a_j/(a_j + kds[i]);
// 	  f[i] = l[i] - x[i];
// 	  m[i] = x[i] + (*bg)*f[i]/(F);
// 	  c[i] = sXc[j*(*n_x) + i];
// 	  C += c[i];
// 	  dmdb[i] = f[i]/F;
// 	  dmdx[i] = (*bg)*f[i]/(F*F);
// 	  delxdela[i] = f[i]/(a_j + kds[i]);
// 	  delXdela += delxdela[i];
// 	  delxdelkds[i] = -x[i]/(a_j + kds[i]);
// 	}
// 	for(i = 0; i < *n_x; i++) {
// 	  dXdkds[i] = delxdelkds[i]/(1 + delXdela);
// 	}
// 	for(i = 0; i < *n_x; i++) {
// 	  dLdm[i] = C/M - c[i]/m[i];
// 	  dgdp[i]          +=  dLdm[i]*dmdxij*delxdelkds[i];
// 	  dgdp[*(n_x)]     +=  dLdm[i]*dmdb[i];
// 	  dgdp[*(n_x) + 1] += A[j]*dLdm[i]/(1 + delXdela)*
// 				  (dmdxij*delxdela[i] + dmdx[i]*delXdela);
// 	  temp       += dLdm[i]*(dmdx[i] - dmdxij*delxdela[i]);
// 	}
// 	for(i = 0; i < *n_x; i++) {
// 	  dgdp[i] += delxdelkds[i]/(1 + delXdela)*temp;
// 	  // dgdp[i] -= dXdkds[i]*temp2;
// 	}
//   }
//   for(i = 0; i < *n_x; i++) {
// 	dgdp[i] = log(10)*kds[i]*dgdp[i];
//   }
//   dgdp[*(n_x) - 1] = 0;
//   dgdp[*(n_x)]     = log(10)*(*bg)*dgdp[*(n_x)];
//   dgdp[*(n_x) + 1] = log(10)*dgdp[*(n_x) + 1];
// }




void
CostEquilMir7All(pars, sXc, L, dils, data_i_start, data_i_end, l_i_p, n_x,
                      cost, cost_each)

double *pars;
double *sXc;
double *L;
double *dils;
int *data_i_start;
int *data_i_end;
int *l_i_p;
int *n_x;
double *cost;
// double *cost2;
// double *l_out;
// double *data_out;
// double *data_out_2;
// double *model_out;
// int *columns_all;
// int *inds_all;
// int *i_all;
// int *j_all;
// int *order_all;
// double *Aj_all;
// double *A_m_all;
// double *bg_m_all;
double *cost_each;
// double *prob_each;
// double *A_j_sum_alt;

{
  int i, j, j_l, j_r, i_m, i_p, i_sXc;
  int i_sample=0, i_column=0;
  int i_added=0;
  double kd;
  double kds[*n_x], l[*n_x], c[*n_x], x[*n_x];
  double YL=0; /* total input counts, for normalizing */
  /* Iterate over the nucleotide positions: */
  double A_m, bg_m, A_j, a_j, F, X;
  for(i_m = 0; i_m < 3; i_m++) {
	A_m = exp(pars[5*(*n_x) + 3 + i_m]);
	bg_m = exp(pars[5*(*n_x) + i_m]);
	for (i_p = 0; i_p < 5; i_p++) {
	  YL = 0;
	  /* Assign the kds and input concentration: */
	  for(i = 0; i < *n_x; i++) {

		kds[i] = pars[i_p*(*n_x) + i];
		l[i] = sXc[(l_i_p[i_p] - 1)*(*n_x) + i] + 1;
		// if(i_m == 0) {    
		// 	l_out[i_p*(*n_x) + i] = l[i];
		// }
		YL += l[i];     
	  }
	  /* Rescale the input l: */
	  for(i = 0; i < *n_x; i++) {
		l[i] = (*L)/YL*l[i];
	  }
	  j_l = data_i_start[i_sample] - 1;
	  j_r = data_i_end[i_sample];
	  for(j = j_l; j < j_r; j++) {
		A_j = A_m*dils[i_column];
		// Aj_all[i_column] = A_j;
		// A_m_all[i_column] = A_m;
		// bg_m_all[i_column] = bg_m;

		a_j = FreeAgo(A_j, l, kds, *n_x);

/* */

		F   = (*L) - A_j + a_j;
		X   = A_j - a_j + bg_m;
		for(i = 0; i < *n_x; i++) {
			kd = exp(kds[i]);
		  // *cost2 -= sXc[(*n_x)*j + i]*log(l[i]/(a_j + kd)/X*(a_j + bg_m*kd/F));
		}
/* */

		// A_j_sum_alt[i_column] = a_j;
		for(i = 0; i < *n_x; i++) {
			kd = exp(kds[i]);
		  c[i] = l[i]*a_j/(a_j + kd);
		  // A_j_sum_alt[i_column] += c[i];
		  x[i] = c[i] + bg_m*(l[i] - c[i])/(*L - A_j + a_j);
		  i_sXc = j*(*n_x) + i;
		  // inds_all[i_column*(*n_x) + i] = i_sXc;
		  *cost -= sXc[j*(*n_x) + i]*log(x[i]/(A_j - a_j + bg_m));
		  cost_each[i_column*(*n_x) + i] = sXc[j*(*n_x) + i]*log(x[i]/(A_j - a_j + bg_m));
		  // prob_each[i_column*(*n_x) + i] = x[i]/(A_j - a_j + bg_m);

		// *cost -= y[(*n_i)*j + i]*log(l[i]/(a_j + kd)/X*
		// 		                                 (a_j + b*kd/F));
		  // data_out[i_column*(*n_x) + i] = sXc[i_sXc];
		  // data_out_2[i_column*(*n_x) + i] = sXc[i_column*(*n_x) + i];
		  // model_out[i_column*(*n_x) + i] = x[i];
		  // i_all[i_column*(*n_x) + i] = i;
		  // j_all[i_column*(*n_x) + i] = j;
		  // order_all[i_column*(*n_x) + i] = i_added;
		  // i_added += 1;
			}
		i_column += 1;
	  }
	  i_sample += 1; 
	}
  }
}



void GradientEquilMir7All(double *pars, double *sXc, double *L, double *dils,
					  int *data_i_start, int *data_i_end, int *l_i_p,
					  int *n_x, double *dgdp)
{
	int i, j, j_l, j_r, i_m, i_p;
	int i_sample=0, i_column=0;
	double kds[*n_x], l[*n_x], c[*n_x], x[*n_x], f[*n_x], y[*n_x];
	double YL; /* total input counts, for normalizing */
	double pdcpda[*n_x], pdCpda, pdcpdkds[*n_x];
	double dmdb[*n_x], dLdx[*n_x], dxdc[*n_x];
	double dxdcij, dxdA[*n_x], dXdkds[*n_x], dgdKd_temp[*n_x];
	double F, X, Y, trace, kd;
	double temp;
	double A_m, bg_m, A_j, a_j;
	for(i_m = 0; i_m < 3; i_m++) {/* START miRNA loop */
		A_m = exp(pars[5*(*n_x) + 3 + i_m]); /* Assign Ago */
		bg_m = exp(pars[5*(*n_x) + i_m]); /* Assign bg */
		for (i_p = 0; i_p < 5; i_p++) {/* START position loop */
			YL=0;
			for(i = 0; i < *n_x; i++) {/* START input count loop */
				kds[i] = pars[i_p*(*n_x) + i];
				l[i] = sXc[(l_i_p[i_p] - 1)*(*n_x) + i] + 1;
				YL += l[i];     
			}/* END input count loop */
			for(i = 0; i < *n_x; i++) {/* START input rescaling loop */
				l[i] = (*L)/YL*l[i];
			}/* END input rescaling loop */
			j_l = data_i_start[i_sample] - 1;
			j_r = data_i_end[i_sample];
			for(j = j_l; j < j_r; j++) {/* START concentration loop */
				A_j = A_m*dils[i_column];
				a_j = FreeAgo(A_j, l, kds, *n_x);
				temp = 0.0;
				Y    = 0.0;
				F    = (*L) - A_j + a_j;
				X    = A_j - a_j + bg_m;
				pdCpda = 0.0;
				for(i = 0; i < *n_x; i++) {/* START site loop 1 */
					kd = exp(kds[i]);
					c[i] = l[i]*a_j/(a_j + kd);
					f[i] = l[i] - c[i];
					x[i] = c[i] + bg_m*f[i]/(F);
					y[i] = sXc[j*(*n_x)+i];
					Y   += y[i];
					pdcpda[i]   = f[i]/(a_j + kd);
					pdCpda     += pdcpda[i];
					pdcpdkds[i] = -c[i]/(a_j + kd);
				}/* END site loop 1 */
				trace = 1.0/(1 + pdCpda);
				for(i = 0; i < *n_x; i++) {/* START site loop 2 */
					dLdx[i]                 = (Y/X - y[i]/x[i])/F;
					dgdp[i_p*(*n_x)  +i  ] += dLdx[i]*(F - bg_m)*pdcpdkds[i];
					dgdp[  5*(*n_x)  +i_m] += dLdx[i]*f[i];
					dgdp[  5*(*n_x)+3+i_m] += A_j*dLdx[i]*trace*((F - bg_m)*pdcpda[i]
																 + bg_m*f[i]/F*pdCpda);
					temp                   += dLdx[i]*(bg_m*f[i]/F -
													   (F - bg_m)*pdcpda[i]);
				}/* END site loop 2 */
				for(i = 0; i < *n_x; i++) {/* START site loop 3 */
					dgdp[i_p*(*n_x)  +i  ] += pdcpdkds[i]*trace*temp;
				}/* END site loop 3 */
				i_column += 1;
			}/* END concentration loop */
			i_sample += 1; 
		}/* END position loop */
	}/* END mirna loop */
	for(i = 0; i < 5*(*n_x); i++) {/* START Kd exponential loop */
		dgdp[i] = exp(pars[i])*dgdp[i];
	}/* END Kd exponential loop */
	// for(i = 1; i < 6; i++) {/* START None_Kd=0 loop */
	// 	dgdp[i*(*n_x)-1] = 0;
	// }/* END None_KD=0 loop */
	for(i = 0; i < 3; i++) {/* START bg exponential loop */
		dgdp[5*(*n_x)+i] = exp(pars[5*(*n_x)+i])*dgdp[5*(*n_x)+i];
	}/* END bg exponential loop */
}




void CostEquilGlobalSingleMirna(double *pars, double *sXc, double *L,
								double *dils, int *data_ls, int *data_rs,
								int *l_is, int *n_x, int *num_exps,
								double *cost, double *kd_check, double *A_check,
								double *b_check, double *data_check,
								double *l_check)
{
	double kds[*n_x], l[*n_x], c[*n_x], x[*n_x];
	double YL=0; /* total input counts, for normalizing */
	/* Iterate over the nucleotide positions: */
	double bg = exp(pars[(*num_exps)*(*n_x)]);
	double A  = exp(pars[(*num_exps)*(*n_x) + 1]);
	int i, j, j_l, j_r;
	int i_column = 0;
	int i_exp;
	double A_j, a_j;
	for(i_exp = 0; i_exp < *num_exps; i_exp++) {
		YL = 0;
		/* Assign the kds and input concentration: */
		for(i = 0; i < *n_x; i++) {
			kds[i] = exp(pars[i_exp*(*n_x) + i]);
			l[i] = sXc[l_is[i_exp]*(*n_x) + i];
			YL += l[i];     
		}
		kd_check[i_exp] = pars[i_exp*(*n_x)];
		/* Rescale the input l: */
		for(i = 0; i < *n_x; i++) {
			l[i] = (*L)/YL*l[i];
		}
		j_l = data_ls[i_exp];
		j_r = data_rs[i_exp];
		for(j = j_l; j < j_r; j++) {
			A_j = A*dils[i_column]/100.0;
			a_j = FreeAgo(A_j, l, kds, *n_x);
			for(i = 0; i < *n_x; i++) {
				c[i] = l[i]*a_j/(a_j + kds[i]);
				x[i] = c[i] + bg*(l[i] - c[i])/(*L - A_j + a_j);
				*cost -= sXc[j*(*n_x) + i]*log(x[i]/(A_j - a_j + bg));
			data_check[i_exp] = sXc[j_l*(*n_x)];
			}
			i_column += 1;
		}
	}
}

void GradEquilGlobalSingleMirna(double *pars, double *sXc, double *L,
								double *dils, int *data_ls, int *data_rs,
								int *l_is, int *n_x, int *num_exps,
								double *dgdp, double *kd_check, double *A_check,
								double *b_check, double *data_check,
								double *l_check)
{
	double kds[*n_x], l[*n_x], c[*n_x], x[*n_x], f[*n_x], y[*n_x];
	double YL=0; /* total input counts, for normalizing */
	/* Iterate over the nucleotide positions: */
	double pdcpda[*n_x], pdCpda, pdcpdkds[*n_x];
	double dmdb[*n_x], dLdx[*n_x], dxdc[*n_x];
	double dxdcij, dxdA[*n_x], dXdkds[*n_x], dgdKd_temp[*n_x];
	double F, X, Y, trace;
	double temp;

	double bg = exp(pars[(*num_exps)*(*n_x)]);
	double A  = exp(pars[(*num_exps)*(*n_x) + 1]);
	int i, j, j_l, j_r;
	int i_column = 0;
	int i_exp;
	double A_j, a_j;
	int ind_bg = (*num_exps)*(*n_x);
	int ind_ago = ind_bg + 1;
	for(i_exp = 0; i_exp < *num_exps; i_exp++) {
		YL=0;
		for(i = 0; i < *n_x; i++) {/* START input count loop */
			kds[i] = exp(pars[i_exp*(*n_x) + i]);
			l[i]   = sXc[l_is[i_exp]*(*n_x) + i];
			YL    += l[i];     
			}/* END input count loop */
		for(i = 0; i < *n_x; i++) {/* START input rescaling loop */
			l[i] = (*L)/YL*l[i];
		}/* END input rescaling loop */
		j_l = data_ls[i_exp];
		j_r = data_rs[i_exp];
		for(j = j_l; j < j_r; j++) {/* START concentration loop */
			A_j = A*dils[i_column]/100.0;
			a_j = FreeAgo(A_j, l, kds, *n_x);
			temp = 0.0;
			Y    = 0.0;
			F    = (*L) - A_j + a_j;
			X    = A_j - a_j + bg;
			pdCpda = 0.0;
			for(i = 0; i < *n_x; i++) {/* START site loop 1 */
				c[i] = l[i]*a_j/(a_j + kds[i]);
				f[i] = l[i] - c[i];
				x[i] = c[i] + bg*f[i]/(F);
				y[i] = sXc[j*(*n_x)+i];
				Y   += y[i];
				pdcpda[i]   = f[i]/(a_j + kds[i]);
				pdCpda     += pdcpda[i];
				pdcpdkds[i] = -c[i]/(a_j + kds[i]);
			}/* END site loop 1 */
			trace = 1.0/(1 + pdCpda);
			for(i = 0; i < *n_x; i++) {/* START site loop 2 */
				dLdx[i]                    = (Y/X - y[i]/x[i])/F;
				dgdp[i_exp*(*n_x)  + i  ] += dLdx[i]*(F - bg)*pdcpdkds[i];
				dgdp[ind_bg] += dLdx[i]*f[i];
				dgdp[ind_ago] += A_j*dLdx[i]*trace*((F - bg)*pdcpda[i]
															 + bg*f[i]/F*pdCpda);
				temp                   += dLdx[i]*(bg*f[i]/F -
												   (F - bg)*pdcpda[i]);
			}/* END site loop 2 */
			for(i = 0; i < *n_x; i++) {/* START site loop 3 */
				dgdp[i_exp*(*n_x) + i] += pdcpdkds[i]*trace*temp;
			}/* END site loop 3 */
			i_column += 1;
		}/* END concentration loop */
	}/* END exp loop */
	for(i = 0; i < (*num_exps)*(*n_x); i++) {/* START Kd exponential loop */
		dgdp[i] = exp(pars[i])*dgdp[i];
	}/* END Kd exponential loop */
	for(i = 0; i < *(num_exps); i++) {/* START None_Kd=0 loop */
		dgdp[(i + 1)*(*n_x) - 1] = 0;
	}/* END None_KD=0 loop */
	dgdp[ind_bg ] = exp(pars[ind_bg ])*dgdp[ind_bg ];
}

void GradEquilGlobal16mersMir7(double *pars, double *sXc, double *L,
								double *dils, int *data_ls, int *data_rs,
								int *l_is, int *n_x_ext, int *n_exp_ext,
								double *dgdp, double *kd_check, double *A_check,
								double *b_check, double *data_check,
								double *l_check)
{
	int n_x = *n_x_ext;
	int n_exp = *n_exp_ext;
	double kds[n_x], l[n_x], c[n_x], x[n_x], f[n_x], y[n_x];
	double YL; /* total input counts, for normalizing */
	/* Iterate over the nucleotide positions: */
	double pdcpda[n_x], pdCpda, pdcpdkds[n_x];
	double dmdb[n_x], dLdx[n_x], dxdc[n_x];
	double dxdcij, dxdA[n_x], dXdkds[n_x], dgdKd_temp[n_x];
	double F, X, Y, trace;
	double temp;
	double A_m, bg_m, A_j, a_j;
	int i, j, j_l, j_r;
	int i_column = 0;
	int i_exp, i_m;
	int bg_ind  = n_exp*n_x;
	int ago_ind = bg_ind + 3;
	for(i_exp = 0; i_exp < n_exp; i_exp++) {
		YL = 0;
		for(i = 0; i < n_x; i++) {/* START input count loop */
			kds[i] = exp(pars[i_exp*n_x + i]);
			l[i]   = sXc[(l_is[i_exp])*n_x + i];
			YL    += l[i];     
			}/* END input count loop */
		for(i = 0; i < n_x; i++) {/* START input rescaling loop */
			l[i] = (*L)/YL*l[i];
		}/* END input rescaling loop */
		for(i_m = 0; i_m < 3; i_m++) {
			A_m  = exp(pars[n_exp*n_x + 3 + i_m]);
			bg_m = exp(pars[n_exp*n_x + i_m]);
			// A_check[3*i_exp + i_m] = A_m;
			j_l = data_ls[3*i_exp + i_m];
			j_r = data_rs[3*i_exp + i_m];
			for(j = j_l; j < j_r; j++) {/* START concentration loop */
				A_j = A_m*dils[i_column]/100.0;
				a_j = FreeAgo(A_j, l, kds, n_x);
				temp = 0.0;
				Y    = 0.0;
				F    = (*L) - A_j + a_j;
				X    = A_j - a_j + bg_m;
				pdCpda = 0.0;
				for(i = 0; i < n_x; i++) {/* START site loop 1 */
					c[i] = l[i]*a_j/(a_j + kds[i]);
					f[i] = l[i] - c[i];
					x[i] = c[i] + bg_m*f[i]/(F);
					y[i] = sXc[j*n_x+i];
					Y   += y[i];
					pdcpda[i]   = f[i]/(a_j + kds[i]);
					pdCpda     += pdcpda[i];
					pdcpdkds[i] = -c[i]/(a_j + kds[i]);
				}/* END site loop 1 */
				trace = 1.0/(1 + pdCpda);
				for(i = 0; i < n_x; i++) {/* START site loop 2 */
					dLdx[i]                = (Y/X - y[i]/x[i])/F;
					dgdp[i_exp*n_x + i  ] += dLdx[i]*(F - bg_m)*pdcpdkds[i];
					dgdp[bg_ind    + i_m] += dLdx[i]*f[i];
					dgdp[ago_ind   + i_m] += A_j*dLdx[i]*trace*(
					  (F - bg_m)*pdcpda[i] + bg_m*f[i]/F*pdCpda
					);
					temp                  += dLdx[i]*(
					  bg_m*f[i]/F - (F - bg_m)*pdcpda[i]
					);
				}/* END site loop 2 */
				for(i = 0; i < n_x; i++) {/* START site loop 3 */
					dgdp[i_exp*n_x + i] += pdcpdkds[i]*trace*temp;
				}/* END site loop 3 */
				i_column += 1;
			}/* END dilution sample loop */
		}/* END mirna loop*/
	}/* END exp loop */
	for(i = 0; i < n_exp*n_x; i++) {/* START Kd exponential loop */
		dgdp[i] = exp(pars[i])*dgdp[i];
	}/* END Kd exponential loop */
	for(i = 0; i < n_exp; i++) {/* START None_Kd=0 loop */
		dgdp[(i + 1)*n_x - 1] = 0;
	}/* END None_KD=0 loop */
	for(i = 0; i < 3; i++) {
		dgdp[n_exp*n_x + i] = exp(pars[n_exp*n_x + i])*dgdp[n_exp*n_x + i];
	}
}



void CostEquilGlobal16mersMir7(double *pars, double *sXc, double *L,
								double *dils, int *data_ls, int *data_rs,
								int *l_is, int *n_x_ext, int *n_exp_ext,
								double *cost, double *kd_check, double *A_check,
								double *b_check, double *data_check,
								double *l_check)
{
	int n_x   = *n_x_ext;
	int n_exp = *n_exp_ext;
	double kds[n_x], l[n_x], c[n_x], x[n_x];
	double YL=0; /* total input counts, for normalizing */
	/* Iterate over the nucleotide positions: */
	int i, j, j_l, j_r;
	int i_column = 0;
	int i_exp, i_m;
	double A_m, bg_m, A_j, a_j;
	for(i_exp = 0; i_exp < n_exp; i_exp++) {
		YL = 0;
		/* Assign the kds and input concentration: */
		for(i = 0; i < n_x; i++) {
			kds[i] = exp(pars[i_exp*n_x + i]);
			l[i]   = sXc[l_is[i_exp]*n_x + i];
			YL    += l[i];     
		}
		kd_check[i_exp] = kds[0];
		/* Rescale the input l: */
		for(i = 0; i < n_x; i++) {
			l[i] = (*L)/YL*l[i];
		}
		for(i_m = 0; i_m < 3; i_m++) {
			cost[3*i_exp + i_m] = 0;
			A_m = exp(pars[n_exp*n_x + 3 + i_m]);
			A_check[3*i_exp + i_m] = A_m;
			bg_m = exp(pars[n_exp*n_x + i_m]);
			l_check[3*i_exp + i_m] = l[0];
			b_check[3*i_exp + i_m] = bg_m;
			j_l = data_ls[3*i_exp + i_m];
			j_r = data_rs[3*i_exp + i_m];
			for(j = j_l; j < j_r; j++) {
				A_j = A_m*dils[i_column]/100.0;
				a_j = FreeAgo(A_j, l, kds, n_x);
				for(i = 0; i < n_x; i++) {
					c[i] = l[i]*a_j/(a_j + kds[i]);
					x[i] = c[i] + bg_m*(l[i] - c[i])/(*L - A_j + a_j);
					cost[3*i_exp + i_m] -= sXc[j*n_x + i]*log(x[i]/(A_j - a_j + bg_m));
				data_check[3*i_exp + i_m] = sXc[j_l*(n_x)];
				}
				i_column += 1;
			}
		}
	}
}

void GradEquilGlobalSingleMirnaParallel(double *pars, double *sXc, double *L,
								double *dils, int *data_ls, int *data_rs,
								int *l_is, int *n_x, int *num_exps,
								double *dgdp, double *kd_check, double *A_check,
								double *b_check, double *data_check,
								double *l_check)
{
	double kds[*n_x], l[*n_x], c[*n_x], x[*n_x], f[*n_x], y[*n_x];
	double YL=0; /* total input counts, for normalizing */
	/* Iterate over the nucleotide positions: */
	double pdcpda[*n_x], pdCpda, pdcpdkds[*n_x];
	double dmdb[*n_x], dLdx[*n_x], dxdc[*n_x];
	double dxdcij, dxdA[*n_x], dXdkds[*n_x], dgdKd_temp[*n_x];
	double F, X, Y, trace;
	double temp;

	double bg = exp(pars[(*num_exps)*(*n_x)]);
	double A  = exp(pars[(*num_exps)*(*n_x) + 1]);
	int i, j, j_l, j_r;
	int i_column = 0;
	int i_exp;
	double A_j, a_j;
	int ind_bg = (*num_exps)*(*n_x);
	int ind_ago = ind_bg + 1;
	for(i_exp = 0; i_exp < *num_exps; i_exp++) {
		YL=0;
		#pragma omp parallel for reduction(+:YL)
		for(i = 0; i < *n_x; i++) {/* START input count loop */
			kds[i] = exp(pars[i_exp*(*n_x) + i]);
			l[i]   = sXc[l_is[i_exp]*(*n_x) + i];
			YL    += l[i];     
		}/* END input count loop */
		#pragma omp parallel for reduction(+:YL)
		for(i = 0; i < *n_x; i++) {/* START input rescaling loop */
			l[i] = (*L)/YL*l[i];
		}/* END input rescaling loop */
		j_l = data_ls[i_exp];
		j_r = data_rs[i_exp];
		for(j = j_l; j < j_r; j++) {/* START concentration loop */
			A_j = A*dils[i_column]/100.0;
			a_j = FreeAgo(A_j, l, kds, *n_x);
			temp = 0.0;
			Y    = 0.0;
			F    = (*L) - A_j + a_j;
			X    = A_j - a_j + bg;
			pdCpda = 0.0;
			#pragma omp parallel for reduction(+:Y, pdCpda)
			for(i = 0; i < *n_x; i++) {/* START site loop 1 */
				c[i] = l[i]*a_j/(a_j + kds[i]);
				f[i] = l[i] - c[i];
				x[i] = c[i] + bg*f[i]/(F);
				y[i] = sXc[j*(*n_x)+i];
				Y   += y[i];
				pdcpda[i]   = f[i]/(a_j + kds[i]);
				pdCpda     += pdcpda[i];
				pdcpdkds[i] = -c[i]/(a_j + kds[i]);
			}/* END site loop 1 */
			trace = 1.0/(1 + pdCpda);
			#pragma omp parallel for reduction(+:dgdp, temp)
			for(i = 0; i < *n_x; i++) {/* START site loop 2 */
				dLdx[i]                    = (Y/X - y[i]/x[i])/F;
				dgdp[i_exp*(*n_x)  + i  ] += dLdx[i]*(F - bg)*pdcpdkds[i];
				dgdp[ind_bg] += dLdx[i]*f[i];
				dgdp[ind_ago] += A_j*dLdx[i]*trace*((F - bg)*pdcpda[i]
															 + bg*f[i]/F*pdCpda);
				temp                   += dLdx[i]*(bg*f[i]/F -
												   (F - bg)*pdcpda[i]);
			}/* END site loop 2 */
			#pragma omp parallel for reduction(+:dgdp)
			for(i = 0; i < *n_x; i++) {/* START site loop 3 */
				dgdp[i_exp*(*n_x) + i] += pdcpdkds[i]*trace*temp;
			}/* END site loop 3 */
			i_column += 1;
		}/* END concentration loop */
	}/* END exp loop */
	for(i = 0; i < (*num_exps)*(*n_x); i++) {/* START Kd exponential loop */
		dgdp[i] = exp(pars[i])*dgdp[i];
	}/* END Kd exponential loop */
	for(i = 0; i < *(num_exps); i++) {/* START None_Kd=0 loop */
		dgdp[(i + 1)*(*n_x) - 1] = 0;
	}/* END None_KD=0 loop */
	dgdp[ind_bg ] = exp(pars[ind_bg ])*dgdp[ind_bg ];
}



void GradientEquilMir7AllTemp(double *pars, double *sXc, double *L, double *dils,
					  int *data_i_start, int *data_i_end, int *l_i_p,
					  int *n_x, double *dgdp, double *kd_check, double *A_check, 
					  double *b_check, double *dils_check, double *data_check,
					  double *l_check)
{
	int i, j, j_l, j_r, i_m, i_p;
	int i_sample=0, i_column=0;
	double kds[*n_x], l[*n_x], c[*n_x], x[*n_x], f[*n_x], y[*n_x];
	double YL; /* total input counts, for normalizing */
	double pdcpda[*n_x], pdCpda, pdcpdkds[*n_x];
	double dmdb[*n_x], dLdx[*n_x], dxdc[*n_x];
	double dxdcij, dxdA[*n_x], dXdkds[*n_x], dgdKd_temp[*n_x];
	double F, X, Y, trace;
	double temp;
	double A_m, bg_m, A_j, a_j;
	for(i_m = 0; i_m < 3; i_m++) {/* START miRNA loop */
		A_m = exp(pars[4*(*n_x) + 3 + i_m]); /* Assign Ago */
		bg_m = exp(pars[4*(*n_x) + i_m]); /* Assign bg */
		A_check[i_column] = A_m;
		b_check[i_column] = bg_m;
		for (i_p = 0; i_p < 4; i_p++) {/* START position loop */
			YL=0;
			for(i = 0; i < *n_x; i++) {/* START input count loop */
				kds[i] = exp(pars[i_p*(*n_x) + i]);
				kd_check[i_column] = kds[0];
				l[i] = sXc[(l_i_p[i_p] - 1)*(*n_x) + i];
				YL += l[i];     
			}/* END input count loop */
			for(i = 0; i < *n_x; i++) {/* START input rescaling loop */
				l[i] = (*L)/YL*l[i];
			}/* END input rescaling loop */
			j_l = data_i_start[i_sample] - 1;
			j_r = data_i_end[i_sample];
			for(j = j_l; j < j_r; j++) {/* START concentration loop */
				A_j = A_m*dils[i_column]/100.0;
				a_j = FreeAgo(A_j, l, kds, *n_x);
				dils_check[i_column] = dils[i_column];
				l_check[i_column] = l[0];
				data_check[i_column] = sXc[j*(*n_x)];
				temp = 0.0;
				Y    = 0.0;
				F    = (*L) - A_j + a_j;
				X    = A_j - a_j + bg_m;
				pdCpda = 0.0;
				for(i = 0; i < *n_x; i++) {/* START site loop 1 */
					c[i] = l[i]*a_j/(a_j + kds[i]);
					f[i] = l[i] - c[i];
					x[i] = c[i] + bg_m*f[i]/(F);
					y[i] = sXc[j*(*n_x)+i];
					Y   += y[i];
					pdcpda[i]   = f[i]/(a_j + kds[i]);
					pdCpda     += pdcpda[i];
					pdcpdkds[i] = -c[i]/(a_j + kds[i]);
				}/* END site loop 1 */
				trace = 1.0/(1 + pdCpda);
				for(i = 0; i < *n_x; i++) {/* START site loop 2 */
					dLdx[i]                 = (Y/X - y[i]/x[i])/F;
					dgdp[i_p*(*n_x)  +i  ] += dLdx[i]*(F - bg_m)*pdcpdkds[i];
					dgdp[  4*(*n_x)  +i_m] += dLdx[i]*f[i];
					dgdp[  4*(*n_x)+3+i_m] += A_j*dLdx[i]*trace*((F - bg_m)*pdcpda[i]
																 + bg_m*f[i]/F*pdCpda);
					temp                   += dLdx[i]*(bg_m*f[i]/F -
													   (F - bg_m)*pdcpda[i]);
				}/* END site loop 2 */
				for(i = 0; i < *n_x; i++) {/* START site loop 3 */
					dgdp[i_p*(*n_x)  +i  ] += pdcpdkds[i]*trace*temp;
				}/* END site loop 3 */
				i_column += 1;
			}/* END concentration loop */
			i_sample += 1; 
		}/* END position loop */
	}/* END mirna loop */
	for(i = 0; i < 4*(*n_x); i++) {/* START Kd exponential loop */
		dgdp[i] = exp(pars[i])*dgdp[i];
	}/* END Kd exponential loop */
	for(i = 1; i < 5; i++) {/* START None_Kd=0 loop */
		dgdp[i*(*n_x)-1] = 0;
	}/* END None_KD=0 loop */
	for(i = 0; i < 3; i++) {/* START bg exponential loop */
		dgdp[4*(*n_x)+i] = exp(pars[4*(*n_x)+i])*dgdp[4*(*n_x)+i];
	}/* END bg exponential loop */
}




void GradientEquilNew(double *pars, double *sXc, double *L, double *dil,
					  int *input_i, int *data_i, int *n_x, int *n_col,
					  double *dgdp)
{
  int i, j;
  double kds[*n_x], c[*n_x], f[*n_x], l[*n_x], x[*n_x], y[*n_x];
  double pdcpda[*n_x], pdCpda, pdcpdkds[*n_x];
  double dmdb[*n_x], dLdx[*n_x], dxdc[*n_x];
  double dxdcij, dxdA[*n_x], dXdkds[*n_x], dgdKd_temp[*n_x];
  double C, F, X, Y, a_j, A_j, trace;
  double temp;
  double bg = exp(pars[(*n_x)]);
  double Y_i; /* total input counts, for normalizing */

  for(i = 0; i < *n_x; i++) {
	kds[i] = exp(pars[i]);
	l[i]   = sXc[(*input_i-1)*(*n_x) + i];
	Y_i   += l[i];
  }
  for(i = 0; i < *n_x; i++) {
	kds[i] = exp(pars[i]);
	l[i] = l[i]*(*L)/Y_i;
  }  


  for(j = 0; j < *n_col; j++) {
	A_j = exp(pars[(*n_x) + 1])*dil[j]/100.0;
	a_j   = FreeAgo(A_j, l, kds, *n_x);
	temp = 0.0;
	Y = 0;
	F = (*L) - A_j + a_j;
	X = A_j - a_j + bg;
	pdCpda = 0;
	for(i = 0; i < *n_x; i++) {
	  c[i] = l[i]*a_j/(a_j + kds[i]);
	  f[i] = l[i] - c[i];
	  x[i] = c[i] + bg*f[i]/(F);
	  y[i] = sXc[(data_i[j]-1)*(*n_x) + i];
	  Y += y[i];
	  pdcpda[i] = f[i]/(a_j + kds[i]);
	  pdCpda += pdcpda[i];
	  pdcpdkds[i] = -c[i]/(a_j + kds[i]);
	}
	trace = 1.0/(1 + pdCpda);
	for(i = 0; i < *n_x; i++) {
	  dLdx[i]           =  (Y/X - y[i]/x[i])/F;
	  dgdp[i]          += dLdx[i]*(F - bg)*pdcpdkds[i];
	  dgdp[*(n_x)]     += dLdx[i]*f[i];
	  dgdp[*(n_x) + 1] += A_j*dLdx[i]*trace*((F - bg)*pdcpda[i] + bg*f[i]/F*pdCpda);
	  temp             += dLdx[i]*(bg*f[i]/F - (F - bg)*pdcpda[i]);
	}
	for(i = 0; i < *n_x; i++) {
	  dgdp[i] += pdcpdkds[i]*trace*temp;
	}
  }
  for(i = 0; i < *n_x; i++) {
	dgdp[i] = kds[i]*dgdp[i];
  }
  dgdp[*(n_x) - 1] = 0;
  dgdp[*(n_x)]     = bg*dgdp[*(n_x)];
}



void ModelEquil(double *kds, double *l, double *L, double *bg, double *dil,
			   double *A, double *a, double *C, int *n_x, int *n_col,
			   double *model)
{
  int i, j;
  double x[*n_x];
  double f[*n_x];
  double m[*n_x];

  double X, F, M, C_j, a_j, A_j;
  for(j = 0; j < *n_col; j++) {
	F = 0;
	M = 0;
	C_j = C[j];
	a_j = a[j];
	A_j = A[j];
	for(i = 0; i < *n_x; i++) {
	  x[i] = l[i]*a_j/(a_j + kds[i]);
	  f[i] = l[i] - x[i];
	  m[i] = x[i] + (*bg)*f[i]/(*L - A_j + a_j);
	  model[j*(*n_x) + i] = m[i]/(A_j - a_j + (*bg))*C_j;
	}
  }
}









/* END file mymod.c */ 