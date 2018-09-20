/******************************************************************************

File        : rlogreg.c

Date        : December 2011

Author      : Jakramate Bootkrajang

Description : MEX implementation of a Bayesian logsitic regression method
with label noise modelling. It is largely based on the training algorithm for 
Shevade and Keerthi's sparse logistic regression procedure [1,2] 
and Bayesian regularisation setting procedure by Cawley and Talbot [3].

References  : [1] Shevade, S. K. and Keerthi, S. S., "A simple and effecient 
algorithm for gene selection using sparse logistic
regression", Technical Report CD-02-22, Control Division,
Department of Mechanical Engineering, National University
of Singapore, Singapore - 117 576, 2002.

[2] Shevade, S. K. and Keerthi, S. S., "A simple and effecient 
algorithm for gene selection using sparse logistic
regression", Bioinformatics, vol. 19, no. 17, pp 2246-2253,
2003.

[3] Cawley, G. C. and Talbot, N. L. C., "Gene selection in
cancer classification using sparse logistic regression with
Bayesian regularisation", Bioinformatics (submitted), 2006.

History     : 02/10/2004 - v1.00
30/03/2006 - v1.10 minor improvements to comments etc.

Copyright   : (c) J. Bootkrajang, December 2011

This program is mxFree software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the mxFree Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the mxFree Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 ******************************************************************************/

#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "mex.h"
#include "matrix.h"

/* symboloc constants */

const int I_z    = 0;
const int I_nz   = 1;
const double LIM = 1000;
const double SN  = 0.0000001;

/* global data */

int     *set;         
double  *w;
double  *xi;
double  *exp_xi;
double   lambda;      /* regularisation parameter        */
int      ntp;         /* number of training patterns     */
int      d;           /* number of input features        */
double  *F;
double  *delta;
double  *tmp1;
double  *tmp2;

double **x;
double  *y;
double   tol;
double   E_w;         /* sum of magnitude of weights     */
double   E_d;         /* the data misfit term            */
int      N;           /* number of active features       */

int      epoch;
int      nm;

double   flip[2][2];

typedef struct {
	double llh;
	double aleph;
} result ;

/* macros */

#define SWAP(a,b,type) { auto type a_b_type = a ; a = b ; b = a_b_type; }

/* utility functions */

double square(double x)
{
	return x*x;
}

double doubleMax3(double a, double b, double c)
{
	return a > b ? (a > c ? a : c) : (b > c ? b : c);
}

/******************************************************************************

Procedure   : objective 

Parameters  : violator, aleph 

Returns     : double - value of objective function

Description : Evaluate the objective function for the current set of
coefficients.

 ******************************************************************************/

double compute_objective(int violator, double aleph)
{
	int i;
	double Q, E_w_tmp, sigmoid;

	/* evaluate the data misfit term */

	E_d    = 0.0;

	for (i = 0; i < ntp; i++)
	{
		if (mxIsInf(exp_xi[i]))
			sigmoid = 0.0;
		else 
			sigmoid = 1.0/(1.0 + exp_xi[i]);
		
		if(y[i] == -1) 
				E_d -= log((flip[1][0]*sigmoid) + (flip[0][0]*(1.0-sigmoid)));
		else
				E_d -= log((flip[1][1]*sigmoid) + (flip[0][1]*(1.0-sigmoid)));
	}

	if (violator == -1)
		E_w_tmp = E_w;
	else
	{
		E_w_tmp = 0.0;

		/* evaluate objective functions based on the new value of w[violator] */
		for (i = 1; i < d; i++)
		{
			if (i == violator)
				E_w_tmp += fabs(aleph);
			else
				E_w_tmp += fabs(w[i]);
		}
	}

	Q = E_d + lambda * E_w_tmp;

	return Q;
}

/******************************************************************************

Procedure   : updateFsAndDeltas

Parameters  : int s - the set of input features to update (I_z or I_nz)

Returns     : void

Description : Update the F and delta statistics for all input features
belong to set I_z or I_nz.

 ******************************************************************************/

void updateFsAndDeltas(int s)
{
	int i, j;
	double delt, f, sn0, sn1, a1, a2, sigmoid;

	for (i = 0; i < d; i++)
	{
		if (set[i] == s) 
		{
			f    = 0.0;
			delt = 0.0;

			sn0   = 0.0;
			sn1   = 0.0;

			for (j = 0; j < ntp; j++)
			{
				if (mxIsInf(exp_xi[j]))
					sigmoid = 0.0;
				else 
					sigmoid = 1.0/(1.0 + exp_xi[j]);

				if (y[j] == -1) {
					sn0 = (flip[1][0] * sigmoid) + (flip[0][0] * (1.0 - sigmoid));

					if (sn0 == 0.0)	sn0 += SN; /* avoiding numerical problem */

					f -= ((flip[0][0] - flip[1][0])/sn0) * sigmoid * (1.0 - sigmoid) * x[i][j];

					a1 = 1.0 - sigmoid - sigmoid ;
					a2 = (flip[0][0]-flip[1][0]) * sigmoid * (1.0 - sigmoid) / sn0;

					delt += (flip[0][0]-flip[1][0]) * sigmoid * (1.0 - sigmoid) * square(x[i][j]) * (a1+a2) / sn0;
				}
				else {
					sn1 = (flip[1][1] * sigmoid) + (flip[0][1] * (1.0 - sigmoid));

					if (sn1 == 0.0)	sn1 += SN; /* avoiding numerical problem */

					f -= ((flip[0][1] - flip[1][1])/sn1) * sigmoid * (1.0 - sigmoid) * x[i][j];

					a1 = 1.0 - sigmoid - sigmoid ;
					a2 = (flip[0][1]-flip[1][1]) * sigmoid * (1.0 - sigmoid) / sn1;

					delt += (flip[0][1]-flip[1][1]) * sigmoid * (1.0 - sigmoid) * square(x[i][j]) * (a1+a2) / sn1;
				}
			}
            
            F[i]     = f == 0 ? f + SN : f ;
            delta[i] = delt == 0 ? delt + SN : delt;

		}
	}
}

/******************************************************************************

Procedure   : updateFAndDelta

Parameters  : int feature

Returns     : void

Description : Update the F and delta statistics for a specified input feature.

 ******************************************************************************/

void updateFAndDelta(int feature)
{
	int i;
	double delt = 0.0, f = 0.0, sn0 = 0.0, sn1 = 0.0, a1, a2, sigmoid; 

	for (i = 0; i < ntp; i++)
	{
		if (mxIsInf(exp_xi[i]))
			sigmoid = 0.0;
		else 
			sigmoid = 1.0/(1.0 + exp_xi[i]);

		if(y[i] == -1) 
		{
			sn0 = (flip[1][0] * sigmoid) + (flip[0][0] * (1.0 - sigmoid));

			if (sn0 == 0.0)	sn0 += SN; /* avoiding numerical problem */

			f -= ((flip[0][0] - flip[1][0])/sn0) * sigmoid * (1.0 - sigmoid) * x[feature][i];

			a1 = 1.0 - sigmoid - sigmoid ;
			a2 = (flip[0][0]-flip[1][0]) * sigmoid * (1.0 - sigmoid) / sn0;

			delt += (flip[0][0]-flip[1][0]) * sigmoid * (1.0 - sigmoid) * square(x[feature][i]) * (a1+a2) / sn0;
		}
		else 
		{
			sn1 = (flip[1][1] * sigmoid) + (flip[0][1] * (1.0 - sigmoid));

			if (sn1 == 0.0)	sn1 += SN; /* avoiding numerical problem */

			f -= ((flip[0][1] - flip[1][1])/sn1) * sigmoid * (1.0 - sigmoid) * x[feature][i];

			a1 = 1.0 - sigmoid - sigmoid ;
			a2 = (flip[0][1]-flip[1][1]) * sigmoid * (1.0 - sigmoid) / sn1;

			delt += (flip[0][1]-flip[1][1]) * sigmoid * (1.0 - sigmoid) * square(x[feature][i]) * (a1+a2) / sn1;
			/* sigmoid could be 1 and then delt will be zero which leads to NaN in Newton's method */
		}
	}

	F[feature]     = f == 0 ? f + SN : f ;
	delta[feature] = delt == 0 ? delt + SN : delt;
   
}

/******************************************************************************

Procedure   : updateXi

Parameters  : int    feature - the feature with the most recently updated
coefficient
double change  - size of the change in the coefficient

Returns     : void

Description : Update the xi and exp_xi statistics for all features, following
a given change in the coefficient for specified feature.

 ******************************************************************************/

void updateXi(int feature, double change)
{
	int i;

	for (i = 0; i < ntp; i++)
	{
		xi[i] += change*x[feature][i];   
		exp_xi[i] = exp(xi[i]);
	}
}

/******************************************************************************

Procedure   : update_flip

Parameters  : none

Returns     : void

Description : Update the label flipping probabilities.

 ******************************************************************************/

void update_flip()
{
	int i;

	double s00 = 0.0, s01 = 0.0, s10 = 0.0, s11 = 0.0;   
	double sn0, sn1, f00, f01, f10, f11, sigmoid; 

	for (i = 0; i < ntp; i++)
	{
		if (mxIsInf(exp_xi[i]))
			sigmoid = 0.0;
		else 
			sigmoid = 1.0/(1.0 + exp_xi[i]);

		sn0 = flip[1][0]*sigmoid + flip[0][0]*(1.0 - sigmoid);
		sn1 = flip[1][1]*sigmoid + flip[0][1]*(1.0 - sigmoid);
		if (y[i] == -1) {
			s00 += (1.0 - sigmoid)/sn0 + SN;
			s10 += sigmoid/sn0 + SN;
		} else {
			s01 += (1.0 - sigmoid)/sn1 + SN;
			s11 += sigmoid/sn1 + SN;
		}
	}

	if (s00 == 0)  s00 += SN;
	if (s01 == 0)  s01 += SN;
	
	/* Do not update flip[0][0] straight away, the following line
	   depends on the old value of flip[0][0]                     */
	f00 = (flip[0][0]*s00) / (flip[0][0]*s00 + flip[0][1]*s01);
	f01 = (flip[0][1]*s01) / (flip[0][0]*s00 + flip[0][1]*s01);

	if (s10 == 0)  s10 += SN;
	if (s11 == 0)  s11 += SN;
	
	f10 = (flip[1][0]*s10) / (flip[1][0]*s10 + flip[1][1]*s11);
	f11 = (flip[1][1]*s11) / (flip[1][0]*s10 + flip[1][1]*s11);	

	flip[0][0] = f00;
	flip[0][1] = f01;
	flip[1][0] = f10;
	flip[1][1] = f11;
}

/******************************************************************************

Procedure   : findMaxViolator

Parameters  : int s - the set of input features to search (I_z or I_nz)

Returns     : int - index of the maximally violating coefficient

Description : Return the index of the coefficient that maximally violates
the optimality conditions from the specified set.  Returns -1
if no feature violates the optimality conditions by more than
the tolerance.

 ******************************************************************************/

int findMaxViolator(int s)
{
	int i, violator;

	double maxViol, viol;

	if (s == I_nz)
	{
		updateFsAndDeltas(I_nz);

		if (fabs(F[0]) > tol)
		{
			maxViol  = fabs(F[0]);
			violator = 0;
		}
		else
		{
			maxViol  = tol;
			violator = -1;
		}


		for (i = 1; i < d; i++)
		{
			if (set[i] == I_nz)
			{
				viol = fabs(w[i] > 0.0 ? lambda - F[i] : lambda + F[i]);

				if (viol > maxViol)
				{
					maxViol  = viol;
					violator = i;
				}
			}
		}
	}
	else
	{
		double scale = N > 0 ? (N+1.0)/N : 1.0;

		updateFsAndDeltas(I_z);

		violator = -1;
		maxViol  = tol;

		for (i = 1; i < d; i++)
		{
			if (set[i] == I_z)
			{
				viol = doubleMax3(F[i] - scale*lambda, -F[i] - scale*lambda, 0.0);

				if (viol > maxViol)
				{
					maxViol  = viol;
					violator = i;
				}
			}
		}
	}

	return violator;
}

/******************************************************************************

Procedure   : newton_method

Parameters  : L = lower bound, H = upper bound, slope = slope at current value

Returns     : result structure

Description : optimise violating parameter via Newton's method.

 ******************************************************************************/

result newton_check (int violator, double aleph, double L, double H, double slope)
{
	int i;
	double a, F_tmp, delta_tmp;

	result r = {DBL_MAX, 0.0};

	/* store original values */ 
	memcpy(tmp1,     xi, ntp*sizeof(double));
	memcpy(tmp2, exp_xi, ntp*sizeof(double));
	F_tmp     = F[violator];
	delta_tmp = delta[violator];

	for (i = 0; i < 100; i++)
	{
		if (fabs(slope) < 0.1*tol)	break;

		if (i == 0)
			a = aleph - 0.5*slope/(delta[violator]);
		else
			a = aleph - slope/(delta[violator]);

		if (a < L || a > H)
		{
			if (mxIsNaN(a))
			{
				r.aleph = 0.0;
				r.llh   = DBL_MAX;

				/* restore old copy */
				SWAP(    xi, tmp1, double*) 
				SWAP(exp_xi, tmp2, double*)
				F[violator]     = F_tmp;
				delta[violator] = delta_tmp;
				return r;
			}
			else{
				a = 0.5*(L+H);
			}
		}

		updateXi(violator, aleph - a);

		updateFAndDelta(violator);

		if (a > 0.0)
			slope = lambda - F[violator];
		else
			slope = -lambda - F[violator];

		if (slope > 0.1*tol)
			H = a;
		else if (slope < -0.1*tol)
			L = a;

		aleph = a;
	}

	/* aleph is the new  w_i, it can't be or close to zero*/
	/* as we already tried the case where it is zero */
	r.aleph = aleph;
	
	if (aleph < tol)
		r.llh = DBL_MAX;
	else
		r.llh = compute_objective(violator, aleph);

	/* restore old copy */
	SWAP(    xi, tmp1, double*) 
	SWAP(exp_xi, tmp2, double*)
	F[violator]     = F_tmp;
	delta[violator] = delta_tmp;
	
	return r;
}

/******************************************************************************

Procedure   : optimiseAlpha

Parameters  : int violator

Returns     : none

Description :

 ******************************************************************************/

void optimiseW(int violator)
{
	int i;

	double L = 0.0, H = 0.0, slope, a, F_tmp, delta_tmp, d_right, d_left;

	double aleph = w[violator];

	result r1, r2;

	/* bracket minimum with interval excluding zero (except at a boundary) */

	if (violator == 0)
	{
		/* case 1 - unregularised bias term */
		L     = -LIM;
		H     = +LIM;
		slope = -F[0];
	}
	else
	{
		/* the violator is a regularised parameter */

		if (aleph == 0.0) 
		{
			/* this section should be OK */
			d_right = +lambda - F[violator];
			d_left  = -lambda - F[violator]; 

			r1 = newton_check (violator, 0.0, 0.0, LIM, d_right);
			r2 = newton_check (violator, 0.0, -LIM, 0.0, d_left);

			if (r1.llh < r2.llh) 
			{
				L 	  = 0.0;
				H 	  = LIM;
				slope = d_right;
			}
			else
			{
				L 	  = -LIM;
				H 	  = +0.0;
				slope = d_left;
			}
		}
		else
		{
			
			slope = (aleph > 0.0 ? +lambda : -lambda) - F[violator];			
			
			/* store a copy of xi, exp_xi, F[violator] and delta[violator] */

			memcpy(tmp1,     xi, ntp*sizeof(double));
			memcpy(tmp2, exp_xi, ntp*sizeof(double));

			F_tmp     = F[violator];
			delta_tmp = delta[violator];

			/* update statistics assuming w[violator] = 0.0 */

			updateXi(violator, aleph);
			updateFAndDelta(violator);

			w[violator]   = 0.0;
			set[violator] = I_z;
			N             = N - 1;
			E_w        	  = E_w - fabs(aleph);

			/* compute left and right derivaltives */

			d_right = +lambda - F[violator];
			d_left  = -lambda - F[violator]; 

			if ((d_right>0.0 && d_left<0.0) || d_right==0.0 || d_left==0.0)
			{
				/* parameter can be safely pruned */

				return;
			}
			else {
				/* restore xi and exp_xi etc. */

				SWAP(    xi, tmp1, double*)
				SWAP(exp_xi, tmp2, double*)

				w[violator]     = aleph;
				set[violator]   = I_nz;
				F[violator]     = F_tmp;
				delta[violator] = delta_tmp;
				N               = N + 1;
				E_w 	        = E_w + fabs(aleph);
			

				if (aleph > 0.0) 
				{
					r1 = newton_check (violator, aleph, 0.0, LIM, slope);
					r2 = newton_check (violator, 0.0, -LIM, 0.0, d_left);

					if (r1.llh < r2.llh) 
					{
						L 	  = 0.0;
						H 	  = LIM;
						slope = slope;
					}
					else
					{
						L 	  = -LIM;
						H 	  = +0.0;
						slope = d_left;
					}
				}
				else 
				{
					r1 = newton_check (violator, 0.0,  0.0, LIM, d_right);
					r2 = newton_check (violator, aleph, -LIM, 0.0, slope);

					if (r1.llh < r2.llh) 
					{
						L 	  = 0.0;
						H 	  = LIM;
						slope = d_right;
					}
					else
					{
						L 	  = -LIM;
						H 	  = +0.0;
						slope = slope;
					}
				}
			}
		}
	}

	/* optimise coefficient for violator via Newton's method */

	for (i = 0; i < 100; i++)
	{
		if (fabs(slope) < 0.1*tol)
		{
			break;
		}

		if (i == 0)
		{
			a = aleph - 0.5*slope/(delta[violator]);
		}
		else
		{
			a = aleph - slope/(delta[violator]);
		}

		if (a < L || a > H)
		{
			a = 0.5*(L+H);
		}
	
		updateXi(violator, aleph - a);

		updateFAndDelta(violator);

		if (violator == 0)
		{
			slope = -F[violator];
		}
		else
		{
			if (a > 0.0)
			{
				slope = lambda - F[violator];
			}
			else
			{
				slope = -lambda - F[violator];
			}
		}

		if (slope > 0.1*tol)
		{
			H = a;
		}
		else if (slope < -0.1*tol)
		{
			L = a;
		}

		aleph = a;

	}

	/* update coefficient vector */

	if (violator != 0)
	{
		N   = w[violator] == 0.0 && aleph != 0.0 ? N+1 : N;
		N   = w[violator] != 0.0 && aleph == 0.0 ? N-1 : N;
		E_w = E_w + fabs(aleph) - fabs(w[violator]);
	}

	w[violator]   = aleph;
	set[violator] = aleph == 0.0 ? I_z : I_nz;
}

/******************************************************************************

Procedure   :

Parameters  :

Returns     :

Description :

 ******************************************************************************/

void rlogreg()
{
	int i, j, violator;

	/* initialisation */

	N       = 0;
	E_w     = 0.0;
	

	for (i = 0; i < 10000; i++)
	{
		/* find maximally violating active parameter */		

		violator = findMaxViolator(I_nz); 

		if (violator != -1)
		{
			optimiseW(violator);
		}
		else
		{	
			if (nm) {
				for (j = 0; j < 3; j++){
					update_flip();
				}
			}
			
			/* find maximally violating inactive parameter */
			violator = findMaxViolator(I_z); 			

			if (violator != -1)
			{
				optimiseW(violator);
			}
		}		

		/* update lambda */
		lambda = E_w == 0.0 ? lambda : N/E_w;

		if (violator == -1) break; 
	}
}

/******************************************************************************

Procedure   : mexFunction

Parameters  : int      nlhs   - number of output arguments
mxArray *plhs[] - array of output arguments
int      nrhs   - number of input arguments
mxArray *prhs[] - array of input arguments

Returns     : void

Description : MEX gateway handling transfer of data to and from MATLAB, the
calling sequence is something like

[W FLIP] = RLOGREG(X, Y, FLIP, TOL)

where X is the design matrix, Y is a vector of target values,
W is the weight vector. LAMBDA is the regularisation parameter, 
TOL is value specific to the algorithm, and FLIP is a matrix 
of label flipping probabilities. 

 ******************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j;

	double *X, *Y, *FLIP_RET;

	/* check number of input arguments */

	if (nrhs != 5)
	{
		mexErrMsgTxt("wrong number of input arguments, see 'rlogreg.c' for usage");
	}

	/* get input patterns */

	ntp  = mxGetM(prhs[0]); 
	d    = mxGetN(prhs[0]); 
	X    = (double*)mxGetPr(prhs[0]);

	/* get target patterns */

	if (mxGetM(prhs[1]) != ntp)
	{
		mexErrMsgTxt("Y must have the same number of rows as X");
	}

	if (mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("Y must be a column vector");
	}

	Y = (double*)mxGetPr(prhs[1]);

	/* parse optional parameters */

	nm        = ((double*)mxGetPr(prhs[3]))[0];
    tol       = ((double*)mxGetPr(prhs[4]))[0];
    
	/* allocate matrix for result */

	plhs[0]   = mxCreateDoubleMatrix(d, 1, mxREAL);	
	plhs[1]   = mxCreateDoubleMatrix(2, 2, mxREAL);
	w   	  = mxGetPr(plhs[0]);	
	FLIP_RET  = (double*)mxGetPr(plhs[1]);


	/* initialisation */

	setvbuf(stdout, NULL, _IONBF, 0L);

	set     = (int*)     mxCalloc(d,   sizeof(int));
	xi      = (double*)  mxCalloc(ntp, sizeof(double));	
	exp_xi  = (double*)  mxCalloc(ntp, sizeof(double));
	F       = (double*)  mxCalloc(d,   sizeof(double));
	delta   = (double*)  mxCalloc(d,   sizeof(double));
	tmp1    = (double*)  mxCalloc(ntp, sizeof(double));
	tmp2    = (double*)  mxCalloc(ntp, sizeof(double));	
	x       = (double**) mxCalloc(d,   sizeof(double*));
	y       = (double*)  mxCalloc(ntp, sizeof(double));

	flip[0][0] = ((double*)mxGetPr(prhs[2]))[0];
	flip[0][1] = ((double*)mxGetPr(prhs[2]))[2];
	flip[1][0] = ((double*)mxGetPr(prhs[2]))[1];
	flip[1][1] = ((double*)mxGetPr(prhs[2]))[3];

	for (i=0; i < d; i++)
	{
		x[i]     = &X[i*ntp];
		set[i]   = I_z;
		w[i]     = 0.0;
	}

	set[0] = I_nz;

	for (i = 0; i < ntp; i++)
	{
		xi[i]     = 0.0;		
		exp_xi[i] = 1.0;
		y[i]      = Y[i];
	}

	rlogreg();
    
    /* Release memory */
    mxFree(set);
    mxFree(xi);
    mxFree(exp_xi);
    mxFree(F);
    mxFree(delta);
    mxFree(tmp1);
    mxFree(tmp2);
    mxFree(x);
    mxFree(y);
    

	FLIP_RET[0] = flip[0][0];
	FLIP_RET[1] = flip[1][0];
	FLIP_RET[2] = flip[0][1];
	FLIP_RET[3] = flip[1][1];
}
