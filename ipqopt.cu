/********** qopt.c
 * a translation of my ipqopt matlab code to CUDA C
 * note: a work-in-progress
 ***********/

#include<stdio.h>
#include<stdlib.h>
#include<cuda_runtime.h>
#include "cublas_v2.h"

int main(void) {
}

struct F{
	const float *c;
	const float *l;
	const float *q;
	const int *n;
}

struct G{
	const float *c;
	const float *l;
	const float *q;
	const int *m;
}

struct H{
	const float *A;
	const float *b;
	const int *p;
}

// using cublas
// pointers must point to device memory
void ipqopt(
	float *x,		// primal optimal vector / host
	float *L,		// dual oprtimal vector / host
	int *I,			// number of iterations / host
	float *fo,		// optimal value of f / host
	F *f,			// struct housing quadratic function f / host
	G *g,			// struct housing quadratic functions g / host
	H *h,			// struct housing linear equality constraints / host
	float *dgap,		// duality gap / host
	float *mu,		// duality gap reduction factor / host
	float *eps,		// duality gap tolerance / host
	float *epsf,		// residual norm tolerance / host
	float *alpha,	// line search step threshold (usually >= 0.01 & <= 0.1) / host
	float *beta) {	// line search step reduction factor (usually >= 0.3 & <= 8) / host

	int *n = *f.n;
	int *m = *g.m;
	int *p = *h.p;

	// declare and malloc all the things
	float *rdua;
	cudaMalloc((void**)&rdua, n*sizeof(float));
	float *rcen;
	cudaMalloc((void**)&rdua, m*sizeof(float));
	float *rpri;
	cudaMalloc((void**)&rpri, p*sizeof(float));
	float *rt;
	cudaMalloc((void**)&rt,(n+m+p)*sizeof(float));
	float *B;
	cudaMalloc((void**)&B, (n+m+p)*(n+m+p)*sizeof(float));
	
	float *t; *t = 0;
	float *invt; *invt = 0;

	float *L;
	cudaMalloc((void**)&L, m*sizeof(float));
	float *v;
	cudaMalloc((void**)&v, p*sizeof(float));
	float *y;
	cudaMalloc((void**)&y, (n+m+p)*sizeof(float));

	float *dx;
	cudaMalloc((void**)&dx, n*sizeof(float));
	float *dL;
	cudaMalloc((void**)&dL, m*sizeof(float));
	float *dv;
	cudaMalloc((void**)&dv, p*sizeof(float));
	float *dy;
	cudaMalloc((void**)&dy, (n+m+p)*sizeof(float));

	float *dxp;
	cudaMalloc((void**)&dxp, n*sizeof(float));
	float *dLp;
	cudaMalloc((void**)&dLp, m*sizeof(float));
	float *dvp;
	cudaMalloc((void**)&dvp, p*sizeof(float));

	float *ge;
	cudaMalloc((void**)&ge, m*sizeof(float));
	float *Dg;
	cudaMalloc((void**)&Dg, (m*m)*sizeof(float));
	float *lgq;
	cudaMalloc((void**)&lgq, (n*n)*sizeof(float));
	float *gqx;
	cudaMalloc((void**)&gqx, m*sizeof(float));
	float *onesv;
	cudaMalloc((void**)&onesv, m*sizeof(float));
	// TODO: instantiate ones vector
	
	float *smax;
	// cudaMalloc((void**)&smax, m*sizeof(float));
	// const float *s; cudaMemCpy(
	// float *oldnorm; *oldnorm = 0;

	int *I; *I = 0;

	// intermediary variables necessitated by CUDA and CUBLAS


	/**** Initialize L *****/
	eval_g(handle, g, x, ge, n, m);
	eval_Dg(g, Dg);

	// TODO: instantiate onesdg vector
	// find inverse of diagonalized ge
	cublasStbsv(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, *n, 1, g.l, *n, ones, 1);

	// constants necessitated by CUDA
	const int *const_one;
	cudaMalloc((void**)&const_one, sizeof(int));
	cudaMemset((void*)const_one, 1, sizeof(int));
	const int *const_neg_one;
	cudaMalloc((void**)&const_neg_one, sizeof(int));
	cudaMemset((void*)const_neg_one, -1, sizeof(int));
	const int *const_zero;
	cudaMalloc((void**)&const_zero, sizeof(int));
	cudaMemset((void*)const_zero, 0, sizeof(int));

	const float *invtones;
	cudaMalloc((void**)&invtones, *m*sizeof(int));

	do {
	/**** I. Determine t ****/
	// MATLAB steps
	//	dgap = -ge'*L;
	//	t = mu*m/dgap;
	//	invt = 1/t;

	// C steps
	// 1. Use cublasSdot: x=ge, y=L
		cublasSdot(handle, *m, ge, 1, L, 1, dgap);
	// 3. write simple kernel that inverts dgap, calculates & inverts t
		*dgap = -*dgap;
		*t = *mu * *g.*m / *dgap;
	// 4. Write kernel that generates invt*ones(m,1) vector (invtones)
		cudaMemset((void*)invtones, invt, *g.*m);


	/**** II. Compute dy ****/
	// MATLAB steps
	//	% calculate rt and build rt
	//	rdua = f.l+Dg'*L;
	//	if (iseq) 
	//		rdua = rdua+h.A'*v;
	//		rpri = h.A*x-h.b;
	//	end
	//	if (isfq) 
	//		rdua = rdua+f.q*x; 
	//	end
	//	rcen = -diag(L)*ge-invt*ones(m,1);
	//	rt = [rdua;rcen;rpri];

	// C steps
	// 1. Set rdua to f.l and call cublasSgemv: y=rdua, A=Dg, x=L, a=1, b=1
		cublasScopy(handle, *n, f.l, 1, rdua, 1);
		cublasSgemv(handle, CUBLAS_OP_T, m, n, const_one, Dg, m, L, 1, const_one, rdua, 1);
	// 2. Call cublasSgemv: y=rdua, A=h.A, x=v, a=1, b=1
		cublasSgemv(handle, CUBLAS_OP_T, p, n, const_one, *h.A, p, v, 1, const_one, rdua, 1);
	// 3. Set rpri to h.b and call cublasSgemv: y=rpri, A=h.A, x=x, a=1, b=-1
		cublasScopy(handle, *p, *h.b, 1, rpri, 1);
		cublasSgemv(handle, CUBLAS_OP_N, p, n, const_one, *h.A, p, x, 1, const_neg_one, rpri, 1);
	// 4. Set rcen to invtones and call cublasSsbmv: y=rcen, A=L, x=ge, a=-1, b=-1, k=1
		cublasScopy(handle, *m, invtones, 1, rcen, 1);
		cublasSsbmv(handle, CUBLAS_FILL_MODE_LOWER, *m, 1, const_neg_one, L, *m, ge, 1, const_neg_one, rcen, 1);
	// 5. Use cublasScopy with clever arithmetic to copy rdua, rcen, and rpri to rt
		cublasScopy(handle, *n, rdua, 1, rt, 1);
		cublasScopy(handle, *m, rcen, 1, rt+*n, 1);
		cublasScopy(handle, *m, rpri, 1, rt+*n+*m, 1);

	// MATLAB steps
	//	% build block matrix and rt
	//	if (isfq) B(1:n,1:n) = (f.q+lgq); end
	//	if (isgq) B(1:n,1:n) = B(1:n,1:n)+lgq; end
	//	B(1:n,n+1:n+m) = Dg';
		
	//	B(n+1:n+m,1:n) = -diag(L)*Dg;
	//	B(n+1:n+m,n+1:n+m) = -diag(ge);
		
	//	if (iseq) 
	//		B(1:n,n+m+1:end) = h.A';
	//		B(n+m+1:end,1:n) = h.A;
	//	end

	// C steps
	// 1. Write kernel that builds the B matrix
	

	// MATLAB steps
	//	% solve system
	//	dy = -(B\rt);
	//	dx = dy(1:n);
	//	dL = dy(n+1:n+m);
	//	dv = dy(n+m+1:n+m+p);

	// C steps
	// 1. Solve system by using LU decomposition, by using someone else's code
	// 2. Use cublas<t>copy() to extract dx, dL, dv

	/**** III. Line search ****/
	// MATLAB steps
	//	% calculate smax
	//	for i=1:m
	//		if dL(i) < 0 && -dL(i) > L(i) % avoid dividing
	//			smax(i) = -L(i)/dL(i);
	//		else
	//			smax(i) = 1;
	//		end
	//	end
	//	s = 0.99*min(smax);
	
	// C steps
	// 1. Write C kernel that carries this out

	// MATLAB steps
	//	% decrease s until g(xp) < 0
	//	while (1)
	//		xp = x+s*dx;
	//		[ge, gqx] = eval_g(g, xp, isgq, ge, gqx);
	//		if (all(ge < zeros(m,1))) break;
	//		else s = beta*s;
	//		end
	//	end

	// C steps
	// 1. Copy x to xp: cublasScopy()
	// 2. Use eval_g
	// 3. Write kernel that returns a boolean and move boolean to host
	// 4. If true, multiply s by beta and move on


	// MATLAB steps
	//	% decrease s until suitable step size is found
	//	oldnorm = norm(rt);
	//	while (1)
	//		xp = x+s*dx;
	//		Lp = L+s*dL;
	//		vp = v+s*dv;

	//		[ge, gqx] = eval_g(g, xp, isgq, ge, gqx);
	//		Dg = eval_Dg(g, xp, gqx, isgq, Dg);

	//		rdua = f.l+Dg'*L;
	//		if (iseq) 
	//			rdua = rdua+h.A'*vp;
	//			rpri = h.A*xp-h.b;
	//		end
	//		if (isfq) 
	//			rdua = rdua+f.q*xp; 
	//		end
	//		rcen = -diag(Lp)*ge-invt*ones(m,1);
	//		rt = [rdua;rcen;rpri];
	//		
	//		if (norm(rt) <= (1-alpha*s)*oldnorm)
	//			break;
	//		else
	//			s = beta*s;
	//		end
	//	end
	//	x = xp;
	//	L = Lp;
	//	v = vp;
	// Repeat steps from above
	// I guess that incorporating code 

	//	I = I+1;

	//	%%%%% IV. Check if x is optimal enough
	//	if (norm(rpri) <= epsf && norm(rdua) <= epsf && dgap <= eps)
	//		break;
	//	end

	} // while (norm(rpri) <= epsf && norm(rdua) <= epsf && dgap <= eps)

	// calculate fo and assign
}

// vectorized evaluations
// TODO: add conditions for quadratic g
/*
void eval_g(
	cublasHandle_t handle,
	G *g,
	const float *x, // x vector
	float *ge,		// g(x)	
	int *n,
	int *m) {
	
	const int *one;
	cudaMalloc((void**)&one, sizeof(int));

	const int *zero; *zero = 0;

		// evaluate g.l'*x
	cublasSgemv(handle, CUBLAS_OP_T, *n, *m, one, *g.l, 1, x, 1, zero, ge, 1);
	// add constants
	cublasSaxpy(handle, *m, one, *g.c, 1, ge, 1);
}

void eval_Dg(
	G *g,
	float *Dg) {

	// might have to copy to a const
	Dg = g.l;
}
*/


// eval_lgq
// disregard for now
