#include<stdio.h>
#include<stdlib.h>
#include "FGH.h"
#include "ipqopt_kernelized.h"
#include<cuda.h>
#include<time.h>

int main(int argc, char *argv[]) {
	// file to parse
	if(argc == 1) {
		fprintf(stderr, "please specify an input file\n");
		return 1;
	}
	else if(argc > 2) {
		fprintf(stderr, "please specify only one input file\n");
		return 1;
	}

	FILE *fp;
	fp = fopen(argv[1], "r");
	
	if(fp == NULL) {
		fprintf(stderr, "file cannot be found\n");
		return 1;
	}

	int check;

	int n;
	int m;
	int p;

	if(fscanf(fp, "%d %d %d", &n, &m, &p) != 3) {
		fprintf(stderr, "n, m, p were not specified validly\n");
		return 2;
	}
	
	int nmp = n+m+p;
	
	float *fq;
	float *fl;
	float *fc;
	fq = (float*)malloc(n*n*sizeof(float));
	fl = (float*)malloc(n*sizeof(float));
	fc = (float*)malloc(sizeof(float));
	
	float *gl;
	float *gc;
	gl = (float*)malloc(m*n*sizeof(float));
	gc = (float*)malloc(m*sizeof(float));

	float *hA;
	float *hb;
	hA = (float*)malloc(p*n*sizeof(float));
	hb = (float*)malloc(p*sizeof(float));

	float *x_init;
	float dgap;
	float mu;
	float eps;
	float epsf;
	float alpha;
	float beta;
	x_init = (float*)malloc(n*sizeof(float));

	int i = 0;
	int fscanf_val;
	int file_in_success = 1;

	fscanf_val = fscanf(fp, "%f %f %f %f %f %f", &dgap, &mu, &eps, &epsf, &alpha, &beta);
	if(fscanf_val != 6 && file_in_success == 1) {
		fprintf(stderr, "execution parameters were not specified validly");
		file_in_success = 0;
	}

	while(i < n && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", x_init+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "x_init has not been initialized validly");
			file_in_success = 0;
		}
		i++;
	}

	i = 0;
	while(i < n*n && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", fq+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "f.q was not specified validly");
			file_in_success = 0;
		}
		i++;
	}

	i = 0;
	while(i < n && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", fl+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "f.l was not specified validly");
			file_in_success = 0;
		}
		i++;
	}

	if(file_in_success == 1 && fscanf(fp, "%f", fc) < 1) {
		fprintf(stderr, "f.c was not specified validly");
		file_in_success = 0;
	}
	
	i = 0;
	while(i < m*n && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", gl+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "g.l was not specified validly");
			file_in_success = 0;
		}
		i++;
	}

	i = 0;
	while(i < m && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", gc+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "g.c was not specified validly");
			file_in_success = 0;
		}
		i++;
	}

	i = 0;
	while(i < p*n && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", hA+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "h.A was not specified correctly");
			file_in_success = 0;
		}
		i++;
	}

	i = 0;
	while(i < p && file_in_success == 1) {
		fscanf_val = fscanf(fp, "%f", hb+i);
		if(fscanf_val != 1) {
			fprintf(stderr, "h.b was not specified correctly");
			file_in_success = 0;
		}
		i++;
	}

	fclose(fp);

	if(file_in_success == 1) {
		/*
		printf("\n(n,m,p) = (%d,%d,%d)\n", n, m, p);
		printf("\ndgap = %f\n", dgap);
		printf("\nmu = %f\n", mu);
		printf("\neps = %f\n", eps);
		printf("\nepsf = %f\n", epsf);
		printf("\nalpha = %f\n", alpha);
		printf("\nbeta = %f\n", beta);
		print_host_matrix(x_init, n, 1, "x_init");
		print_host_matrix(fq, n, n, "f.q");
		print_host_matrix(fl, n, 1, "f.l");
		print_host_matrix(fc, 1, 1, "f.c");
		print_host_matrix(gl, n, m, "g.l");
		print_host_matrix(gc, 1, m, "g.c");
		print_host_matrix(hA, p, n, "h.A");
		print_host_matrix(hb, p, 1, "h.b");
		*/

		ipqopt(fq, fl, fc, n, gl, gc, m, hA, hb, p, x_init, dgap, mu, eps, epsf, alpha, beta, 0);
	}

	free((void*)fq);
	free((void*)fl);
	free((void*)fc);

	free((void*)gl);
	free((void*)gc);;

	free((void*)hA);
	free((void*)hb);

	free((void*)x_init);

	//float x_init[] = {0.8, 1, 1};
	//float h_dgap = 10;
	//float h_mu = 10;
	//float h_eps = 1e-3;
	//float h_epsf = 1e-3;
	//float h_alpha = 0.01;
	//float h_beta = 0.8;

	//float fq[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	//float fq[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	//float fl[] = {5, 3, 8};
	//float fc[] = {0};
	//int n = 3;
	
	//float gl[] = {-1, 0, 0, 0, -1, 0, 0, 0, -1};
	//float gc[] = {0, 0, 0};
	//int m = 3;
	
	//float hA[] = {1, 1, 2};
	//float hb[] = {4};
	//int p = 1;

	//int nmp = n+m+p;

	// ipqopt(fq, fl, fc, n, gl, gc, m, hA, hb, p, x_init, h_dgap, h_mu, h_eps, h_epsf, h_alpha, h_beta);
}

void ipqopt(float *fq, float *fl, float *fc, int n, float *gl, float *gc, int m, float *hA, float *hb, int p, float *x_init, float dgap, float mu, float eps, float epsf, float alpha, float beta, int time_execution) {
	// structs on host
	F *h_f;
	G *h_g;
	H *h_h;

	// structs on device
	F *f;
	G *g;
	H *h;

	// pointers to contents of device content
	float *d_fq;
	float *d_fl;
	float *d_fc;
	int *d_fn;
	
	float *d_gl;
	float *d_gc;
	int *d_gm;
	
	float *d_hA;
	float *d_hb;
	int *d_hp;

	// struct allocation
	cudaMallocHost((void**)&h_f, sizeof(F));
	cudaMallocHost((void**)&h_g, sizeof(G));
	cudaMallocHost((void**)&h_h, sizeof(H));

	cudaMalloc((void**)&f, sizeof(F));
	cudaMalloc((void**)&g, sizeof(G));
	cudaMalloc((void**)&h, sizeof(H));

	cudaMalloc((void**)&d_fq, n*n*sizeof(float));
	cudaMalloc((void**)&d_fl, n*sizeof(float));
	cudaMalloc((void**)&d_fc, sizeof(float));
	cudaMalloc((void**)&d_fn, sizeof(int));

	cudaMalloc((void**)&d_gl, m*n*sizeof(float));
	cudaMalloc((void**)&d_gc, m*sizeof(float));
	cudaMalloc((void**)&d_gm, sizeof(int));

	cudaMalloc((void**)&d_hA, p*n*sizeof(float));
	cudaMalloc((void**)&d_hb, p*sizeof(float));
	cudaMalloc((void**)&d_hp, sizeof(int));

	cudaMemcpy((void*)d_fq, (void*)fq, n*n*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_fl, (void*)fl, n*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_fc, (void*)fc, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_fn, (void*)&n, sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy((void*)d_gl, (void*)gl, m*n*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_gc, (void*)gc, m*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_gm, (void*)&m, sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy((void*)d_hA, (void*)hA, m*p*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_hb, (void*)hb, p*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)d_hp, (void*)&p, sizeof(int), cudaMemcpyHostToDevice);

	h_f->q = d_fq;
	h_f->l = d_fl;
	h_f->c = d_fc;
	h_f->n = d_fn;

	h_g->l = d_gl;
	h_g->c = d_gc;
	h_g->m = d_gm;

	h_h->A = d_hA;
	h_h->b = d_hb;
	h_h->p = d_hp;
	
	cudaMemcpy((void*)f, (void*)h_f, sizeof(F), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)g, (void*)h_g, sizeof(G), cudaMemcpyHostToDevice);
	cudaMemcpy((void*)h, (void*)h_h, sizeof(H), cudaMemcpyHostToDevice);

	// vector and constant allocation
	float *x;
	float *L;
	float *v;
	float *ge;
	
	float *xp;
	float *Lp;
	float *vp;

	float *d_dgap;
	float *d_mu;
	float *d_eps;
	float *d_epsf;
	float *d_alpha;
	float *d_beta;
	float *d_invt;
	int *d_m;
	int *d_n;
	int *d_p;
	int *d_nmp;

	float *h_dgap;
	float *h_mu;
	float *h_eps;
	float *h_epsf;
	float *h_alpha;
	float *h_beta;
	float *h_t;
	float *h_invt;
	int *h_n;
	int *h_m;
	int *h_p;
	int *h_nmp;

	cudaMalloc((void**)&x, n*sizeof(float));
	cudaMalloc((void**)&L, m*sizeof(float));
	cudaMalloc((void**)&v, p*sizeof(float));
	cudaMalloc((void**)&ge, m*sizeof(float));
	cudaMalloc((void**)&xp, n*sizeof(float));
	cudaMalloc((void**)&Lp, m*sizeof(float));
	cudaMalloc((void**)&vp, p*sizeof(float));
	// x is set in timer loop
	// v is set in timer loop

	// host allocation of constants and parameters
	cudaMallocHost((void**)&h_dgap, sizeof(float));
	cudaMallocHost((void**)&h_mu, sizeof(float));
	cudaMallocHost((void**)&h_eps, sizeof(float));
	cudaMallocHost((void**)&h_epsf, sizeof(float));
	cudaMallocHost((void**)&h_alpha, sizeof(float));
	cudaMallocHost((void**)&h_beta, sizeof(float));
	cudaMallocHost((void**)&h_invt, sizeof(float));
	cudaMallocHost((void**)&h_t, sizeof(float));
	cudaMallocHost((void**)&h_n, sizeof(int));
	cudaMallocHost((void**)&h_m, sizeof(int));
	cudaMallocHost((void**)&h_p, sizeof(int));
	cudaMallocHost((void**)&h_nmp, sizeof(int));
	
	*h_dgap = dgap;
	*h_mu = mu;
	*h_eps = eps;
	*h_epsf = epsf;
	*h_alpha = alpha;
	*h_beta = beta;
	*h_invt = 0;
	*h_t = 0;
	*h_n = n;
	*h_m = m;
	*h_p = p;
	*h_nmp = n+m+p;
	
	// device allocation of constants and parameters
	cudaMalloc((void**)&d_dgap, sizeof(float));
	cudaMalloc((void**)&d_mu, sizeof(float));
	cudaMalloc((void**)&d_eps, sizeof(float));
	cudaMalloc((void**)&d_epsf, sizeof(float));
	cudaMalloc((void**)&d_alpha, sizeof(float));
	cudaMalloc((void**)&d_beta, sizeof(float));
	cudaMalloc((void**)&d_invt, sizeof(float));
	cudaMalloc((void**)&d_n, sizeof(int));
	cudaMalloc((void**)&d_m, sizeof(int));
	cudaMalloc((void**)&d_p, sizeof(int));
	cudaMalloc((void**)&d_nmp, sizeof(int));
	
	// system
	float *rt;
	float *B;

	cudaMalloc((void**)&rt, *h_nmp*sizeof(float));
	cudaMalloc((void**)&B, *h_nmp**h_nmp*sizeof(float));
	cudaMemset((void*)rt, 0, *h_nmp*sizeof(float));
	cudaMemset((void*)B, 0, *h_nmp**h_nmp*sizeof(float));

	// solution of system
	float *B_gauss;
	float *dy;
	cudaMalloc((void**)&B_gauss, *h_nmp**h_nmp*sizeof(float));
	cudaMalloc((void**)&dy, *h_nmp*sizeof(float));
	cudaMemset((void*)B_gauss, 0, *h_nmp**h_nmp*sizeof(float));
	cudaMemset((void*)dy, 0, *h_nmp*sizeof(float));

	// line search parameters
	float *s;
	int *ge_neg;
	int *h_ge_neg;
	float *oldnormsq;
	cudaMalloc((void**)&s, sizeof(float));
	cudaMalloc((void**)&ge_neg, sizeof(int));
	cudaMallocHost((void**)&h_ge_neg, sizeof(int));
	cudaMalloc((void**)&oldnormsq, sizeof(float));

	// DEBUG
	
	//	print_dev_matrix(h_f->q, n, n, "f.q");
	//	print_dev_matrix(h_f->l, n, 1, "f.l");
	//	print_dev_matrix(h_f->c, 1, 1, "f.c");

	//	print_dev_matrix(h_g->l, n, m, "g.l");
	//	print_dev_matrix(h_g->c, 1, m, "g.c");

	//	print_dev_matrix(h_h->A, p, n, "h.A");
	//	print_dev_matrix(h_h->b, p, 1, "h.b");
	
	
	// used in forward gaussian elimination
	float *max_val;
	int *index;
	int *loop_index;
	cudaMalloc((void**)&max_val, sizeof(float));
	cudaMalloc((void**)&index, sizeof(int));
	cudaMalloc((void**)&loop_index, sizeof(int));

	float *h_max_val;
	int *h_index;
	cudaMallocHost((void**)&h_max_val, sizeof(float));
	cudaMallocHost((void**)&h_index, sizeof(int));

	// line search
	int *h_keep_searching;
	int *keep_searching;
	cudaMallocHost((void**)&h_keep_searching, sizeof(int));
	cudaMalloc((void**)&keep_searching, sizeof(int));

	// norms used to determine exit conditions
	float *h_normrprisq;
	float *h_normrduasq;
	float *normrprisq;
	float *normrduasq;
	cudaMallocHost((void**)&h_normrprisq, sizeof(float));
	cudaMallocHost((void**)&h_normrduasq, sizeof(float));
	cudaMalloc((void**)&normrprisq, sizeof(float));
	cudaMalloc((void**)&normrduasq, sizeof(float));

	// counter
	int I;

	// grids used for the execution of kernels	
	const dim3 grid1(1,1,1);
	const dim3 block1(m32(*h_nmp), *h_nmp, 1);
	
	const dim3 grid2(*h_m, 1, 1);
	const dim3 block2(1, m32(*h_n), 1);
		
	const dim3 grid3(*h_nmp, 1, 1);
	const dim3 block3(1, m32(*h_nmp), 5);

	const dim3 grid4(1,1,1);
	const dim3 block4(m32(*h_nmp),1,1);
		
	// grid  - (1, nmp+1, 1)
	// block - (nmp, 1, 1)
	const dim3 grid5(1, *h_nmp+1, 1);
	const dim3 block5(m32(*h_nmp), 1, 1);
		
	// grid  - (1,  1, 1)
	// block - (n,  1, 1)
	const dim3 grid6(1, 1, 1);
	const dim3 block6(m32(*h_nmp), 1, 1);

	// grid  - (1,     1, 1)
	// block - (max(nmp,2*m), 1, 1)
	const dim3 grid7(1, 1, 1);
	const dim3 block7(m32(*h_nmp > 2**h_m ? *h_nmp : 2**h_m), 1, 1);

	// blocks  - (m, 1, 1)
	// threads - (1, n, 1)
	const dim3 grid8(*h_m, 1, 1);
	const dim3 block8(1, m32(*h_n), 1);

	// grid - (1, 1, 1)
	// block - (m, 1, 1)
	const dim3 grid9(1, 1, 1);
	const dim3 block9(m32(*h_m), 1, 1);
	
	// grid  - (max(n,m,p), 1, 1)
	// block - (1, max(n,m,p), 6)
	const dim3 grid10(*h_nmp, 1, 1);
	const dim3 block10(1, m32(*h_nmp), 6);

	// grid - (1, 1, 1)
	// block - (nmp, 1, 1)
	const dim3 grid11(1, 1, 1);
	const dim3 block11(m32(*h_nmp), 1, 1);

	// grid - (1, 1, 1)
	// block - (m, 1, 1)
	const dim3 grid12(1, 1, 1);
	const dim3 block12(m32(*h_m), 1, 1);

	// streams
	const int stream_count = 14;
	cudaStream_t stream[stream_count];
	for(int i = 0;i < stream_count;i++) {
		cudaStreamCreate(&stream[i]);
	}

	//+++++ Initialize all parameters
	*h_dgap = dgap;
	*h_mu = mu;
	*h_eps = eps;
	*h_epsf = epsf;
	*h_alpha = alpha;
	*h_beta = beta;
	*h_invt = 0;
	*h_t = 0;
	*h_n = n;
	*h_m = m;
	*h_p = p;
	*h_nmp = n+m+p;

	cudaMemcpyAsync((void*)d_dgap, (void*)h_dgap, sizeof(float), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync((void*)d_mu, (void*)h_mu, sizeof(float), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync((void*)d_eps, (void*)h_eps, sizeof(float), cudaMemcpyHostToDevice, stream[2]);
	cudaMemcpyAsync((void*)d_epsf, (void*)h_epsf, sizeof(float), cudaMemcpyHostToDevice, stream[3]);
	cudaMemcpyAsync((void*)d_alpha, (void*)h_alpha, sizeof(float), cudaMemcpyHostToDevice, stream[4]);
	cudaMemcpyAsync((void*)d_beta, (void*)h_beta, sizeof(float), cudaMemcpyHostToDevice, stream[6]);
	cudaMemcpyAsync((void*)d_invt, (void*)h_invt, sizeof(float), cudaMemcpyHostToDevice, stream[7]);
	cudaMemcpyAsync((void*)d_n, (void*)h_n, sizeof(int), cudaMemcpyHostToDevice, stream[8]);
	cudaMemcpyAsync((void*)d_m, (void*)h_m, sizeof(int), cudaMemcpyHostToDevice, stream[9]);
	cudaMemcpyAsync((void*)d_p, (void*)h_p, sizeof(int), cudaMemcpyHostToDevice, stream[10]);
	cudaMemcpyAsync((void*)d_nmp, (void*)h_nmp, sizeof(int), cudaMemcpyHostToDevice, stream[11]);
	
	cudaMemcpyAsync((void*)x, (void*)x_init, n*sizeof(float), cudaMemcpyHostToDevice, stream[12]);
	cudaMemsetAsync((void*)v, 0, p*sizeof(float), stream[13]);
	cudaDeviceSynchronize();

	I = 0;

	//+++++ Initialize B, ge, and L
	init_B<<<grid1, block1>>>(f, g, h, d_n, d_m, d_p, d_nmp, B);
	init_ge_and_L<<<grid2, block2>>>(g, x, L, ge, d_dgap, d_m, d_n);

	// DEBUG
		
	//	print_dev_matrix(B, *h_nmp, *h_nmp, "B_init");
	//	print_dev_matrix(x, *h_n, 1, "x_init");
	//	print_dev_matrix(L, *h_m, 1, "L_init");
	//	print_dev_matrix(v, *h_p, 1, "v_init");
	//	print_dev_matrix(ge, *h_m, 1, "ge_init");

	do {
		//+++++ I. Determine t
		*h_t = *h_mu**h_m/(*h_dgap);
		*h_invt = 1/(*h_t);
		
		cudaMemcpy((void*)d_invt, (void*)h_invt, sizeof(float), cudaMemcpyHostToDevice);

		//+++++ II. Compute dy
		// 1. Calculate and build rt and B
		build_rt_and_B<<<grid3, block3>>>(g, f, h, x, L, v, ge, d_invt, d_n, d_m, d_p, d_nmp, rt, B);
		
		// DEBUG
		//	print_dev_matrix(rt, nmp, 1, "rt");
		//	print_dev_matrix(B, nmp, nmp, "B");

		// 2. Ax=b
		cudaMemcpyAsync((void*)B_gauss, (void*)B, *h_nmp**h_nmp*sizeof(float), cudaMemcpyDeviceToDevice, stream[0]);
		cudaMemcpyAsync((void*)dy, (void*)rt, *h_nmp*sizeof(float), cudaMemcpyDeviceToDevice, stream[1]);
		cudaMemsetAsync((void*)loop_index, 0, sizeof(int), stream[2]);
		cudaDeviceSynchronize();

		find_max_pivot<<<grid4, block4>>>(B_gauss, d_nmp, loop_index, index, max_val);
		// forward gauss
		for(int i = 0;i < *h_nmp-1;i++) {
			gauss_forward<<<grid5, block5>>>(B_gauss, dy, d_nmp, loop_index, index, max_val);
		}

		// takes results from gauss_forward and solves the system
		gauss_backward<<<grid6, block6>>>(B_gauss, dy, d_nmp);

		// DEBUG
		
		//	print_dev_matrix(B_gauss, nmp, nmp, "B_gauss");
		//	print_dev_matrix(dy, nmp, 1, "dy");

		//+++++ III. Line search
		find_s_oldnormsq<<<grid7, block7>>>(dy, L, rt, s, oldnormsq, d_n, d_m, d_p, d_nmp);

		// DEBUG
		
		//	print_dev_matrix(s, 1, 1, "s");
		//	print_dev_matrix(oldnormsq, 1, 1, "oldnormsq");
		

		*h_ge_neg = 0;
		cudaMemset((void*)ge_neg, 0, sizeof(int));
		while(*h_ge_neg == 0) {
			eval_g_sdx<<<grid8, block8>>>(g, x, dy, s, ge, d_m, d_n);
			// check_ge_neg
			check_ge_neg<<<grid9, block9>>>(ge, ge_neg, s, d_beta, d_m);
			cudaMemcpy((void*)h_ge_neg, (void*)ge_neg, sizeof(int), cudaMemcpyDeviceToHost);
		}

		// DEBUG
			// print_dev_matrix(ge, m, 1, "ge_neg");

		*h_keep_searching = 1;
		cudaMemset((void*)keep_searching, 0, sizeof(int));

		do {
			linesearch<<<grid10, block10>>>(f, g, h, x, L, v, xp, Lp, vp, dy, ge, rt, s, d_invt, d_n, d_m, d_p, d_beta);
			compare_norms<<<grid11, block11>>>(rt, d_n, d_m, d_p, oldnormsq, s, d_alpha, d_beta, keep_searching, normrprisq, normrduasq);
			cudaMemcpy((void*)h_keep_searching, (void*)keep_searching, sizeof(int), cudaMemcpyDeviceToHost);
		} while(*h_keep_searching == 1);

		// DEBUG
		
		//	print_dev_matrix(xp, n, 1, "xp");
		//	print_dev_matrix(Lp, m, 1, "Lp");
		//	print_dev_matrix(vp, p, 1, "vp");
		//	print_dev_matrix(rt, nmp, 1, "rtp");
		//	print_dev_matrix(ge, m, 1, "ge_line");

		//	print_dev_matrix(s, 1, 1, "s");
		//	printf("\nkeep_searching: %d\n", *h_keep_searching);
		//	print_dev_matrix(normrprisq, 1, 1, "norm rpri sq");
		//	print_dev_matrix(normrduasq, 1, 1, "norm rdua sq");
		

		// copy normrduasq, normrprisq
		// copy xp, Lp vp into x, L, v
		cudaMemcpyAsync((void*)x, (void*)xp, *h_n*sizeof(float), cudaMemcpyDeviceToDevice, stream[0]);
		cudaMemcpyAsync((void*)L, (void*)Lp, *h_m*sizeof(float), cudaMemcpyDeviceToDevice, stream[1]);
		cudaMemcpyAsync((void*)v, (void*)vp, *h_p*sizeof(float), cudaMemcpyDeviceToDevice, stream[2]);
		cudaMemcpyAsync((void*)h_normrprisq, (void*)normrprisq, sizeof(float), cudaMemcpyDeviceToHost, stream[3]);
		cudaMemcpyAsync((void*)h_normrduasq, (void*)normrduasq, sizeof(float), cudaMemcpyDeviceToHost, stream[4]);
		cudaDeviceSynchronize();

		I += 1;

		//++++++ IV. Evaluate dgap and check if x is optimal enough
		evaluate_dgap<<<grid12, block12>>>(ge, L, d_m, d_dgap);
		cudaMemcpy((void*)h_dgap, (void*)d_dgap, sizeof(float), cudaMemcpyDeviceToHost);
	//+++++ END LOOP
	} while((sqrtf(*h_normrprisq) > (*h_epsf) || sqrtf(*h_normrduasq) > (*h_epsf) || *h_dgap > *h_eps) && I < 200);	

	// print out results
	
	print_dev_matrix(x, n, 1, "optimal x");
	print_dev_matrix(L, m, 1, "optimal L");
	print_dev_matrix(v, p, 1, "optimal v");
	printf("\nIterations: %d\n", I);
	
	// free up memory
	// TODO: Check that everything is being freed
	cudaFreeHost((void*)h_f);
	cudaFreeHost((void*)h_g);
	cudaFreeHost((void*)h_h);

	cudaFree((void*)f);
	cudaFree((void*)g);
	cudaFree((void*)h);

	cudaFree((void*)d_fq);
	cudaFree((void*)d_fl);
	cudaFree((void*)d_fc);
	cudaFree((void*)d_fn);

	cudaFree((void*)d_gl);
	cudaFree((void*)d_gc);
	cudaFree((void*)d_gm);

	cudaFree((void*)d_hA);
	cudaFree((void*)d_hb);
	cudaFree((void*)d_hp);
	
	cudaFree((void*)x);
	cudaFree((void*)L);
	cudaFree((void*)v);
	cudaFree((void*)ge);
	cudaFree((void*)xp);
	cudaFree((void*)Lp);
	cudaFree((void*)vp);

	cudaFreeHost((void*)h_dgap);
	cudaFreeHost((void*)h_mu);
	cudaFreeHost((void*)h_eps);
	cudaFreeHost((void*)h_epsf);
	cudaFreeHost((void*)h_alpha);
	cudaFreeHost((void*)h_beta);
	cudaFreeHost((void*)h_invt);
	cudaFreeHost((void*)h_t);
	cudaFreeHost((void*)h_n);
	cudaFreeHost((void*)h_m);
	cudaFreeHost((void*)h_p);
	cudaFreeHost((void*)h_nmp);

	cudaFree((void*)d_dgap);
	cudaFree((void*)d_mu);
	cudaFree((void*)d_eps);
	cudaFree((void*)d_epsf);
	cudaFree((void*)d_alpha);
	cudaFree((void*)d_beta);
	cudaFree((void*)d_invt);
	cudaFree((void*)d_n);
	cudaFree((void*)d_m);
	cudaFree((void*)d_p);
	cudaFree((void*)d_nmp);

	cudaFree((void*)rt);
	cudaFree((void*)B);

	cudaFree((void*)B_gauss);
	cudaFree((void*)dy);

	cudaFree((void*)s);
	cudaFree((void*)ge_neg);
	cudaFreeHost((void*)h_ge_neg);
	cudaFree((void*)oldnormsq);

	cudaFree((void*)max_val);
	cudaFree((void*)index);
	cudaFree((void*)loop_index);

	cudaFreeHost((void*)h_max_val);
	cudaFreeHost((void*)h_index);

	cudaFreeHost((void*)h_keep_searching);
	cudaFree((void*)keep_searching);

	cudaFreeHost((void*)h_normrprisq);
	cudaFreeHost((void*)h_normrduasq);
	cudaFree((void*)normrprisq);
	cudaFree((void*)normrduasq);
	return;
}

// grid  - (1, 1, 1)
// block - (m, 1, 1)
__global__ void evaluate_dgap(float *ge, float *L, const int *d_m, float *dgap) {
	const int N = 1000;
	__shared__ float scratch[N];

	int col = threadIdx.x;
	int m = *d_m;

	if(col < m) {
		scratch[col] = -ge[col]*L[col];
	}
	__syncthreads();
	if(col < m) {
		log2add_vector(scratch, m, col);
	}
	__syncthreads();
	*dgap = scratch[0];
}

// tested
__global__ void init_B(const F *f, const G *g, const H *h, const int *d_n, const int *d_m, const int *d_p, const int *d_nmp, float *B) {
	// (1, 1) - f.q
	// (1, 2) - g.l
	// (1, 3) - h.A'
	// (3, 1) - h.A

	int row = threadIdx.x + blockDim.x * blockIdx.x;
	int col = threadIdx.y + blockDim.y * blockIdx.y;
	int n = *d_n;
	int m = *d_m;
	int p = *d_p;
	int nmp = *d_nmp;

	if(row < n) {
		if(col < n) {
			B[col*nmp+row] = f->q[col*n+row];
		}
		else if(col < n+m) {
			B[col*nmp+row] = g->l[(col-n)*n+row];
		}
		else if(col < nmp) {
			B[col*nmp+row] = h->A[row*p+(col-(n+m))];
		}
	}
	else if(row >= (n+m) && row < nmp && col < n) {
		B[col*nmp+row] = h->A[col*p+row-(n+m)];
	}
}

// tested
// blocks  - (m, 1, 1)
// threads - (1, n, 1)
__global__ void init_ge_and_L(const G *g, const float *x, float *L, float *ge, const float *dgap, const int *d_m, const int *d_n) {
	const int N = 1000;
	__shared__ float ge_scratch[N];
	__shared__ float gl_scratch[N];
	__shared__ float gc_scratch;
	
	// dims
	const int row = blockIdx.x;
	const int col = threadIdx.y;
	const int m = *d_m;
	const int n = *d_n;

	if(col < n && row < m) {
		// store x in ge_scratch, for fewer memory acceses
		ge_scratch[col] = x[col];
		// store column 'row' of g->l
		gl_scratch[col] = g->l[row*n+col];
	}
	if(col == 0 && row < m) {
		gc_scratch = g->c[row];
	}
	__syncthreads();

	// g.l'diag(x)
	if(col < n) {
		ge_scratch[col] = gl_scratch[col]*ge_scratch[col];
	}
	__syncthreads();

	// sum the columns
	log2add_vector(ge_scratch, n, col);
	__syncthreads();

	// add g.c and return results
	if(row < m && col == 0) {
		ge_scratch[0] = ge_scratch[0]+gc_scratch;
		ge[row] = ge_scratch[0];
		L[row] = -(*dgap/ge_scratch[0]);
	}
	//else if(row < m && col == 1) {
	//	L[row] = -(*dgap/ge_scratch[0]);
	//}
}

// tested
// blocks  - (nmp, 1, 1)
// threads - (1, nmp, 5)
__global__ void build_rt_and_B(const G *g, const F *f, const H *h, const float *x, const float *L, const float *v, const float *ge, float *d_invt, int *d_n, int *d_m, int *d_p, int *d_nmp, float *rt, float *B) {
	const int N = 1000;
	__shared__ float rdua_scratch[3][N];
	__shared__ float rcen_scratch;
	__shared__ float rpri_scratch[N];

	const int row = blockIdx.x;
	const int col = threadIdx.y;
	const int k = threadIdx.z; // index for computation
	const float invt = *d_invt;
	const int n = *d_n;
	const int m = *d_m;
	const int p = *d_p;
	const int nmp = *d_nmp;

	// prefetch
	//switch(k) {
	//	case 0: 

	// perform matrix mults
	// rdua0 = Dg'*L or g.l*L
	switch(k) {
		// rdua0 = g.l*L
		case 0: AdiagX(g->l, L, rdua_scratch[0], n, m, row, col); break;
		// rdua1 = f.q*x
		case 1: AdiagX(f->q, x, rdua_scratch[1], n, n, row, col); break;
		// rdua2 = h.A'*v
		case 2: ATdiagX(h->A, v, rdua_scratch[2], n, p, row, col); break;
		// rcen = -L.*ge
		case 3: if(col == 0 && row < m) rcen_scratch = -L[row]*ge[row]; break;
		// rpri = h.A*x
		case 4: AdiagX(h->A, x, rpri_scratch, p, n, row, col); break;
	}
	__syncthreads();

	// log2 sum of vectors
	switch(k) {
		case 0: log2add_vector(rdua_scratch[0], m, col); break;
		case 1: log2add_vector(rdua_scratch[1], n, col); break;
		case 2: log2add_vector(rdua_scratch[2], p, col); break;
		case 4: log2add_vector(rpri_scratch, n, col); break;
	}
	__syncthreads();

	// add or subtract vectors
	if(col == 0) {
		switch(k) {
			// rdua0 += f.l+rdua1+rdua2
			case 0: if(row < n) rdua_scratch[0][0] += (f->l[row]+rdua_scratch[1][0]+rdua_scratch[2][0]); break;
			// rdua1 += rdua2
			// case 1: if(row < n) rdua_scratch[1][0] += rdua_scratch[2][0]; break;
			// rcen -= invt*ones
			case 3: if(row < m) rcen_scratch -= invt; break;
			// pri -= h.b
			case 4: if(row < p) rpri_scratch[0] -= h->b[row]; break;
		}
	}
	__syncthreads();

	// rdua, rcen, rpri have been built
	
	// time to update block matrix
	// B should already have the following values
	// (1, 1) - f.q
	// (1, 2) - g.l
	// (1, 3) - h.A'
	// (3, 1) - h.A
	// 
	// (2, 1) will be set to -diag(L)*Dg
	// (2, 2) will be set to -diag(ge)
	// rdua_scratch[1] will be reused for calculating -diag(L)*g.l'

	if(k == 0) {
		alphadiagXAT(-1, L, g->l, rdua_scratch[1], m, n, row, col); 
	}
	__syncthreads();

	// send everything back to global memory
	// rt
	if(k < 3 && col == 0) {
		switch(k) {
			// rdua
			case 0: if(row < n) rt[row] = rdua_scratch[0][0]; break;
			// rcen
			case 1: if(row < m) rt[n+row] = rcen_scratch; break;
			// rpri
			case 2: if(row < p) rt[n+m+row] = rpri_scratch[0]; break;
		}
	}
	// B(2,1) = -diag(L)*Dg, which is stored in rdua_scratch[1]
	else if(k == 3 && row < m && col < n) {
		B[col*nmp+n+row] = rdua_scratch[1][col];
	}
	// B(2,2) = -diag(ge)
	else if(k == 4 && col == row && col < m) {
		B[(col+n)*nmp+n+row] = -ge[col];
	}
	__syncthreads();
}

// tested
// blocks  - (1,   1, 1)
// threads - (nmp, 1, 1)
// log2(n) execution time
__global__ void find_max_pivot(const float *A, const int *d_n, const int *d_k, int *index, float *value) {
	const int N = 1000;
	__shared__ float Ak[N];
	__shared__ int indices[N];

	int row2 = 0;
	int i;
	int i_prev;

	float a = 0;
	float b = 0;
	const int n = *d_n;
	const int k = *d_k;
	const int offset = k;

	const int row = threadIdx.x;
	const int row_off = row-offset;

	if(row < n && row_off >= 0) {
		Ak[row_off] = A[k*n+row];
		indices[row_off] = row;
	}
	i_prev = (n-offset);
	i = (i_prev+1)/2;
	__syncthreads();

	while(row_off >= 0 && row_off < i) {
		row2 = row_off+i;
		if(row2 < i_prev) {
			a = Ak[row_off];
			b = Ak[row2];
			if(a < 0) {
				a = -a;
			}
			if(b < 0) {
				b = -b;
			}
			if(b > a) {
				Ak[row_off] = Ak[row2];
				indices[row_off] = indices[row2];
			}
		}
		i_prev = i;
		i = ((i == 1) ? 0 : (i+1)/2);
	}
	__syncthreads();

	if(row == 0) {
		*index = indices[0];
		*value = Ak[0];
	}
}

// grid  - (1, nmp+1, 1)
// block - (nmp, 1, 1)
__global__ void gauss_forward(float *A, float *x, const int *d_n, int *d_k, int *pivot_index, float *pivot_value) {
	const int N = 1000;

	__shared__ float Ac[N];
	__shared__ float Ajk[N];

	const int row = threadIdx.x;
	const int col = blockIdx.y;
	const int n = *d_n;
	const int k = *d_k;

	// BEGIN NEW ADDITIONS
	__shared__ float Apiv[N];
	__shared__ int indices[N];
	int row2 = 0;
	int i;
	int i_prev;
	float a = 0;
	float b = 0;
	const int offset = k+1;
	const int row_off = row-offset;
	// END NEW ADDITIONS

	// load slice of A
	if(row < n && col < n) {
		// pivoting and storage of block column
		if(row == k) Ac[row] = A[col*n+*pivot_index];
		else if(row == *pivot_index) Ac[row] = A[col*n+k];
		else Ac[row] = A[col*n+row];
		
		// storage of k column
		if(row > k && row == *pivot_index) {
			Ajk[row] = A[k*n+k];
		}
		else if(row > k) {
			Ajk[row] = A[k*n+row];
		}
	}
	else if(row < n && col == n) {
		// pivoting and storage
		if(row == k) Ac[row] = x[*pivot_index];
		else if(row == *pivot_index) Ac[row] = x[k];
		else Ac[row] = x[row];
		
		if(row > k && row == *pivot_index) {
			Ajk[row] = A[k*n+k];
		}
		else if(row > k) {
			Ajk[row] = A[k*n+row];
		}
	}
	__syncthreads();

	// subtract row
	if(row > k && row < n && col <= n) {
		Ac[row] -= Ajk[row]/(*pivot_value)*Ac[k];
	}

	// BEGIN NEW ADDITIONS 2
	if(col == offset && row < n && row_off >= 0) {
		Apiv[row_off] = Ac[row];
		indices[row_off] = row;
		i_prev = (n-offset);
		i = (i_prev+1)/2;
	
		while(row_off < i) {
			row2 = row_off+i;
			if(row2 < i_prev) {
				a = Apiv[row_off];
				b = Apiv[row2];
				if(a < 0) {
					a = -a;
				}
				if(b < 0) {
					b = -b;
				}
				if(b > a) {
					Apiv[row_off] = Apiv[row2];
					indices[row_off] = indices[row2];
				}
			}
			i_prev = i;
			i = ((i == 1) ? 0 : (i+1)/2);
		}
	}
	__syncthreads();

	// store back
	if(col == offset && row == 0) {
		*pivot_index = indices[0];
		*pivot_value = Apiv[0];
	}
	if(row > k && row < n && col == k) {
		A[col*n+row] = 0;
	}
	else if(row < n && col < n) {
		A[col*n+row] = Ac[row];
	}
	else if(row < n && col == n) {
		x[row] = Ac[row];
	}
	if(row == 0 && col == 0) {
		*d_k = k+1;
	}
}

// takes results from gauss_forward and solves the system
// grid  - (1,  1, 1)
// block - (n,  1, 1)
__global__ void gauss_backward(float *A, float *b, const int *d_n) {
	const int N = 1000;
	__shared__ float Ac[N];
	__shared__ float bc[N];

	const int row = threadIdx.x;
	const int col = blockIdx.y;
	const int n = *d_n;

	for(int k = n-1;k >= 0;k--) {
		// load column
		if(row <= k && col == 0) {
			Ac[row] = A[k*n+row];
			bc[row] = b[row];
		}
		__syncthreads();

		// normalize
		if(row == k && col == 0) {
			bc[row] /= Ac[row];
			Ac[row] = 1;
		}
		__syncthreads();

		// backwards row elimination
		if(row < k && col == 0) {
			bc[row] -= Ac[row]*bc[k];
			Ac[row] = 0;
		}
		__syncthreads();

		// write back
		if(row <= k && col == 0) {
			A[k*n+row] = Ac[row];
			b[row] = bc[row];
		}
		__syncthreads();
	}

	if(row <= n && col == 0) {
		b[row] = -b[row];
	}
}

// grid  - (1,     1, 1)
// block - (max(nmp,2*m), 1, 1)
// dx - 0 to n-1
// dL - n to n+m-1
// dv - n+m to n+m+p-1
__global__ void find_s_oldnormsq(float *dy, float *L, float *rt, float *s, float *oldnormsq, const int *d_n, const int *d_m, const int *d_p, const int *d_nmp) {
	const int N = 1000;

	__shared__ float scratch[N];
	__shared__ float scratch2[N];
	__shared__ float scratch3[N];

	int row = threadIdx.x;
	int n = *d_n;
	int m = *d_m;
	int p = *d_p;
	int nmp = *d_nmp;

	// log2 indices
	__shared__ int i;
	__shared__ int i_prev;
	int row2 = 0;
	
	// calculate norm of rt
	if(row < nmp) {
		scratch[row] = rt[row];
		scratch[row] *= scratch[row];
	}
	__syncthreads();

	log2add_vector(scratch, nmp, row);
	__syncthreads();

	if(row == 0) {
		*oldnormsq = scratch[0];
	}
	__syncthreads();
	
	// load dL and L
	if(row < m) {
		scratch[row] = dy[n+row];
		scratch2[row] = L[row];
	}
	__syncthreads();

	// calculate s max
	// dL in scratch
	// L in scratch2
	// store smax in scratch3
	if(row < m) {
		if(scratch[row] < 0 && (-scratch[row]) > scratch2[row]) {
			scratch3[row] = -scratch2[row]/scratch[row];
		}
		else {
			scratch3[row] = 1;
		}
	}

	if(row == 0) {
		i = (m+1)/2;
		i_prev = m;
	}
	__syncthreads();

	float a;
	float b;

	// find min smax
	while(i > 0) {
		if(row < i) {
			row2 = row+i;
			if(row2 < i_prev) {
				a = scratch3[row];
				b = scratch3[row2];
				if(a > b) {
					scratch3[row] = b;
				}
			}
		}
		if(row == 0) {
			i_prev = i;
			i = (i == 1) ? 0 : (i+1)/2;
		}
		__syncthreads();
	}
	if(row == 0) {
		*s = 0.99f*scratch3[0];
	}
}

// block  - (m, 1, 1)
// thread - (1, n, 1)
__global__ void eval_g_sdx(const G *g, const float *x, float *dy, float *s, float *ge, const int *d_m, const int *d_n) {
	const int N = 1000;
	__shared__ float ge_scratch[N];
	__shared__ float x_scratch[N];
	
	// dims
	int row = blockDim.x * blockIdx.x;
	int col = threadIdx.y;
	int m = *d_m;
	int n = *d_n;

	if(col < n) {
		x_scratch[col] = x[col]+*s*dy[col];
	}

	// g.l'*x
	ATdiagX(g->l, x_scratch, ge_scratch, m, n, row, col);
	__syncthreads();

	// sum the columns
	log2add_vector(ge_scratch, n, col);
	__syncthreads();

	if(col == 0 && row < m) {
		ge_scratch[0] += g->c[row];
	}
	__syncthreads();


	// return results
	if(row < m && col == 0) {
		ge[row] = ge_scratch[0];
	}
	__syncthreads();
}

// block - (1, 1, 1)
// thread - (m, 1, 1)
__global__ void check_ge_neg(const float *ge, int* ge_neg, float *s, float *beta, const int *d_m) {
	const int N = 1000;
	__shared__ int ge_neg_elem[N];
	__shared__ int i;
	__shared__ int i_prev;
	
	int row = threadIdx.x;
	int row2;
	int m = *d_m;
	
	if(row < m) {
		ge_neg_elem[row] = (ge[row] < 0);
	}
	if(row == 0) {
		i = (m+1)/2;
		i_prev = m;
	}

	// log2 verification

	while(i > 0) {
		if(row < i) {
			row2 = row+i;
			if(row2 < i_prev) {
				ge_neg_elem[row] = (ge_neg_elem[row] == 1) && (ge_neg_elem[row2] == 1);
			}
		}
		if(row == 0) {
			i_prev = i;
			i = (i == 1) ? 0 : (i+1)/2;
		}
		__syncthreads();
	}

	if(row == 0) {
		if(ge_neg_elem[0] == 0) {
			*s *= *beta;
		}
		*ge_neg = ge_neg_elem[0];
	}
}

// grid  - (max(n,m,p), 1, 1)
// block - (1, max(n,m,p), 6)
__global__ void linesearch(const F *f, const G *g, const H *h, float *x, float *L, float *v, float *xp, float *Lp, float *vp, float *dy, float *ge, float *rt, float *s, float *invt, const int *d_n, const int *d_m, const int *d_p, float *beta) {
	const int N = 790;
	__shared__ float rdua_scratch[3][N];
	__shared__ float rcen_scratch;
	__shared__ float rpri_scratch[N];
	__shared__ float ge_scratch[N];

	int row = blockIdx.x;
	int col = threadIdx.y;
	int k = threadIdx.z; // index for computation
	int n = *d_n;
	int m = *d_m;
	int p = *d_p;

	// calculate and store xp, Lp, vp
	switch(k) {
		case 0: if(col < m) rdua_scratch[0][col] = L[col]+*s*dy[n+col]; break;
		case 1: if(col < n) rdua_scratch[1][col] = x[col]+*s*dy[col]; break;
		case 2: if(col < p) rdua_scratch[2][col] = v[col]+*s*dy[n+m+col]; break;
		case 3: if(col == 0 && row < m) rcen_scratch = L[row]+*s*dy[n+row]; break;
		case 4: if(col < n) rpri_scratch[col] = x[col]+*s*dy[col]; break;
		case 5: if(col < n) ge_scratch[col] = x[col]+*beta**s*dy[col]; break;
	}
	__syncthreads();

	// send xp, Lp, vp back
	if(col == 0) {
		switch(k) {
			case 0: if(row < m) Lp[row] = rdua_scratch[0][row]; break;
			case 1: if(row < n) xp[row] = rdua_scratch[1][row]; break;
			case 2: if(row < p) vp[row] = rdua_scratch[2][row]; break;
		}
	}
	__syncthreads();

	// perform matrix mults
	switch(k) {
		// rdua0 = g.l*L
		case 0: AdiagX(g->l, rdua_scratch[0], rdua_scratch[0], n, m, row, col); break;
		// rdua1 = f.q*x
		case 1: AdiagX(f->q, rdua_scratch[1], rdua_scratch[1], n, n, row, col); break;
		// rdua2 = h.A'*v
		case 2: ATdiagX(h->A, rdua_scratch[2], rdua_scratch[2], n, p, row, col); break;
		// rcen = -L.*ge
		case 3: if(col == 0 && row < m) rcen_scratch = -rcen_scratch*ge[row]; break;
		// rpri = h.A*x
		case 4: AdiagX(h->A, rpri_scratch, rpri_scratch, p, n, row, col); break;
		// ge = g.l'*x
		case 5: ATdiagX(g->l, ge_scratch, ge_scratch, m, n, row, col); break;
	}
	__syncthreads();

	// log2 sum of vectors
	switch(k) {
		case 0: log2add_vector(rdua_scratch[0], m, col); break;
		case 1: log2add_vector(rdua_scratch[1], n, col); break;
		case 2: log2add_vector(rdua_scratch[2], p, col); break;
		case 4: log2add_vector(rpri_scratch, n, col); break;
		case 5: log2add_vector(ge_scratch, n, col); break;
	}
	__syncthreads();

	// add or subtract vectors
	if(col == 0) {
		switch(k) {
			// rdua0 += f.l+rdua1+rdua2
			case 0: if(row < n) rdua_scratch[0][0] += (f->l[row]+rdua_scratch[1][0]+rdua_scratch[2][0]); break;
			// rdua1 += rdua2
			// case 1: if(row < n) rdua_scratch[1][0] += rdua_scratch[2][0]; break;
			// rcen -= invt*ones
			case 3: if(row < m) rcen_scratch -= *invt; break;
			// pri -= h.b
			case 4: if(row < p) rpri_scratch[0] -= h->b[row]; break;
			// ge += g.c
			case 5: if(row < m) ge_scratch[0] += g->c[row]; break;
		}
	}
	__syncthreads();

	// send everything back to global memory
	// rt
	if(col == 0) {
		switch(k) {
			// rdua
			case 0: if(row < n) rt[row] = rdua_scratch[0][0]; break;
			// rcen
			case 1: if(row < m) rt[n+row] = rcen_scratch; break;
			// rpri
			case 2: if(row < p) rt[n+m+row] = rpri_scratch[0]; break;
			// ge
			case 5: if(row < m) ge[row] = ge_scratch[0]; break;
		}
	}
	__syncthreads();	
}

// compare_norms(rt, nmp, oldnorm, s, alpha, beta);
// grid - (1, 1, 1)
// block - (nmp, 1, 1)
__global__ void compare_norms(float *rt, const int *d_n, const int *d_m, const int *d_p, float *oldnorm, float *s, float *alpha, float *beta, int *keep_searching, float *normrprisq, float *normrduasq) {
	const int N = 1000;
	
	__shared__ float scratch[3][N];
	__shared__ float newnormsq;

	int row = threadIdx.x;
	int n = *d_n;
	int m = *d_m;
	int p = *d_p;

	// store square of rt
	if(row < n) {
		scratch[0][row] = rt[row];
		scratch[0][row] *= scratch[0][row];
	}
	else if(row < n+m) {
		scratch[1][row-n] = rt[row];
		scratch[1][row-n] *= scratch[1][row-n];
	}
	else if(row < n+m+p) {
		scratch[2][row-n-m] = rt[row];
		scratch[2][row-n-m] *= scratch[2][row-n-m];
	}
	__syncthreads();

	if(row < n) {
		log2add_vector(scratch[0], n, row);
	}
	else if(row < n+m) {
		log2add_vector(scratch[1], m, row-n);
	}
	else if(row < n+m+p) {
		log2add_vector(scratch[2], p, row-n-m); 
	}
	__syncthreads();

	if(row == 0) {
		newnormsq = scratch[0][0] + scratch[1][0] + scratch[2][0];
		if(newnormsq <= (1-*alpha**s)*((1-*alpha**s)**oldnorm)) {
			*keep_searching = 0;
			*normrprisq = scratch[2][0];
			*normrduasq = scratch[0][0];
		}
		else {
			*s *= *beta;
		}
	}
}

// sums elements of vector A of length n
// k represents the sum index
__device__ void log2add_vector(float *A, const int n, const int k) {
	int i = (n+1)/2;
	int i_prev = n;
	int k2 = 0;
	while(k < i) {
		k2 = k+i;
		if(k2 < i_prev) {
			A[k] += A[k2];
		}
		i_prev = i;
		i = (i == 1 ? 0 : (i+1)/2);
	}
}

// assings to B the row at index row of A*diag(x)
// A - (m, n)
// B - (1, n)
// x - (m, 1)
__device__ void AdiagX(const float *A, const float *x, float *B, const int m, const int n, const int row, const int col) {
	if(row < m && col < n) {
		B[col] = A[col*m+row]*x[col];
	}
}

// assigns to B the row at index row of A'*diag(x)
// A - (n, m)
// B - (1, n)
// x - (n, 1)
__device__ void ATdiagX(const float *A, const float *x, float *B, const int m, const int n, const int row, const int col) {
	if(row < m && col < n) {
		B[col] = A[row*n+col]*x[col];
	}
}

// equivalent to MATLAB's .* syntax
// ie, z = alpha*x.*y
// all vectors are of dim(n, 1)
// k is the index of the thread
__device__ void xdotmulty(const float alpha, const float *x, const float *y, float *z, int n, int k) {
	if(k < n) {
		z[k] = alpha*x[k]*y[k];
	}
}

__device__ void alphadiagXA(const float alpha, const float *x, const float *A, float *B, const int m, const int n, const int row, const int col) {
	if(row < m && col < n) {
		B[col] = alpha*x[row]*A[col*m+row];
	}
}

__device__ void alphadiagXAT(const float alpha, const float *x, const float *A, float *B, const int m, const int n, const int row, const int col) {
	if(row < m && col < n) {
		B[col] = alpha*x[row]*A[row*n+col];
	}
}

void print_host_matrix(float *A, int m, int n, char *ident) {
	printf("\n%s:\n", ident);
	for (int r = 0;r < m;r++) {
		for (int c = 0;c < n;c++) {
			if(A[c*m+r] >= 0) printf(" %f ", A[c*m+r]);
			else printf("%f ", A[c*m+r]);
		}
		printf("\n");
	}
}

void print_dev_matrix(float *A, int m, int n, char *ident) {
	float *h_A;
	cudaMallocHost((void**)&h_A, m*n*sizeof(float));
	cudaMemcpy(h_A, A, m*n*sizeof(float), cudaMemcpyDeviceToHost);
	print_host_matrix(h_A, m, n, ident);
	cudaFreeHost((void*)h_A);
}

int m32(int n) {
	if(n <= 0) {
		return n;
	}

	return 32*((n+31)/32);
}
