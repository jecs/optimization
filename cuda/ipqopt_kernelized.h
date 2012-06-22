void ipqopt(float *fq, float *fl, float *fc, int n, float *gl, float *gc, int m, float *hA, float *hb, int p, float *x_init, float dgap, float mu, float eps, float epsf, float alpha, float beta, int time_execution);

__global__ void evaluate_dgap(float *ge, float *L, const int *d_m, float *dgap);

__global__ void init_B(const F *f, const G *g, const H *h, const int *d_n, const int *d_m, const int *d_p, const int *d_nmp, float *B);

__global__ void init_ge_and_L(const G *g, const float *x, float *L, float *ge, const float *dgap, const int *d_m, const int *d_n);

__global__ void build_rt_and_B(const G *g, const F *f, const H *h, const float *x, const float *L, const float *v, const float *ge, float *d_invt, int *d_n, int *d_m, int *d_p, int *d_nmp, float *rt, float *B);

__global__ void find_max_pivot(const float *A, const int *d_n, const int *d_k, int *index, float *value);

__global__ void gauss_forward(float *A, float *x, const int *d_n, int *d_k, int *pivot_index, float *pivot_value);

__global__ void gauss_backward(float *A, float *b, const int *d_n);

__global__ void find_s_oldnormsq(float *dy, float *L, float *rt, float *s, float *oldnormsq, const int *d_n, const int *d_m, const int *d_p, const int *d_nmp);

__global__ void eval_g_sdx(const G *g, const float *x, float *dy, float *s, float *ge, const int *d_m, const int *d_n);

__global__ void check_ge_neg(const float *ge, int* ge_neg, float *s, float *beta, const int *d_m);

__global__ void linesearch(const F *f, const G *g, const H *h, float *x, float *L, float *v, float *xp, float *Lp, float *vp, float *dy, float *ge, float *rt, float *s, float *invt, const int *d_n, const int *d_m, const int *d_p, float *beta);

__global__ void compare_norms(float *rt, const int *d_n, const int *d_m, const int *d_p, float *oldnorm, float *s, float *alpha, float *beta, int *keep_searching, float *normrprisq, float *normrduasq);

__device__ void log2add_vector(float *A, const int n, const int k);

__device__ void AdiagX(const float *A, const float *x, float *B, const int m, const int n, const int row, const int col);

__device__ void ATdiagX(const float *A, const float *x, float *B, const int m, const int n, const int row, const int col);

__device__ void xdotmulty(const float alpha, const float *x, const float *y, float *z, int n, int k);

__device__ void alphadiagXA(const float alpha, const float *x, const float *A, float *B, const int m, const int n, const int row, const int col);

__device__ void alphadiagXAT(const float alpha, const float *x, const float *A, float *B, const int m, const int n, const int row, const int col);

void print_host_matrix(float *A, int m, int n, char *ident);

void print_dev_matrix(float *A, int m, int n, char *ident);

int m32(int n);
