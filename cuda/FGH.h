// declaration of structs that are used in ipqopt.cu

typedef struct F {
	float *c;
	float *l;
	float *q;
	int *n;
} F;

typedef struct G{
	float *c;
	float *l;
	int *m;
} G;

typedef struct H{
	float *A;
	float *b;
	int *p;
} H;
