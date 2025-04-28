
#define PI 3.141592654


/* the following are for use with GEMINI */


typedef double r_real;
typedef int s_intg;

/*
typedef s_intg *itertype;
typedef r_real *vector;
typedef vector  *matrix;
*/


#ifndef MAXN
#define MAXN            50   /* maximum number of iterated parameters for gemini */
#endif

typedef s_intg itertype [MAXN];
typedef r_real vector [MAXN];
typedef vector matrix [MAXN];


