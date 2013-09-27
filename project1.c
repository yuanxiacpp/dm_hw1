#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double BAR = 1.0E-13;

int UPDATE_FIRST_MATRIX = 1;
int UPDATE_SECOND_MATRIX = 2;


void printMatrix(double *a, int n) {
  printf("***************** Matrix %d x %d *********************\n", n, n);
  int i = 0;
  for (i = 0; i < n * n; ++i) {
    printf("%8.3f ", a[i]);
    if ((i+1) % n == 0)
      printf("\n");
  }
  return;
  
}



void findMax(double *a, int n, double *max, int *mi, int *mj) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i == j)
	continue;
      if (fabs(a[i*n + j]) > *max) {
	*max = fabs(a[i*n + j]);
	*mi = i;
	*mj = j;
      }
    }
  }
  return;
}

//check when to stop processing the matrix
int checkZeros(double *a, int n) {
  int i,j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i != j && a[i*n + j] > BAR)
	return 1;
    }
  }
  return 0;
}
int sgn(double x) {
  if (x == 0)
    return 0;
  else if (x > 0)
    return 1;
  return -1;
}
void assignIdentity(double *a, int n) {
  int r, c;
  for (r = 0; r < n; ++r) {
    for (c = 0; c < n; ++c) {
      if (r == c)
	a[r*n+c] = 1;
      else
	a[r*n+c] = 0;
    }
  }
  return;
}
void initializeUtVt(double *s, double *u, double *v, int n, int mi, int mj) {
  assignIdentity(u, n);
  assignIdentity(v, n);


  int i = mi;
  int j = mj;
  double kk = s[i*n+i];
  double kl = s[i*n+j];
  double ll = s[j*n+j];
  

  double beta = (ll - kk)/(2 * kl);
  double t = sgn(beta) / (fabs(beta) + sqrt(beta*beta + 1));
  
  double cos = 1 / sqrt(t*t+1);
  double sin = cos * t;
  printf("cos=%f, sin=%f\n", cos, sin);
  u[i*n+i] = cos;
  u[i*n+j] = sin;
  u[j*n+i] = -sin;
  u[j*n+j] = cos;

  v[i*n+i] = cos;
  v[i*n+j] = -sin;
  v[j*n+i] = sin;
  v[j*n+j] = cos;

  return;
  
}

//flag determines which input(a or b) to free. T to free a, F to free b. 
//the return will malloc new mem for either a or b
void multiply(double *a, double *b, int n, int perseve_flag) {
  int k, i, j;
  double *result = (double*)malloc(n * n * sizeof(double));
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      double sum = 0;
      for (k = 0; k < n; ++k) 
	sum += a[i*n + k] * b[k*n + j];
      
      result[i*n + j] = sum;
    }
  }
  if (perseve_flag == 1) {
    memcpy(a, result, n*n*sizeof(double));
  }
  else {
    memcpy(b, result, n*n*sizeof(double));
  }

  free(result);
  return;
}

void transpose(double *a, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      double tmp = a[i*n + j];
      a[i*n + j] = a[j*n + i];
      a[j*n + i] = tmp;
    }
  }
  return;
}

void massage(double *u, double *s, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if (s[i*n+i] < 0) {
      s[i*n+i] = -s[i*n+i];
      int j;
      for (j = 0; j < n; j++)
	u[j*n+i] = -u[j*n+i];
    }
  }
  return;
}

void swapColumn(double *a, int n, int i, int j) {
  int row;
  for (row = 0; row < n; ++row) {
    double tmp = a[row*n + i];
    a[row*n + i] = a[row*n + j];
    a[row*n + j] = tmp;
  }
  return;
}
void swapRow(double *a, int n, int i, int j) {
  int col;
  for (col = 0; col < n; ++col) {
    double tmp = a[i*n + col];
    a[i*n + col] = a[j*n + col];
    a[j*n + col] = tmp;
  }
  return;
}
void adjustMatrix(double *s, double *u, double *v, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n - i - 1; ++j) {
      if (s[j*n + j] > s[(j+1)*n + (j+1)]) {
	swapColumn(s, n, j, j+1);
	swapColumn(u, n, j, j+1);
	swapRow(v, n, j, j+1);
      }
    }
  }
  return;
}

void jacobi(double *a, int n, double *s, double *u, double *v) {
  double *max = (double*)malloc(sizeof(double));
  int *mi = (int*)malloc(sizeof(int));
  int *mj = (int*)malloc(sizeof(int));
  *max = 0;
  *mi = 0;
  *mj = 0;

  //intialize U, V, S
  //memcpy(&u, &a, sizeof(a));
  assignIdentity(u, n);
  memcpy(s, a, n*n*sizeof(double));
  assignIdentity(v, n);
  
  int count = 0;

  while (checkZeros(s, n) != 0) {
    count++;
    printf("Processing %dth round...\n", count);
    findMax(s, n, max, mi, mj);
    //printf("max=%f, mi=%d, mj=%d\n", *max, *mi, *mj);
    double *ut = (double*)malloc(n*n*sizeof(double));
    double *vt = (double*)malloc(n*n*sizeof(double));
    initializeUtVt(s, ut, vt, n, *mi, *mj);
    
    printf("Matrix UT VT S\n");
    printMatrix(ut, n);
    printMatrix(vt, n);
    printMatrix(s, n);

    multiply(ut, s, n, UPDATE_SECOND_MATRIX);
    multiply(s, vt, n, UPDATE_FIRST_MATRIX);

    multiply(ut, u, n, UPDATE_SECOND_MATRIX);
    multiply(v, vt, n, UPDATE_FIRST_MATRIX);
    
    printf("Matrix U, S, V\n");
    printMatrix(u, n);
    printMatrix(s, n);
    printMatrix(v, n);
    getchar();
    //if (count == 3)
    //  break;
  }
  //transpose(u, n);
  //massage(u, s, n);
  //transpose(v, n);
  //adjustMatrix(s, u, v, n);
  //transpose(v, n);
  return;
  
  
}

int main() {
  int i = 0;
  int j = 0;
  //generate test case for b)
  int len1[3] = {10, 20, 40};
  double **a = (double**)malloc(3 * sizeof(double*));
  int k = 0;
  for (k = 0; k < 3; ++k) {
    int n = len1[k];
    double *tmp = (double*)malloc(n * n * sizeof(double));
    for (i = 0; i < n; ++i) {
      for (j = 0; j < n; ++j) {
	tmp[i*n + j] = sqrt(i * i + j * j);
      }
    }
    a[k] = tmp;
  }


  //generate test case for c)
  int len2 = 3;
  double *b = (double*)malloc(len2 * len2 * sizeof(double));
  for (i = 0; i < len2; ++i) {
    for (j = 0; j < len2; ++j) {
      b[i*len2 + j] = (i+1) * (i+1) + (j+1) * (j+1);
    }
  }
  //for (k = 0; k < 3; ++k) 
  //  printMatrix(a[k], len1[k]);
  int n = len2;
  printMatrix(b, n);
  //transpose(b, n);
  //printMatrix(b, n);
  double *s = (double *)malloc(n*n*sizeof(double));
  double *v = (double *)malloc(n*n*sizeof(double));
  double *u = (double *)malloc(n*n*sizeof(double));
  jacobi(b, n, s, u, v);
  //printMatrix(s, n);
  //printMatrix(u, n);
  //printMatrix(v, n);
  return 0;
}
