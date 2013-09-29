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
    printf("%20.15f ", a[i]);
    if ((i+1) % n == 0)
      printf("\n");
  }
  return;
  
}



void findMax(double *a, int n, double *max, int *mi, int *mj) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i != j && fabs(a[i*n + j]) > *max) {
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
void initializeUtVt(double *s, double *u, double *v, int n, int i, int j) {
  assignIdentity(u, n);
  assignIdentity(v, n);


  double kk;
  double kl;
  double ll;
  if (i < j) {
    kk = s[i*n+i];
    kl = s[i*n+j];
    ll = s[j*n+j];
  }
  else {
    kk = s[j*n+j];
    kl = s[j*n+i];
    ll = s[i*n+i];
  }
  double beta = (ll - kk)/(2.0 * kl);
  double t = sgn(beta) / (fabs(beta) + sqrt(beta*beta + 1));
  
  double cosa = 1 / sqrt(t*t+1);
  double sina = cosa * t;
  printf("cos=%f, sin=%f\n", cosa, sina);
  u[i*n+i] = cosa;
  u[i*n+j] = sina;
  u[j*n+i] = -sina;
  u[j*n+j] = cosa;

  v[i*n+i] = cosa;
  v[i*n+j] = -sina;
  v[j*n+i] = sina;
  v[j*n+j] = cosa;

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
	sum += a[i*n+k] * b[k*n+j];
      result[i*n+j] = sum;
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
void sortMatrix(double *s, double *u, double *v, int n) {
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
  //intialize U, V, S
  assignIdentity(u, n);
  memcpy(s, a, n*n*sizeof(double));
  assignIdentity(v, n);
  
  int count = 0;

  while (checkZeros(s, n) != 0) {
    count++;
    printf("Processing %dth round...\n", count);
    double max = 0;
    int mi = -1;
    int mj = -1;
    findMax(s, n, &max, &mi, &mj);
    printf("mi=%d, mj=%d\n", mi, mj);
    double *ut = (double*)malloc(n*n*sizeof(double));
    double *vt = (double*)malloc(n*n*sizeof(double));
    initializeUtVt(s, ut, vt, n, mi, mj);
    
    //printf("Matrix UT VT\n");
    //printMatrix(ut, n);
    //printMatrix(vt, n);

    multiply(ut, s, n, UPDATE_SECOND_MATRIX);
    multiply(s, vt, n, UPDATE_FIRST_MATRIX);

    multiply(ut, u, n, UPDATE_SECOND_MATRIX);
    multiply(v, vt, n, UPDATE_FIRST_MATRIX);
    
    printf("Matrix U, V, S\n");
    printMatrix(u, n);
    printMatrix(v, n);
    printMatrix(s, n);
    getchar();
  }


  transpose(u, n);
  massage(u, s, n);
  transpose(v, n);
  sortMatrix(s, u, v, n);
  transpose(v, n);


  printf("Final U, S, V: \n");
  printMatrix(u, n);
  printMatrix(s, n);
  printMatrix(v, n);

  return;
  
  
}

void problemB(int n) {
  double *a = (double*)malloc(n*n*sizeof(double));
  int i, j;
  for (i = 0; i < n; ++i) 
    for (j = 0; j < n; ++j) 
      a[i*n+j] = sqrt((i+1)*(i+1) + (j+1)*(j+1));
  
  printMatrix(a, n);

  double *s = (double *)malloc(n*n*sizeof(double));
  double *v = (double *)malloc(n*n*sizeof(double));
  double *u = (double *)malloc(n*n*sizeof(double));
  jacobi(a, n, s, u, v);
  return;
  
}

void problemC(int n) {
  double *a = (double*)malloc(n*n*sizeof(double));
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      a[i*n+j] = (i+1)*(i+1) + (j+1)*(j+1);
    }
  }
  printMatrix(a, n);


  double *s = (double *)malloc(n*n*sizeof(double));
  double *v = (double *)malloc(n*n*sizeof(double));
  double *u = (double *)malloc(n*n*sizeof(double));
  jacobi(a, n, s, u, v);
  return;
}

int main() {
  problemC(3);
  return 0;
}
