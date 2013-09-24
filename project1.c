#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int LEFT[4] = {1, 2, 3, 4};
int RIGHT[4] = {4, 3, 2, 1};
double BAR = 1.0E-13;

int FREE_FIRST_MATRIX = 1;
int FREE_SECOND_MATRIX = 2;


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
      if (a[i*n + j] > BAR)
	return 1;
    }
  }
  return 0;
}

void initializeUV(double *u, double *v, int n, int *mi, int *mj) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i == j) {
	u[i * n + j] = 1;
	v[i * n + j] = 1;
      }
      else {
	u[i * n + j] = 0;
	v[i * n + j] = 0;
      }
    }
  }
  
  u[i * n + i] = LEFT[0];
  u[i * n + j] = LEFT[1];
  u[j * n + i] = LEFT[2];
  u[j * n + j] = LEFT[3];

  v[i * n + i] = RIGHT[0];
  v[i * n + j] = RIGHT[1];
  v[j * n + i] = RIGHT[2];
  v[j * n + j] = RIGHT[3];
  return;
  
}

//flag determines which input(a or b) to free. T to free a, F to free b. 
//the return will malloc new mem for either a or b
double* multiply(double *a, double *b, int n, int flag) {
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
  if (flag == 1)
    free(a);
  else
    free(b);
  return result;
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
void sortMatrix(double *a, int n) {
  int i, j;
  for (i = 0; i < n; ++i) {
    for (j = i; j < n - 1; ++j) {
      if (a[j*n + j] > a[(j+1)*n + (j+1)])
	swapColumn(a, n, j, j+1);
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
  memcpy(&s, &a, sizeof(a));
  
  int count = 0;

  while (checkZeros(s, n) != 0) {
    count++;
    printf("Processing %dth round...\n", count);
    findMax(s, n, max, mi, mj);
    double *ut = (double*)malloc(n*n*sizeof(double));
    double *vt = (double*)malloc(n*n*sizeof(double));
    initializeUV(ut, vt, n, mi, mj);

    s = multiply(ut, s, n, FREE_SECOND_MATRIX);
    s = multiply(s, vt, n, FREE_FIRST_MATRIX);

    u = multiply(ut, u, n, FREE_SECOND_MATRIX);
    v = multiply(v, vt, n, FREE_FIRST_MATRIX);
  }

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
  int len2 = 10;
  double *b = (double*)malloc(len2 * len2 * sizeof(double));
  for (i = 0; i < len2; ++i) {
    for (j = 0; j < len2; ++j) {
      b[i*len2 + j] = i * i + j * j;
    }
  }
  //for (k = 0; k < 3; ++k) 
  //  printMatrix(a[k], len1[k]);
  printMatrix(b, len2);
  int n = len2;
  double *s = (double *)malloc(n*n*sizeof(double));
  double *v = (double *)malloc(n*n*sizeof(double));
  double *u = (double *)malloc(n*n*sizeof(double));
  jacobi(b, n, s, u, v);
  printMatrix(s, n);
  printMatrix(u, n);
  printMatrix(v, n);
  return 0;
}
