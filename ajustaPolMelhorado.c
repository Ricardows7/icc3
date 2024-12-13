#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <fenv.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <likwid.h>

#include "utils.h"

#define BLOCK_SIZE 64

/////////////////////////////////////////////////////////////////////////////////////
//   AJUSTE DE CURVAS
/////////////////////////////////////////////////////////////////////////////////////

void zera_vetor (double *restrict array, int init){
	array[init] = 0.0;
	array[init + 1] = 0.0;
	array[init + 2] = 0.0;
	array[init + 3] = 0.0;

	return;
}

double atribui_matriz (double ** restrict matrix, int init, double *holster){
	matrix[init][0] = holster[0];
	matrix[init+1][0] = holster[1];
	matrix[init+2][0] = holster[2];
	matrix[init+3][0] = holster[3];

	return holster[3];
}

double atribui_vetor (double * restrict array, double mult1, double mult2, int init){
	double a = pow(mult1, init);
	double b = pow(mult1, init + 1);
	double c = pow(mult1, init + 2);
	double d = pow(mult1, init + 3);

	//printf ("Na linha %d os values sao: %f %f %f %f\n", init, holster[0], holster[1], holster[2], holster[3]);

	array[init] += a * mult2;
	array[init+1] += b * mult2;
	array[init+2] += c * mult2;
	array[init+3] += d * mult2;	//?????????????

	return d;
}

void seta_matriz (double ** restrict matrix, double value, int init, int col, double tax){
	matrix[init][col] += value;
	matrix[init+1][col] += value * tax;
	matrix[init+2][col] += value * tax * tax;
	matrix[init+3][col] += value * tax * tax * tax;
	return;
}

void zera_matriz (double ** restrict matriz, int init, int col){
	matriz[init][col] = 0.0;
	matriz[init+1][col] = 0.0;
	matriz[init+2][col] = 0.0;
	matriz[init+3][col] = 0.0;
	return;
}

void montaSL(double **A, double *b, int n, long long int p, double *x, double *y) {
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      A[i][j] = 0.0;
      for (long long int k = 0; k < p; ++k) {
	A[i][j] += pow(x[k], i+j);
      }
    }

  for (int i = 0; i < n; ++i) {
    b[i] = 0.0;
    for (long long int k = 0; k < p; ++k)
      b[i] += pow(x[k],i) * y[k];
  }
}

void printaMatriz (double **A, double *b, int n){
	printf("COMECANDO A PRINTAR O VETOR!:"); 
	for (int i = 0; i < n; i++)
		printf ("%f ", b[i]);
	printf ("\n\n");
	printf ("COMECANDO A PRINTAR A MATRIZ!:\n");
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++)
			printf ("%f ", A[i][j]);
		printf ("\n");
	}
	printf ("\n\n");
	return;
}

void montaSL_V2(double ** restrict A, double * restrict b, int n, long long int p, double * restrict x, double * restrict y, double * restrict safe, double * restrict imp){
	double aux = 0.0, helper = 0.0;

	for (long long int k = 0; k < p; k++){
		safe[k] = 1;	//inicializa x[k] ^ i com i = 0!
	}

	for (int i = 0; i < n-(n%4); i += 4){
		zera_vetor (b, i);		//zera as primeiras 4 posições
		zera_matriz (A, i, 0);
		for (long long int k = 0; k < p; ++k){	//logica da primeira coluna e vetor!
			aux = x[k];
			imp[k] = safe[k];		//vai guardar x[k] ^ i para o calculo da matriz   
			safe[k] = atribui_vetor (b, aux, y[k], i) * aux;	//seta os valores das i primeiras linhas da matriz e do vetor
		}
		
		//printf ("linha atual %d!\n\n", i);
		for (int j = 1; j < n; j++){
			zera_matriz (A, i, j); 
			for (long long int k = 0; k < p; ++k){
				aux = x[k];
				imp[k] *= aux;	//calcula x[k] ^i + j!	
				seta_matriz (A, imp[k], i, j, aux);	//to passando pras proximas linhas!
			}//tem que multiplicar por maix aux aqui por causa do imp!
		}
	}

	for (int i = n-(n%4); i < n; ++i){
		b[i] = 0.0;
		for (int j = 0; j < n; j++)
			A[i][j] = 0.0;
		for (long long int k = 0; k < p; k++){
			aux = x[k];
			helper = pow (aux, i);
			b[i] += helper * y[k];
			for (int j = 0; j < n; j++){
				A[i][j] += helper;
				helper *= aux;
			}
		}
	}	//ACHO QUE O PROBLEMA ESTA NO VETOR!!!
	for (int i = 0; i < n-(n%4); i++){
		A[i][0] = 0.0;
		for (long long int k = 0; k < p; k++){
			aux = x[k];
			helper = pow(aux, i);
			A[i][0] += helper;
		}
	}
	return;
}

double pivo (double **A, double *b, int ori, int max){
	double *tmp, aux;

	tmp = A[ori];
	A[ori] = A[max];
	A[max] = tmp;

	aux = b[ori];
	b[ori] = b[max];
	b[max] = aux;
	
	return b[ori];
}

void eliminacaoGauss(double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    int iMax = i;
    for (int k = i+1; k < n; ++k)
      if (A[k][i] > A[iMax][i])
	iMax = k;
    if (iMax != i) {
      double *tmp, aux;
      tmp = A[i];
      A[i] = A[iMax];
      A[iMax] = tmp;

      aux = b[i];
      b[i] = b[iMax];
      b[iMax] = aux;
    }

    for (int k = i+1; k < n; ++k) {
      double m = A[k][i] / A[i][i];
      A[k][i]  = 0.0;	//!
      for (int j = i+1; j < n; ++j)
	A[k][j] -= A[i][j]*m;
      b[k] -= b[i]*m;
    }
  }
}

void eliminacaoGauss_V2 (double **A, double *b, int n){
  for (int i = 0; i < n; ++i) {
    // Encontra o pivô máximo na coluna
    int iMax = i;
    double maxVal = fabs(A[i][i]);
    for (int k = i + 1; k < n; ++k) {
      double val = fabs(A[k][i]);
      iMax = (val > maxVal) ? k : iMax;
      maxVal = (val > maxVal) ? val : maxVal;
    }

    // Troca as linhas
    pivo(A, b, i, iMax);

    // Calcula o inverso do pivô para evitar divisão dentro dos loops
    double invPivo = 1.0 / A[i][i];

    // Eliminação por blocos
    for (int k = i + 1; k < n; k += BLOCK_SIZE) {
      int kEnd = k + BLOCK_SIZE;
      if (kEnd > n) kEnd = n;

      for (int kk = k; kk < kEnd; ++kk) {
        double m = A[kk][i] * invPivo; // Multiplicação mais eficiente que divisão
        b[kk] -= b[i] * m;

        // Processa em blocos de 8 elementos com loop unrolling
        int j = i + 1;
        for (; j + 7 < n; j += 8) {
          A[kk][j] -= A[i][j] * m;
          A[kk][j + 1] -= A[i][j + 1] * m;
          A[kk][j + 2] -= A[i][j + 2] * m;
          A[kk][j + 3] -= A[i][j + 3] * m;
          A[kk][j + 4] -= A[i][j + 4] * m;
          A[kk][j + 5] -= A[i][j + 5] * m;
          A[kk][j + 6] -= A[i][j + 6] * m;
          A[kk][j + 7] -= A[i][j + 7] * m;
        }

        // Processa os elementos restantes
        for (; j < n; ++j) {
          A[kk][j] -= A[i][j] * m;
        }
      }
    }
  }
}


void retrossubs(double **A, double *b, double *x, int n) {
  for (int i = n-1; i >= 0; --i) {
    x[i] = b[i];
    for (int j = i+1; j < n; ++j)
      x[i] -= A[i][j]*x[j];
    x[i] /= A[i][i];
  }
}

double P(double x, int N, double *alpha) {
  double Px = alpha[0];
  for (int i = 1; i <= N; ++i)
    Px += alpha[i]*pow(x,i);
  
  return Px;
}

void imprime (double* alpha, double* y, double* x, int n, int N, long long int p, long long int K, double tSL, double tEG){
  // Imprime coeficientes
  for (int i = 0; i < n; ++i)
    printf("%1.15e ", alpha[i]);
  puts("");

  // Imprime resíduos
  for (long long int i = 0; i < p; ++i)
    printf("%1.15e ", fabs(y[i] - P(x[i],N,alpha)) );
  puts("");

  // Imprime os tempos
  printf("%lld %1.10e %1.10e\n", K, tSL, tEG);
}

void resultado (double** A, double *B, int n){
  for (int i = 0; i < n; i++){
          for (int j = 0; j < n; j++){
                  printf ("%lf ", A[i][j]);
          }
          printf ("| %lf\n", B[i]);
  }
  return;
}

int main() {

  int N, n;
  long long int K, p;

  scanf("%d %lld", &N, &K);
  p = K;   // quantidade de pontos
  n = N+1; // tamanho do SL (grau N + 1)

  double *safe = (double *) malloc (p * sizeof(double));
  double *imp = (double *) malloc (p * sizeof(double));

  double *x = (double *) malloc(sizeof(double)*p);
  double *y = (double *) malloc(sizeof(double)*p);

  // ler numeros
  for (long long int i = 0; i < p; ++i)
    scanf("%lf %lf", x+i, y+i);

  double **A = (double **) malloc(sizeof(double *)*n);
  for (int i = 0; i < n; ++i)
    A[i] = (double *) malloc(sizeof(double)*n);

  double *b = (double *) malloc(sizeof(double)*n);

  double *alpha = (double *) malloc(sizeof(double)*n); // coeficientes ajuste

  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START ("SL");
  // (A) Gera SL
  double tSL = timestamp();
  montaSL_V2(A, b, n, p, x, y, safe, imp);
  tSL = timestamp() - tSL;
  LIKWID_MARKER_STOP("SL");

//printf("chegamos!\n");
/*
  LIKWID_MARKER_START("SL2");
  double tSL2 = timestamp();
  montaSL_V2(C, d, n, p, x, y, safe, imp);
  tSL2 = timestamp() - tSL2;
  LIKWID_MARKER_STOP("SL2");
*/
  //printaMatriz (A, b, n);
  //printaMatriz (C,d,n);
  // (B) Resolve SL

  LIKWID_MARKER_START("EG");
  double tEG = timestamp();
  eliminacaoGauss_V2(A, b, n); 
  retrossubs(A, b, alpha, n); 
  tEG = timestamp() - tEG;
  LIKWID_MARKER_STOP("EG");

/*
 LIKWID_MARKER_START("EG2");
  double tEG2 = timestamp();
  eliminacaoGauss_V2(C,d,n);
  retrossubs(C,d,beta,n);
  tEG2 = timestamp() - tEG2;
  LIKWID_MARKER_STOP("EG2");

  bool certo = true;
 int coor_x = 0, coor_y = 0;
  for (int i = 0; i < n; i++){
	  for (int j = 0; j < n; j++)
		  if (A[i][j] != C[i][j]){
			  certo = false;
			  coor_x = i;
			  coor_y = j;
          }
	  if (b[i] != d[i])
		  certo = false;
  }
  if (certo)
	  printf ("DEU BOAAA!\n");
  else{
	  printf ("NAO deu boa :(\n");
		printf ("X: %d / Y: %d\n", coor_x, coor_y);
  }*/

  printf("%lld %1.10e %1.10e\n", K, tSL, tEG);  

  for (int i = 0; i < n; i++){
	free(A[i]);
  }
	
	free(A);
	free(b);
	free(alpha);
	free (safe);
	free (imp);
	free (x);
	free (y);

  LIKWID_MARKER_CLOSE;
  return 0;
}
