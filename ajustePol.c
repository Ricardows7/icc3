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

void atribui (double ** restrict matrix, double * restrict array, double value, double value1, double value2, double value3, double mult, int init){
//esse nao sei se otimizou tanto :(
	matrix[init][0] = value;
	matrix[init+1][0] = value1;
	matrix[init+2][0] = value2;
	matrix[init+3][0] = value3;

	array[init] += value * mult;
	array[init+1] += value1 * mult;
	array[init+2] += value2 * mult;
	array[init+3] += value3 * mult;	//?????????????/

	return;
}

void seta_matriz (double ** restrict matrix, double value, int init, int col, double tax, double tax2, double tax3){
	matrix[init][col] += value;
	matrix[init+1][col] += value * tax;
	matrix[init+2][col] += value * tax2;
	matrix[init+3][col] += value * tax3;
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

void montaSL_V2(double ** restrict A, double * restrict b, int n, long long int p, double * restrict x, double * restrict y, double * restrict safe, double * restrict imp){
	double aux = 0.0, helper = 0.0;

	for (long long int k = 0; k < p; k++)
		safe[k] = 1;	//inicializa x[k] ^ i com i = 0!

	for (int i = 0; i < n-(n%4); i += 4){
		zera_vetor (b, i);
		zera_matriz (A, i, 0);
		for (long long int k = 0; k < p; ++k){	//logica da primeira coluna e vetor!
			aux = x[k];
			atribui (A, b, safe[k], aux, aux*aux, aux*aux*aux, y[k], i);	//seta os valores das i primeiras linhas da matriz e do vetor
			imp[k] = safe[k];	//vai guardar x[k] ^ i para o calculo da matriz
			safe[k] *= aux*aux*aux*aux;		//ja atualiza as postencias de x[k] para o proximo i! (JA QUE SAO 4 LINHAS POR VEZ NAO SERIA 3 * AUX?
		}
		for (int j = 1; j < n; j++){
			zera_matriz (A, i, j); 
			for (long long int k = 0; k < p; ++k){
				aux = x[k];
				imp[k] *= aux;	//calcula x[k] ^i + j!	
				seta_matriz (A, imp[k], i, j, aux, aux*aux, aux*aux*aux);	//to passando pras proximas linhas!
			}//tem que multiplicar por maix aux aqui por causa do imp!
		}
	}

	for (int i = n-n%4; i < n; ++i){
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
	return;
}

long pivo (double **A, double *b, int ori, int max){
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

void eliminacaoGauss_V2 (double **A, double *b, int n) {
  for (int i = 0; i < n; ++i) {
    // Encontra o pivô máximo na coluna
    int iMax = i;
    for (int k = i + 1; k < n; ++k) {
      if (fabs(A[k][i]) > fabs(A[iMax][i]))
        iMax = k;
    }

    // Troca as linhas, se necessário
    if (iMax != i) {
      pivo(A,b,i,iMax);
    }

    // Eliminação com blocos
    for (int k = i + 1; k < n; k += BLOCK_SIZE) {
      int kEnd = (k + BLOCK_SIZE > n) ? n : k + BLOCK_SIZE;
      
      // Loop unrolling para processar linhas em blocos
      for (int kk = k; kk < kEnd; ++kk) {
        double m = A[kk][i] / A[i][i];
        b[kk] -= b[i] * m;

        // Unrolling em j
        for (int j = i + 1; j + 3 < n; j += 4) {
          A[kk][j] -= A[i][j] * m;
          A[kk][j + 1] -= A[i][j + 1] * m;
          A[kk][j + 2] -= A[i][j + 2] * m;
          A[kk][j + 3] -= A[i][j + 3] * m;
        }

        // Restante (caso n não seja múltiplo de 4)
        for (int j = n - (n % 4); j < n; ++j) {
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
  double **C = (double **) malloc (sizeof(double *)*n);
  for (int i = 0; i < n; ++i)
    C[i] = (double *) malloc(sizeof(double)*n);

  double *b = (double *) malloc(sizeof(double)*n);
  double *d = (double *) malloc(sizeof(double)*n);

  double *alpha = (double *) malloc(sizeof(double)*n); // coeficientes ajuste
  double *beta = (double *) malloc(sizeof(double)*n);

  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START ("SL");
  // (A) Gera SL
  double tSL = timestamp();
  montaSL(A, b, n, p, x, y);
  tSL = timestamp() - tSL;
  LIKWID_MARKER_STOP("SL");

//printf("chegamos!\n");
  LIKWID_MARKER_START("SL2");
  double tSL2 = timestamp();
  montaSL_V2(C, d, n, p, x, y, safe, imp);
  tSL2 = timestamp() - tSL2;
  LIKWID_MARKER_STOP("SL2");

  // (B) Resolve SL
  LIKWID_MARKER_START("EG");
  double tEG = timestamp();
  eliminacaoGauss(A, b, n); 
  retrossubs(A, b, alpha, n); 
  tEG = timestamp() - tEG;
  LIKWID_MARKER_STOP("EG");

 LIKWID_MARKER_START("EG2");
  double tEG2 = timestamp();
  eliminacaoGauss(C,d,n);
  retrossubs(C,d,beta,n);
  tEG2 = timestamp() - tEG2;
  LIKWID_MARKER_STOP("EG2");

  imprime (alpha, y, x, n, N, p, K, tSL, tEG);
printf ("\n\n\n");
  imprime (beta, y, x, n, N, p, K, tSL2, tEG2);

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
  }

  LIKWID_MARKER_CLOSE;
  return 0;
}
