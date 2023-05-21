#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

int MAX = 500;
int TAMconjuntoN = 0;

struct grafoD{
	int **MatrizDeAdjacencias;
	int *conjuntoN;
	int *vetorD;
	int *vetorS;
	int ordemDaMatriz;
	int X;
	int Y;
};

typedef struct grafoD GrafoD;

struct grafoBF{
	int **MatrizDeAdjacencias;
	int *vetorD;
	int *vetorS;
	int X;
	int I;
	int ordemDaMatriz;
};

typedef struct grafoBF GrafoBF;

struct grafoF{
	int **MatrizDeAdjacencias;
	int ordemDaMatriz;
};

typedef struct grafoF GrafoF;

static int **CriarMatriz(int tam);
static void AlgoritmoDijkstra(GrafoD *Dijkstra);
static void InserirNoConjuntoN(int *conjuntoN, int v, int pos);
static bool naoEstaNoConjuntoN(int v, int *conjuntoN, int tam);
static int PegarVerticeDeDistanciaMinima(int X, int *vetorD, int *conjuntoN, int tam);
static int MenorDistancia(int v, int w);
static void AlgoritmoBellmanFord(GrafoBF *BellmanFord);
static void AlgoritmoFloyd(GrafoF *Floyd);

void main(){
	int i = 0,j, k, tam, verticeX, verticeY;
	FILE *arquivo;
	int matriz[MAX];


	arquivo = fopen("matriz_de_adjacencias_caminho_minimo.txt", "r");

	while(!feof(arquivo) && i <= MAX){
		fscanf(arquivo, "%d", &matriz[i]);
		i++;
	}

	MAX = i;
	fclose(arquivo);
	tam = sqrt(i);

	printf("O grafo tem %d vértices. De 0 a %d\nInsira os vértices X e Y para calcular o caminho:\n",tam, tam-1);
	printf("\nInsira o vértice X:\n");
	scanf("%d",&verticeX);

	if(verticeX < 0 || verticeX > (tam - 1))
		printf("\nvertice inserido não existe\n");
	else{
		printf("Insira o vértice Y:\n");
		scanf("%d", &verticeY);
		if(verticeY < 0 || verticeY > (tam - 1))
			printf("\nvertice inserido não existe\n");
		else{
			printf("Vértices selecionados, X: %d, Y: %d\n", verticeX, verticeY);

			GrafoD Dijkstra;
			Dijkstra.ordemDaMatriz = tam;
			Dijkstra.MatrizDeAdjacencias = CriarMatriz(tam);
			Dijkstra.X = verticeX;
			Dijkstra.Y = verticeY;

			GrafoBF BellmanFord;
			BellmanFord.MatrizDeAdjacencias = CriarMatriz(tam);
			BellmanFord.X = verticeX;
			BellmanFord.ordemDaMatriz = tam;

			GrafoF Floyd;
			Floyd.MatrizDeAdjacencias = CriarMatriz(tam);
			Floyd.ordemDaMatriz = tam;


			k = 0;
			for(i = 0; i< tam; ++i){
				for(j = 0; j< tam; ++j){
					if(i == j){
						BellmanFord.MatrizDeAdjacencias[i][j] = 0;
						Floyd.MatrizDeAdjacencias[i][j] = 0;
					}else{
						BellmanFord.MatrizDeAdjacencias[i][j] = matriz[k];
						Floyd.MatrizDeAdjacencias[i][j] = matriz[k];
					}

					Dijkstra.MatrizDeAdjacencias[i][j] = matriz[k];
					k++;

				}
			}

			AlgoritmoDijkstra(&Dijkstra);
			AlgoritmoBellmanFord(&BellmanFord);
			AlgoritmoFloyd(&Floyd);

		}
	}

	printf("\nFim da execução do programa\n");

}

static int **CriarMatriz(int tam){
	int **m = malloc(tam * sizeof(int*));
	
	for(int i = 0;i<tam; ++i)
		m[i] = malloc(tam * sizeof(int));

	return m;
}

static void AlgoritmoDijkstra(GrafoD *Dijkstra){
	int i, pos = 0, DistanciaAnterior;
	printf("\nAlgoritmo de Dijkstra\n");
	Dijkstra->vetorD = calloc(Dijkstra->ordemDaMatriz, sizeof(int));
	Dijkstra->vetorS = calloc(Dijkstra->ordemDaMatriz,sizeof(int));
	Dijkstra->conjuntoN = calloc(Dijkstra->ordemDaMatriz,sizeof(int));

	InserirNoConjuntoN(Dijkstra->conjuntoN, Dijkstra->X, pos);
	TAMconjuntoN = 1;
	pos++;

	Dijkstra->vetorD[Dijkstra->X] = 0;

	int Z, P;
	for(Z = 0; Z < Dijkstra->ordemDaMatriz; ++Z){
		if(naoEstaNoConjuntoN(Z, Dijkstra->conjuntoN, Dijkstra->ordemDaMatriz)){
			Dijkstra->vetorD[Z] = Dijkstra->MatrizDeAdjacencias[Dijkstra->X][Z];
			Dijkstra->vetorS[Z] = Dijkstra->X;
		}else{
			Dijkstra->vetorD[Z] = 0;
		}
	}

	while(naoEstaNoConjuntoN(Dijkstra->Y, Dijkstra->conjuntoN, Dijkstra->ordemDaMatriz)){
		P = PegarVerticeDeDistanciaMinima(Dijkstra->X, Dijkstra->vetorD, Dijkstra->conjuntoN, Dijkstra->ordemDaMatriz);
		InserirNoConjuntoN(Dijkstra->conjuntoN, P, pos);
		pos++;
		for( Z = 0; Z < Dijkstra->ordemDaMatriz; ++Z){
			if(naoEstaNoConjuntoN(Z, Dijkstra->conjuntoN, Dijkstra->ordemDaMatriz)){
				DistanciaAnterior = Dijkstra->vetorD[Z];
				Dijkstra->vetorD[Z] = MenorDistancia(Dijkstra->vetorD[Z], (Dijkstra->vetorD[P] + Dijkstra->MatrizDeAdjacencias[P][Z]));
				if(Dijkstra->vetorD[Z] != DistanciaAnterior){
					Dijkstra->vetorS[Z] = P;
				}
			}
		}
	}

	printf("\n    Em ordem inversa o caminho mínimo é: ");
	printf("%d ", Dijkstra->Y);
	Z = Dijkstra->Y;
	do{
		printf("%d ", Dijkstra->vetorS[Z]);
		Z = Dijkstra->vetorS[Z];
	}while(Z != Dijkstra->X);

	printf("\n    O peso do caminho mínimo é: %d\n", Dijkstra->vetorD[Dijkstra->Y]);

} 

static void InserirNoConjuntoN(int *conjuntoN, int v, int pos){
	conjuntoN[pos] = v;
	TAMconjuntoN++;	
}

static bool naoEstaNoConjuntoN(int v, int *conjuntoN, int tam){
	for(int i = 0; i < TAMconjuntoN; ++i){
		if(conjuntoN[i] == v){
			return false;
		}
	}

	return true;
}

static int PegarVerticeDeDistanciaMinima(int X, int *vetorD, int *conjuntoN, int tam){
	int i, minimo = 1000, j = 0, v;
	for(i = 0; i < tam; ++i){
		if(i != X && naoEstaNoConjuntoN(i, conjuntoN, tam)){
			if(vetorD[i] <= minimo){
				minimo = vetorD[i];
				v = i;
			}
		}
	}

	return v;
}

static int MenorDistancia(int v, int w){
	if(v >= 1000  && w >= 1000){
		return 1000;
	}else{
		if(v < w)
			return v;
		else
			return w;
	}
}

static void AlgoritmoBellmanFord(GrafoBF *BellmanFord){
	printf("\nAlgoritmo de Bellman-Ford\n");
	int k,p,pos,z, i,minimo, soma;
	BellmanFord->vetorD = calloc(BellmanFord->ordemDaMatriz, sizeof(int));
	BellmanFord->vetorS = calloc(BellmanFord->ordemDaMatriz, sizeof(int));

	BellmanFord->vetorD[BellmanFord->X] = 0;
	int t[BellmanFord->ordemDaMatriz];

	for(z =0; z < BellmanFord->ordemDaMatriz; ++z){
		if(z != BellmanFord->X){
			BellmanFord->vetorD[z] = BellmanFord->MatrizDeAdjacencias[BellmanFord->X][z];
			BellmanFord->vetorS[z] = BellmanFord->X;
		}
	}

	for(i = 2; i < BellmanFord->ordemDaMatriz-1;++i){
		
		for(k = 0; k < BellmanFord->ordemDaMatriz; ++k){
			t[k] = BellmanFord->vetorD[k];
		}

		for(z=0;z<BellmanFord->ordemDaMatriz;++z){
			minimo = 1000;

			if(z != BellmanFord->X){
				p = 0;
				while(p <= BellmanFord->ordemDaMatriz - 1){
					if(p != BellmanFord->X){
						soma = BellmanFord->vetorD[p]+BellmanFord->MatrizDeAdjacencias[p][z];
						
						if(soma < minimo){
							minimo = soma;
							t[z] = minimo;
							pos = p;
						}

					}
					p++;

				}

				if( pos != z)
					BellmanFord->vetorS[z] = pos;
				
			}


		}

		for(k = 0;k < BellmanFord->ordemDaMatriz; ++k){
			BellmanFord->vetorD[k] = t[k];
		}
	}
	int y;

	for( y = 0; y < BellmanFord->ordemDaMatriz; ++y){
		printf("    O peso do caminho minimo de %d até o vértice %d é: %d\n",BellmanFord->X, y, BellmanFord->vetorD[y]);
		printf("    O caminho mínimo em ordem inversa é: %d ", y);
		i = y;
		while(i != BellmanFord->X){
			printf("%d ", BellmanFord->vetorS[i]);
			i = BellmanFord->vetorS[i];
		};
		printf("\n\n");
		
	}

}

static void AlgoritmoFloyd(GrafoF *Floyd){
	printf("\nAlgoritmo de Floyd\n");
	int i,j,k;
	int n = Floyd->ordemDaMatriz;
	int A[n][n];
	
	for(i =0;i<n;++i)
		for(j=0;j<n;++j)
			A[i][j] = Floyd->MatrizDeAdjacencias[i][j];

		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				for(k = 0; k < n; k++){
					if(A[i][k] + A[k][j] < A[i][j]){
						A[i][j] = A[i][k] + A[k][j];
					}
				}
			}		
		}

		for(i = 0; i < n;++i){
			for(j = 0; j < n;++j){
				printf("%4d",A[i][j]);
			}
			printf("\n");
		}
}