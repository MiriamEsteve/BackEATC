
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define __TIEMPOS__

#include "00lib/distrib.h"
#include "00lib/data.h"
#include "00lib/vtree.h"
#include "eat.h"

void Datos(DFrameNum *df, int semilla)
{
	int row, col;
	double arr[10000];

	int numX = df->nCols-1;
	double alpha = 0.5/numX;
	double lambda = 3.0;

	srand(semilla);

	for(col=0 ; col<numX ; col++)
	{
		RandUniform(10, 50, arr, df->nRows);
		dfNumSetColumn(df, col, arr);
	}
	
	// ESCENARIO COBB DOUGLAS para calcular la 'Y'
	for(row=0 ; row<df->nRows ; row++)
	{
		// Función de 'dios'
		df->data[row][numX]=1;
		for(col=0 ; col<numX ; col++)
			df->data[row][numX] *= pow( df->data[row][col] , alpha);

		// Ineficiencia
		df->data[row][numX] *=  exp( log(1-(rand()/(1.0+RAND_MAX)))/lambda );
	}
}



int main(int argc, char *argv[])
{
	time_t t;

	int numNoFinalLeaves;

	int nRows = 15;
	int nCols = 8;
	char *header[]={"X1", "X2", "X3", "X4", "X5", "X6", "X7", "Y1"};
	int X[] = { 0, 1, 2, 3, 4, 5, 6};
	int nX=nCols-1;
	int nY = 1;
	int Y[]={nX};

	int semilla = (int)clock();
	int random = 0;
	int numStop1 = 2;
	int numStop2 = 9;
	int folder = 5;
	double errH, errBest;

	clock_t t_ini, t_fin;


	VLIST *treeAlphaList = listNew(TreeAlphaPrint, TreeAlphaNew, TreeAlphaDel, TreeAlphaCpy);
	VLIST *leaves = listNew(NEATPrintVTN, NEATNewVTN, NEATDelVTN, NEATCpyVTN);
	VTREE *tree = treeNew(NEATPrint, NEATNew, NEATDel, NEATCpy);

	DFrameNum *df = dfNumNew(nRows, nCols, header);
	Datos(df, semilla); // COBB DUGLAS
	double border = 0;
	double noise = 0;

	// 1 - Deep Tree Heuristic
	numNoFinalLeaves=1;
	printf("1- EATStart\n");
	EATStart(df, X, nX, Y, nY, tree, leaves, treeAlphaList);
	
	/*
	printf("2- EATHeuristic\n");
	t_ini = clock();
	EATHeuristic(df, X, nX, Y, nY, numStop1, random, tree, leaves, numNoFinalLeaves, treeAlphaList); // 0-> Heuristico por la izquierda, -1->Heurístico sofisticado
	t_fin = clock();

	errH = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;

	printf("End heuristic....\n=======\n");
	printf("nRows: %d, nCols: %d, numStop: %d, size: %d, error: %f, time : %f\n", nRows, nCols, numStop1, tree->size, errH, (double)(t_fin - t_ini) / CLOCKS_PER_SEC);
	//treePrint(tree);

	//1-2 - Update leaves with numStop2
	numNoFinalLeaves = CheckLeaves(leaves, numStop2);
	/**/

	// 2 - SEATBackTracking
	printf("3- EATBacktracking (numStop=%d)\n", numStop1);
	printf("nRows: %d, nCols: %d, numStop: %d\n", nRows, nCols, numStop2);

	t_ini = clock();
	EATBackTracking(df, X, nX, Y, nY, numStop2, numStop1, random, tree, leaves, numNoFinalLeaves, treeAlphaList);
	t_fin = clock();

	errBest = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;

	printf("End backtraking....\n=======\n");
	printf("nRows: %d, nCols: %d, numStop: %d, size: %d, error: %f, time : %f\n", nRows, nCols, numStop1, tree->size, errBest, (double)(t_fin - t_ini) / CLOCKS_PER_SEC);
	//treePrint(tree);
    
    if(CheckEAT(tree))
        printf("tree OK\n");
    else
        printf("tree KK\n");

	dfNumFree(df);
	listFree(treeAlphaList);
	listFree(leaves);
	treeFree(tree);

	getchar();
	return 0;
}

