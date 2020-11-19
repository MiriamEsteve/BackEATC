#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include "eat.h"

#pragma warning(disable: 4996)

double MiROUND(double val, int div)
{
	return round(val * div) / div;
}

/****************************************************************************************
 ###################### Nodo para el arbol EAT #####################
******************************************************************************************/
// *********** Mensajes de error EAT ********************
void EATErr(char* err)
{
	fprintf(stderr, "  ## EAT-ERROR >> %s\n", err);
	getchar();
	exit(0);
}

void  NEATPrint(void* d)
{
	int i;
	NEAT* n = (NEAT*)d;

	if (!d)
	{
		printf("[NEAT NULL]\n");
		return;
	}

	printf("[ME:%d (F:%d), xi: %d, s: %f, R: %f(%f), y->[%.4f ", n->id, n->idF, n->xi, n->s, n->R, n->errMin, n->y[0]);
	for (i = 1; i < n->nY; i++)
		printf(",%.4f", n->y[i]);

	printf("], [index:%d -> (%d", n->nIndex, n->index[0]);
	if (n->nIndex < 11)
	{
		for (i = 1; i < n->nIndex; i++)
			printf(",%d", n->index[i]);
	}
	else
	{
		for (i = 1; i < 10; i++)
			printf(",%d", n->index[i]);
		printf(", ...,%d", n->index[(n->nIndex - 1)]);
	}
	printf(")]");

	/*******************************************************************
	// a
	printf("		   a:[%.5f", n->a[0]);
	for(i=1 ; i<n->nX ; i++)
	   printf(", %.5f", n->a[i]);
	printf("]\n");

	// b
	printf("		   b:[%.5f", n->b[0]);
	for(i=1 ; i<n->nX ; i++)
	   printf(", %.5f", n->b[i]);
	printf("]\n");

	// varInfo
	printf("		   vInfo:[[SL:%.4f, SR:%.4f, S:%.4f]", n->varInfo[0].SL, n->varInfo[0].SR, n->varInfo[0].S);
	for(i=1 ; i<n->nX ; i++)
		printf(" [SL:%.4f, SR:%.4f, S:%.4f]", n->varInfo[i].SL, n->varInfo[i].SR, n->varInfo[i].S);
	/*******************************************************************/

	printf("]\n");

}

void* NEATNew(void* d)
{
	char aux[200];
	int i;
	NEAT* newD;

	if (!(newD = (NEAT*)malloc(sizeof(NEAT))))
		return NULL;

	((NEAT*)newD)->nIndex = ((NEAT*)d)->nIndex;
	if (!(((NEAT*)newD)->index = (int*)malloc(sizeof(int) * ((NEAT*)newD)->nIndex)))
		EATErr("NEATNew [memory fault-1]");
	for (i = 0; i < ((NEAT*)newD)->nIndex; i++)
		((NEAT*)newD)->index[i] = ((NEAT*)d)->index[i];

	newD->nX = ((NEAT*)d)->nX;
	newD->nY = ((NEAT*)d)->nY;

	((NEAT*)newD)->varInfo = (VARINFO*)malloc(sizeof(VARINFO) * newD->nX);
	((NEAT*)newD)->y = (double*)malloc(sizeof(double) * newD->nY);
	if (!((NEAT*)newD)->varInfo || !((NEAT*)newD)->y)
	{
		sprintf(aux, "NEATNew [memory fault-2] varInfo=%p, y=%p [nX=%d]", ((NEAT*)newD)->varInfo, ((NEAT*)newD)->y, newD->nX);
		EATErr(aux);
	}

	for (i = 0; i < newD->nX; i++)
	{
		((NEAT*)newD)->varInfo[i].SL = ((NEAT*)d)->varInfo[i].SL;
		((NEAT*)newD)->varInfo[i].SR = ((NEAT*)d)->varInfo[i].SR;
		((NEAT*)newD)->varInfo[i].S = ((NEAT*)d)->varInfo[i].S;
	}

	((NEAT*)newD)->id = ((NEAT*)d)->id;
	((NEAT*)newD)->idF = ((NEAT*)d)->idF;
	((NEAT*)newD)->R = ((NEAT*)d)->R;
	((NEAT*)newD)->errMin = ((NEAT*)d)->errMin;
	((NEAT*)newD)->xi = ((NEAT*)d)->xi;
	((NEAT*)newD)->s = ((NEAT*)d)->s;
	((NEAT*)newD)->sMax = ((NEAT*)d)->sMax;
	for (i = 0; i < newD->nY; i++)
		((NEAT*)newD)->y[i] = ((NEAT*)d)->y[i];

	((NEAT*)newD)->a = (double*)malloc(sizeof(double) * newD->nX);
	((NEAT*)newD)->b = (double*)malloc(sizeof(double) * newD->nX);
	if (!((NEAT*)newD)->a || !((NEAT*)newD)->b)
		EATErr("NEATNew [memory fault-3]");

	for (i = 0; i < newD->nX; i++)
	{
		((NEAT*)newD)->a[i] = ((NEAT*)d)->a[i];
		((NEAT*)newD)->b[i] = ((NEAT*)d)->b[i];
	}

	return (void*)newD;
}

void  NEATDel(void* d)
{
	if (d)
	{
		free(((NEAT*)d)->index);
		free(((NEAT*)d)->varInfo);
		free(((NEAT*)d)->y);
		free(((NEAT*)d)->a);
		free(((NEAT*)d)->b);
		free(((NEAT*)d));
	}
}

void  NEATCpy(void* d1, void* d2)
{
	int i;

	((NEAT*)d1)->nIndex = ((NEAT*)d2)->nIndex;
	if (((NEAT*)d1)->nIndex == 0)
		((NEAT*)d1)->index = NULL;
	else
	{
		if (!(((NEAT*)d1)->index = (int*)malloc(sizeof(int) * ((NEAT*)d1)->nIndex)))
			EATErr("NEATCpy [memory fault-1]");
		for (i = 0; i < ((NEAT*)d1)->nIndex; i++)
			((NEAT*)d1)->index[i] = ((NEAT*)d2)->index[i];
	}

	((NEAT*)d1)->nX = ((NEAT*)d2)->nX;
	((NEAT*)d1)->nY = ((NEAT*)d2)->nY;
	((NEAT*)d1)->varInfo = (VARINFO*)malloc(sizeof(VARINFO) * ((NEAT*)d1)->nX);
	((NEAT*)d1)->y = (double*)malloc(sizeof(double) * ((NEAT*)d1)->nY);
	if (!((NEAT*)d1)->varInfo || !((NEAT*)d1)->y)
		EATErr("NEATCpy [memory fault-2]");

	for (i = 0; i < ((NEAT*)d1)->nX; i++)
	{
		((NEAT*)d1)->varInfo[i].SL = ((NEAT*)d2)->varInfo[i].SL;
		((NEAT*)d1)->varInfo[i].SR = ((NEAT*)d2)->varInfo[i].SR;
		((NEAT*)d1)->varInfo[i].S = ((NEAT*)d2)->varInfo[i].S;
	}

	((NEAT*)d1)->id = ((NEAT*)d2)->id;
	((NEAT*)d1)->idF = ((NEAT*)d2)->idF;
	((NEAT*)d1)->R = ((NEAT*)d2)->R;
	((NEAT*)d1)->errMin = ((NEAT*)d2)->errMin;
	((NEAT*)d1)->xi = ((NEAT*)d2)->xi;
	((NEAT*)d1)->s = ((NEAT*)d2)->s;
	((NEAT*)d1)->sMax = ((NEAT*)d2)->sMax;
	for (i = 0; i < ((NEAT*)d1)->nY; i++)
		((NEAT*)d1)->y[i] = ((NEAT*)d2)->y[i];

	((NEAT*)d1)->a = (double*)malloc(sizeof(double) * ((NEAT*)d1)->nX);
	((NEAT*)d1)->b = (double*)malloc(sizeof(double) * ((NEAT*)d1)->nX);
	if (!((NEAT*)d1)->a || !((NEAT*)d1)->b)
		EATErr("NEATCpy [memory fault-3]");
	for (i = 0; i < ((NEAT*)d1)->nX; i++)
	{
		((NEAT*)d1)->a[i] = ((NEAT*)d2)->a[i];
		((NEAT*)d1)->b[i] = ((NEAT*)d2)->b[i];
	}
}

// copia info para el hijo (izq. o der.)
void  NEATCpySon(void* d1, void* d2)
{
	int i;

	((NEAT*)d1)->s = ((NEAT*)d1)->xi = -1;

	for (i = 0; i < ((NEAT*)d2)->nX; i++)
	{
		((NEAT*)d1)->a[i] = ((NEAT*)d2)->a[i];
		((NEAT*)d1)->b[i] = ((NEAT*)d2)->b[i];
	}
}

void  NEATSet(void* d, int nX, int nY, int N)
{
	int i;

	((NEAT*)d)->nIndex = N;
	if (N == 0)
		((NEAT*)d)->index = NULL;
	else
	{
		if (!(((NEAT*)d)->index = (int*)malloc(sizeof(int) * N)))
			EATErr("NEATSet [memory fault-1]");
		for (i = 0; i < N; i++)
			((NEAT*)d)->index[i] = i;
	}

	((NEAT*)d)->nX = nX;
	((NEAT*)d)->nY = nY;

	((NEAT*)d)->varInfo = (VARINFO*)malloc(sizeof(VARINFO) * ((NEAT*)d)->nX);
	((NEAT*)d)->y = (double*)malloc(sizeof(double) * ((NEAT*)d)->nY);
	if (!((NEAT*)d)->varInfo || !((NEAT*)d)->y)
		EATErr("NEATSet [memory fault-2]");
	for (i = 0; i < ((NEAT*)d)->nX; i++)
		((NEAT*)d)->varInfo[i].SL = ((NEAT*)d)->varInfo[i].SR = ((NEAT*)d)->varInfo[i].S = 0;

	((NEAT*)d)->R = -1;
	((NEAT*)d)->errMin = -1;
	((NEAT*)d)->xi = -1;
	((NEAT*)d)->s = -1;
	((NEAT*)d)->sMax = -1;
	for (i = 0; i < ((NEAT*)d)->nY; i++)
		((NEAT*)d)->y[i] = -1;

	((NEAT*)d)->a = (double*)malloc(sizeof(double) * ((NEAT*)d)->nX);
	((NEAT*)d)->b = (double*)malloc(sizeof(double) * ((NEAT*)d)->nX);
	if (!((NEAT*)d)->a || !((NEAT*)d)->b)
		EATErr("NEATSet [memory fault-2]");

	for (i = 0; i < ((NEAT*)d)->nX; i++)
	{
		((NEAT*)d)->a[i] = 0;
		((NEAT*)d)->b[i] = INFINITY;
	}
}

//-----------------
void  NEATPrintVTN(void* d) { NEATPrint(((VTN*)d)->data); }
void* NEATNewVTN(void* d) { return d; }
void  NEATDelVTN(void* d) {  }
void  NEATCpyVTN(void* d1, void* d2) { d1 = d2; }


/****************************************************************************************
 ########################## Pruning EAT ##########################
******************************************************************************************/
void  TreeAlphaPrint(void* ta) // 'ta' -> treeAlpha
{
	printf(" [alpha: %.15f, score: %f, errFold: %f, size: %d]\n",
		((TREEALPHA*)ta)->alpha, ((TREEALPHA*)ta)->score, ((TREEALPHA*)ta)->errFolders, ((TREEALPHA*)ta)->tree->size);
	//listPrint(((TREEALPHA*)ta)->tree);
}

void* TreeAlphaNew(void* ta)
{
	TREEALPHA* newTA;

	if (!ta)
		EATErr("TreeAlphaNew [memory fault - 1]");

	if (!(newTA = (TREEALPHA*)malloc(sizeof(TREEALPHA))))
		EATErr("TreeAlphaNew [memory fault - 2]");

	((TREEALPHA*)newTA)->tree = treeNew(NEATPrint, NEATNew, NEATDel, NEATCpy);
	treeCopy(((TREEALPHA*)newTA)->tree, ((TREEALPHA*)ta)->tree);
	((TREEALPHA*)newTA)->alpha = ((TREEALPHA*)ta)->alpha;
	((TREEALPHA*)newTA)->score = ((TREEALPHA*)ta)->score;
	((TREEALPHA*)newTA)->errFolders = ((TREEALPHA*)ta)->errFolders;

	if (newTA->tree == NULL)
		EATErr("TreeAlphaNew [memory fault - 3]");

	return (void*)newTA;
}

void  TreeAlphaDel(void* ta)
{
	if ((TREEALPHA*)ta)
	{
		treeFree(((TREEALPHA*)ta)->tree);
		free((TREEALPHA*)ta);
	}
}

void  TreeAlphaCpy(void* ta1, void* ta2)
{
	((TREEALPHA*)ta1)->alpha = ((TREEALPHA*)ta2)->alpha;
	((TREEALPHA*)ta1)->score = ((TREEALPHA*)ta2)->score;
	((TREEALPHA*)ta1)->errFolders = ((TREEALPHA*)ta2)->errFolders;
	treeCopy(((TREEALPHA*)ta1)->tree, ((TREEALPHA*)ta2)->tree);
}


/****************************************************************************************
 ############################## EAT ############################
******************************************************************************************/
void ArrayPrint(double* array, int nArray)
{
	int i;
	for (i = 0; i < nArray; i++)
		printf("[%3d] -> %.5f\n", i, array[i]);
}

void ArrayPrintInt(int* array, int nArray)
{
	int i;
	for (i = 0; i < nArray; i++)
		printf("[%3d] -> %3d\n", i, array[i]);
}

double MaxDF(DFrameNum* df, int* index, int nIndex, int y)
{
	int i;
	double max;

	if (y < 0 || y >= df->nCols)
		EATErr("MaxDF [y-col error]");

	max = df->data[index[0]][y];
	for (i = 1; i < nIndex; i++)
		if (max < df->data[index[i]][y])
			max = df->data[index[i]][y];
	return max;
}

double QuadErr(DFrameNum* df, int* index, int nIndex, int* Y, double* yEst, int nY)
{
	int iCol, row;
	double R = 0, aux;


	for (iCol = 0; iCol < nY; iCol++)
	{
		if (Y[iCol] < 0 || Y[iCol] >= df->nCols)
			EATErr("SqrErr [y-col error]");

		for (row = 0; row < nIndex; row++)
		{
			aux = df->data[index[row]][Y[iCol]] - yEst[iCol];
			R += aux * aux;
		}
	}
	return MiROUND((R / (nY * df->nRows)), 10000);
}


int ComparePareto(const void* a, const void* b)
{
	int x, cont1 = 0, cont2 = 0;
	NEAT* t1 = (NEAT*)a;
	NEAT* t2 = (NEAT*)b;

	for (x = 0; x < t1->nX; x++)
	{
		if (t1->a[x] < t2->b[x])
			cont1++;
		if (t2->a[x] < t1->b[x])
			cont2++;
	}
	if (cont1 == t1->nX)
		return -1;
	if (cont2 == t1->nX)
		return 1;
	return 0;
}

int CompareDouble(const void* a, const void* b)
{
	double fa = *(const double*)a;
	double fb = *(const double*)b;
	return (fa > fb) - (fa < fb);
}

int IsValueInArray(double val, double* array, int nArray)
{
	int i;
	for (i = 0; i < nArray; i++)
		if (array[i] == val)
			return i;
	return -1;
}

int GetSortArrayFormIndex(double* array, DFrameNum* df, int column, int* index, int nIndex)
{
	int row, nArray = 0;

	array[nArray] = df->data[index[0]][column];
	nArray++;

	for (row = 1; row < nIndex; row++)
		if (IsValueInArray(df->data[index[row]][column], array, nArray) == -1)
		{
			array[nArray] = df->data[index[row]][column];
			nArray++;
		}
	qsort(array, nArray, sizeof(double), CompareDouble);
	return nArray;
}

int* NumValsLeft(DFrameNum* df, int* tIndex, int nIndex, int xi, double s, int* nIndexL)
{
	int i;
	int* indL = NULL;
	(*nIndexL) = 0;
	char aux[200];

	for (i = 0; i < nIndex; i++)
		if (df->data[tIndex[i]][xi] < s)
			(*nIndexL)++;

	if ((*nIndexL) == 0 || !(indL = (int*)malloc(sizeof(int) * (*nIndexL))))
	{
		sprintf(aux, "NumValsLeft [memory fault] (n=%d)", (*nIndexL));
		EATErr(aux);
	}

	(*nIndexL) = 0;
	for (i = 0; i < nIndex; i++)
		if (df->data[tIndex[i]][xi] < s)
		{
			indL[(*nIndexL)] = tIndex[i];
			(*nIndexL)++;
		}

	return indL;
}

int* NumValsRight(DFrameNum* df, int* tIndex, int nIndex, int xi, double s, int* nIndexR)
{
	int i;
	int* indR = NULL;
	(*nIndexR) = 0;

	for (i = 0; i < nIndex; i++)
		if (df->data[tIndex[i]][xi] >= s)
			(*nIndexR)++;

	indR = (int*)malloc(sizeof(int) * (*nIndexR));

	if ((*nIndexR) == 0 || !indR)
		EATErr("NumValsRight [memory fault]");

	(*nIndexR) = 0;
	for (i = 0; i < nIndex; i++)
		if (df->data[tIndex[i]][xi] >= s) {
			indR[(*nIndexR)] = tIndex[i];
			(*nIndexR)++;
		}

	return indR;
}

int IsFinalNode(NEAT* t, DFrameNum* df, int* X, int numStop)
{
	int row, xi;

	if (t->nIndex <= numStop)
		return 1;

	for (row = 1; row < t->nIndex; row++)
	{
		for (xi = 0; xi < t->nX; xi++)
		{
			if (df->data[t->index[0]][X[xi]] != df->data[t->index[row]][X[xi]])
				return 0;
		}
	}
	return 1;
}

void EstimEAT(NEAT* tData, DFrameNum* df, VLIST* leaves, int* Y, int* indexL, int nIndexL, double* estL)
{
	int i;
	double maxL;
	double yInfLeft[100];

	VListIter iter;
	VTN* n;

	if (tData->nY > 100)
		EATErr("EstimEAT [memory fault for 'yInfLeft']");

	for (i = 0; i < tData->nY; i++)
		yInfLeft[i] = 0;

	if (nIndexL == 0)
		for (i = 0; i < tData->nY; i++)
			estL[i] = INFINITY;
	else
	{
		listRefMoveFirst(leaves, &iter);
		while ((n = (VTN*)listRefMoveNext(leaves, &iter)))
		{
			if (ComparePareto(tData, ((NEAT*)n->data)) == 1)
				for (i = 0; i < tData->nY; i++)
					if (yInfLeft[i] < ((NEAT*)n->data)->y[i])
						yInfLeft[i] = ((NEAT*)n->data)->y[i];
		}
		for (i = 0; i < tData->nY; i++)
		{
			maxL = MaxDF(df, indexL, nIndexL, Y[i]);
			estL[i] = (maxL >= yInfLeft[i]) ? maxL : yInfLeft[i];
		}
	}
}

void ErrEAT(NEAT* tData, DFrameNum* df, VLIST* leaves, int* Y, int xCol, double s,
	double* errL, double* errR, double* estL)
{
	double bAux = tData->b[xCol];
	int nIndexL, nIndexR;
	int* indexL = NumValsLeft(df, tData->index, tData->nIndex, xCol, s, &nIndexL);
	int* indexR = NumValsRight(df, tData->index, tData->nIndex, xCol, s, &nIndexR);

	tData->b[xCol] = s; // tData simula su hijo izquierdo 'TL'
	EstimEAT(tData, df, leaves, Y, indexL, nIndexL, estL);
	tData->b[xCol] = bAux; // tData vuelve a ser 'father'

	(*errL) = QuadErr(df, indexL, nIndexL, Y, estL, tData->nY);
	(*errR) = QuadErr(df, indexR, nIndexR, Y, tData->y, tData->nY); // estimR = Father

	free(indexL);
	free(indexR);
}

void Split(NEAT* tData, VLIST* leaves, DFrameNum* df, int* X, int* Y, NEAT* TL, NEAT* TR, double* array)
{
	double err, errMin = INFINITY;
	int i, i2, row, xi, nArray;

	int xiRes = 0;
	double sRes = 0.0;

	double errL, errR;
	double* yMaxLeft, * yEstimL;

	yEstimL = (double*)malloc(sizeof(double) * tData->nY);
	yMaxLeft = (double*)malloc(sizeof(double) * tData->nY);
	if (!yEstimL || !yMaxLeft)
		EATErr("Split [memory fault]");
	for (i = 0; i < tData->nY; i++)
		yMaxLeft[i] = 0.0;

	// Recorre todas las variables (columnas)
	for (xi = 0; xi < tData->nX; xi++)
	{
		nArray = GetSortArrayFormIndex(array, df, X[xi], tData->index, tData->nIndex);
		tData->varInfo[xi].SL = tData->varInfo[xi].SR = INFINITY;
		tData->varInfo[xi].S = -1;

		// El siguiente bucle empieza en '1' porque la rama izquierda
		// ( < array[i]) debe tener almenos 1 elemento (si empezara en
		// 0, en la primera iteracion la rama izquierda estaria vacia)
		for (i = 1; i < nArray; i++)
		{
			ErrEAT(tData, df, leaves, Y, X[xi], array[i], &errL, &errR, yEstimL);
			err = errL + errR;

			// R y s menor alternativo para cad xi
			if ((tData->varInfo[xi].SL + tData->varInfo[xi].SR) > err)
			{
				tData->varInfo[xi].SL = errL;
				tData->varInfo[xi].SR = errR;
				tData->varInfo[xi].S = array[i];
			}

			// Si el Error del split s en xi es menor se realiza el split con esa 's'
			if (err < errMin)
			{
				errMin = err;
				xiRes = xi;
				sRes = array[i];
				for (i2 = 0; i2 < tData->nY; i2++)
					yMaxLeft[i2] = yEstimL[i2];
			}
		}
	}

	NEATCpySon(TL, tData);
	NEATCpySon(TR, tData);

	for (i = 0; i < tData->nY; i++)
	{
		TL->y[i] = yMaxLeft[i];
		TR->y[i] = tData->y[i];
	}

	TL->index = NumValsLeft(df, tData->index, tData->nIndex, X[xiRes], sRes, &(TL->nIndex));
	TL->b[xiRes] = sRes;
	TL->R = QuadErr(df, TL->index, TL->nIndex, Y, TL->y, tData->nY);

	TR->index = NumValsRight(df, tData->index, tData->nIndex, X[xiRes], sRes, &(TR->nIndex));
	TR->a[xiRes] = sRes;
	TR->R = QuadErr(df, TR->index, TR->nIndex, Y, TR->y, tData->nY);

	TL->idF = TR->idF = tData->id;

	tData->errMin = errMin;
	tData->xi = X[xiRes];
	tData->s = sRes;

	free(yEstimL);
	free(yMaxLeft);
}

// ******************************* ALPHA y EAT **********************************
int IsContent(NEAT* t1, NEAT* t2)
{
	int i;
	for (i = 0; i < t1->nX; i++)
	{
		if (t1->a[i] > t2->a[i])
			return 0;
		if (t1->b[i] < t2->b[i])
			return 0;
	}
	return 1;
}

int IsPareto(VTREE* tree, NEAT* i)
{
	VTN* t;
	int order;

	VTreeIter iter;
	treeIterSet(tree, &iter, LEVELS, START);

	while ((t = treeIterNextNode(&iter)))
	{
		if (((NEAT*)t->data) == i || t->left != NULL ||
			IsContent(i, ((NEAT*)t->data))) // t->data contenido dentro de i
			continue;

		order = ComparePareto(i, ((NEAT*)t->data));
		if (order == -1 && (i->y <= ((NEAT*)t->data)->y))
			return 1;
		if (order == 1 && (i->y >= ((NEAT*)t->data)->y))
			return 1;
	}
	return 0;
}

void RBranch(VTN* t, VTREE* tree, double* errBranch, int* numLeaves)
{
	double errLeft, errRight;
	int numLeft, numRight;

	if (t->left == NULL)
	{
		*errBranch = ((NEAT*)t->data)->R;
		*numLeaves = 1;
		return; // fin recursividad
	}
	RBranch(t->left, tree, &errLeft, &numLeft);
	RBranch(t->right, tree, &errRight, &numRight);

	*errBranch = (errLeft + errRight);
	*numLeaves = (numLeft + numRight);
}

double Alpha(VTREE* tree)
{
	VTN* t, * tMin = NULL;
	double alpha;
	double errBranch;
	double alphaMin = INFINITY;
	int numLeaves;

	VTreeIter iter;
	treeIterSet(tree, &iter, LEVELS, START);

	while ((t = treeIterNextNode(&iter)))
	{
		if ((NEAT*)t->left == NULL) // descartar nodos hoja
			continue;
		//if(!IsPareto(tree, ((NEAT*)t->data)))
		//	continue;

		RBranch(t, tree, &errBranch, &numLeaves);
		alpha = (((NEAT*)t->data)->R - errBranch) / (numLeaves - 1);
		if (alpha < alphaMin)
		{
			alphaMin = alpha;
			tMin = t;
		}
	}
	return alphaMin;
}

VTN* SelectFatherHeuristic(VLIST* leaves, DFrameNum* df, int* X, int numStop)
{
	VTN *t1, *t2, *tSup;
	NEAT *tData1, *tData2;
	VListIter iter;

	int pos = 0;
	int conInf = 0; //Nodos inferiores
	int conClose = 0; //Nodos inferiores cerrados
	int conMaxInf = -1; //Nodo con más superiores
	int auxConMaxInf = -1;

	tSup = NULL;

	if (leaves->size == 1)
	{
		t1 = (VTN*)listGetRef(leaves, 0);
		listNulData(leaves, t1);
		listDel(leaves, 0);

		return t1;
	}

	listRefMoveFirst(leaves, &iter);

	while (t1 = (VTN*)listRefMoveNext(leaves, &iter))
	{
		pos++;

		tData1 = (NEAT*)t1->data;

		if (tData1->s != -1 || IsFinalNode(tData1, df, X, numStop))
			continue;

		conInf = 0;
		conClose = 0;
		auxConMaxInf = -1;

		while (t2 = (VTN*)listRefMoveNext(leaves, &iter))
		{
			tData2 = (NEAT*)t2->data;
			//Nodos inferiores
			if (ComparePareto(tData1, tData2) == 1)
			{
				conInf++;

				if (IsFinalNode(tData2, df, X, numStop))
					conClose++;
			}

			//Nodos superiores
			if (ComparePareto(tData1, tData2) == -1)
			{
				auxConMaxInf++;
			}
		}
		//Devuelve el VTN con cero inferiores
		if (conInf == 0)
		{
			listNulData(leaves, t1);
			listDel(leaves, 0);
			printf("\t\tNodo con 0 inferiores\n");

			return t1;
		}

		//Devuelve el VTN con todos sus inferiores cerrados
		if (conInf == conClose)
		{
			listNulData(leaves, t1);
			listDel(leaves, 0);
			printf("\t\tNodo con inferiores cerrados\n");

			return t1;
		}

		//Devuelve el VTN con más superiores
		if (auxConMaxInf > conMaxInf)
		{
			conMaxInf = auxConMaxInf;
			NEATCpyVTN(tSup, t1);
		}
	}
	printf("\t\tNodo con más superiores\n");
	return tSup;
}

VTN* SelectFather(VLIST* leaves, int num, int random)
{
	int pos = 0;
	VTN* t;

	// Random
	if (random)
		pos = rand() % num;

	t = (VTN*)listGetRef(leaves, pos);
	listNulData(leaves, t);
	listDel(leaves, pos);

	return t;
}

void EATStart(DFrameNum* df, int* X, int nX, int* Y, int nY, VTREE* tree, VLIST* leaves, VLIST* treeAlphaList)
{
	int pos, iCol;
	NEAT* t0;
	TREEALPHA* TAlpha = (TREEALPHA*)malloc(sizeof(TREEALPHA));

	if (!TAlpha)
		EATErr("EATStart [memory fault (TAlpha)]");

	// INICIALIZANDO NODO RAIZ ****************************************************
	if (!(t0 = (NEAT*)malloc(sizeof(NEAT))))
		EATErr("EAT [memory fault-2]");
	t0->id = 0;
	t0->idF = -1;
	NEATSet(t0, nX, nY, df->nRows);

	// inicializar 'a' con los valores manimos...
	for (iCol = 0; iCol < nX; iCol++)
		t0->a[iCol] = dfNumColumMinVal(df, X[iCol], &pos);

	// Calcular la estimacion (max) para el nodo raiz
	for (iCol = 0; iCol < nY; iCol++)
		t0->y[iCol] = dfNumColumMaxVal(df, Y[iCol], &pos);

	t0->R = QuadErr(df, t0->index, t0->nIndex, Y, t0->y, nY);
	treeAppendRoot(tree, t0);

	// TreeAlphList, primer arbol = Nodo Raiz (size=1)
	TAlpha->alpha = Alpha(tree);
	TAlpha->score = INFINITY;
	TAlpha->errFolders = 0.0;
	TAlpha->tree = treeNew(NEATPrint, NEATNew, NEATDel, NEATCpy);
	treeCopy(TAlpha->tree, tree);
	listAppend(treeAlphaList, TAlpha);

	// Inicializar lista de hojas "leaves" con el nodo raiz
	listAppend(leaves, tree->first);
}

int CheckLeaves(VLIST* leaves, int numStop)
{
	int n, numNoFinalLeaves = 0;
	VListIter iter;
	void* data;

	listRefMoveLast(leaves, &iter);
	for (n = 0; n < leaves->size; n++)
	{
		if (((NEAT*)(((VTN*)(iter->data))->data))->nIndex > numStop)
		{
			numNoFinalLeaves++;
			data = iter->data;
			iter->data = NULL;
			listDelIter(leaves, &iter, -1);
			listAdd(leaves, 0, data);
		}
		else
			listRefMovePrev(leaves, &iter);
	}
	return numNoFinalLeaves;
}

void EATHeuristic(DFrameNum* df, int* X, int nX, int* Y, int nY, int numStop, int random, VTREE* tree, VLIST* leaves, int numNoFinalLeaves, VLIST* treeAlphaList)
{
	NEAT* TL, * TR;
	VTN* t;
	NEAT* tData;

	TREEALPHA TAlpha;
	TAlpha.tree = tree;
	TAlpha.score = INFINITY;
	TAlpha.alpha = ((TREEALPHA*)(treeAlphaList->first->data))->alpha;
	TAlpha.errFolders = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;

	double* arrayForSplit;
	if (!(arrayForSplit = (double*)malloc(sizeof(double) * df->nRows)))
		EATErr("EAT [memory fault-0 (arrayForSplit)]");

	while (numNoFinalLeaves)
	{
		// SELECT FATHER .....
		t = SelectFather(leaves, numNoFinalLeaves, random);
		tData = (NEAT*)t->data;
		numNoFinalLeaves--;

		if (!(TR = (NEAT*)malloc(sizeof(NEAT))))
			EATErr("EAT [memory fault-3]");
		if (!(TL = (NEAT*)malloc(sizeof(NEAT))))
			EATErr("EAT [memory fault-4]");
		NEATSet(TR, nX, nY, 0);
		NEATSet(TL, nX, nY, 0);

		Split(tData, leaves, df, X, Y, TL, TR, arrayForSplit);
		TL->id = tree->size;
		TR->id = tree->size + 1;

		treeAppendSons(t, TL, TR, tree);

		if (IsFinalNode(TR, df, X, numStop))
		{
			listAppend(leaves, t->right); // Nodo final: Se pone al final de 'leaves'
			TAlpha.errFolders += TR->R;
		}
		else
		{
			listAdd(leaves, 0, t->right); // Nodo no-final: Se pone al principio de 'leaves'
			numNoFinalLeaves++;
		}
		if (IsFinalNode(TL, df, X, numStop))
		{
			listAppend(leaves, t->left); // Nodo final: Se pone al final de 'leaves'
			TAlpha.errFolders += TL->R;
		}
		else
		{
			listAdd(leaves, 0, t->left); // Nodo no-final: Se pone al principio de 'leaves'
			numNoFinalLeaves++;
		}

		TAlpha.alpha = Alpha(tree);
		listAddCopy(treeAlphaList, 0, &TAlpha);
	}
	free(arrayForSplit);
}

/********* BackTracking *****************************************************************************/

void NextSplit(NEAT* tData, VLIST* leaves, DFrameNum* df, int* X, int* Y, NEAT* TL, NEAT* TR, double* array)
{
	int i, nArray = 0, idS = 1;

	array = (double*)malloc(sizeof(double) * tData->nIndex);
	if (!array)
		EATErr("NextSplit [memory fault]");


	if (tData->xi == -1)
	{
		tData->xi = 0;
		nArray = GetSortArrayFormIndex(array, df, X[tData->xi], tData->index, tData->nIndex);
	}
	else if (tData->s < tData->sMax)
	{
		nArray = GetSortArrayFormIndex(array, df, X[tData->xi], tData->index, tData->nIndex);
		while (array[idS] <= tData->s)
			idS++;
	}
	else if (tData->xi < (tData->nX - 1))
	{
		tData->xi++;
		nArray = GetSortArrayFormIndex(array, df, X[tData->xi], tData->index, tData->nIndex);
	}
	else
		EATErr("NextSplit [Nooooooooooooooo :-( !!!!!]");

	tData->s = array[idS];
	tData->sMax = array[nArray - 1];

	NEATCpySon(TL, tData); // Copia 'a' y 'b' del padre a los hijos
	NEATCpySon(TR, tData);
	TL->b[X[tData->xi]] = tData->s;
	TR->a[X[tData->xi]] = tData->s;

	TL->index = NumValsLeft(df, tData->index, tData->nIndex, X[tData->xi], tData->s, &(TL->nIndex));
	TR->index = NumValsRight(df, tData->index, tData->nIndex, X[tData->xi], tData->s, &(TR->nIndex));

	EstimEAT(TL, df, leaves, Y, TL->index, TL->nIndex, TL->y);

	for (i = 0; i < tData->nY; i++)
		TR->y[i] = tData->y[i];

	TL->R = QuadErr(df, TL->index, TL->nIndex, Y, TL->y, tData->nY);
	TR->R = QuadErr(df, TR->index, TR->nIndex, Y, TR->y, tData->nY);

	TL->idF = TR->idF = tData->id;
}

VTN* FindSon(VLIST* list, VTN* t, VListIter* iterLeaves)
{
	static int nSons;
	VTN* aux;

	if ((*iterLeaves) == NULL)
	{
		nSons = 0;
		listRefMoveLast(list, iterLeaves);
	}

	while (nSons < 2 && (aux = (VTN*)listRefMovePrev(list, iterLeaves)))
	{
		if (((NEAT*)(aux->data))->idF == ((NEAT*)(t->data))->id)
		{
			nSons++;
			if ((*iterLeaves) == NULL)
				listRefMoveLast(list, iterLeaves);
			else
				listRefMoveNext(list, iterLeaves);
			return aux;
		}
	}
	return NULL;
}

int IsSplittable(NEAT* t)
{
	if ((t->xi == (t->nX - 1) && t->s == t->sMax) || t->nIndex < 2)
		return 0;
	return 1;
}

void EATBackTracking(DFrameNum* df, int* X, int nX, int* Y, int nY, int numStop, int numStopHeur,
	int random, VTREE* tree, VLIST* leaves, int numNoFinalLeaves, VLIST* treeAlphaList)
{
	int countBackTrack = 0;
	int numNoFinalLeavesAux;
	double bestErr = INFINITY;

	NEAT* TL, * TR;
	VTN* t, * son;
	NEAT* tData;
	VListIter iterLeaves;

	VLIST* inners = listNew(NEATPrintVTN, NEATNewVTN, NEATDelVTN, NEATCpyVTN);
	VLIST* leavesAux = listNew(NEATPrintVTN, NEATNewVTN, NEATDelVTN, NEATCpyVTN);
	VLIST* treeAlphaListAux = listNew(TreeAlphaPrint, TreeAlphaNew, TreeAlphaDel, TreeAlphaCpy);
	VTREE* treeAux = treeNew(NEATPrint, NEATNew, NEATDel, NEATCpy);
	VLIST* treeAlphaListBest = listNew(TreeAlphaPrint, TreeAlphaNew, TreeAlphaDel, TreeAlphaCpy);
	VTREE* treeBest = treeNew(NEATPrint, NEATNew, NEATDel, NEATCpy);

	TREEALPHA TAlpha;
	TAlpha.tree = tree;
	TAlpha.score = INFINITY;
	TAlpha.alpha = ((TREEALPHA*)(treeAlphaList->first->data))->alpha;
	TAlpha.errFolders = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;

	double* arrayForSplit;
	if (!(arrayForSplit = (double*)malloc(sizeof(double) * df->nRows)))
		EATErr("EATBackTracking [memory fault-0 (arrayForSplit)]");

	// *****************************************************************************

	while (numNoFinalLeaves)
	{
		countBackTrack++;
		if (countBackTrack > 100000) {
			printf("countBackTrack > 100000\n");
			break;
		}

		while (numNoFinalLeaves)
		{	
			// SELECT FATHER .....
			t = SelectFather(leaves, numNoFinalLeaves, 0); // 0-> Tomar siempre el primer nodo de la izquierda
			tData = (NEAT*)t->data;
			numNoFinalLeaves--;
			if (!IsSplittable(tData))
				EATErr("EstimEAT [tData NO Splittable]");

			listAdd(inners, 0, t); // al hacer split de 't', se convierte en nodo interno

			if (!(TR = (NEAT*)malloc(sizeof(NEAT))))
				EATErr("EAT [memory fault-3]");
			if (!(TL = (NEAT*)malloc(sizeof(NEAT))))
				EATErr("EAT [memory fault-4]");
			NEATSet(TR, nX, nY, 0);
			NEATSet(TL, nX, nY, 0);

			NextSplit(tData, leaves, df, X, Y, TL, TR, arrayForSplit);
			TL->id = tree->size;
			TR->id = tree->size + 1;

			treeAppendSons(t, TL, TR, tree);

			if (IsFinalNode(TR, df, X, numStop))
			{
				listAppend(leaves, t->right); // Nodo final: Se pone al final de 'leaves'
				TAlpha.errFolders += TR->R;
			}
			else
			{
				listAdd(leaves, 0, t->right); // Nodo no-final: Se pone al principio de 'leaves'
				numNoFinalLeaves++;
			}
			if (IsFinalNode(TL, df, X, numStop))
			{
				listAppend(leaves, t->left); // Nodo final: Se pone al final de 'leaves'
				TAlpha.errFolders += TL->R;
			}
			else
			{
				listAdd(leaves, 0, t->left); // Nodo no-final: Se pone al principio de 'leaves'
				numNoFinalLeaves++;
			}
			TAlpha.alpha = Alpha(tree);
			listAddCopy(treeAlphaList, 0, &TAlpha);
		}

		// Llamar al HEURISTICO para terminar el arbol --------------------------------------------------------
		if (numStopHeur >= 1 && numStopHeur < numStop)
		{
			//printf("Heuristic working...\n");
			listCopy(leavesAux, leaves);
			listCopy(treeAlphaListAux, treeAlphaList);
			treeCopy(treeAux, tree);
			numNoFinalLeavesAux = CheckLeaves(leavesAux, numStopHeur);
			EATHeuristic(df, X, nX, Y, nY, numStopHeur, random, treeAux, leavesAux, numNoFinalLeaves, treeAlphaListAux);
			if (((TREEALPHA*)(treeAlphaListAux->first->data))->errFolders < bestErr)
			{
				bestErr = ((TREEALPHA*)(treeAlphaListAux->first->data))->errFolders;
				listCopy(treeAlphaListBest, treeAlphaListAux);
				treeCopy(treeBest, treeAux);
				printf("\n\nmejora: err= %f, size: %d (countBackTrack=%d)\n", bestErr, treeBest->size, countBackTrack);
				//treePrint(treeBest);
			}
		}
		else if (((TREEALPHA*)(treeAlphaList->first->data))->errFolders < bestErr)
		{
			bestErr = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;
			listCopy(treeAlphaListBest, treeAlphaList);
			treeCopy(treeBest, tree);
			printf("\n\nmejora: err= %f, size: %d (countBackTrack=%d)\n", bestErr, treeBest->size, countBackTrack);
			//treePrint(treeBest);
		}

		while (inners->size > 0)
		{
			// 1- saca el ultimo nodo interior 't' al que se le hizo split
			t = SelectFather(inners, 0, 0);
			tData = (NEAT*)t->data;

			// 2- Pasar 't' a leaves y eliminar sus hijos de (1)leales y de (2)tree EN ESE ORDEN
			listAdd(leaves, 0, t);

			iterLeaves = NULL;
			while ((son = FindSon(leaves, t, &iterLeaves)))
			{
				if (!(IsFinalNode((NEAT*)(son->data), df, X, numStop)))
					numNoFinalLeaves--;
				else
					TAlpha.errFolders -= ((NEAT*)(son->data))->R;
				listDelIter(leaves, &iterLeaves, -1);
			}
			//TAlpha.errFolders += tData->R;
			treeDelBranches(tree, t);
			numNoFinalLeaves++;
			if (numNoFinalLeaves > leaves->size)
				numNoFinalLeaves = leaves->size;

			// 3- borrar el ultimo arbol de treeAlphaList
			listDelLast(treeAlphaList);

			// 4- Verificar si 't' es "splitable"
			if (IsSplittable(tData))
				break;
			else
				tData->xi = tData->sMax = -1.0;
		}
		// SI YA NO HAY NODOS INTERNOS (inners==0) SI ADEMAS EL UNICO
		// NODO HOJA (leaves==1) ESTA INICIALIZADO (x1==-1) ES PORQUE
		// DICHO NODO HOJA ES EL RAIZ Y YA HIZO TODOS SUS SPLITS, POR
		// TANTO FINALIZAMOS EL BACTRACKING ==> BREAK
		if (inners->size == 0 && leaves->size == 1 && ((NEAT*)(leaves->first->data))->xi == -1)
			break;
	}

	listEmpty(treeAlphaList);
	treeAlphaList->first = treeAlphaListBest->first;
	treeAlphaList->last = treeAlphaListBest->last;
	treeAlphaList->size = treeAlphaListBest->size;

	treeEmpty(tree);
	tree->first = treeBest->first;
	tree->size = treeBest->size;

	listFree(inners);
	listFree(leavesAux);
	listFree(treeAlphaListAux);
	treeFree(treeAux);
	free(arrayForSplit);

	printf("countBackTrack: %d\n", countBackTrack);
}


VLIST* EATBackHeur(DFrameNum* df, int* X, int nX, int* Y, int nY, int numStop, double errH, VTREE* treeBest)
{
	int pos, iCol;
	int N = df->nRows;
	int numNoFinalLeaves, countBackTrack = 0;
	double errHeuristic = errH;

	NEAT* t0, * TL, * TR;
	VTN* t, * son;
	NEAT* tData;
	VListIter iterLeaves;
	VLIST* leaves, * inners, * treeAlphaList, * treeBestAlphaList = NULL;
	VTREE* tree;
	TREEALPHA TAlpha;

	double* arrayForSplit;
	if (!(arrayForSplit = (double*)malloc(sizeof(double) * df->nRows)))
		EATErr("EAT [memory fault-0 (arrayForSplit)]");

	// Crea e inicializa 'leaves' (lista de nodos hoja) y 'treeAlphaList'
	leaves = (VLIST*)malloc(sizeof(VLIST));
	inners = (VLIST*)malloc(sizeof(VLIST));
	treeAlphaList = (VLIST*)malloc(sizeof(VLIST));
	treeBestAlphaList = (VLIST*)malloc(sizeof(VLIST));
	tree = (VTREE*)malloc(sizeof(VTREE));
	if (!leaves || !inners || !treeAlphaList || !treeBestAlphaList || !tree)
		EATErr("EAT [memory fault-1]");
	listIniAll(leaves, NEATPrintVTN, NEATNewVTN, NEATDelVTN, NEATCpyVTN);
	listIniAll(inners, NEATPrintVTN, NEATNewVTN, NEATDelVTN, NEATCpyVTN);
	listIniAll(treeAlphaList, TreeAlphaPrint, TreeAlphaNew, TreeAlphaDel, TreeAlphaCpy);
	listIniAll(treeBestAlphaList, TreeAlphaPrint, TreeAlphaNew, TreeAlphaDel, TreeAlphaCpy);
	treeIniAll(tree, NEATPrint, NEATNew, NEATDel, NEATCpy);

	// INICIALIZANDO NODO RAIZ ****************************************************
	if (!(t0 = (NEAT*)malloc(sizeof(NEAT))))
		EATErr("EAT [memory fault-2]");
	t0->id = 0;
	t0->idF = -1;
	NEATSet(t0, nX, nY, N);

	// inicializar 'a' con los valores manimos...
	for (iCol = 0; iCol < nX; iCol++)
		t0->a[iCol] = dfNumColumMinVal(df, X[iCol], &pos);

	// Calcular la estimacion (max) para el nod raiz
	for (iCol = 0; iCol < nY; iCol++)
		t0->y[iCol] = dfNumColumMaxVal(df, Y[iCol], &pos);

	t0->R = QuadErr(df, t0->index, t0->nIndex, Y, t0->y, nY);
	treeAppendRoot(tree, t0);
	listAppend(leaves, tree->first); // tree->first <=> t0
	numNoFinalLeaves = 1;

	// TreeAlphList, primer arbol = Nodo Raiz (size=1)
	TAlpha.alpha = Alpha(tree);
	TAlpha.score = INFINITY;
	TAlpha.errFolders = 0.0; //t0->R;
	TAlpha.tree = tree;
	listAppendCopy(treeAlphaList, &TAlpha);

	// 1- BACKTRAKING **********************************************************
	while (numNoFinalLeaves)
	{
		countBackTrack++;

		while (numNoFinalLeaves)
		{
			// SELECT FATHER .....
			t = SelectFather(leaves, numNoFinalLeaves, 0); // 0-> Tomar siempre el primer nodo de la izquierda
			tData = (NEAT*)t->data;
			numNoFinalLeaves--;
			if (!IsSplittable(tData))
				EATErr("EstimEAT [tData NO Splittable]");

			listAdd(inners, 0, t); // al hacer spli de 't', se convierte en nodo interno

			if (!(TR = (NEAT*)malloc(sizeof(NEAT))))
				EATErr("EAT [memory fault-3]");
			if (!(TL = (NEAT*)malloc(sizeof(NEAT))))
				EATErr("EAT [memory fault-4]");
			NEATSet(TR, nX, nY, 0);
			NEATSet(TL, nX, nY, 0);

			NextSplit(tData, leaves, df, X, Y, TL, TR, arrayForSplit);
			TL->id = tree->size;
			TR->id = tree->size + 1;

			treeAppendSons(t, TL, TR, tree);
			TAlpha.alpha = Alpha(tree);
			//TAlpha.errFolders += TL->R + TR->R - tData->R;
			TAlpha.tree = tree;
			listAddCopy(treeAlphaList, 0, &TAlpha);

			if (IsFinalNode(TR, df, X, numStop))
			{
				listAppend(leaves, t->right); // Nodo final: Se pone al final de 'leaves'
				TAlpha.errFolders += TR->R;
			}
			else
			{
				listAdd(leaves, 0, t->right); // Nodo no-final: Se pone al principio de 'leaves'
				numNoFinalLeaves++;
			}
			if (IsFinalNode(TL, df, X, numStop))
			{
				listAppend(leaves, t->left); // Nodo final: Se pone al final de 'leaves'
				TAlpha.errFolders += TL->R;
			}
			else
			{
				listAdd(leaves, 0, t->left); // Nodo no-final: Se pone al principio de 'leaves'
				numNoFinalLeaves++;
			}
		}

		if (((TREEALPHA*)(treeAlphaList->first->data))->errFolders < errH)
		{
			errH = ((TREEALPHA*)(treeAlphaList->first->data))->errFolders;
			listFree(treeBestAlphaList);
			listCopy(treeBestAlphaList, treeAlphaList);
			treeFree(treeBest);
			treeCopy(treeBest, tree);
		}

		while (inners->size > 0)
		{
			// 1- saca el ultimo nodo interior 't' al que se le hizo split
			t = SelectFather(inners, 0, 0);
			tData = (NEAT*)t->data;

			// 2- Pasar 't' a leaves y eliminar sus hijos de (1)leales y de (2)tree EN ESE ORDEN
			listAdd(leaves, 0, t);

			iterLeaves = NULL;
			while ((son = FindSon(leaves, t, &iterLeaves)))
			{
				if (!(IsFinalNode((NEAT*)(son->data), df, X, numStop)))
					numNoFinalLeaves--;
				else
					TAlpha.errFolders -= ((NEAT*)(son->data))->R;
				listDelIter(leaves, &iterLeaves, -1);
			}
			//TAlpha.errFolders += tData->R;
			treeDelBranches(tree, t);
			numNoFinalLeaves++;
			if (numNoFinalLeaves > leaves->size)
				numNoFinalLeaves = leaves->size;

			// 3- borrar el ultimo arbol de treeAlphaList
			listDelLast(treeAlphaList);

			// 4- Verificar si 't' es "splitable"
			if (IsSplittable(tData))
				break;
			else
				tData->xi = tData->sMax = -1.0;
		}
		// SI YA NO HAY NODOS INTERNOS (inners==0) SI ADEMAS EL UNICO
		// NODO HOJA (leaves==1) ESTA INICIALIZADO (x1==-1) ES PORQUE
		// DICHO NODO HOJA ES EL RAIZ Y YA HIZO TODOS SUS SPLITS, POR
		// TANTO FINALIZAMOS EL BACTRACKING ==> BREAK
		if (inners->size == 0 && leaves->size == 1 && ((NEAT*)(leaves->first->data))->xi == -1)
			break;
	}

	// 2- HEURISTICO ***********************************************************


	free(leaves);
	listFree(treeAlphaList);
	free(treeAlphaList);
	treeFree(tree);
	free(tree);
	free(arrayForSplit);

	printf("countBackTrack: %d\n", countBackTrack);
	return treeBestAlphaList;
}


int CheckEAT(VTREE* tree)
{
	VTreeIter it1, it2;
	VTN* t1, * t2;
	NEAT* n1, * n2;
	int resComp, check = 1;
	int i, nY = ((NEAT*)tree->first->data)->nY;

	treeIterSet(tree, &it1, LEVELS, START);
	while ((t1 = treeIterNextNode(&it1)))
	{
		if (t1->left == NULL) // es nodo hoja
		{
			n1 = (NEAT*)t1->data;
			treeIterSet(tree, &it2, LEVELS, START);
			while ((t2 = treeIterNextNode(&it2)))
			{
				if (t2->left == NULL) // es nodo hoja
				{
					n2 = (NEAT*)t2->data;
					resComp = ComparePareto(n1, n2);
					for (i = 0; i < nY; i++)
					{
						if (resComp == -1)
						{
							if (n1->y[i] <= n2->y[i])
							{
								//printf("OK");
							}
							else
							{
								printf("ERR (%d)\n", resComp);
								NEATPrint(n1);
								NEATPrint(n2);
								NEATPrint((NEAT*)(t2->father->data));
								printf("-----------------------------------------\n\n");
								check = 0;
								break;
							}
						}
						else if (resComp == 1)
						{
							if (n1->y[i] >= n2->y[i])
							{
								//printf("OK");
							}
							else
							{
								printf("ERR (%d)\n", resComp);
								NEATPrint(n1);
								NEATPrint(n2);
								printf("-----------------------------------------\n\n");
								check = 0;
								break;
							}
						}
						else
						{
							//printf("ni fu no fa...");
						}
					}
				}
			}
		}
	}
	return check;
}

/****************************************************************************************
 ########################## Pruning EAT ##########################
*****************************************************************************************/
double* Predictor(VTREE* tree, double* reg)
{
	VTN* ti = tree->first; // Nodo raiz del arbol

	while (ti->left != NULL) // Mientras 'ti' no sea nodo hoja
	{
		if (reg[((NEAT*)ti->data)->xi] < ((NEAT*)ti->data)->s)
			ti = ti->left;
		else
			ti = ti->right;
	}

	return ((NEAT*)ti->data)->y;
}

void TreesForRCV(DFrameNum* df, int* X, int nX, int* Y, int nY, int folder, int numStop, VLIST** TAiv)
{
	int v, N = df->nRows;
	int numRowsFolder = N / folder;
	DFrameNum notLv; // datos para Trainig (sub 'v')

	VTREE Tiv;
	treeIniAll(&Tiv, NEATPrint, NEATNew, NEATDel, NEATCpy);

	dfNumSetPartition(&notLv, df, (N - numRowsFolder));
	for (v = 0; v < folder; v++)
	{
		dfNumPartitionLess(&notLv, df, v, numRowsFolder);
		//TAiv[v]=EATHeuristic(&notLv, X, nX, Y, nY, numStop, 0, &Tiv);
		printf("Podar noooooooooooo !!!!! (lin 1365 mas o menos)\n\n");
		exit(0);
		treeFree(&Tiv);
	}
	free(notLv.data);
}

double RCV(DFrameNum* Lv, int* Y, int nY, double alphaIprim, int folder,
	VLIST** TAiv, TREEALPHA** BestTivs)
{
	int reg, v, i;
	double* predictor, Rcv = 0.0;
	VListIter iter;

	TREEALPHA* Tiv, * TivAux;

	for (v = 0; v < folder; v++)
	{
		Tiv = (TREEALPHA*)TAiv[v]->first->data;
		listRefMoveFirst(TAiv[v], &iter);
		while ((TivAux = listRefMoveNext(TAiv[v], &iter)) != NULL && TivAux->tree->size > 1)
		{
			if (TivAux->alpha <= alphaIprim)
				Tiv = TivAux;
		}

		BestTivs[v] = Tiv;
		//TEST +++++++++++++++++++++++++++++
		for (reg = 0; reg < Lv[v].nRows; reg++)
		{
			predictor = Predictor(Tiv->tree, Lv[v].data[reg]);

			for (i = 0; i < nY; i++)
				Rcv += (Lv[v].data[reg][Y[i]] - predictor[i])
				* (Lv[v].data[reg][Y[i]] - predictor[i]);
		}
	}
	return Rcv;
}

// *****************************  Pruning EAT *********************************
TREEALPHA* EATPruning(VLIST* treeAlphaList, DFrameNum* df,
	int* X, int nX, int* Y, int nY, int folder, int numStop)
{
	int v, N = df->nRows;
	double alphaIprim;
	int numRowsFolder = N / folder;

	DFrameNum* Lv = (DFrameNum*)malloc(sizeof(DFrameNum) * folder);
	if (!Lv)
		EATErr("EATPruning [memory fault-1]");
	for (v = 0; v < folder; v++)
	{
		dfNumSetPartition(&Lv[v], df, numRowsFolder);
		dfNumPartition(&Lv[v], df, v, numRowsFolder);
	}

	//******************************************************************************
	// PRUNING
	TREEALPHA* ta1, * ta2, * Tk;
	VListIter iter;
	VLIST** TAiv = (VLIST**)malloc(sizeof(VLIST*) * folder);
	TREEALPHA** BestTivs = (TREEALPHA**)malloc(sizeof(TREEALPHA*) * folder);
	TREEALPHA** BestTivsAux = (TREEALPHA**)malloc(sizeof(TREEALPHA*) * folder);

	if (TAiv == NULL || BestTivs == NULL || BestTivsAux == NULL)
		EATErr("TAiv / BestTivs / BestTivsAux is NULL");

	TreesForRCV(df, X, nX, Y, nY, folder, numStop, TAiv);

	listRefMoveFirst(treeAlphaList, &iter);
	Tk = ta1 = listRefMoveNext(treeAlphaList, &iter);

	BestTivs[0] = NULL;
	while (ta1 != NULL && (ta2 = listRefMoveNext(treeAlphaList, &iter)) != NULL && ta2->tree->size > 1)
	{
		alphaIprim = sqrt(ta1->alpha * ta2->alpha);
		ta1->score = RCV(Lv, Y, nY, alphaIprim, folder, TAiv, BestTivsAux) / N / nY;

		if (Tk->score > ta1->score || BestTivs[0] == NULL)
		{
			Tk = ta1; // TREEALPHA con el menor score
			for (v = 0; v < folder; v++)
				BestTivs[v] = BestTivsAux[v];
		}
		ta1 = ta2;
	}

	//******************************************************************************
	// SERule (vuelve a hacer Cross-Validation... solo Validation)
	int reg, i;
	double SE, s2 = 0.0; // 's' al cuadrado
	double* predictor;
	double dif1, dif2;

	for (v = 0; v < folder; v++)
	{

		for (reg = 0; reg < numRowsFolder; reg++)
		{
			predictor = Predictor(BestTivs[v]->tree, Lv[v].data[reg]);
			for (i = 0; i < nY; i++)
			{
				dif1 = (Lv[v].data[reg][Y[i]] - predictor[i]);
				dif2 = dif1 * dif1 - Tk->score;
				s2 += dif2 * dif2;
			}
		}
	}

	for (v = 0; v < folder; v++)
		free(Lv[v].data);
	free(Lv);
	SE = sqrt(s2 / N / N / nY / nY);

	//******************************************************************************
	// Seleccionar arbol definitivo: el de menor tamanno con un score con holgura SE
	double margin = Tk->score + SE;

	printf("Tk->score=%f (size:%d)\n", Tk->score, Tk->tree->size);
	printf("SE = %f\n", SE);
	printf("margin = %f\n", margin);

	listRefMoveLast(treeAlphaList, &iter);
	while ((Tk = listRefMovePrev(treeAlphaList, &iter)) != NULL)
		if (Tk->score <= margin)
			break;

	listNulData(treeAlphaList, Tk);
	listFree(treeAlphaList);
	for (v = 0; v < folder; v++)
		listFree(TAiv[v]);
	free(BestTivsAux);
	free(BestTivs);
	free(TAiv);
	return Tk;
}
/****************************************************************************************/

