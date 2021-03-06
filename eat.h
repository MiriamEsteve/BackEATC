/*
Estructura de nodo para algoritmo EAT
*/


#ifndef _DATA_DEA_
#define _DATA_DEA_

//#define INF (1.0/0.0)

#include "00lib/vlist.h"
#include "00lib/data.h"
#include "00lib/vtree.h"


// int NCOLS; // Numero de columnas del fichero de datos
// Indicar su valor antes de usar 'newNEAT', ...

/************************************************************************************************************************
 #################################### NEAT ###########################################
**************************************************************************************************************************/
// ************************ Estructura Var Info *****************************
typedef struct
{
	double SL, SR, S;
} VARINFO;

// **************************** Estructura NEAT *****************************
typedef struct NodeEAT
{
	int id, idF; // ID de nodo / ID del padre
	int nX, nY;  // nº de Xs (inputs) / nº de Ys (outputs)
	int *index, nIndex; // array datos compatibles con el nodo (soporte)
	VARINFO *varInfo; // tamanyo: nX (num. variables independientes)
	double R, errMin;
	int xi;
	double s, sMax; // 'sMax' es para el backtracking
	double *y;
	double *a, *b; // tamanyo: nX (num. variables independientes)
} NEAT;

// **************************** Funciones NEAT *****************************
void  EATErr(char *err); // Mensajes de error EAT
void  NEATPrint(void *d);
void* NEATNew(void *d);
void  NEATDel(void *d);
void  NEATCpy(void *d1, void *d2);
void  NEATCpySon(void *d1, void *d2); // copia info para el hijo (izq. o der.)
void  NEATSet(void *d, int nX, int nY, int N);

void  NEATPrintVTN(void *d);
void* NEATNewVTN(void *d);
void  NEATDelVTN(void *d);
void  NEATCpyVTN(void *d1, void *d2);

/************************************************************************************************************************
 ################################## Pruning EAT ########################################
**************************************************************************************************************************/
// ************************* Estructura Pruning EAT **************************
typedef struct
{
    double alpha;
    double score;
	double errFolders;
    VTREE *tree;
}TREEALPHA;

// ************************* Funciones Pruning EAT ***************************
void  TreeAlphaPrint(void *ta); // 'ta' -> treeAlpha
void* TreeAlphaNew(void *ta);
void  TreeAlphaDel(void *ta);
void  TreeAlphaCpy(void *ta1, void *ta2);

/************************************************************************************************************************
 #################################### EAT ###########################################
**************************************************************************************************************************/
void   ArrayPrint(double *array, int nArray);
void   ArrayPrintInt(int *array, int nArray);

double MaxDF(DFrameNum *df, int *index, int nIndex, int y);
double QuadErr(DFrameNum *df, int *index, int nIndex, int *Y, double *yEst, int nY);

int    ComparePareto(const void *a, const void *b);
int    CompareDouble(const void *a, const void *b);
int    IsValueInArray(double val, double *array, int nArray); // find duplicates
int    GetSortArrayFormIndex(double *array, DFrameNum *df, int column, int *index, int nIndex);

int*   NumValsLeft(DFrameNum *df, int *tIndex, int nIndex, int xi, double s, int *nIndexL); //Datos que quedan a la izq.
int*   NumValsRight(DFrameNum *df, int *tIndex, int nIndex, int xi, double s, int *nIndexR); //Datos que quedan a la der.

int    IsFinalNode(NEAT *t, DFrameNum *df, int *X, int numStop);

void   EstimEAT(NEAT *tData, DFrameNum *df, VLIST *leaves, int *Y, int *indexL, int nIndexL, double *estL);
void   ErrEAT(NEAT *tData, DFrameNum *df, VLIST *leaves, int *Y, int xCol, double s,
              double *errL, double *errR, double *estL);

void   Split(NEAT *tData, VLIST *leaves, DFrameNum *df, int *X, int *Y, NEAT *TL, NEAT *TR, double *array);

// ******************************* ALPHA y EAT **********************************
int  IsContent(NEAT *t1, NEAT *t2);
int  IsPareto(VTREE *tree, NEAT *i);
void RBranch(VTN *t, VTREE *tree, double *errBranch, int *numLeaves);
double Alpha(VTREE *tree);

VTN* SelectFatherHeuristic(VLIST *leaves, DFrameNum* df, int* X, int numStop);
VTN* SelectFather(VLIST *leaves, int num, int random);
void EATStart(DFrameNum *df, int *X, int nX, int *Y, int nY, VTREE *tree, VLIST *leaves, VLIST *treeAlphaList);

void EATHeuristic(DFrameNum *df, int *X, int nX, int *Y, int nY, int numStop, int random, VTREE *tree, VLIST *leaves, int numNoFinalLeaves, VLIST *treeAlphaList); //random indica

int CheckLeaves(VLIST* leaves, int numNoFinalNodes);
void NextSplit(NEAT *tData, VLIST *leaves, DFrameNum *df, int *X, int *Y, NEAT *TL, NEAT *TR, double *array);
int	 IsSplittable(NEAT *t);
void EATBackTracking(DFrameNum *df, int *X, int nX, int *Y, int nY, int numStop, int numStopHeur, 
					 int random, VTREE* tree, VLIST* leaves, int numNoFinalLeaves, VLIST* treeAlphaList);

int CheckEAT(VTREE *tree);


/************************************************************************************************************************
 ################################## EAT Prunning #######################################
**************************************************************************************************************************/
double* Predictor(VTREE *tree, double *reg);
void TreesForRCV(DFrameNum *df, int *X, int nX, int *Y, int nY, int folder, int numStop, VLIST **TAiv);
double RCV(DFrameNum *Lv, int *Y, int nY, double alphaIprim, int folder, VLIST **TAiv, TREEALPHA **BestTivs);

// *****************************  Pruning EAT *********************************
TREEALPHA* EATPruning(VLIST *treeAlphaList, DFrameNum *df, int *X, int nX, int *Y, int nY, int folder, int numStop);

#endif







