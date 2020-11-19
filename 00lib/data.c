#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "data.h"

#pragma warning(disable: 4996)


/****************************************************************************************
 ############### Lectura de un fichero -> ###########################
 ############### matriz de datos de tipo DOUBLE ####################
******************************************************************************************/
// *********** Mensajes de error DATAFRAME ********************
void DATAErr(char *err)
{
	fprintf(stderr, "  ## DFRAME-ERROR >> %s\n", err);
	getchar();
	exit(0);
}

void dfNumIni(DFrameNum *df)
{
	df->nRows = df->nCols = 0;
	df->header = NULL;
	df->data = NULL;
}

DFrameNum* dfNumNew(int rows, int cols, char **head)
{
	DFrameNum *df = (DFrameNum*)malloc(sizeof(DFrameNum));
	int i;

	if(!df)
		DATAErr("dfNumNew [memory fault-0]");

	df->nRows = rows;
	df->nCols = cols;

	if( !(df->header = (char**)malloc(sizeof(char*) * cols)) )
		DATAErr("dfNumNew [memory fault-1]");
	for(i=0 ; i<cols ; i++)
		if(!(df->header[i] = _strdup(head[i])))
			DATAErr("dfNumNew [memory fault-2]");

	if( !(df->data = (double**)malloc(sizeof(double*) * rows)) )
		DATAErr("dfNumNew [memory fault-3]");
	for(i=0 ; i<rows ; i++)
		if(!(df->data[i] = (double*)malloc(sizeof(double) * cols)))
			DATAErr("dfNumNew [memory fault-4]");
	return df;
}

DFrameNum* dfNumLoad(char *fileName, char *sep)
{
	DFrameNum *df;
	char cad[1000], cad2[1000];
	char *tokens[100];
	int rows=0, cols=0;
	FILE *pf = fopen(fileName, "r");

	if(!pf)
		DATAErr("dfLoad [File not open]");

	fgets(cad, 1000, pf);
	tokens[cols] = strtok(cad, sep);
	cols++;
	while( (tokens[cols] = strtok(NULL, sep)) )
		cols++;
	tokens[cols-1][strlen(tokens[cols-1])-1]='\0';

	while(!feof(pf))
		if(fgets(cad2, 1000, pf))
			rows++;

	df = dfNumNew(rows, cols, tokens);
	fseek(pf, SEEK_SET, 0);
	fgets(cad, 1000, pf);
	rows=0;
	while(!feof(pf))
		if(fgets(cad, 1000, pf))
		{
			cols=0;
			tokens[cols] = strtok(cad, sep);
			df->data[rows][cols] = atof(tokens[cols]);
			cols++;
			while( (tokens[cols] = strtok(NULL, sep)) ){
				df->data[rows][cols] = atof(tokens[cols]);
				cols++;
			}
			rows++;
		}
	fclose(pf);
	return df;
}

void dfNumFree(DFrameNum *df)
{
	int i;
	
	if(!df)
		return;

	for(i=0 ; i < df->nCols ; i++)
		free(df->header[i]);
	free(df->header);

	for(i=0 ; i < df->nRows ; i++)
		free(df->data[i]);
	free(df->data);
	free(df);
}

void dfNumPrint(DFrameNum *df)
{
	int rows, cols;

	printf("rows: %d, cols: %d\n", df->nRows, df->nCols);

	for(cols=0 ; cols<df->nCols ; cols++)
		printf("\t%s", df->header[cols]);
	printf("\n");

	for(rows=0 ; rows<df->nRows ; rows++)
	{
		printf("[%3d]->", rows);
		for(cols=0 ; cols<df->nCols ; cols++)
			printf("%12.6f", df->data[rows][cols]);
		printf("\n");
	}
}

void dfNumPrintFile(DFrameNum *df, char *nameFile)
{
	int rows, cols;
	FILE *f=fopen(nameFile, "w");
	
	if(!f)
		DATAErr("dfPrintFile [file open???]");

	fprintf(f, "%s", df->header[0]);
	for(cols=1 ; cols<df->nCols ; cols++)
		fprintf(f, ";%s", df->header[cols]);
	fprintf(f, "\n");

	for(rows=0 ; rows<df->nRows ; rows++)
	{
		fprintf(f, "%12.9f", df->data[rows][0]);
		for(cols=1 ; cols<df->nCols ; cols++)
			fprintf(f, ";%12.9f", df->data[rows][cols]);
		fprintf(f, "\n");
	}
	fclose(f);
}

void dfNumRandomize(DFrameNum *df)
{
	int iRow, iRowRand;
	double *rowAux;

	for(iRow=0 ; iRow<df->nRows ; iRow++)
	{
		iRowRand = rand() % df->nRows;

		rowAux = df->data[iRow];
		df->data[iRow] = df->data[iRowRand];
		df->data[iRowRand] = rowAux;
	}
}

void dfNumSetPartition(DFrameNum *df1, DFrameNum *df2, int nRowsPart)
{
	df1->header = df2->header;
	df1->nRows = nRowsPart;
	df1->nCols = df2->nCols;
	if( !(df1->data = (double**)malloc(sizeof(double*) * nRowsPart)) )
		DATAErr("dfNumSetPartition [data]");
}

void dfNumPartition(DFrameNum *df1, DFrameNum *df2, int v, int length)
{
	int i, ini = v*length;
	
	for(i=0 ; i<length ; i++)
		df1->data[i] = df2->data[(ini+i)];
}

void dfNumPartitionLess(DFrameNum *df1, DFrameNum *df2, int v, int length)
{
	int i, ini = v*length;
	for(i=0 ; i<df1->nRows ; i++)
		if(i<ini)
			df1->data[i] = df2->data[i];
		else
			df1->data[i] = df2->data[i+length];
}

void dfNumBag(DFrameNum *dfBag, DFrameNum *df)
{
	int i;

	for(i=0 ; i<df->nRows ; i++)
		dfBag->data[i] = df->data[ (int)( (rand()/(double)RAND_MAX)*df->nRows ) ];
}


///////////////////////////////////////////////////////////////////////////////
double dfNumColumMaxVal(DFrameNum *df, int col, int *pos)
{
	int row;
	double maxVal;
	
	if(col<0 || col>=df->nCols)
		DATAErr("## ERR 'dfColumMaxVal' (col)");
	
	maxVal = df->data[0][col];
	for(row=1 ; row<df->nRows ; row++)
	{
		if(maxVal < df->data[row][col])
		{
			maxVal = df->data[row][col];
			(*pos) = row;
		}
	}
	return maxVal;
}

double dfNumColumMinVal(DFrameNum *df, int col, int *pos)
{
	int row;
	double minVal;
	
	if(col<0 || col>=df->nCols)
		DATAErr("## ERR 'dfColumMinVal' (col)");
	
	minVal = df->data[0][col];
	for(row=1 ; row<df->nRows ; row++)
	{
		if(minVal > df->data[row][col])
		{
			minVal = df->data[row][col];
			(*pos) = row;
		}
	}
	return minVal;
}

///////////////////////////////////////////////////////////////////////////////

void dfNumSetColumn(DFrameNum *df, int col, double *arr)
{
	int i;
	
	for(i=0 ; i<df->nRows ; i++)
		df->data[i][col] = arr[i];
}

void dfNumAddColumn(DFrameNum *df, int pos, char *name, double val)
{
	int r, c; // rows, cols

	if(pos<0 || pos>df->nCols)
		DATAErr("## ERR 'dfAddColumn' (pos)");

	// HEADER
	df->header = (char**)realloc(df->header, sizeof(char*)*(df->nCols+1));
	if(!df->header)
		DATAErr("## ERR 'dfAddColumn' (header)");

	for(c=df->nCols ; c>0 ; c--)
		if(c>pos)
			df->header[c] = df->header[c-1];
		else
			break;
	df->header[pos] = strdup(name);

	// DATA
	for(r=0 ; r<df->nRows ; r++)
	{
		df->data[r] = (double*)realloc(df->data[r], sizeof(double)*(df->nCols+1));
		if(!df->data[r])
			DATAErr("## ERR 'dfAddColumn' (data)");
		for(c=df->nCols ; c>0 ; c--)
			if(c>pos)
				df->data[r][c] = df->data[r][c-1];
			else
				break;
		df->data[r][pos] = val;
	}
	df->nCols++;
}

void dfNumCpyColumn(DFrameNum *df, DFrameNum *df2, int pos, int pos2)
{
	int i;

	if(pos<0 || pos>=df->nCols || pos2<0 || pos2>=df2->nCols)
		DATAErr("## ERR 'dfCpyColumn' (column pos)");
	if(df->nRows != df2->nRows)
		DATAErr("## ERR 'dfCpyColumn' (num rows)");

	for(i=0 ; i<df->nRows ; i++)
		df->data[i][pos] = df2->data[i][pos2];
}

