/* ###########################################################################################################################
## Organization         : The University of Arizona
##                      :
## File name            : GaB.c
## Language             : C (ANSI)
## Short description    : Gallager-B Hard decision Bit-Flipping algorithm
##                      :
##                      :
##                      :
## History              : Modified 19/01/2016, Created by Burak UNAL
##                      :
## COPYRIGHT            : burak@email.arizona.edu
## ######################################################################################################################## */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

//#####################################################################################################
__global__ void DataPassGB(int *VtoC,int *CtoV,int *Receivedword,int *InterResult,int *Interleaver,int *ColumnDegree,int N,int NbBranch)
{
	int t,numB,n,buf;
	int Global;
	numB=0;
	//for (n=0;n<N;n++)
	//{
        n = threadIdx.x + blockIdx.x*blockDim.x;
        numB = ColumnDegree[n] * n;
		//Global=(Amplitude)*(1-2*ReceivedSymbol[n]);
		Global=(1-2*Receivedword[n]); 
		//Global=(1-2*(Decide[n] + Receivedword[n])); //Decide[n]^Receivedword[n];
		for (t=0;t<ColumnDegree[n];t++) Global+=(-2)*CtoV[Interleaver[numB+t]]+1;

		for (t=0;t<ColumnDegree[n];t++)
		{
		  buf=Global-((-2)*CtoV[Interleaver[numB+t]]+1);
		  if (buf<0)  VtoC[Interleaver[numB+t]]= 1; //else VtoC[Interleaver[numB+t]]= 1;
		  else if (buf>0) VtoC[Interleaver[numB+t]]= 0; //else VtoC[Interleaver[numB+t]]= 1;
		  else  VtoC[Interleaver[numB+t]]=Receivedword[n];
		}
		numB=numB+ColumnDegree[n];
	//}
}
//#####################################################################################################
//#####################################################################################################
__global__ void DataPassGBIter0(int *Dev_VtoC,int *Dev_CtoV,int *Dev_Receivedword,int *Dev_Interleaver,int *Dev_ColumnDegree,int N,int NbBranch)
{
	int t,numB,buf;
    int n = threadIdx.x + blockIdx.x*blockDim.x;
     numB = Dev_ColumnDegree[n] * n;
    // if(n > 1280) printf("\n %d and %d ", n , numB);
     if (n < N) {
 	   for (t=0;t<Dev_ColumnDegree[n];t++) {    
        Dev_VtoC[Dev_Interleaver[numB+t]]=Dev_Receivedword[n];
        }
     }
}

//##################################################################################################
__global__ void CheckPassGB(int *CtoV,int *VtoC,int M,int NbBranch,int *RowDegree)
{
   int t,numB=0,m,signe;
   m = threadIdx.x + blockIdx.x*blockDim.x;
   numB= RowDegree[m] * m;
     if (m < M) {
		signe=0;for (t=0;t<RowDegree[m];t++) signe^=VtoC[numB+t];
	    for (t=0;t<RowDegree[m];t++) 	CtoV[numB+t]=signe^VtoC[numB+t];
    }

}
//#####################################################################################################
__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int *ColumnDegree,int N,int M,int NbBranch)
{
   	int t,numB,n,buf;
	int Global;
    n = threadIdx.x + blockIdx.x*blockDim.x;
	numB=ColumnDegree[n] * n;
    
//	for (n=0;n<N;n++)
//	{
    if (n < N) {
		Global=(1-2*Receivedword[n]);
		for (t=0;t<ColumnDegree[n];t++) Global+=(-2)*CtoV[Interleaver[numB+t]]+1;
        if(Global>0) Decide[n]= 0;
        else if (Global<0) Decide[n]= 1;
        else  Decide[n]=Receivedword[n];
    }
//		numB=numB+ColumnDegree[n];
//	}
}
//#####################################################################################################
int ComputeSyndrome(int *Decide,int **Mat,int *RowDegree,int M)
{
	int Synd,k,l;

	for (k=0;k<M;k++)
	{
		Synd=0;
		for (l=0;l<RowDegree[k];l++) Synd=Synd^Decide[Mat[k][l]];
		if (Synd==1) break;
	}
	return(1-Synd);
}

