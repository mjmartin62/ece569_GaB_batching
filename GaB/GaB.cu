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
__global__ void DataPassGB(int *VtoC, int *CtoV, int *Receivedword, int *Interleaver,int *ColumnDegree,int N,int NbBranch, int iter)
{
	int t,numB,n,buf;
	int Global;
	numB=0;
	



        n = threadIdx.x + blockIdx.x*blockDim.x;

    //debug code ####################################################################################
    /*
    if (n == 0) {
        printf("VtoC[0] is: %d \n",VtoC[0]);
        printf("VtoC[1] is: %d \n",VtoC[1]);
        printf("VtoC[2] is: %d \n",VtoC[2]);
        
        printf("Receivedword[0] is: %d \n",Receivedword[0]);
        printf("Receivedword[1] is: %d \n",Receivedword[1]);
        printf("Receivedword[2] is: %d \n",Receivedword[2]);
       
        
        }
        */


        numB = ColumnDegree[n] * n;
    if (n < N) {
        if (iter == 0) {
              for (t=0;t<ColumnDegree[n];t++)     
               VtoC[Interleaver[numB+t]]=Receivedword[n];





        } else {
        
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
           }
     }
	
}
//#####################################################################################################
//#####################################################################################################
/*
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
*/

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

__global__ void CheckPassGB(int *CtoV,int *VtoC,int M,int NbBranch,int *RowDegree, *numWords)
{
    int t,numB=0,i,signe;
    // Calc relative global memory index where i spans multiple arrays for multiple words
    i = threadIdx.x + blockIdx.x*blockDim.x;
    // Calculate strided position for message arrays
    numB = RowDegree * i;

    // Conditional is boundary check
    if (i < M*numWords) {
        signe=0;
        for (t=0;t<RowDegree;t++) 
            signe^=VtoC[numB+t];
        for (t=0;t<RowDegree;t++)   
            CtoV[numB+t]=signe^VtoC[numB+t];
    }

    

}

//#####################################################################################################
__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int *ColumnDegree,int N,int M,int NbBranch)
{
   	int t,numB,n,buf;
	int Global;
    n = threadIdx.x + blockIdx.x*blockDim.x;
	numB=ColumnDegree[n] * n;
    

    if (n < N) {
		Global=(1-2*Receivedword[n]);
		for (t=0;t<ColumnDegree[n];t++) Global+=(-2)*CtoV[Interleaver[numB+t]]+1;
        if(Global>0) Decide[n]= 0;
        else if (Global<0) Decide[n]= 1;
        else  Decide[n]=Receivedword[n];
    }

}
//#####################################################################################################
__global__ void ComputeSyndrome(int *Decide,int *Mat,int *RowDegree,int M, int *Dev_Syndrome)
{
	int Synd,k,l;
    //This needs reduction function 
    __shared__ int sh_Synd[648];
    
     int n = threadIdx.x + blockIdx.x*blockDim.x;
     int thd_id = threadIdx.x;

     if(n ==0 ) *Dev_Syndrome = 1;
     
     for (l=0;l<RowDegree[n];l++)Synd=Synd^Decide[Mat[n*8 + l]];    

     if (n < M) sh_Synd[thd_id] = Synd; 
     __syncthreads();
     
    //Reduce to a single value 
    for(int stride = blockDim.x/2 ; stride > 0; stride = stride/2) {
     sh_Synd[thd_id] = sh_Synd[thd_id] | sh_Synd[thd_id + stride];
     __syncthreads();
     }
    
     if (thd_id == 0 ) atomicMin(Dev_Syndrome, (1 - sh_Synd[0])); 

}

