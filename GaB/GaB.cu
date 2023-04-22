/* ########################################################################################################################
## Organization         : The University of Arizona
##                      :
## File name            : GaB.cu
## Language             : C (ANSI)
## Short description    : Gallager-B Hard decision Bit-Flipping algorithm
##                      :
##                      :
##                      :
## History              : Modified 12/04/2023
## ########################################################################################################################*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

//#####################################################################################################
// Update VN to CN message array [Includes Multi Codeword Processing]
__global__ void DataPassGB(int *VtoC, int *CtoV, int *Receivedword, int *Interleaver,int ColumnDegree,int N,int NbBranch, int iter, int numWords)
{
	int t,numB,n,buf;
	int Global;
	// Calc relative global memory index where n spans multiple concatenated arrays for multiple words
    n = threadIdx.x + blockIdx.x*blockDim.x;
    // Spaced position in interleaver matrix where modulo operation allows for multi word wrapping
    numB = (ColumnDegree * n) % NbBranch;
    // Find which CW in concatenated array the thread is associated and calculate the offset for the concatenated array
    int CW_offset = (n / N) * NbBranch;

                    // DEBUG CODE
            //if (n < N)
             //  printf("CW offset = %d \n",CW_offset);



    // Conditional is boundary check
    if (n < N*numWords) {
        if (iter == 0) {
            for (t=0;t<ColumnDegree;t++)     
               VtoC[Interleaver[numB+t] + CW_offset]=Receivedword[n];
        }
        else {
		    //Global=(Amplitude)*(1-2*ReceivedSymbol[n]);
		    Global=(1-2*Receivedword[n]); 
		    //Global=(1-2*(Decide[n] + Receivedword[n])); //Decide[n]^Receivedword[n];
		    for (t=0;t<ColumnDegree;t++) 
                Global+=(-2)*CtoV[Interleaver[numB+t] + CW_offset]+1;

		    for (t=0;t<ColumnDegree;t++) {
		        buf=Global-((-2)*CtoV[Interleaver[numB+t] + CW_offset]+1);
		        if (buf<0)  
                    VtoC[Interleaver[numB+t] + CW_offset]= 1; //else VtoC[Interleaver[numB+t]]= 1;
		        else if (buf>0) 
                    VtoC[Interleaver[numB+t] + CW_offset]= 0; //else VtoC[Interleaver[numB+t]]= 1;
		        else  
                    VtoC[Interleaver[numB+t] + CW_offset]=Receivedword[n];
		    }
        }
    }
}

//##################################################################################################
// Update the CN to VN message array [Includes Multi Codeword Processing]
/*
__global__ void CheckPassGB(int *CtoV,int *VtoC,int M,int NbBranch,int RowDegree, int numWords)
{
    int t,numB=0,m,signe;
    // Calc relative global memory index where m spans multiple concatenated arrays for multiple words
    m = threadIdx.x + blockIdx.x*blockDim.x;
    // Calculate strided position for message arrays
    numB = (RowDegree * m) % NbBranch;
    // Find which CW in concatenated array the thread is associated and calculate the offset for the concatenated array
    int CW_offset = (m / M) * NbBranch;


    // Conditional is boundary check
    if (m < M*numWords) {
        signe=0;
        for (t=0;t<RowDegree;t++) {
            signe^=VtoC[numB+t + CW_offset];
            

        }

        for (t=0;t<RowDegree;t++) {     
            CtoV[numB+t + CW_offset]=signe^VtoC[numB+t + CW_offset];
        }
            
    }
}
*/

//##################################################################################################
// Update the CN to VN message array [Includes Multi Codeword Processing]

__global__ void CheckPassGB(int *CtoV,int *VtoC,int M,int NbBranch,int RowDegree, int numWords)
{
    // Calc relative global memory index where m spans multiple concatenated message arrays for multiple words
    int m = threadIdx.x + blockIdx.x*blockDim.x;
    int tid = threadIdx.x;

    // shared memory declaration for a copy of the VtoC message array (blocksize*row degree)
    __shared__ int sh_VtoC[128];
    // shared memory declaration for signe (blocksize* 1/2 * row degree)
    __shared__ int sh_signe[128];

    // Find which CW in concatenated array the thread is associated with
    //int CW_offset = m % NbBranch;

    // Conditional is boundary check
    if (m < NbBranch*numWords) {

        // pull in from global to shared memory and sync threads before subsequent computations occur
        sh_VtoC[tid] = VtoC[m];
        __syncthreads();
        // make copy and sync threads
        sh_signe[tid] = sh_VtoC[tid];
        __syncthreads();

        // Reduction to single signe value for each CN
        // Reduction loop set up to limit thread control divergence
        for (int boundary = blockDim.x; boundary > blockDim.x/RowDegree; boundary = boundary/2) {
            if (tid < boundary/2) {
                if (m == 686)
                    printf("in this loop boundary = %d\n", boundary);
                int tmp = sh_signe[tid*2] ^ sh_signe[tid*2 + 1];
                __syncthreads();
                sh_signe[tid] = tmp;
                __syncthreads();

            }
        }
        // Sync threads before writing to global memory then construct
        __syncthreads();
        CtoV[m] = sh_signe[tid / RowDegree] ^ sh_VtoC[tid];
    
    }
}


//#####################################################################################################
//  Update the VN's [Includes Multi Codeword Processing]
__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int ColumnDegree,int N,int M,int NbBranch, int numWords)
{
   	int t,numB,n,buf;
	int Global;
    // Calc relative global memory index where n spans multiple concatenated arrays for multiple words
    n = threadIdx.x + blockIdx.x*blockDim.x;
	// Spaced position in interleaver matrix where modulo operation allows for multi word wrapping
    numB = (ColumnDegree * n) % NbBranch;
    // Find which CW in concatenated array the thread is associated and calculate the offset for the concatenated array
    int CW_offset = (n / N) * NbBranch;



    
    // Conditional is boundary check
    if (n < N*numWords) {
		Global=(1-2*Receivedword[n]);



		for (t=0;t<ColumnDegree;t++) 
            Global+=(-2)*CtoV[Interleaver[numB+t] + CW_offset]+1;

            // DEBUG CODE
            //if (n > N)
             //   printf("Global value 2nd loc = %d \n",Global);


        if(Global>0) 
            Decide[n]= 0;
        else if (Global<0) 
            Decide[n]= 1;
        else  
            Decide[n]=Receivedword[n];
    }


}

//#####################################################################################################
// Calculate Syndrome; determine if corrected word is valid codeword

// This kernel is from the Single CW per kernel prototype 
/*
__global__ void ComputeSyndrome(int *Decide,int *Mat,int *RowDegree,int M, int *Dev_Syndrome, int numWords)
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
*/

// This is a temp kernel w/o optimzation in mind
__global__ void ComputeSyndrome(int *Decide,int *Mat,int RowDegree,int M, int *Dev_Syndrome, int numWords)
{
	int Synd,k,l,i;

    // Single thread per CW Syndrome calculation
    i = threadIdx.x * 1296;
    
	for (k=0;k<M;k++) {
		Synd=0;
		for (l=0;l<RowDegree;l++) {
            Synd=Synd^Decide[Mat[k*RowDegree+l] + i];
            //printf("Kernel Internal Synd =  %d  \n",Synd);

        }
        
        if (Synd == 1)
            break;

    }

    // Update Syndrome tracker array; each entry in array is assigned to single CW syndrome result
    Dev_Syndrome[threadIdx.x] = 1-Synd;


}