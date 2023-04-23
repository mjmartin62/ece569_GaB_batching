__global__ void DataPassGB(int *VtoC,int *CtoV,int *Receivedword, int *Interleaver,int ColumnDegree,int N,int NbBranch, int inter, int numWords);

__global__ void DataPassGBIter0(int *Dev_VtoC,int *Dev_CtoV,int *Dev_Receivedword, int *Interleaver,int ColumnDegree,int N,int NbBranch, int numWords);

__global__ void CheckPassGB(int *DCtoV, int *VtoC,int M,int NbBranch,int RowDegree, int numWords);

__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int ColumnDegree,int N,int M,int NbBranch, int numWords);

__global__ void ComputeSyndrome(int *Decide,int *Mat,int RowDegree,int M, int *Dev_Syndrome, int numWords);
