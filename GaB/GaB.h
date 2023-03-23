__global__ void DataPassGB(int *VtoC,int *CtoV,int *Receivedword, int *Interleaver,int *ColumnDegree,int N,int NbBranch, int inter);

__global__ void DataPassGBIter0(int *Dev_VtoC,int *Dev_CtoV,int *Dev_Receivedword, int *Interleaver,int *ColumnDegree,int N,int NbBranch);
void DataPassGBIter0_old(int *Dev_VtoC,int *Dev_CtoV,int *Dev_Receivedword, int *Interleaver,int *ColumnDegree,int N,int NbBranch);
__global__ void CheckPassGB(int *DCtoV, int *VtoC,int M,int NbBranch,int *RowDegree);

__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int *ColumnDegree,int N,int M,int NbBranch);
__global__ void ComputeSyndrome(int *Decide,int *Mat,int *RowDegree,int M, int *Dev_Syndrome);
