#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>
// Define here as constant for easy change
#define REAL double

void printCheck ( REAL V[], int N )
{
  int x;

  REAL S=0;
  for (x=0; x<=N+1; x++)
    S = S + V[x];

  printf("\nCheckSum = %1.10e\nSome values: ", S);

  for (x=0; x<10; x++)
    printf("(%d)=%1.10f, ", x*N/10, V[x*N/10]);

  printf("(%d)=%1.10f\n", x*N/10, V[x*N/10]);
}
  


void SimulationStep ( REAL *In, REAL *Out, REAL L, int N )
{
  for (int x=1; x<=N; x++)
    Out[x] = In[x] - 2.0f*L*In[x] + L*In[x+1] + L*In[x-1];
}

  
void CopyVector ( REAL *In, REAL *Out, int N )
{
  for (int x=1; x<=N; x++)
    Out[x] = In[x];
} 
 

int main(int argc, char **argv)
{  
  int  x, t, N= 10000000, T=1000;
  REAL L= 0.123456, LK=1-2*L;  
  REAL * __restrict U1;
  REAL *U2,*U3;

  if (argc>1) { T = atoi(argv[1]); } // get  first command line parameter
  if (argc>2) { N = atoi(argv[2]); } // get second command line parameter
  if (argc>3) { L = atof(argv[3]); } // get  third command line parameter
 
  if (N < 1 || T < 1 || L >= 0.5) {
    printf("arguments: T N L (T: steps, N: vector size, L < 0.5)\n");
    return 1;
  }

  U1 = (REAL *) malloc ( sizeof(REAL)*(N+2) );
  U2=U1;
  if (!U1) { printf("Cannot allocate vectors\n"); exit(1); }

  // initialize temperatures at time t=0  
  for (x=1; x<=N; x++)
    U1[x] = x*3.1416;
 
  // initialize fixed boundary conditions on U1
  {
    U1[0]  = 1.2345678e+12;
    U1[N+1]= -1.2345678e+16;
  } 
  REAL pxz;  

  printf("Challenge #3: Simulate %d steps on 1-D vector of %d elements with L=%1.10e\n", T, N, L);
  U3 = (REAL *) malloc ( sizeof(REAL)*(N+2) );
  
    REAL px=U1[0]*L;
    REAL px2;
    int r=1;     
  if(N%(16*20)==0){   
  for(int p=0; p<20;p++){ 
     
    
   for (r=0; r<250;r++){ 
        _mm_prefetch((const char*)&U1[(p*N/20)+r*8],_MM_HINT_NTA);
       } 
   r=250;    
       
   pxz=U1[((p+1)*N/20)]*L;
          
   for (t=1; t<=T; t++)
  {         
    for (int x=1+p*N/20; x<=((p+1)*N/20); x+=16){ 
        _mm_prefetch((REAL*)&U1[x+r*8],_MM_HINT_NTA);
        _mm_prefetch((REAL*)&U1[x+8+r*8],_MM_HINT_NTA); 
       for (int i=0; i<16; i+=2){ 
            px2=*(U2+x+i)*L;
            U1[x+i] = (*(U2+x+i)*LK) + L*(*(U2+x+1+i)) + px;
            px=px2; 
            
            px2=*(U2+x+i+1)*L;
            U1[x+i+1] = (*(U2+x+1+i)*LK) + L*(*(U2+x+1+i+1)) + px;
            px=px2;
            } 
            }   
       
       px=L*U2[p*N/20];
    }
   px=pxz;
  }
}






else{


 
for (t=1; t<=T; t++)
  {  // loop on time
     
      for (int x=1; x<=N; x++)
        
        U3[x] = U1[x] - 2.0f*L*U1[x] + L*U1[x+1] + L*U1[x-1];
          CopyVector( U3, U1, N );
         
  }
  
  


}
  printCheck(U1,N);
  free(U1); 
  
  
  
}



