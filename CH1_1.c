
#include <stdio.h>
#include <stdlib.h>
 

// Variable used to generate pseudo-random numbers
unsigned int seed;

// Function to generate pseudo-random numbers
inline int myRandom() {
	seed = (214013 * seed + 2531011);
	return (seed >> 13);
}
struct A {


	unsigned char x, y, z;
};
struct pixel { unsigned char x; unsigned char y; unsigned char z; };

void InitKernel(struct pixel *k, unsigned Klen)
{
	unsigned  val;//o int o unsigned uno de los dos va un pco mas rapido
	for (int i = 0; i < Klen; i++) {
		val = myRandom() & 7;
		/*k[i].x = val & 4;
		k[i].y = val & 2;
		k[i].z = val & 1;*/

		k[i].x = (val & 4) ? 255 : 0;
		k[i].y = (val & 2) ? 255 : 0;
		k[i].z = (val & 1) ? 255 : 0;


	}
}

 

void InitImage(struct pixel *I, unsigned N)
{
	for (int j = 0; j < N; j++)
	{
		I[j].x = myRandom() % 31;
		I[j].y = myRandom() % 31;
		I[j].z = myRandom() % 31;
	}
}

void TransfImageMenor(struct pixel * I, int N, struct pixel *K, int Klen)
{
	struct pixel *T;
	// copy I to T in order to prevent data dependencies
	T = I;
	unsigned char vx = 0, vy = 0, vz = 0; unsigned sum;
	for (int i = 0; i < N; i++) {
		vx = 0; vy = 0; vz = 0;

		for (int k = 0; k < N - i; k++)
		{
			 
			vx += (K[k].x & T[k].x);
			vy += (K[k].y & T[k].y);
			vz += (K[k].z & T[k].z);
		}
		T++;
   
		sum = vx + vy + vz + 1;
		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;
   
   
	} 
} 


void TransfImage(struct pixel * __restrict  I, int N, struct pixel * __restrict K, int Klen)
{
	struct pixel * __restrict T ;
	// copy I to T in order to prevent data dependencies
	T = I;
  
	unsigned char vx = 0, vy = 0, vz = 0; unsigned sum;

  
	for (int i = 0; i < N-Klen ; i++) {
 
   // int xr=(i>=(N-Klen))?Klen-(i%(N-Klen)):Klen; 
    vx = 0; vy = 0; vz = 0;
    
		for (int k = 0; k < Klen ; k++)
		{
			vx += (K[k].x & T[k].x) ;
			vy += (K[k].y & T[k].y) ;
			vz += (K[k].z & T[k].z) ;
		}
		T++; 
		sum = vx + vy + vz + 1;
		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;
  }
  
	for (int i = N - Klen; i < N; i++) {
		vx = 0; vy = 0; vz = 0;

		int k;
		for (k = 0; k < N - i; k++)
		{
			/*vx = (K[k].x ) ? vx + I[i + k].x : vx;
			vy = (K[k].y ) ? vy + I[i + k].y : vy;
			vz = (K[k].z ) ? vz + I[i + k].z : vz;*/

			vx += (K[k].x & T[k].x);
			vy += (K[k].y & T[k].y);
			vz += (K[k].z & T[k].z);

		}
		T++;
		sum = vx + vy + vz + 1;
		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;

	}
}




int main(int argc, char **argv)
{
	int i, N = 10000, Klen = 2000, Rep = 1000;

	seed = 12345;

	// obtain parameters at run time
	if (argc > 1) { N = atoi(argv[1]); }
	if (argc > 2) { Klen = atoi(argv[2]); }
	if (argc > 3) { Rep = atoi(argv[3]); }

	printf("Challenge #1: Vector size is %d. Kernel size is %d. Repeat %d times\n", N, Klen, Rep);

	// Create Image & Kernel
	struct pixel *image = (struct pixel*) malloc(N * sizeof(struct pixel));
	struct pixel *K = (struct pixel *) malloc(Klen * sizeof(struct pixel));

	InitImage(image, N);
	InitKernel(K, Klen);


	if (N >= Klen) {
		for (i = 0; i < Rep; i++)
		{
			TransfImage(image, N, K, Klen);
			int ii;
			ii = myRandom() % N; image[ii].x = myRandom() % 31;
			ii = myRandom() % N; image[ii].y = myRandom() % 31;
			ii = myRandom() % N; image[ii].z = myRandom() % 31;
		}
	}
	else {
		for (i = 0; i < Rep; i++)
		{
			TransfImageMenor(image, N, K, Klen);
			int ii;
			ii = myRandom() % N; image[ii].x = myRandom() % 31;
			ii = myRandom() % N; image[ii].y = myRandom() % 31;
			ii = myRandom() % N; image[ii].z = myRandom() % 31;
		}
	}
 
	int sumX = 0, sumY = 0, sumZ = 0;
	for (i = 0; i < N; i++)
	{
		sumX += image[i].x;
		sumY += image[i].y;
		sumZ += image[i].z;
	}

	printf("Result: sumX= %d, sumY= %d, sumZ= %d\n", sumX, sumY, sumZ);

	free(image); free(K);

	return 0;
}


