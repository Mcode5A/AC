
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


	char x; char y; char z;
};

void InitKernel(struct A *k, unsigned Klen)
{
	int val ;
	for (int i = 0; i < Klen; i++) {
		val= myRandom() & 7;
		/*k[i].x = val & 4;
		k[i].y = val & 2;
		k[i].z = val & 1;*/
		 
		k[i].x = (val & 4) ? 255 : 0;
		k[i].y = (val & 2) ? 255 : 0;
		k[i].z = (val & 1) ? 255 : 0;


	}
}

struct pixel { unsigned char x; unsigned char y; unsigned char z; };

void InitImage(struct pixel *I, unsigned N)
{
	for (int j = 0; j < N; j++)
	{
		I[j].x = myRandom() % 31;
		I[j].y = myRandom() % 31;
		I[j].z = myRandom() % 31;
	}
}

void TransfImage(struct pixel * I, unsigned N, const struct A *K, unsigned Klen)
{


	// copy I to T in order to prevent data dependencies




	for (int i = 0; i < N; i++) {
		unsigned vx = 0, vy = 0, vz = 0, sum;


    int k=0;
    int image_i = i + k;

    while (( image_i < N)&&(  k < Klen))
		{
			
		 
				/*vx = (K[k].x ) ? vx + I[i + k].x : vx;
				vy = (K[k].y ) ? vy + I[i + k].y : vy;
				vz = (K[k].z ) ? vz + I[i + k].z : vz;*/

				vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);
        image_i++;
        k++;
        
        if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
        
        
        
        
        
        
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         
         
          
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;  
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
         if (image_i < N){
        vx += (K[k].x & I[image_i].x);
				vy += (K[k].y & I[image_i].y);
				vz += (K[k].z & I[image_i].z);}
        image_i++;
        k++;
        
		}










		sum = (vx & 255) + (vy & 255) + (vz & 255) + 1;
		I[i].x = (vx & 255) * 31 / sum;
		I[i].y = (vy & 255) * 31 / sum;
		I[i].z = (vz & 255) * 31 / sum;

	}










}

int main(int argc, char **argv)
{
	int i, sumX, sumY, sumZ, N = 10000, Klen = 2000, Rep = 1000;

	seed = 12345;

	// obtain parameters at run time
	if (argc > 1) { N = atoi(argv[1]); }
	if (argc > 2) { Klen = atoi(argv[2]); }
	if (argc > 3) { Rep = atoi(argv[3]); }

	printf("Challenge #1: Vector size is %d. Kernel size is %d. Repeat %d times\n", N, Klen, Rep);

	// Create Image & Kernel
	struct pixel *image = (struct pixel*) malloc(N * sizeof(struct pixel));
	struct A *K = (struct A *) malloc(Klen * sizeof(struct A));

	InitImage(image, N);
	InitKernel(K, Klen);

	// Repeat
	for (i = 0; i < Rep; i++)
	{
		TransfImage(image, N, K, Klen);
		int ii;
		ii = myRandom() % N;  image[ii].x = myRandom() % 31;
		ii = myRandom() % N;  image[ii].y = myRandom() % 31;
		ii = myRandom() % N;  image[ii].z = myRandom() % 31;
	}

	for (i = 0, sumX = sumY = sumZ = 0; i < N; i++)
	{
		sumX += image[i].x;
		sumY += image[i].y;
		sumZ += image[i].z;
	}

	printf("Result: sumX= %d, sumY= %d, sumZ= %d\n", sumX, sumY, sumZ);

	free(image); free(K);

	return 0;
}



