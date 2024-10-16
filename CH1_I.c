
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>


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

void InitKernel(struct A *k, unsigned Klen)
{
	unsigned  val;//o int o unsigned uno de los dos va un pco mas rapido
	int i;
	for (i = 0; i < Klen; i++) {
		val = myRandom() & 7;
		/*k[i].x = val & 4;
		k[i].y = val & 2;
		k[i].z = val & 1;*/

		k[i].x = (val & 4) ? 0xff : 0;
		k[i].y = (val & 2) ? 0xff : 0;
		k[i].z = (val & 1) ? 0xff : 0;

	}
}


void InitImage(struct pixel *I, unsigned N)
{
	int j;
	for (j = 0; j < N; j++)
	{
		I[j].x = myRandom() % 31;
		I[j].y = myRandom() % 31;
		I[j].z = myRandom() % 31;
	}
}
void TransfImageMenor(struct pixel * I, unsigned N, struct A *K, unsigned Klen)
{
	struct pixel *T;
	// copy I to T in order to prevent data dependencies
	int i;
	for (i = 0; i < N; i++) {
		unsigned char vx = 0, vy = 0, vz = 0; unsigned sum;
		T = &I[i];
		int k;
		for (k = 0; k < N - i; k++)
		{

			vx += (K[k].x & T[k].x);
			vy += (K[k].y & T[k].y);
			vz += (K[k].z & T[k].z);

		}
		sum = vx + vy + vz + 1;
		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;
	}
}
void TransfImage(struct pixel * I, unsigned N, struct A *K, unsigned Klen)
{
	struct pixel *T;
	// copy I to T in order to prevent data dependencies
	int i;
	for (i = 0; i < N - Klen; i++) {
		unsigned char vx = 0, vy = 0, vz = 0; unsigned sum;



		__m256i t0, t1, t2, tmask_v1, tmask_v2, tmask_v3, tr1, tr2, tr3;

		tr1 = _mm256_setzero_si256();
		tr2 = _mm256_setzero_si256();
		tr3 = _mm256_setzero_si256();




		int k;
		for (k = 0; k < Klen; k += 32)
		{
			int image_i = i + k;
			//bucle +32
			tmask_v1 = _mm256_load_si256((__m256i *)(&K[k].x));///epi8
			tmask_v2 = _mm256_load_si256((__m256i *)(&K[k + 10].z));///epi8
			tmask_v3 = _mm256_load_si256((__m256i *)(&K[k + 21].y));///epi8


			t0 = _mm256_load_si256((__m256i *)(&I[image_i].x));///epi8
			t1 = _mm256_load_si256((__m256i *)(&I[image_i + 10].z));///epi8
			t2 = _mm256_load_si256((__m256i *)(&I[image_i + 21].y));///epi8

			t0 = _mm256_and_si256(t0, tmask_v1);
			t1 = _mm256_and_si256(t1, tmask_v2);
			t2 = _mm256_and_si256(t2, tmask_v3);


			tr1 = _mm256_add_epi8(tr1, t0);
			tr2 = _mm256_add_epi8(tr2, t1);  //opcional inicializar trx
			tr3 = _mm256_add_epi8(tr3, t2);
			//fin bucle

		}

		// creacion mascaras para separar variables xyz

		tmask_v1 = _mm256_set_epi8(0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00);
		tmask_v2 = _mm256_set_epi8(0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff);
		tmask_v3 = _mm256_set_epi8(0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00, 0xff, 0x00, 0x00);


		//fase separacion de variables x y z
		t0 = _mm256_add_epi8(_mm256_and_si256(tr1, tmask_v1), _mm256_and_si256(tr2, tmask_v2));
		t0 = _mm256_add_epi8(_mm256_and_si256(tr3, tmask_v3), t0);



		t1 = _mm256_add_epi8(_mm256_and_si256(tr1, tmask_v2), _mm256_and_si256(tr2, tmask_v3));
		t1 = _mm256_add_epi8(_mm256_and_si256(tr3, tmask_v1), t1);



		t2 = _mm256_add_epi8(_mm256_and_si256(tr1, tmask_v3), _mm256_and_si256(tr2, tmask_v1));
		t2 = _mm256_add_epi8(_mm256_and_si256(tr3, tmask_v2), t2);


		///fase suma byte a byte
		tr1 = _mm256_setzero_si256();

		t0 = _mm256_sad_epu8(tr1, t0);
		t1 = _mm256_sad_epu8(tr1, t1);
		t2 = _mm256_sad_epu8(tr1, t2);


		/// fase guardar los ultimos 8 bits de cada registro (mejorable)           
		char a[32];
		_mm256_storeu_si256((__m256i *)(&a[0]), t0);
		vx = a[31];
		_mm256_storeu_si256((__m256i *)(&a[0]), t1);
		vy = a[31];
		_mm256_storeu_si256((__m256i *)(&a[0]), t2);
		vz = a[31];


		sum = vx + vy + vz + 1;
		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;

	}


	for (i = N - Klen; i < N; i++) {
		unsigned char vx = 0, vy = 0, vz = 0; unsigned sum;

		int k = 0;

		T = &I[i];
		while (i + k < N)
		{
			/*vx = (K[k].x ) ? vx + I[i + k].x : vx;
			vy = (K[k].y ) ? vy + I[i + k].y : vy;
			vz = (K[k].z ) ? vz + I[i + k].z : vz;*/

			vx += (K[k].x & T[k].x);
			vy += (K[k].y & T[k].y);
			vz += (K[k].z & T[k].z);

			k++;
		}
		sum = vx + vy + vz + 1;

		I[i].x = vx * 31 / sum;
		I[i].y = vy * 31 / sum;
		I[i].z = vz * 31 / sum;

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
	if (N < Klen) {
		for (i = 0; i < Rep; i++)
		{
			TransfImageMenor(image, N, K, Klen);
			int ii;
			ii = myRandom() % N; image[ii].x = myRandom() % 31;
			ii = myRandom() % N; image[ii].y = myRandom() % 31;
			ii = myRandom() % N; image[ii].z = myRandom() % 31;
		}
	}
	else {
		for (i = 0; i < Rep; i++)
		{
			TransfImage(image, N, K, Klen);
			int ii;
			ii = myRandom() % N; image[ii].x = myRandom() % 31;
			ii = myRandom() % N; image[ii].y = myRandom() % 31;
			ii = myRandom() % N; image[ii].z = myRandom() % 31;


		}
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

/*

__forceinline __m128 _mm_hadd4_ps(__m128 i)
{ __m128 t; t = _mm_movehl_ps(t, i);
i = _mm_add_ps(i, t);
t = _mm_shuffle_ps(i, i, 0x55);
i = _mm_add_ps(i, t); return i; }*/