
/**
 * description:
 * 	Genetic algorithm for finding heuristic solution of Traveller Salesman Problem
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


 //  Define Constants: Size of chromosome population, number of cities and size of elite
#define POPSIZE	    40000
#define NUM_CITIES  250
#define ELITSIZE    ((int)(POPSIZE*0.1f))

// Variable used to generate pseudo-random numbers
unsigned int seed;


// Function to generate pseudo-random numbers
inline int myRandom() {
    seed = (214013 * seed + 2531011);
    return (seed >> 13);
}

// Type definition of a point
typedef struct
{
    int x, y;
} point;


// Structure base of chromosome. A chromosome is interpreted as a path
typedef struct
{
    float distance;
    unsigned char * gen;
} chrom;


// Array containing the positions of the cities
point Cities[NUM_CITIES];


//Matrix of all gens 
unsigned char gen[POPSIZE][NUM_CITIES];



//Arrays of chroms, refers to matrix of gens
chrom Population[POPSIZE];



// Helping structure to account for cities that has been already visited
unsigned char MASK[NUM_CITIES];


//Matrix of distances to avoid square
float distances[NUM_CITIES][NUM_CITIES];

//Array to save orderen random numbers
unsigned short Random[2][POPSIZE];
unsigned char Random1[3][POPSIZE];

// Generate random positions for the cities
void Random_Cities()
{
    int i;
    for (i = 0; i < NUM_CITIES; i++)
    {
        Cities[i].x = myRandom() % 4096;
        Cities[i].y = myRandom() % 4096;
    }
}


// Distance between two points
float distance(const point P1, const point P2)
{
    return sqrtf((P1.x - P2.x) * (P1.x - P2.x) +
        (P1.y - P2.y) * (P1.y - P2.y));
}


void InitCityDistance() {

    for (int i = 0; i < NUM_CITIES; i++)
        for (int j = 0; j < NUM_CITIES; j++)
            distances[i][j] = distance(Cities[i], Cities[j]);


}

// Initialize a population of N chromosomes
void init_population(chrom* P, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        P[i].distance = 0;
        P[i].gen = gen[i];
        for (j = 0; j < NUM_CITIES; j++)
            gen[i][j] = j;
    }
}




// A path is a permutation of cities
// Calculate the total distance of a path
float compute_path_distance(unsigned char* path)
{
    int i;
    float dist = 0.0f;

    for (i = 1; i < NUM_CITIES; i++)
        dist = dist + distances[path[i - 1]][path[i]];

    return dist + distances[path[NUM_CITIES - 1]][path[0]];
}


//  Calculate each individual fitness in population. Fitness is path distance
void compute_fitness(chrom* P, int N)
{

    for (int i = 0; i < N; i++)
        P[i].distance = compute_path_distance(P[i].gen);

}


// Display path into screen
void print_path(unsigned char* path)
{
    int i;
    float dist = 0.0f;
    for (i = 1; i < NUM_CITIES; i++)
    {
        printf("%d,", path[i - 1]);
        dist = dist + distance(Cities[path[i - 1]], Cities[path[i]]);
    }
    printf("%d\nTotal Distance: %f\n", path[i - 1],
        dist + distance(Cities[path[i - 1]], Cities[path[0]]));
}


// Merge two sorted arrays, A with szA chromosomes and B with szB chromosomes, into a sorted array C
void merge(chrom* A, int szA, chrom* B, int szB, chrom* C)
{
    int i = 0, j = 0;
    while (i + j < szA + szB)
    {
        if (j == szB || ((i < szA) && (A[i].distance <= B[j].distance)))
        {
            C[i + j] = A[i];  // copy full struct
            i++;
        }
        else
        {
            C[i + j] = B[j];  // copy full struct
            j++;
        }
    }
}


// Sort array A with n chromosomes using recursive merge-sort algorithm
void merge_sort(chrom* A, int n)
{
    int i, n1, n2;
    chrom* A1, * A2;

    if (n < 2)
        return;   // the array is sorted when n=1
    if (n == 2)
    {
        if (A[0].distance > A[1].distance)
        {  // swap values of A[0] and A[1]
            chrom T = A[1];
            A[1] = A[0];
            A[0] = T;
        }
        return; // elements sorted
    }

    // divide A into two arrays, A1 and A2
    n1 = n / 2;   // number of elements in A1
    n2 = n - n1;  // number of elements in A2
    A1 = (chrom*)malloc(sizeof(chrom) * n1);
    A2 = (chrom*)malloc(sizeof(chrom) * n2);

    // move first n/2 elements to A1
    for (i = 0; i < n1; i++) {
        A1[i] = A[i]; // copy full entry
    }
    // move the rest to A2
    for (i = 0; i < n2; i++) {
        A2[i] = A[i + n1]; // copy full entry
    }

    // recursive calls
    merge_sort(A1, n1);
    merge_sort(A2, n2);

    // merge
    merge(A1, n1, A2, n2, A);

    // free allocated memory
    free(A1);
    free(A2);
}


// copy input population to output population   
void copy_population(chrom* P_in, chrom* P_out, int N)
{
    int i, j;
    for (i = 0; i < N; i++) {
        P_out[i].distance = P_in[i].distance;
        for (j = 0; j < NUM_CITIES; j++)
            P_out[i].gen[j] = P_in[i].gen[j];
    }

}


// Checks is a path is valid: does not contain repeats
int check_valid(unsigned char* IN)
{
    int i;

    // clear mask
    for (i = 0; i < NUM_CITIES; i++)
        MASK[i] = 0;

    // check if city has been already visited, otherwise insert city in mask
    for (i = 0; i < NUM_CITIES; i++)
        if (MASK[IN[i]] == 0)
            MASK[IN[i]] = 1;
        else
            return 0;

    return 1;
}


// mate randomly the elite population in P_in into P_out
void mate(chrom* P_in, int Nin, chrom* P_out, int Nout)
{


    // mate the elite population to generate new genes

    #pragma omp parallel for  num_threads(4) proc_bind(close) shared(P_out, P_in)
    for (int m = 0; m < Nout; m++)
    {
        char MASK1[NUM_CITIES] = { 0 };
        // Create new gene in Output population by mating to genes from the elite input population
        // select two random genes from elite population and mate them at random position pos

        int i1, i2, pos;
        //#pragma omp critical

        i1 = Random[0][m];
        i2 = Random[1][m];
        pos = Random1[0][m];

        unsigned char city;

         for (int i = 0; i < pos; i++){
             P_out[m].gen[i] = P_in[i1].gen[i];

               city = P_out[m].gen[i]; 
              MASK1[city] = 1;}



          
        // copy cities in input gene i2 to last part of output gene, 
        //    maintaining the ordering in gene i2
        // copy those cities that are not in the first part of gene i1

        int j = 0; // points to the consecutive positions in gen i2


        for (int i = pos; i < NUM_CITIES; i++)
        {
            unsigned char city;
            do { // skip cities in gen i2 already visited
                city = P_in[i2].gen[j];
                j++;
            } while (MASK1[city] == 1);

            MASK1[city] = 1;         // mark city as seen
            P_out[m].gen[i] = city; // copy city to output gene
        }

        int apos = Random1[1][m];
        int bpos = Random1[2][m];
        unsigned char CityA = P_out[m].gen[apos];
        P_out[m].gen[apos] = P_out[m].gen[bpos];
        P_out[m].gen[bpos] = CityA;

        //compute_path_distance(P_out[m].gen);

      /*  float dist = 0.0f;
        
        unsigned char * __restrict ptr1 = P_out[m].gen;
        
        #pragma acc kernels
        for (int i = 1; i < NUM_CITIES; i++)
            dist = dist + distances[ptr1[i - 1]][ptr1[i]];

        P_out[m].distance = dist + distances[P_out[m].gen[NUM_CITIES - 1]][P_out[m].gen[0]];*/

    }
}


// mutate population: swap cities from two random positions in chromosome
void mutate(chrom* Pop, int N)
{
    int m;

    for (m = 0; m < N; m++)
    { // generate 2 random positions to swap
        int apos = myRandom() % NUM_CITIES;
        int bpos = myRandom() % NUM_CITIES;
        int CityA = Pop[m].gen[apos];
        Pop[m].gen[apos] = Pop[m].gen[bpos];
        Pop[m].gen[bpos] = CityA;
    }
}
void init_random_array(int N) {


    for (int i = 0; i < N; i++) {
        Random[0][i] = myRandom() % ELITSIZE;
        Random[1][i] = myRandom() % ELITSIZE;
        Random1[0][i] = myRandom() % NUM_CITIES;
    }

    for (int i = 0; i < N; i++) {
        Random1[1][i] = myRandom() % NUM_CITIES;
        Random1[2][i] = myRandom() % NUM_CITIES;

    }

}


int main(int argc, char** argv)
{
    int i, EPOCHS = 5;
    seed = 12345;

    // obtain parameters at run time
    if (argc > 1) { EPOCHS = atoi(argv[1]); }
    if (argc > 2) { seed = atoi(argv[2]); }

    printf("Find shortest path for %d cities. %d Epochs. Population Size: %d\n",
        NUM_CITIES, EPOCHS, POPSIZE);


    Random_Cities();
    InitCityDistance();
    init_population(Population, POPSIZE);
   


    // generate random mutations into initial population
    for (i = 0; i < 10; i++)
        mutate(Population, POPSIZE);


    compute_fitness(Population, POPSIZE);
    merge_sort(Population, POPSIZE);


    // generate new populations from initial population
    #pragma acc data copyin(distances[0:NUM_CITIES][0:NUM_CITIES])
    {
    for (i = 0; i < EPOCHS; i++)
    {
        init_random_array(POPSIZE - ELITSIZE);
        mate(Population, ELITSIZE, Population + ELITSIZE, POPSIZE - ELITSIZE); // LOOP FUSION MATE, MUTATE
        
       
        // COMPUTE_PATH_DISTANCE
    
    
        //  #pragma acc kernels
         for (int i = 0; i < POPSIZE; i++){
            float dist = 0.0f;
           #pragma acc kernels
           for (int j = 1; j < NUM_CITIES; j++)  
              dist = dist + distances[Population[i].gen[j - 1]][Population[i].gen[j]];
           Population[i].distance = dist+distances[Population[i].gen[248]][Population[i].gen[249]];}
        
        // END COMPUTE PATH DISTANCE

        merge_sort(Population, POPSIZE);  // sort population by lower fitness, to generate new elite


        // display progress
        if (i % 50 == 1)
        {
            // print current best individual
            printf("Fitness: %f\n", Population[0].distance);

            // sanity check
            if (!check_valid(Population[0].gen))
            {
                printf("ERROR: gen is not a valid permutation of Cities");
                 exit(1);
            }
        }
    }
}
    // print final result
    print_path(Population[0].gen);
    // sanity check
    if (!check_valid(Population[0].gen))
    {
        printf("ERROR: gen is not a valid permutation of Cities");
        exit(1);
    }

    return 0;
}


