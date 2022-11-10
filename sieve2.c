/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

// A prettier way of identifying which processor prints what, can be defined to print nothing when final compilation is needed
// #define pprintf(id, str, ...) printf("Proc %i: " str, id, __VA_ARGS__)
#define pprintf(id, str, ...)

/*
   module load mpich-3.2.1/gcc-4.8.5

   cd build && make -j && cd ../script && ./submit.sh && cd ..
   watch -n 1 squeue
   cd script && ./data.sh && cd ..
*/

/**
 * We want the primes up to and including max_num.
 * ex. if max_num = 11, we expect to get 3, 5, 7, 11 (2 isnt included)
 * and num_primes set to 4
 */
unsigned long int * primes_up_to(int max_num, unsigned long int *num_primes) {
   unsigned long int first;
   unsigned long int size;
   unsigned long int prime;
   unsigned long int i;
   unsigned long int index;

   *num_primes = 0;
   size = (max_num - 1)/2; // eg. max_num = 12, allocate room for: 3,5,7,9,11
   char *marked = (char *) malloc(sizeof(char) * size);
   unsigned long int *out = malloc(sizeof(unsigned long int) * size); // This is technically overprovisioned
   for(i=0; i<size; i++) {
      marked[i] = 0;
   }

   prime = 3;
   index = 0;

   do {
      // pprintf(-1, "Prime: %i\n", prime);
      first = (prime*prime - 3) / 2;

      for(i=first; i<size; i+=prime)
         marked[i] = 1;

      while(marked[++index]);
      prime = index*2 + 3;
   } while(prime * prime <= max_num);

   int c=0;
   for(i=0; i<size; i++) {
      if (c+1 == *num_primes) break;
      if (!marked[i]) {
         (*num_primes)++;
         out[c] = i*2+3;
         c++;
      }
   }

   pprintf(-1, "Num primes: %i, First: %i, Last: %i\n", (*num_primes), out[0], out[(*num_primes)-1]);

   free(marked);

   return out;
}

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   unsigned long int    j;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */


   /* My Code Start */

   low_value = 3 + id * (n - 2) / p;
   if (low_value % 2 == 0) low_value++; // If the low number is an even number, bump it up to make it odd
   high_value = 2 + (id + 1) * (n - 2) / p;
   if (high_value % 2 == 0) high_value--; // If the high number is even, bump it down to make it odd.
   size = (high_value - low_value + 1)/2 + 1; // We divide by 2 to remove the even elements

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n - 1) / p;

   if ((2 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }

   /* Allocate this process's share of the array. */

   marked = (char *) malloc(size);

   if (marked == NULL) {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < size; i++)
      marked[i] = 0;

   if (!id) index = 0;

   unsigned long int num_primes;
   unsigned long int *sprimes = primes_up_to(sqrt(n), &num_primes); // num_primes will be set to the number of primes found, starting at 3

   for(j=0; j<num_primes; j++) {
      prime = sprimes[j];

      if (prime*prime > low_value) { // If we need to start indexing from the middle of the array (as opposed to from around the front)
         first = (prime*prime - low_value) / 2; // Get index mid-array of first prime^2
      } else {
         if (!(low_value % prime)) { // If the low_value is a product of the prime, then were good! Easy case!
            first = 0;
         } else {
            int temp = prime - (low_value % prime); // p, p-1, p-2, ..., 1
            if ((low_value + temp) % 2 == 1) { // If we would end up on our multiple, and it is already odd, then good!
               first = temp/2;
            } else { // Otherwise, if the next multiple and its even, then we need the next one
               first = (temp + prime)/2;
            }
         }
      }

      // These lines are more-or-less the same as the original,
      // outside of playing with the indexing when finding the next prime.
      for (i = first; i < size; i += prime) // Increment by prime
         marked[i] = 1;
   }

   free(sprimes);

   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i])
         count++;

   if (p > 1)
      MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

   global_count++; // We need to account for 2! Poor guy... lost but not forgotten.

   /* My Code End */


   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

