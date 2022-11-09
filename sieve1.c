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

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */


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

   // Some debugging stuff
   // int num_odd = 0;
   // for(i=low_value; i<=high_value; i++) {
   //    if (i%2 == 1) num_odd++; // Count the number of odd numbers
   // }
   // pprintf(id, "low: %i, high: %i, size: %i, num odd: %i\n", low_value, high_value, size, num_odd);

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

   prime = 3; // Multiples of 2 are not included, start with 3

   do {
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
      pprintf(id, "Prime: %i, First: %i\n", prime, first);

      // These lines are more-or-less the same as the original,
      // outside of playing with the indexing when finding the next prime.
      for (i = first; i < size; i += prime) // Increment by prime
         marked[i] = 1;

      if (!id) {
         while (marked[++index]);
         prime = index*2 + 3;
      }

      if (p > 1)
         MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);

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

