#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

   int      p_id;                      // Processor ID
   int      comm_size;                 // Number of processors
   double   elapsed_time;              // Execution time
   long     n;                         // find prime from 2 to n
   long     low_index;                 // start index of this process
   long     high_index;                // end index of this process
   long     size;                      // number of elements that a process holds
   long     proc0_size;                // array size of process 0
   long     *primes;                   // hold primes number in the end
   long     i;
   long     prime;                     // prime found so far, this is global variable
   long     current_index;             // index of current prime
   long     first;                     // first index to search
   long     local_min_prime;           // hold the minimum local prime
   long     global_max_distance = 0;   // hold the final result


   MPI_Init(&argc, &argv);

   // Start program
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
   MPI_Barrier(MPI_COMM_WORLD);

   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (p_id == 0) printf("Command line: %s <number_of_processors>\n", argv[0]);
      MPI_Finalize();
      exit(1);
   }

   // get user input
   n = atoi(argv[1]);

   // find a way to partition
   low_index = 2 + p_id * (n - 1)/comm_size;
   high_index = 1 + (p_id + 1) * (n  - 1)/comm_size;
   size = high_index - low_index + 1;

   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = (n-1)/comm_size;

   if ((2 + proc0_size) < (long) sqrt((double) n)) {
      if (p_id == 0) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   // printf("(%d) low_index: %ld\n", p_id, low_index);
   // printf("(%d) high_index: %ld\n", p_id, high_index);
   // printf("(%d) size: %ld\n", p_id, size);

   // allocate memory for primes
   primes = (long *) malloc((size + 1) * sizeof(long));

   // // give primes array the number of its process
   for (i = 0; i < size; i++) primes[i] = low_index + i;
   if (p_id == 0) current_index = 0;
   prime = 2;
   do {
      if (prime * prime > low_index)
         first = prime * prime - low_index;
      else {
         if (!(low_index % prime)) first = 0;
         else first = prime - (low_index % prime);
      }
      for (i = first; i < size; i += prime) {
         primes[i] = 0;
      }
      if (p_id == 0) {
         while (!primes[++current_index]);
         prime = current_index + 2;
      }
      // if (comm_size > 1) 
         MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);

   

   // find local_min 
   local_min_prime = 0; // if a process doesn't have any prime, it will send 0
   for (i = 0; i < size; i++) {
      if (primes[i]) {
         local_min_prime = low_index + i;
         break;
      }
   }

   // printf("%d local min prime %d\n", p_id, local_min_prime);

   
   
   // send local_min prime to p_id - 1 process
   if (p_id != 0) 
      MPI_Send(&local_min_prime, 1, MPI_INT, p_id - 1, 0, MPI_COMM_WORLD);
   // receive local_min prime from p_id + 1 process
   if (p_id != (comm_size - 1)) 
      MPI_Recv(&primes[size], 1, MPI_INT, p_id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

   // mark out the primes[size] of last process
   if (p_id == comm_size - 1) primes[size] = 0;

   // printf("%d last prime: %d\n", p_id, primes[size]);

   // find local_max distance
   long local_max_distance = -1;       // hold the max distance in local process
   long last_prime = -1;
   long current_prime = -1;
   for (i = 0; i <= size; i++) {
      if (primes[i]) {
         if (current_prime == -1) {
            current_prime = primes[i];
         } else {
            current_prime = primes[i];
            long new_distance = current_prime - last_prime;
            local_max_distance = new_distance > local_max_distance ? new_distance : local_max_distance;
         }
         last_prime = current_prime;
         // printf("%d prime: %d\n", p_id, primes[i]);
      }
   }


   MPI_Reduce(&local_max_distance, &global_max_distance, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);


   elapsed_time += MPI_Wtime();
   if (p_id == 0) {
      printf("Max distance between two consecutive prime numbers under %ld is %ld\n", n, global_max_distance);
      printf("SIEVE(%d) %10.6f\n", comm_size, elapsed_time);
   }

   // deallocate primes
   free(primes);

   // terminate MPI
   MPI_Finalize();

   return 0;
}