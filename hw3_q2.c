#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main() {

	int p_id;				// Processor ID
	int comm_size;			// Number of processors
	int local_i;			
	int global_sum;			// Parallel sum		
	double elapsed_time;	// Execution time
	int correct_sum;		// Sequential sum

	MPI_Init(NULL, NULL);

	// Start program
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
	MPI_Barrier(MPI_COMM_WORLD);

	elapsed_time = -MPI_Wtime();

	// if (argc != 2) {
	// 	if (p_id == 0) printf("Command line: %s <number_of_processors>\n", argv[0]);
	// 	MPI_Finalize();
	// 	exit(1);
	// }


	local_i = p_id + 1;

	MPI_Reduce(&local_i, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


	elapsed_time += MPI_Wtime();
	if (p_id == 0) {
		printf("Sum from 1 to %d: %d\n", comm_size, global_sum);
		correct_sum = (comm_size * (comm_size + 1)) / 2;
		printf("Correct sum: %d\n", correct_sum);
		printf("SUM_SERIES(%d) %10.6f\n", comm_size, elapsed_time);
	}
	MPI_Finalize();
}