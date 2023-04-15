#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <omp.h>

int mpi_size, mpi_rank;

void print_array(int *data, int length) {
    if (mpi_rank == 0) {
        for (int i = 0; i < length; i++) {
            printf("%d ", data[i]);
        }
        if (mpi_size > 1) {
            int *values = malloc(sizeof(int) * length);
            for (int rank = 1; rank < mpi_size; rank++) {
                MPI_Recv(values, length, MPI_INT, rank, 1, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                for (int i = 0; i < length; i++) {
                    printf("%d ", values[i]);
                }
            }
        }
        printf("\n");
    } else {
        MPI_Send(data, length, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
}

/*
Gets the number of digits of an integer
https://stackoverflow.com/a/3068420/2521647
*/
int int_length(int value) {
    int nDigits = 1;
    if (value != 0)
        nDigits = floor(log10(abs(value))) + 1;
    return nDigits;
}

int index_max(int *data, int length) {
    int max_value = data[0];
    int index = 1;
    for (int i = 1; i < length; i++) {
        if (data[i] > max_value) {
            max_value = data[i];
            index = i;
        }
    }
    return index;
}

int get_digit_at_position(int num, int position) {
    return (int)(num / (pow(10.0, position))) % 10;
}

int *counting_sort(int *data, int length, int position) {
    int counting[10] = {0};
    int *digits = malloc(length * sizeof(int));
    // Will be written to before first read
    // memset(digits, 0, length * sizeof(int));
    int *output = malloc(length * sizeof(int));
    memset(output, 0, length * sizeof(int));

#pragma omp parallel for reduction(+ : counting[ : 10])
    for (int i = 0; i < length; i++) {
        int digit = get_digit_at_position(data[i], position);
        digits[i] = digit;
        ++counting[digit];
    }

    // Update to contain cumulative counts of elements
    for (int i = 1; i < 10; i++) {
        counting[i] += counting[i - 1];
        // printf("Cumulative: %d\n", counting[i]);
    }

    // Not possible: #pragma omp parallel forreduction(- : counting[:10])
    // Start from the back, otherwise would reposition sorted/short elements
    for (int i = length - 1; i >= 0; i--) {
        int new_position = counting[digits[i]] - 1;
        // Reduce count of occurrences of this digit
        --counting[digits[i]];
        output[new_position] = data[i];
    }

    // Copy over to input
    /*
    for (int i = 0; i < length; i++) {
        data[i] = output[i];
    }
    */

    // free(output);
    // TODO: Check if works correctly, check for leaks

    free(data);
    free(digits);
    return output;
}

int *counting_sort_mpi(int *data, int length, int position) {
    int counting[10] = {0};
    int counting_cumulative[10] = {0};
    int counting_cumulative_global[10] = {0};
    int counting_global[10] = {0};
    int counting_global_up_to_rank[10] = {0};

    int *shared_sorting;
    MPI_Win win;

    int *digits = malloc(length * sizeof(int));
    // Will be written to before first read
    // memset(digits, 0, length * sizeof(int));
    int *output = malloc(length * sizeof(int));
    memset(output, 0, length * sizeof(int));

#pragma omp parallel for reduction(+ : counting[ : 10])
    for (int i = 0; i < length; i++) {
        int digit = get_digit_at_position(data[i], position);
        digits[i] = digit;
        ++counting_cumulative[digit];
    }

    // Update to contain cumulative counts of elements
    counting[0] = counting_cumulative[0];
    for (int i = 1; i < 10; i++) {
        counting[i] = counting_cumulative[i];
        counting_cumulative[i] += counting_cumulative[i - 1];
        // printf("Cumulative: %d\n", counting[i]);
    }

    // Not possible: #pragma omp parallel forreduction(- : counting[:10])
    // Start from the back, otherwise would reposition sorted/short elements
    for (int i = length - 1; i >= 0; i--) {
        int new_position = counting_cumulative[digits[i]] - 1;
        // Reduce count of occurrences of this digit
        --counting_cumulative[digits[i]];
        output[new_position] = data[i];
    }

    /*
     * At this point we need to determine how we should redistribute our values
     * to each process
     */
    int *counting_per_rank = malloc(sizeof(int) * 10 * mpi_size);
    MPI_Allgather(counting, 10, MPI_INT, counting_per_rank, 10, MPI_INT,
                  MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        for (int i = 0; i < (10 * mpi_size); i++) {
            printf("[%d]: %d\n", i, counting_per_rank[i]);
        }
    }

    // Count globally
    for (int i = 0; i < 10; i++) {
        for (int rank = 0; rank < mpi_size; rank++) {
            counting_global[i] += counting_per_rank[i + 10 * rank];
        }
        counting_cumulative_global[i] = counting_global[i];
    }

    // Count globally up to rank
    for (int rank = 0; rank <= mpi_rank; rank++) {
        for (int i = 0; i < 10; i++) {
            counting_global_up_to_rank[i] += counting_per_rank[i + 10 * rank];
        }
    }

    // Count globally cumulative
    for (int i = 1; i < 10; i++) {
        counting_cumulative_global[i] += counting_cumulative_global[i - 1];
    }

    /*
    // Redistribute data
    MPI_Win_allocate(length * sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &shared_sorting, &win);
    for (int i = 0; i < length; i++) {
        int new_position = counting_cumulative_global[digits[i]] - 1;
        --counting_cumulative[digits[i]];
        output[new_position] = data[i];
    }
    */

    /*
     * Moves data to a new position, which might be in the memory space of
     * another process
     *
     * Position is calculated like in the following:
     * Index according to global cumulative position - (position of digit of up
     * to globally up to rank - local position) - 1
     */

    MPI_Win_allocate(length * sizeof(int), sizeof(int), MPI_INFO_NULL,
                     MPI_COMM_WORLD, &shared_sorting, &win);
    MPI_Win_fence(0, win);

    for (int i = 0; i < length; i++) {
        int new_position =
            counting_cumulative_global[digits[i]] -
            (counting_global_up_to_rank[digits[i]] - counting[digits[i]]) - 1;
        int move_to_rank = new_position / mpi_size;
        int rank_local_position = new_position % length;
        printf(
            "Rank: %d: [%d]: %d, New Rank: %d, Local position: %d, Value: %d\n",
            mpi_rank, i, new_position, move_to_rank, rank_local_position,
            digits[i]);

        if (move_to_rank == mpi_rank)
            shared_sorting[rank_local_position] = output[i];
        else {
            MPI_Put(&output[i], 1, MPI_INT, move_to_rank, rank_local_position,
                    1, MPI_INT, win);
        }
        // Reduce count of occurrences of this digit
        --counting_cumulative_global[digits[i]];
        // output[new_position] = data[i];

        // MPI_Barrier(MPI_COMM_WORLD);
        // print_array(shared_sorting, length);
    }
    MPI_Win_fence(0, win);
    MPI_Barrier(MPI_COMM_WORLD);

    print_array(shared_sorting, length);

    MPI_Win_free(&win);

    free(counting_per_rank);
    free(data);
    free(digits);
    return output;
}

int *radix_sort_basic(int *data, int length) {
    // int length = 5;
    // int data[] = {2, 40, 1099, 99, 782};

    /*
    for (int i = 0; i < length; i++) {
        printf("Input: %d, Digits: %d\n", data[i], int_length(data[i]));
    }
    */
    int index_max_element = index_max(data, length);
    int max_digits = int_length(data[index_max_element]);
    // printf("Max index: %d, Value: %d, max_digits: %d\n",
    // index_max_element,
    //       data[index_max_element], max_digits);

    for (int digit = 0; digit < max_digits; digit++) {
        data = counting_sort(data, length, digit);
    }
    return data;
}

int *radix_sort_mpi(int *data, int length) {
    // Local max
    int index_max_element = index_max(data, length);
    int max_digits = int_length(data[index_max_element]);

    // Global max
    int global_max_digits;
    MPI_Allreduce(&max_digits, &global_max_digits, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (mpi_rank == 1)
        printf("Max global: %d\n", global_max_digits);

    for (int digit = 0; digit < global_max_digits; digit++) {
        data = counting_sort_mpi(data, length, digit);
    }
    return data;
}

void random_data(int *data, int lower, int upper, int count, int startIdx) {
    for (int i = 0; i < count; i++) {
        data[startIdx + i] = (rand() % (upper - lower + 1)) + lower;
    }
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    srand(mpi_rank);
    double itime, ftime, exec_time = 0.0;

    if (argc < 4) {
        printf("Missing parameters, show help\n");
        return 1;
    }
    // 0: Prints, 1: No debug prints
    int debug = atoi(argv[1]);
    // Number of values to sort
    int N = atoi(argv[2]);
    /*
     * Variant
     * 0: no parallelization, 1: OpenMP, 2: MPI, 3: Hybrid
     */
    int variant = atoi(argv[3]);
    // Overwrite independet of environment variable
    if (variant == 0 || variant == 2) {
        omp_set_num_threads(1);
    }
    if (mpi_rank == 0) {
        printf("Running with %d OpenMP threads and %d MPI processes\n",
               omp_get_max_threads(), mpi_size);
    }

    int *data = malloc(N * sizeof(int));
    if (data == NULL) {
        printf("Memory not allocated.\n");
        exit(0);
    }
    random_data(data, 0, N, N, 0);

    if (debug) {
        if (mpi_rank == 0)
            printf("Unsorted:\n");
        print_array(data, N);
    }

    // Non-MPI version
    if (mpi_size == 1) {
        itime = omp_get_wtime();
        data = radix_sort_basic(data, N);
        ftime = omp_get_wtime();
        exec_time = ftime - itime;
    } else { // MPI version
             // TODO: Time measurement
        data = radix_sort_mpi(data, N);
    }

    printf("Elapsed time: %f\n", exec_time);

    if (debug) {
        if (mpi_rank == 0)
            printf("Sorted:\n");
        print_array(data, N);
    }

    free(data);
    MPI_Finalize();
    return 0;
}
