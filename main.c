#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#define ROOT 0

void Print_Matrix(int *a, int n, int m) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            printf("%3i", a[i * m + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {

    int N1 = atoi(argv[1]);
    int N2 = atoi(argv[2]);
    int BLOCK_SIZE = atoi(argv[3]);

    int PROCESS_NUM;
    int PROCESS_RANK;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PROCESS_NUM);
    MPI_Comm_rank(MPI_COMM_WORLD, &PROCESS_RANK);

    MPI_Status status;
    MPI_Request request;

    int ROWS_IN_TAPE = N1 / PROCESS_NUM;
    int ELEMENTS_IN_TAPE = ROWS_IN_TAPE * N2;
    int BLOCKS_IN_TAPE = N2 / BLOCK_SIZE;

    int NEXT_PROC = PROCESS_RANK + 1;
    int PREV_PROC = PROCESS_RANK - 1;

    int *TAPE = (int *) calloc((size_t) ELEMENTS_IN_TAPE, sizeof(int));

    int PREV_TAPE_LAST_ROW[BLOCK_SIZE];
    for (int i = 0; i < BLOCK_SIZE; i++) {
        PREV_TAPE_LAST_ROW[i] = 0;
    }

    for (int block_i = 0; block_i < BLOCKS_IN_TAPE; block_i++) {

        if (PREV_PROC >= ROOT) {
            MPI_Recv(PREV_TAPE_LAST_ROW, BLOCK_SIZE, MPI_INT, PREV_PROC, block_i, MPI_COMM_WORLD, &status);
        }

        for (int i = 0; i < ROWS_IN_TAPE; i++) {
            for (int j = block_i * BLOCK_SIZE, k = 0; j < (block_i + 1) * BLOCK_SIZE; j++, k++) {
                if (i == 0) {
                    TAPE[i * N2 + j] = PREV_TAPE_LAST_ROW[k] + 1;
                } else {
                    TAPE[i * N2 + j] = TAPE[(i - 1) * N2 + j] + 1;
                }
            }
        }

        if (NEXT_PROC < PROCESS_NUM) {
            int *SEND_BLOCK_P = &(TAPE[(ROWS_IN_TAPE - 1) * N2]) + block_i * BLOCK_SIZE;
            MPI_Isend(
                    SEND_BLOCK_P,
                    BLOCK_SIZE,
                    MPI_INT,
                    NEXT_PROC,
                    block_i,
                    MPI_COMM_WORLD,
                    &request
            );
        }
    }


    int *RESULT_MATRIX;
    if (PROCESS_RANK == ROOT) {
        RESULT_MATRIX = (int *) malloc(N1 * N2 * sizeof(int));
    }

    MPI_Gather(TAPE, ELEMENTS_IN_TAPE, MPI_INT, RESULT_MATRIX, ELEMENTS_IN_TAPE, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (PROCESS_RANK == ROOT) {
        Print_Matrix(RESULT_MATRIX, N1, N2);
    }

    MPI_Finalize();
    return 0;
}