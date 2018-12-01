#include <stdio.h>
#include <math.h>
#include <mpi.h>
//#include <omp.h>

// We expect N <= M
#define N 100
#define M 100
#define EPS 0.0001

double matrix[N][M];

int num_no_zero;
int ind_no_zero;
int ind_pos;
int lines[N];
int ans = 0;

int rank, ranksize;

void matrix_init() {
    int ind = 1;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            //scanf("%lf", &matrix[i][j]);
            matrix[i][j] = (i*j*j*i) + 8;
            ind++;
        }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void matrix_print() {
    printf("\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}
void find_no_zero_all() {
    int num_no_zero_new = 0;
    for (int i = 0; i < num_no_zero; ++i) {
        int flag = 0;
        if (lines[i] == ind_no_zero) {
            continue;
        }
        for (int j = 0; j < M; ++j) {
            if (fabs(matrix[lines[i]][j]) > EPS) {
                flag = 1;
                break;
            }
        }
        if (flag) {
            //printf("ind %d\n", lines[i]);
            lines[num_no_zero_new++] = lines[i];
        }
    }
    num_no_zero = num_no_zero_new;
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("find_no_zero_all %d\n", num_no_zero);
}

void find_no_zero_first() {
    int pos = ind_pos + 1;
    for (; pos < M; ++pos) {
        for (int i = 0; i < num_no_zero; ++i) {
            if (fabs(matrix[lines[i]][pos]) > EPS) {
                ind_no_zero = lines[i];
                ind_pos = pos;
                ans++;
                //MPI_Barrier(MPI_COMM_WORLD);
                return;
            }
        }
    }
    ind_no_zero = -1;
    ind_pos = -1;
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void transform_matrix(int ind, int pos) {
    //printf("Transform ind %d pos %d\n", ind, pos);
    for (int i = 0; i < num_no_zero; ++i) {
        if (lines[i] == ind) {
            continue;
        }
        if (matrix[lines[i]][pos] == 0) {
            continue;
        }
        double coef = matrix[lines[i]][pos] / matrix[ind][pos];
        for (int j = pos; j < M; ++j) {
            matrix[lines[i]][j] -= matrix[ind][j] * coef;
        }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
    MPI_Barrier(MPI_COMM_WORLD);

    matrix_init();
    if (rank == 0) {
        //matrix_print();
    }
    num_no_zero = N;
    ind_no_zero = -1;
    ind_pos = -1;

    for (int i = 0; i < N; ++i) {
        lines[i] = i;
    }

    find_no_zero_first();

    if (ind_pos < 0) {
        if (rank == 0) {
            printf("ANS 0\n");
        }
    } else {
        int old_num_no_zero = num_no_zero;

        transform_matrix(ind_no_zero, ind_pos);
        if (rank == 0) {
            //matrix_print();
        }
        find_no_zero_all();
        while (old_num_no_zero >= num_no_zero) {
            old_num_no_zero = num_no_zero;
            find_no_zero_first();
            if (ind_pos < 0) {
                break;
            }
            transform_matrix(ind_no_zero, ind_pos);
            if (rank == 0) {
                //matrix_print();
            }
            find_no_zero_all();
        }
        if (rank == 0) {
            printf("ANS %d\n", ans);
        }
    }
    //MPI_Finalize();
    return 0;
}
