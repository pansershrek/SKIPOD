#include <stdio.h>
#include <math.h>
//#include <mpi.h>
//#include <omp.h>

// We expect N <= M
#define N 3
#define M 3
#define EPS 0.0001

double matrix[N][M];

int num_no_zero;
int ind_no_zero;
int ind_pos;
int lines[N];
int ans = 0;

void matrix_init() {
    int ind = 1;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            //scanf("%lf", &matrix[i][j]);
            matrix[i][j] = i*j;
            ind++;
        }
    }
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
                return;
            }
        }
    }
    ind_no_zero = -1;
    ind_pos = -1;
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
    return;
}

int main(int argc, char *argv[]) {
    /*int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);

    #pragma omp parallel
    printf("(%04d, %d): Hello, world!\n", rank, omp_get_thread_num());

    printf("!!!!!!!!!!!!!!!!!\n");
    MPI_Finalize(); */

    matrix_init();
    matrix_print();
    num_no_zero = N;
    ind_no_zero = -1;
    ind_pos = -1;

    for (int i = 0; i < N; ++i) {
        lines[i] = i;
    }

    find_no_zero_first();

    if (ind_pos < 0) {
        printf("ANS 0\n");
        return 0;
    }

    int old_num_no_zero = num_no_zero;

    transform_matrix(ind_no_zero, ind_pos);
    matrix_print();
    find_no_zero_all();
    while (old_num_no_zero >= num_no_zero) {
        old_num_no_zero = num_no_zero;
        find_no_zero_first();
        if (ind_pos < 0) {
            break;
        }
        transform_matrix(ind_no_zero, ind_pos);
        matrix_print();
        find_no_zero_all();
    }

    printf("ANS %d\n", ans);
    return 0;
}
