#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

// We expect N <= M
#define N 1000
#define M 1000
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
            matrix[i][j] = 1;
            if (i==j)
                matrix[i][j] = 2;
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
    int i,j,flag;
    int lines_num[M] = {0};
    #pragma omp parallel for private(i, j, flag)
    for (i = 0; i < num_no_zero; ++i) {
        flag = 0;
        if (lines[i] == ind_no_zero) {
            continue;
        }
        for (j = 0; j < M; ++j) {
            if (fabs(matrix[lines[i]][j]) > EPS) {
                flag = 1;
                break;
            }
        }
        if (flag) {
            lines_num[i] = 1;
        }
    }

    for ( i = 0; i < M; ++i) {
        if (lines_num[i]) {
            lines[num_no_zero_new++] = lines[i];
        }
    }
    num_no_zero = num_no_zero_new;
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
    double coef;
    int i,j;
    #pragma omp parallel for private(i, j, coef) shared(matrix)
    for (i = 0; i < num_no_zero; ++i) {
        if (lines[i] == ind) {
            continue;
        }
        if (matrix[lines[i]][pos] == 0) {
            continue;
        }
        coef = matrix[lines[i]][pos] / matrix[ind][pos];
        for (j = pos; j < M; ++j) {
            matrix[lines[i]][j] -= matrix[ind][j] * coef;
        }
    }
}


int main(int argc, char *argv[]) {
    matrix_init();
    num_no_zero = N;
    ind_no_zero = -1;
    ind_pos = -1;

    for (int i = 0; i < N; ++i) {
        lines[i] = i;
    }
    double timerOpenMp = omp_get_wtime();
    find_no_zero_first();

    if (ind_pos < 0) {
        printf("ANS 0\n");
    } else {
        int old_num_no_zero = num_no_zero;

        transform_matrix(ind_no_zero, ind_pos);

        find_no_zero_all();

        while (old_num_no_zero >= num_no_zero) {
            old_num_no_zero = num_no_zero;
            find_no_zero_first();
            if (ind_pos < 0) {
                break;
            }

            transform_matrix(ind_no_zero, ind_pos);
            find_no_zero_all();
        }
        printf("ANS %d\n", ans);
    }
    return 0;
}
