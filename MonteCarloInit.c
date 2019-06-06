#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define STATES 20
#define GRID_SIZE 920
#define ITERATIONS 1000000
#define CELL_SIZE 2

typedef long long int LLint;

void printMatrix(LLint a[GRID_SIZE][GRID_SIZE]) {
    LLint i, j;
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            printf("%lld ", a[i][j]);
        }
        printf("\n");
    }
}

void exportMatrixToFile(LLint a[GRID_SIZE][GRID_SIZE]) {
    FILE *fp;
    fp = fopen("monteCarlo_init_matrix.csv", "w+");
    int i, j;
    fprintf(fp, "\n");
    for (i = 0; i < GRID_SIZE; i++) {
        for(j = 0; j < GRID_SIZE; j++) {
            if (j == GRID_SIZE - 1) {
                fprintf(fp, "%lld", a[i][j]);
                continue;
            }
            fprintf(fp, "%lld,", a[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    printf("\n");
    printf("\nMatrix Values written to file monteCarlo_init_matrix.csv!\n");
}

LLint** monteCarloInit(LLint mat[GRID_SIZE][GRID_SIZE]) {
    LLint i, j;
    LLint* arr = malloc(sizeof(LLint) * STATES);
    LLint index = 0;
    for (j = 0; j < STATES; j++) {
        arr[j] = 0;
    }
    for (i = 0; i < ITERATIONS; i++) {
        LLint x = rand() % GRID_SIZE;
        LLint y = rand() % GRID_SIZE;
        if ( x > 0 && y > 0)        arr[(LLint)mat[x - 1][y - 1] - 1] += 1;
        if (x > 0)                     arr[(LLint)mat[x - 1][y] - 1] += 1;
        if(x > 0 && y < GRID_SIZE - 1)   arr[(LLint)mat[x - 1][y + 1] - 1] += 1;
        if(y > 0)                     arr[(LLint)mat[x][y - 1] - 1] += 1;
        if(y < GRID_SIZE - 1)                   arr[(LLint)mat[x][y + 1] - 1] += 1;
        if(x < GRID_SIZE - 1 && y > 0)           arr[(LLint)mat[x + 1][y - 1] - 1] += 1;
        if(x < GRID_SIZE - 1)                   arr[(LLint)mat[x + 1][y] - 1] += 1;
        if(x < GRID_SIZE - 1 && y < GRID_SIZE - 1)         arr[(LLint)mat[x + 1][y + 1] - 1] += 1;
        LLint temp_max = 10000;
        for (j = 0; j < STATES; j++) {
            if (arr[j] > temp_max) {
                temp_max = arr[j];
                index = j + 1;
            }
        }
        mat[x][y] = index;
    }

    return mat;
}

int main(int argc, char* argv[]) {

    srand(time(0));

    LLint matrix[GRID_SIZE][GRID_SIZE];

    LLint i, j;

    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {
            matrix[i][j] = rand() % STATES + 1;
        }
    }
    printMatrix(matrix);
    monteCarloInit(matrix);
    printMatrix(matrix);
    exportMatrixToFile(matrix);

    return 0;
}
