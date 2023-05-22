#include "hw2_output.h"
#include <string.h> 
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/* GLOBAL VARIABLES */
int **matrix1; 
int **matrix2; 
int **matrix3; 
int **matrix4; 
int **resultMatrix1; 
int **resultMatrix2; 
int **finalResultMatrix;
int sumTidCount, mulTidCount, threadCnt; 
pthread_t *sumTids;  
pthread_t *mulTids;
sem_t *semSum1;
sem_t *semSum2;  

/* STDIN reading elements */ 
int readLine = 255; 
char line[255]; 
char *token ; 
int matrixSizes[2][2]; 

/* FUNCTIONS */
void parseFile(int* row, int* column); 
int **allocateMatrix(int width, int height); 
void printMatrix(int **matrix, int width, int height); 
void fillMatrix(int **matrix, int width, int height); 
void matrixRowSum1(int width, int rowNum); 
void matrixRowSum2(int width, int rowNum);
void multiplyRow(int rowNum);
void readInput(); 
void sumMatrixes();
void freeAll(); 
void *mulThread(void *vargp); 
void *sumThread(void *vargp); 
void mainOperation();

int main(){ 
    if (stdin == NULL) {
        printf("enter a standart input file!\n"); 
        exit(EXIT_FAILURE); 
    }
    hw2_init_output();
    readInput();
    mainOperation(); 
    exit(EXIT_SUCCESS); 
    return 0; 
}

void mainOperation(){ 
    int row1, column1, row2, column2; 
    row1    = matrixSizes[0][0]; 
    column1 = matrixSizes[0][1]; 
    row2    = matrixSizes[1][0]; 
    column2 = matrixSizes[1][1]; 

    int mat2Size = row2*column2; 
    resultMatrix1 = allocateMatrix(column1, row1); 
    resultMatrix2 = allocateMatrix(column2, row2); 
    finalResultMatrix = allocateMatrix(column2, row1); 

    sumTidCount = threadCnt - matrixSizes[0][0]; 
    mulTidCount = matrixSizes[0][0];
    sumTids = malloc(sumTidCount*sizeof(pthread_t)); 
    mulTids = malloc(mulTidCount*sizeof(pthread_t)); 
    semSum1 = malloc(sumTidCount*sizeof(sem_t)); 
    semSum2 = malloc(mat2Size*sizeof(sem_t));

    for (int t=0; t<mat2Size; t++){ 
        sem_init(&semSum2[t], 0, 0);
    }
    for (int t=0; t<sumTidCount; t++){ 
        sem_init(&semSum1[t], 0, 0);
    }

    for (int t=0; t<sumTidCount; t++){ 
        pthread_create(&sumTids[t], NULL, sumThread, (void *)&sumTids[t]);
    }

    for (int t=0; t<mulTidCount; t++){ 
        pthread_create(&mulTids[t], NULL, mulThread, (void *)&mulTids[t]);
    }

    for (int t=0; t<sumTidCount; t++){ 
        pthread_join(sumTids[t], NULL);
    }

    for (int t=0; t<mulTidCount; t++){ 
        pthread_join(mulTids[t], NULL);
    }

    printMatrix(finalResultMatrix, column2, row1);
    pthread_exit(NULL); 

}

void *sumThread(void *vargp){
    int *threadId = (int *)vargp; 
    int index; 
    for (int i=0; i<sumTidCount; i++){ 
        int *sumTidID = (int *)sumTids[i];
        if (*threadId == *sumTidID) {
            index = i; 
            break; 
        }
    }
    if (index < matrixSizes[0][0]) { 
        matrixRowSum1(matrixSizes[0][1], index); 
        sem_post(&semSum1[index]); 
    }
    else { 
        index = index - matrixSizes[0][0]; 
        matrixRowSum2(matrixSizes[1][1], index); 
        for (int i=0; i<matrixSizes[1][1]; i++){
            int postIndex = index*matrixSizes[1][1]+i;
            sem_post(&semSum2[postIndex]);
        }
    }
}

void *mulThread(void *vargp){ 
    int *threadId = (int *)vargp; 
    int index; 
    for (int i=0; i<mulTidCount; i++){ 
        int *mulTidID = (int *)mulTids[i];
        if (*threadId == *mulTidID) {
            index = i; 
            break; 
        }
    }
    multiplyRow(index); 
}

void readInput(){
    int row, column;
    threadCnt = 0; 
    parseFile(&row, &column); 
    matrix1 = allocateMatrix(column, row); 
    fillMatrix(matrix1, column, row); 
    matrixSizes[0][0] = row; 
    matrixSizes[0][1] = column; 
    threadCnt += row; 

    parseFile(&row, &column); 
    matrix2 = allocateMatrix(column, row); 
    fillMatrix(matrix2, column, row); 
    threadCnt += column; 

    parseFile(&row, &column); 
    matrix3 = allocateMatrix(column, row);
    fillMatrix(matrix3, column, row); 
    matrixSizes[1][0] = row; 
    matrixSizes[1][1] = column; 
    threadCnt += row; 

    parseFile(&row, &column); 
    matrix4 = allocateMatrix(column, row); 
    fillMatrix(matrix4, column, row); 
}

void parseFile(int* row, int* column){ 
    fgets(line, readLine, stdin); 
    token = strtok(line, " "); 
    *row = atoi(token); 
    token = strtok(NULL, " "); 
    *column = atoi(token); 

}

void matrixRowSum1(int width, int rowNum){ 
    for(int i=0; i<width; i++){
        resultMatrix1[rowNum][i] = matrix1[rowNum][i] + matrix2[rowNum][i]; 
        hw2_write_output(0, rowNum+1, i+1,  resultMatrix1[rowNum][i]); 
    }
}

void matrixRowSum2(int width, int rowNum){ 
    for(int i=0; i<width; i++){
        resultMatrix2[rowNum][i] = matrix3[rowNum][i] + matrix4[rowNum][i]; 
        hw2_write_output(1, rowNum+1, i+1,  resultMatrix2[rowNum][i]); 
    }
}

void sumMatrixes(){
    int row1, column1, row2, column2; 
    row1    = matrixSizes[0][0]; 
    column1 = matrixSizes[0][1]; 
    row2    = matrixSizes[1][0]; 
    column2 = matrixSizes[1][1]; 

    resultMatrix1 = allocateMatrix(column1, row1); 
    resultMatrix2 = allocateMatrix(column2, row2); 

    for(int i=0; i<row1; i++){ 
        matrixRowSum1(column1, i);
    }
    for(int i=0; i<row2; i++){ 
        matrixRowSum2(column2, i);
    }
    printMatrix(resultMatrix1, column1, row1);
    printMatrix(resultMatrix2, column2, row2);
}

int **allocateMatrix(int width, int height){ 
    int **mMatrix = malloc(height*sizeof(int*)); 
    for (int m=0; m<height; ++m){ 
        mMatrix[m] = malloc(width*sizeof(int)); 
    }
    return mMatrix; 
    } 

void fillMatrix(int **matrix, int width, int height){ 
    for (int a=0; a<height; a++){ 
        fgets(line, readLine, stdin); 
        token = strtok(line, " "); 
        int value = atoi(token); 
        matrix[a][0] = value;
        for(int b=1; b<width; b++){
            token = strtok(NULL, " "); 
            value = atoi(token);
            matrix[a][b] = value; 
        }
    }
}


void multiplyRow(int rowNum){
    int column1 = matrixSizes[0][1]; 
    int column2 = matrixSizes[1][1];
    sem_wait(&semSum1[rowNum]);
    for (int i=0; i<column2; i++) { 
        for (int j=0; j<column1; j++){
            int waitIndex = j*column2 + i ; 
            sem_wait(&semSum2[waitIndex]);
            finalResultMatrix[rowNum][i] += resultMatrix1[rowNum][j] * resultMatrix2[j][i];
            sem_post(&semSum2[waitIndex]);
        }
    hw2_write_output(2, rowNum+1, i+1, finalResultMatrix[rowNum][i]);
    }
}

void printMatrix(int **matrix, int width, int height){ 
    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){ 
            printf("%d ", matrix[i][j]); 
        }
        printf("\n"); 
    }
}