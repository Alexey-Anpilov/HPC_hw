#include <omp.h>
#include <iostream>
#include <string.h>

// вывод матрицы
void print_matrix(double* A, int n){
    for(int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cout << A[i*n + j] << " ";
        }
        std::cout << std::endl;
    }
}


void LU_Decomposition(double *A, double* L, double* U, int n, int N){
    double a_11 = *A;
    L[0] = 1;
    U[0] = a_11;
    
    for(int i = 1; i < n; i++){
        U[i] = A[i];
        L[i * N] = A[i * N]/a_11;   
    }

    #pragma omp parallel for collapse(2)
    for(int i = 1; i < n; i++){
        for(int j = 1; j < n; j++){
            A[i*N +j] = A[i*N + j] - L[i*N] * U[j];
        }
    }
    
    if (n == 1){
        return;
    } else{
        LU_Decomposition(A + N + 1, L + N + 1, U + N + 1, n - 1, N);
    }
}

// считывание матрицы из файла
void get_matrix_from_file(char* file_name, double* A, int n){
    FILE* file;
    file=fopen(file_name, "rb");
    fread(A, sizeof(double), n*n, file);
    fclose(file);
} 

int main(int argc, char* argv[]) {
    // обработка аргументов
    int n = atoi(argv[1]);
    char* file_name;
    file_name = argv[2];

    double* A = new double[n*n];
    double* L = new double[n*n];
    double* U = new double[n*n];
    get_matrix_from_file(file_name, A, n);

    double t_start, t_stop;
    double* T = new double[n*n];
    FILE *file = fopen("performance.txt", "w");

    for (int i = 1; i < 30; i++){
        memcpy(T, A, sizeof(double)*n*n);
        omp_set_num_threads(i);
        t_start = omp_get_wtime();
        LU_Decomposition(T, L, U, n, n);
        t_stop = omp_get_wtime();
        fprintf(file, "%d %f\n", i, t_stop-t_start);
    }
    fclose(file);
    return 0;
}