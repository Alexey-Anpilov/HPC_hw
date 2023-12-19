#include <sw/redis++/redis++.h>
#include <iostream>
#include <string>

using namespace sw::redis;


void get_matrix_from_file(char* file_name, double* A, int n){
    FILE* file;
    file=fopen(file_name, "rb");
    fread(A, sizeof(double), n*n, file);
    fclose(file);
} 


int main(int argc, char* argv[]) {
    auto redis = Redis("tcp://127.0.0.1:6380");
    int n = atoi(argv[1]);
    char* file_name = argv[2];
    double* A = new double[n*n];
    get_matrix_from_file(file_name, A, n);
    std::string key, t;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            key =  "A_" + std::to_string(i) + "_" + std::to_string(j);
            t = std::to_string(A[i*n + j]);
            redis.set(key, t);
        }
    }
}