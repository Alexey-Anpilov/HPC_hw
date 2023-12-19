#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <cmath>

/*
    Представлена реализация алгоритма ортогонализации Грамма-Шмидта.
    На вход подается n векторов размерности n. Расчеты разбиваются между m процессами.
    Каждому процессу достается набор из n/m векторов. Последний процесс получает дополнительно n%m векторов.
    Вектора обрабатываются поочередно, после нормировки очередного вектора, он передается остальным процессам,
    которые вычитают из своих векторов, проекции на полученный вектор. Групповая передача/приём производится с помощью функции Bcast. 

*/

double scalar_product(double* a, double*b, int n){
    // расчет скалярно произведения двух векторов
    double res = 0;
    for(int i = 0; i < n; i++){
        res += a[i] * b[i];
    }
    return res;
}

void gram_schmidt(double *M, int n, int group_size)
{
    // реализация процесса Грама-Шмидта
    int proc_amt, proc_num;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_amt);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
    double *buff = new double[n];
    double t;
    for (int i = 0; i < n; ++i)
    {
        // определим каким процессом будет обрабатываться вектор i 
        int group_num = i / (n / proc_amt);
        if (group_num >= proc_amt)
        {
            group_num = proc_amt - 1;
        }
        int loc = i - group_num * (n / proc_amt);

        // нормировка вектора
        if (group_num == proc_num)
        {
            t = std::sqrt(scalar_product(M + loc*n, M + loc*n, n));
            for (int i = 0; i < n; i++)
            {
                M[loc * n + i] /= t;
                buff[i] = M[loc * n + i];
            }
        }

        // передаем вектор остальным процессам
        MPI_Bcast(buff, n, MPI_DOUBLE, group_num, MPI_COMM_WORLD);

        
        if (group_num > proc_num)
        {
            continue;
        }
        if (group_num < proc_num)
        {
            loc = -1;
        }

        // вычитаем проекции
        for (int j = loc + 1; j < group_size; ++j)
        {
            t = scalar_product(buff, M + j * n, n);
            for (int k = 0; k < n; ++k)
            {
                M[j * n + k] -= t * buff[k];
            }
        }
    }
    delete[] buff;
}


void get_matrix_from_file(char* filename, double* M, int n, int proc_num){
    // Заполняет матрицу данными из файла
    MPI_File file;
    MPI_Status status;
    MPI_Offset low_border;
    int count;
    double *buf;
    int proc_amt;
    
    // подготовим отступ для каждого процесса и размер буффера
    MPI_Comm_size(MPI_COMM_WORLD, &proc_amt);
    int bufsize = n / proc_amt * n;
    low_border = bufsize * proc_num;
    if (proc_num == proc_amt-1){
        bufsize += n%proc_amt * n;
    }
    buf = new double[bufsize];

    // считаем данные из файла
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &file);
    MPI_File_read_at(file, low_border * sizeof(double),
                              buf, bufsize,
                              MPI_DOUBLE, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    printf("process %d read %d floats\n", proc_num, count);
    MPI_File_close(&file);
    
    // заполним матрицу считанными в буфер значениями
    for (int i = 0; i < bufsize; ++i){
        M[i] = buf[i];
    }
    delete[] buf;
}


void print_matrix(double* M, int n, int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            std::cout << M[i * m + j] << " ";
        }
        std::cout << std::endl;
    }
}


int main(int argc, char *argv[])
{
    // инициализируем MPI
    MPI_Init(&argc, &argv);
    int proc_amt, proc_num;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_amt);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);

    // обработаем аргументы командной строки
    if (argc < 3){
        std::cout << "usage: GS_mpi <кол-во векторов> <входной файл>";
        MPI_Finalize();
        return -1;
    }
    int n = atoi(argv[1]);    
    char * file_name;
    file_name = argv[2];

    // задаем размер набора векторов для каждого процесса
    int group_size = n / proc_amt;
    if (proc_num == proc_amt - 1)
    {
        group_size += n % proc_amt;
    }
    
    // выделяем память для матрицы и заполняем из файла/случайными значениями
    double *M = new double[group_size * n];
    get_matrix_from_file(file_name, M, n, proc_num);
    //print_matrix(M, group_size, n);

    // замеряем время исполнения
    double t_start; 
    double t_stop;
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();
    gram_schmidt(M, n, group_size);
    MPI_Barrier(MPI_COMM_WORLD);
    t_stop = MPI_Wtime();
    
    
    // Вывод получившейся матрицы
    for(int i = 0; i < proc_amt; i++){
        if (proc_num == i){
            print_matrix(M, group_size, n);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // записываем замеры в файл
    if (proc_num == 0)
    {
        FILE *file = fopen("performance.txt", "a");
        if (file != NULL) {
        fprintf(file, "%d %f\n", proc_amt, t_stop-t_start);
        fclose(file);
        }
    }
    delete[] M;
    MPI_Finalize();
    return 0;
}