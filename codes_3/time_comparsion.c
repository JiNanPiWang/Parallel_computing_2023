//
// Created by Administrator on 2023/7/14.
//
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define N 100000000
int a[N];

int main() {
    int i;
    double sum = 0.0;
    clock_t start, end;
    double cpu_time_used;

    // Serial loop
    start = clock();
    for (i = 0; i < N; i++) {
         a[i] = i;
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Serial loop time: %f seconds\n", cpu_time_used);

    // Parallel loop using OpenMP
    sum = 0.0; // Reset sum
    start = clock();
    #pragma omp parallel for
    for (i = 0; i < N; i++) {
         a[i] = i;
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Parallel loop time: %f seconds\n", cpu_time_used);

    return 0;
}
