#include <omp.h>
#include <stdio.h>

int main()
{
    int tid;
    int num_threads;

    // 设置使用的线程数为4
    omp_set_num_threads(4);

    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();

        // 只有主线程打印线程数
        #pragma omp master
        {
            num_threads = omp_get_num_threads();
            printf("Number of threads: %d\n", num_threads);
        }

        printf("Hello World from thread = %d\n", tid);
    }

    return 0;
}
