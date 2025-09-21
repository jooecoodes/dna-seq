#include <iostream>
#include <omp.h>

int main() {
    omp_set_num_threads(4);  // Set number of threads BEFORE parallel region
    
    int sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < 100; i++) {
        std::cout << i;
        int thread_id = omp_get_thread_num();
        std::cout << "Thread " << thread_id << " processes iteration " << i << std::endl;
    }
    
    return 0;
}