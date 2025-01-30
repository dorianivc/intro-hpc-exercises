#include <chrono>
#include <iostream>
#include <math.h>
#include <mpi.h>

#define PI25DT 3.141592653589793238462643

int main(int argc, char** argv){
    //MPI
    // Initialize MPI communicator
    MPI_Init(&argc, &argv);

    // Initialize VTK MPI handler
    vtkNew<vtkMPIController> mpicontr;
    mpicontr->Initialize(&argc, &argv, 1); 
    // 

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // end MPI


    int N, recv_data;

    if (argc < 2) {
      std::cerr << "Usage: " << argv[0]
              << " <N> \n";
      exit(-1);
    } else {
      N = std::atoi(argv[1]);
    }

    double pi, h, sum, x;

    h = 1.0 / (double) N;
    sum = 0.0; 

    std::chrono::duration<double> total_time;
    auto start{std::chrono::steady_clock::now()};

    for (int i = 0; i <= N; i++) {
        x = h*((double)i-0.5);
        sum += 4.0 / (1.0 + x*x);  
    }
    MPI_Reduce(&sum, &recv_data, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    pi = h*sum;
    total_time += std::chrono::steady_clock::now() - start;

    printf("Total time: %.5f seconds\n", total_time.count());    
    printf("pi is approximately %.16f, Error is %.16f\n",
      pi, fabs(pi - PI25DT));

    return 0;
}
