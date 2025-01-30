#include <chrono>
#include <iostream>
#include <math.h>
#include <mpi.h>

#define PI25DT 3.141592653589793238462643

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    double pi, h, sum, x;
    int N;
    double recv_data; // Fixed type

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <N>\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    } else {
        N = std::atoi(argv[1]);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    h = 1.0 / (double) N;
    sum = 0.0;
    int inicio = (N / size) * rank;
    int finalizacion = (N / size) * (rank + 1); // Exclude last value

    std::chrono::duration<double> total_time;
    auto start{std::chrono::steady_clock::now()};

    for (int i = inicio; i < finalizacion; i++) { // Fixed loop condition
        x = h * ((double)i + 0.5); // Correct midpoint
        sum += 4.0 / (1.0 + x * x);
    }
    
    MPI_Reduce(&sum, &recv_data, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    total_time += std::chrono::steady_clock::now() - start;

    if (rank == 0) {
        pi = h * recv_data; // Compute final pi only in rank 0
        printf("Total time: %.5f seconds\n", total_time.count());
        printf("pi is approximately %.16f, Error is %.16f\n",
               pi, fabs(pi - PI25DT));
    }

    MPI_Finalize();
    return 0;
}
