#include <iostream>
#include <chrono>
#include <mpi.h>

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  int Nrows; 
  int Ncols;
   int rank, size;
   int* recv_data;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
            << " <Nrows> <Ncols>\n";
    exit(-1);
  } else {
    Nrows = atoi(argv[1]);
    Ncols = atoi(argv[2]);
  }

  int *A = new int[Nrows*Ncols];
  int *x = new int[Ncols];
  int *b = new int[Ncols];
MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for (int i = 0; i < Ncols; i++) {
      x[i] = i + 1;
  }    

  for (int j = 0; j < Nrows; j++) {
      b[j] = 0;
      for (int i = 0; i < Ncols; i++) {
          A[j*Ncols+i] = 1;
      } 
  }

  // compute A*x
  std::chrono::duration<double> total_time;
  auto start{std::chrono::steady_clock::now()};
  for (int j = 0; j < Nrows; j++) {
      for (int i = 0; i < Ncols; i++){
          b[j] += A[j*Ncols+i]*x[i];
      }
  }
MPI_Gather(&x, 1, MPI_INT, recv_data, 1, MPI_INT, 0, MPI_COMM_WORLD);
total_time += std::chrono::steady_clock::now() - start;
if (rank == 0) {
        std::cout << "Datos recolectados en el proceso raíz: ";
        for (int i = 0; i < size; i++) {
            std::cout << recv_data[i] << " ";
        }
        std::cout << std::endl;
        delete[] recv_data; // Liberar memoria
    }



  printf("Total time: %.8f seconds\n", total_time.count());    
  //printf("Result:\n");
  //for (int i = 0; i < Nrows; i++) {
  //    printf("%d ", b[i]);
  //}
  //printf("\n");
  MPI_Finalize();
  return 0; 
}