#include <chrono>
#include <iostream>
#include <math.h>
#include <Kokkos_Core.hpp>
#define PI25DT 3.141592653589793238462643

using namespace std;
int main(int argc, char** argv){
    
    int N;

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
    Kokkos::initialize(argc, argv);
    {
    Kokkos::View<double *, Kokkos::HostSpace> data_host("Interval Data", N);
    Kokkos::parallel_reduce(
      "Calculo de PI", N,
      KOKKOS_LAMBDA(const int N, double &sum) { 
        int x=h*((double)N-0.5);
        sum += 4.0 / (1.0 + x*x);
        data_host(N) = sum;
        }, sum);
    
    total_time += std::chrono::steady_clock::now() - start;
   
    
   cout <<"Resultado paso a paso "<< endl;
    for(int i=0; i<N,i++;){
    cout << "Para el numero: "<< i<< "el valor es: "<< data_host(i)<<endl;
    }
    }
    pi = h*sum;

    cout << "Pi = "<< pi << endl;    
    
    
    Kokkos::finalize();
    printf("Total time: %.5f seconds\n", total_time.count());    
    printf("pi is approximately %.16f, Error is %.16f\n",
      pi, fabs(pi - PI25DT));
    return 0;
}