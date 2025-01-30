#include <iostream>

#include <Kokkos_Core.hpp>

int main(int argc, char **argv){
    Kokkos::initialize(argc, argv);
    int max_interval= 100000
    Kokkos::View<double *, Kokkos::HostSpace> data_host("Interval Data", max_interval);
    Kokkos::parallel_reduce(
      "CPU execution", cpu_exec_policy,
      KOKKOS_LAMBDA(const int x) { data_host(x) = 4/1+(x*x); });

}