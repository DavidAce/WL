//
// Created by david on 2016-08-11.
//

#ifndef CLASS_WL_PARALLEL_ALGORITHMS_H
#define CLASS_WL_PARALLEL_ALGORITHMS_H
#include <thread>
#include <chrono>
#include "class_tic_toc.h"
#include <Eigen/Dense>
#include <Eigen/Core>


#include "class_WL_worker.h"
//class class_worker;
namespace mpi {
    extern void swap                           (class_worker &) ;
    extern void merge                          (class_worker &, bool broadcast, bool setNaN) ;
    extern void broadcast_merger               (class_worker &) ;
    extern void divide_global_range_dos_area   (class_worker &) ;
    extern void divide_global_range_dos_volume (class_worker &) ;
    extern void check_help                     (class_worker &)__attribute__((hot)) ;
    extern void take_help                      (class_worker &)__attribute__((hot)) ;
    extern void take_help2                      (class_worker &)__attribute__((hot)) ;
    extern void setup_help                     (class_worker &,class_backup &) ;


    template <typename Derived, typename mpitype>
    void __attribute__((hot)) bcast_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE , int master_id){
        int rows = (int) arr.rows();
        int cols = (int) arr.cols();
        MPI_Bcast(&rows, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        MPI_Bcast(&cols, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        arr.derived().resize(rows, cols);
        MPI_Bcast(arr.derived().data(), cols*rows, MPI_TYPE, master_id, MPI_COMM_WORLD);
    }
    template <typename Derived, typename mpitype, typename mpicommtype>
    void __attribute__((hot)) bcast_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE , int master_id, mpicommtype MPI_COMM){
        int rows = (int) arr.rows();
        int cols = (int) arr.cols();
        MPI_Bcast(&rows, 1, MPI_INT, master_id, MPI_COMM);
        MPI_Bcast(&cols, 1, MPI_INT, master_id, MPI_COMM);
        arr.derived().resize(rows, cols);
        MPI_Bcast(arr.derived().data(), cols*rows, MPI_TYPE, master_id, MPI_COMM);
    }



    template <typename Derived, typename mpitype>
    void __attribute__((hot)) send_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE, int dest_id){
        int rows = (int) arr.rows();
        int cols = (int) arr.cols();
        MPI_Send(&rows, 1, MPI_INT, dest_id, dest_id + 1, MPI_COMM_WORLD);
        MPI_Send(&cols, 1, MPI_INT, dest_id, dest_id + 2, MPI_COMM_WORLD);
        MPI_Send(arr.derived().data(), cols*rows, MPI_TYPE, dest_id, dest_id + 3, MPI_COMM_WORLD);
    }

    template <typename Derived, typename mpitype>
    void __attribute__((hot)) recv_dynamic (ArrayBase<Derived> &arr, mpitype MPI_TYPE , int src_id ){
        int rows = (int) arr.rows();
        int cols = (int) arr.cols();
        int id;
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Recv(&rows, 1, MPI_INT, src_id, id + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&cols, 1, MPI_INT, src_id, id + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        arr.derived().resize(rows, cols);
        MPI_Recv(arr.derived().data(),  cols*rows, MPI_TYPE, src_id, id + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    template <typename T>
    void debug_print        (class_worker &worker, T input){
        MPI_Barrier(MPI_COMM_WORLD);
        if (worker.world_ID == 0) {
            cout << input;
            cout.flush();
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


};


#endif //CLASS_WL_PARALLEL_ALGORITHMS_H
