//
// Created by david on 2016-10-18.
//

#ifndef WL_NMSPC_MPI_EXTENSIONS_H
#define WL_NMSPC_MPI_EXTENSIONS_H
#include <Eigen/Core>
#include <mpi.h>

namespace mpi {

    template<typename Derived, typename mpitype>
    void __attribute__((hot)) bcast_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int master_id) {
        auto rows = static_cast<int> (arr.rows());
        auto cols = static_cast<int> (arr.cols());
        MPI_Bcast(&rows, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        MPI_Bcast(&cols, 1, MPI_INT, master_id, MPI_COMM_WORLD);
        arr.derived().resize(rows, cols);
        MPI_Bcast(arr.derived().data(), cols * rows, MPI_TYPE, master_id, MPI_COMM_WORLD);
    }
    template<typename Derived, typename mpitype, typename mpicommtype>
    void __attribute__((hot)) bcast_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int master_id, mpicommtype MPI_COMM) {
        auto rows = static_cast<int>(arr.rows());
        auto cols = static_cast<int>(arr.cols());
        MPI_Bcast(&rows, 1, MPI_INT, master_id, MPI_COMM);
        MPI_Bcast(&cols, 1, MPI_INT, master_id, MPI_COMM);
        arr.derived().resize(rows, cols);
        MPI_Bcast(arr.derived().data(), cols * rows, MPI_TYPE, master_id, MPI_COMM);
    }

    template<typename Derived, typename mpitype>
    void __attribute__((hot)) send_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int dest_id) {
        auto rows = static_cast<int>(arr.rows());
        auto cols = static_cast<int>(arr.cols());
        MPI_Send(&rows, 1, MPI_INT, dest_id, dest_id + 1, MPI_COMM_WORLD);
        MPI_Send(&cols, 1, MPI_INT, dest_id, dest_id + 2, MPI_COMM_WORLD);
        MPI_Send(arr.derived().data(), cols * rows, MPI_TYPE, dest_id, dest_id + 3, MPI_COMM_WORLD);
    }

    template<typename Derived, typename mpitype, typename mpicommtype>
    void __attribute__((hot)) send_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int dest_id, mpicommtype MPI_COMM) {
        auto rows = static_cast<int>(arr.rows());
        auto cols = static_cast<int>(arr.cols());
        MPI_Send(&rows, 1, MPI_INT, dest_id, dest_id + 1, MPI_COMM);
        MPI_Send(&cols, 1, MPI_INT, dest_id, dest_id + 2, MPI_COMM);
        MPI_Send(arr.derived().data(), cols * rows, MPI_TYPE, dest_id, dest_id + 3, MPI_COMM);
    }

    template<typename Derived, typename mpitype>
    void __attribute__((hot)) recv_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int src_id) {
        auto rows = static_cast<int>(arr.rows());
        auto cols = static_cast<int>(arr.cols());
        int id;
        MPI_Comm_rank(MPI_COMM_WORLD, &id);
        MPI_Recv(&rows, 1, MPI_INT, src_id, id + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&cols, 1, MPI_INT, src_id, id + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        arr.derived().resize(rows, cols);
        MPI_Recv(arr.derived().data(), cols * rows, MPI_TYPE, src_id, id + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    template<typename Derived, typename mpitype, typename mpicommtype>
    void __attribute__((hot)) recv_dynamic(Eigen::ArrayBase<Derived> &arr, mpitype MPI_TYPE, int src_id, mpicommtype MPI_COMM) {
        auto rows = static_cast<int>(arr.rows());
        auto cols = static_cast<int>(arr.cols());
        int id;
        MPI_Comm_rank(MPI_COMM, &id);
        MPI_Recv(&rows, 1, MPI_INT, src_id, id + 1, MPI_COMM, MPI_STATUS_IGNORE);
        MPI_Recv(&cols, 1, MPI_INT, src_id, id + 2, MPI_COMM, MPI_STATUS_IGNORE);
        arr.derived().resize(rows, cols);
        MPI_Recv(arr.derived().data(), cols * rows, MPI_TYPE, src_id, id + 3, MPI_COMM, MPI_STATUS_IGNORE);
    }
}

#endif // WL_NMSPC_MPI_EXTENSIONS_H
