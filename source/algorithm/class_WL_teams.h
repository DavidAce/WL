//
// Created by david on 2018-05-26.
//

#ifndef WL_CLASS_WL_TEAMS_H
#define WL_CLASS_WL_TEAMS_H


#include <mpi.h>
#include <memory>

class class_worker;

class class_WL_teams {
private:
    int  num_teams                  = -1;
    bool team_commander             = false;
    bool team_leader                = false;
    bool active                     = false;
    bool team_finished              = false;
    int  team_id                    = -1;
    int  team_rank                  = -1;
    int  team_size                  = -1;

    int  team_leader_rank           = -1;
    int  team_leader_size           = -1;
    int  team_commander_world_id    = -1;


    MPI_Comm MPI_COMM_TEAM  = MPI_COMM_NULL; //Communicator among team members
    MPI_Comm MPI_COMM_LEAD  = MPI_COMM_NULL; //Communicator among team leaders
    class_worker &worker;    //A ref to parent class
public:
    class_WL_teams(class_worker &worker_ref, int num_teams_);
    void sync_teams();
    void setup_teams();
    void setup_comms();
    void set_defaults();

    int get_team_id() const{ return team_id;}
    int get_team_rank() const{ return team_rank;}
    int get_team_size() const{ return team_size;}
    int get_team_commander_world_id() const{ return team_commander_world_id;}

    const MPI_Comm & get_MPI_COMM_TEAM() const {return MPI_COMM_TEAM;}
    const MPI_Comm & get_MPI_COMM_LEAD() const {return MPI_COMM_LEAD;}


    bool is_leader()    const{return team_leader;}
    bool is_commander() const{return team_commander;}
    bool is_finished()  const{return team_finished;}
    bool is_active()    const{return active;}

    template <int on = 0, typename T>
    void debug_print_team_leader (T input){
        if constexpr(on == 1) {
            MPI_Barrier(MPI_COMM_TEAM);
            for (int i = 0; i < num_teams; i++) {
                if (i == team_id and team_leader) {
                    cout << input;
                    cout.flush();
                    std::this_thread::sleep_for(std::chrono::microseconds(100));
                }
                MPI_Barrier(MPI_COMM_TEAM);
            }
        }
    }

    template <int on = 0, typename T>
    void debug_print_team_commander (T input){
        if constexpr(on == 1) {
            MPI_Barrier(MPI_COMM_TEAM);
            if (team_commander) {
                cout << input;
                cout.flush();
                std::this_thread::sleep_for(std::chrono::microseconds(100));
            }
            MPI_Barrier(MPI_COMM_TEAM);
        }
    }

//        ArrayXXi histogram_recv; //Receive histogram from helpers
};


#endif //WL_CLASS_WL_TEAMS_H
