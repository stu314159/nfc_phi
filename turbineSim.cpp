#include <mpi.h>
#include <omp.h>

#include <iostream>
#include "TurbineChannel3D.h"


#ifdef CRAYPAT
#include "pat_api.h"
#endif

using namespace std;

int main(int argc, char * argv[])
{
     #ifdef CRAYPAT
     PAT_record(PAT_STATE_OFF);
     #endif
    const char * paramFN;
    const char * obstFN;
    int rank, size;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
          
    TurbineChannel3D pp(rank, size);
    
    if(rank==0){
        cout<< "Problem initialized!" << endl;
    }

    // write out bc data
    pp.write_bc_arrays(MPI_COMM_WORLD);    
    
    if(rank==0){
        cout<< "BC arrays written" << endl;
    }

    #ifdef CRAYPAT
    PAT_record(PAT_STATE_ON);
    #endif
    
    int* inl = pp.inl;
    int* onl = pp.onl;
    int* snl = pp.snl;
    float* u_bc = pp.u_bc;
    int nnodes = pp.nnodes;
    
    int * Mspeeds = pp.Mspeeds;
    int * Pspeeds = pp.Pspeeds;
    
    int numPspeeds = pp.numPspeeds;
    int numMspeeds = pp.numMspeeds;
    
    int numSpd = pp.numSpd;
    float * fEven = pp.fEven;
    float * fOdd = pp.fOdd;
    
    if(rank==0){
        cout<< "Writing initial data" << endl;
    }
        // write initial data
        // (data processing script will be expecting it)
        pp.write_data(MPI_COMM_WORLD,true);
        

        if (rank == 0){
          cout << "Entering time stepping loop..." << endl;
        }
        double time_start, time_end, ex_time, LPU_sec, gNumLP;
        time_start = MPI_Wtime();
        for(int ts = 0; ts<pp.Num_ts;ts++){
            // say something comforting about the problem progress
            if((ts+1)%(pp.ts_rep_freq)==0){
                if(rank==0){
                    cout << "Executing time step number " << ts+1 << endl;
                }
                
            }
            pp.take_lbm_timestep(ts%2==0,MPI_COMM_WORLD); // weird function call sig.
            MPI_Barrier(MPI_COMM_WORLD); // make sure all processes complete this time step.
            if((ts+1)%(pp.plot_freq)==0){
                // write data at requested intervals.
                pp.write_data(MPI_COMM_WORLD,ts%2==0);
            }
        }
        time_end = MPI_Wtime();
        
        if(rank==0){
            ex_time = time_end - time_start;
            gNumLP = pp.Nx*pp.Ny*pp.Nz;
            LPU_sec = ((double)gNumLP*(double)pp.Num_ts)/ex_time;
            cout << "Estimiated LPU/sec = " << LPU_sec << endl;
        }
        
  
    
    #ifdef CRAYPAT
    PAT_record(PAT_STATE_OFF);
    #endif
    
    MPI_Finalize();
    return 0;
}
