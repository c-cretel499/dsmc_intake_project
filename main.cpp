
#include <iostream>
#include <string>
#include <vector>
#include "Vec.h"
#include "World.h"
#include "Output.h"
#include "Source.h"
#include "Collisions.h"
#include "Species.h"
#include <thread>

using namespace std;

//create the rnd function
Rnd rnd;

//Simulation using for loop for parallelization
//void perform_simulation(ColdBeamSource source, Species neutrals, World world, DSMC_MEX dsmc_mex, int num_cols, int i_start, int i_end ) {
//
//	for (int i = i_start; i < i_end; i++) {
//		source.sample();		//injects particles
//
//		neutrals.advance();		//integrate velocity and position
//
//		//periodically do collisions, every col_skip time steps
//		num_cols = dsmc_mex.perform();
//
//
//		// periodically save results
//		if (world.getTs() % 100 == 0) {
//			//cout << "ts=" << world->getTs() << ", np=" << neutrals->getNp() << ", num_cols_inst=" << num_cols << endl;
//			//Output::saveVTI(world, neutrals);
//		}
//	}
//
//}


int main(int argc, char* argv[]) {
	// set some simulation parameters
	// CC: define at least 2 of these with commandline arguments
	const double nd0   = 1e21;		// desired injection number density
	const double vel0  = 10000;		// injection velocity,
	const int N_sample = 100;		// number of particles to create

	World world{41,41,81};		// 41x41x81 mesh
	world.setBoundingBox({0,0,0}, {0.4,0.4,1.2});   //set origin and diagonal corner
	// there should be a physical box added here
	world.setNumThreads(8);

	// set time step such that particles move 0.1*dz at vel0
	double dz = world.getDh()[2];	// get cell spacing in z
	world.setTime(0.1*dz/vel0, 3000);			// time step and max number of time steps

	cout<<"Running for "<<world.getMaxTs()<<" "<<world.getDt()<<" time steps"<<endl;

	//world.addSphere({0.2, 0.2, 0.15}, 0.05);
	

	// Lx * Ly
	double3 xd = world.getXd();
	double3 x0 = world.getX0();
	double A_source = (xd[0]-x0[0])*(xd[1]-x0[1]); // source area

	// set macroparticle weight such that we create 100 particles per time step
	double n_dot  = nd0*vel0*A_source;		// number of real particles created per second
	double n_real = n_dot*world.getDt();	// number of real particles created per time step
	double mpw    = n_real/N_sample;		// set mpw such that N_create particles are generated

	// create species container for N2 neutrals
	Species neutrals("N2", 14*Const::AMU, mpw, world);  // container for neutrals
	//Species oxygen("O2", 32 * Const::AMU, mpw, world);

	// create a cold beam source on the z=0 face, samples N_sample particles with velocity vel0
	ColdBeamSource source(neutrals,world,N_sample,vel0); // where is the ColdBeamSource

	// initialize DSMC,
	/* collisions will be performed only once every 10 time steps
	 * this "subcycling" helps reduce run time but is sometimes also helpful when colision rate is too low*/
	DSMC_MEX dsmc_mex(neutrals, world, 2);
	int num_cols = 0;

	// main loop, using a while loop instead of for
	while(world.advanceTime()) {
		//The main function perform simulation goes here 
		source.sample();		//injects particles

		neutrals.advance();		//integrate velocity and position

		//periodically do collisions, every col_skip time steps
		num_cols = dsmc_mex.perform();


		// periodically save results
		if (world.getTs() % 100 == 0) {
			cout << "ts=" << world.getTs() << ", np=" << neutrals.getNp() << ", num_cols_inst=" << num_cols << endl;
			Output::saveVTI(world, neutrals);
		}
	
	}



	//Parallelizing the simulation code
	cout << "Max number of supported cores is " << thread::hardware_concurrency() << endl;
	int max_num_threads = thread::hardware_concurrency();
	int num_threads = 0;
	cout << "Enter the number of threads to run with = " << endl;
	cin >> num_threads;
	cout << "Running with " << num_threads <<" threads" <<endl;

	//vector<double> thread_time;		//To store time for each thread from 1 to Num_Threads

	//vector<thread> threads;

	//for (int thread_count = 0; thread_count < max_num_threads; thread_count++) {
	//	int i_start = thread_count * world.getMaxTs() / max_num_threads;
	//	int i_end = (thread_count+1) * world.getMaxTs() / max_num_threads;

	//	if (i_end > world.getMaxTs()) i_end = world.getMaxTs();

	//	//threads.emplace_back(perform_simulation, source, neutrals, world, dsmc_mex, num_cols, i_start, i_end);
	//	threads.emplace_back(perform_simulation, source ,neutrals, world, dsmc_mex, num_cols, i_start, i_end);
	//}
	//
	////wait for the threads to join
	//for (thread& th : threads) {
	//	th.join();
	//}

	/* grab ending time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;
}
