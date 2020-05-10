
#include <iostream>
#include <string>
#include <vector>
#include "Vec.h"
#include "World.h"
#include "Output.h"
#include "Source.h"
#include "Collisions.h"
#include "Species.h"
/*
In Main I have access to double3 so far. 
*/


using namespace std;

//create the rnd function
Rnd rnd;

int main() {
	// set some simulation parameters
	const double nd0   = 1.6e15;	// desired injection number density
	const double vel0  = 7850;		// injection velocity,
	const int N_sample = 100;		// number of particles to create
	
	// print out simulation parameters
	cout << "inlet particle density = " << nd0 << endl;
	cout << "particle beam velocity = " << vel0 << endl;
	cout << "N_sample = " << N_sample << endl;
	
	// define the world boundaries
	double AR = 10; // aspect ratio
	double h = 0.1; // length of the x, y sides
	World world{41,41,81};		// 41x41x81 mesh
	world.setBoundingBox({0,0,0}, {h,h,h*AR});   //set origin and diagonal corner
	// there should be a physical box added here
	
	// set time step such that particles move 0.1*dz at vel0
	double3 dh = world.getDh();
	world.setTime(0.1*dh[2]/vel0, 7000); // can't access max_ts turns out
	
	cout<<"Running for "<<world.getMaxTs()<<" "<<world.getDt()<<" time steps"<<endl;
	
	// Lx * Ly
	double3 xd = world.getXd();
	double3 x0 = world.getX0();
	double A_source = (xd[0]-x0[0])*(xd[1]-x0[1]); // source area
	cout << "Source area: " << A_source << endl;
	
	// set particle weight such that we create 100 particles per time step
	double n_dot = nd0 * vel0 * A_source;	// number of real particles created per second 
	double n_real = n_dot * world.getDt();	// number of real particles created per time step
	double mpw = n_real / N_sample; 		// set mpw such that N_create particles are generated
	
	// create species container for N2 neutrals 
	//Species neutrals("N2", 28 * Const::AMU, mpw, world);
	Species neutrals("O2", 32 * Const::AMU, mpw, world);
	
	// WarmBeamSource is from Source.h 
	//ColdBeamSource source(neutrals, world, N_sample, vel0);
	WarmBeamSource source(neutrals, world, N_sample, vel0, 12000);
	
	DSMC_MEX dsmc_mex(neutrals, world, 2); // 2 is instances to skip
	
	int num_cols = 0; // what is this... it gets updated with dsmc_mex.perform
	
	//main loop, using a while loop instead of a for loop
	while(world.advanceTime()){
		// sample particles, injects particles, only for the first 6000 timesteps
		source.sample();
		
		// integrate velocity and position
		neutrals.advance(); 
		
		// periodically save results
		if (world.getTs()%100==0) {
			cout << "ts=" <<world.getTs()<<", np=" << neutrals.getNp()<< ", num_cols_inst=" << num_cols << endl;
			Output::saveVTI(world, neutrals);
		}	
	}
	// grab ending time
	cout << " Simulation took " << world.getWallTime()<<"seconds" <<endl; // the simulation is not getting to this point, its ending early

	
	return 0;
}