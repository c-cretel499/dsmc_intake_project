#ifndef _SOURCE_H
#define _SOURCE_H

#include "Species.h"
#include "Vec.h"

// samples N_sample particles on the zmin face with uniform velocity
class ColdBeamSource {
public:
// update vel0 to be a double3 not a double, should include thermal velocity
	ColdBeamSource(Species &sp, World &world, int N_sample, double vel0) :
		sp{sp}, world{world}, N_sample{N_sample}, vel0{vel0} {
			x0 = world.getX0()[0];
			y0 = world.getX0()[1];
			Lx = world.getXd()[0]-x0;
			Ly = world.getXd()[1]-y0;
		}

	// samples N_sample randomly positioned particles
	void sample(double T) {
		for (int p = 0; p < N_sample; p++) {
			double3 pos { x0 + rnd() * Lx,  y0 + rnd() * Ly,  0};  // random x, y, z=0
			double3 vel = {0, 0, vel0}; // no problems
			sp.addParticle(pos,vel);
		}
	}
protected:
	Species &sp;
	World &world;
	int N_sample;
	double vel0;
	double T;
	double x0, y0;		// origin
	double Lx, Ly;		// domain length
};
/* WarmBeamSource accounts for thermal velocities of the particles */
class WarmBeamSource{
public:
	WarmBeamSource(Species &sp, World &world, int N_sample, double vel0) :
		sp{sp}, world{world}, N_sample{N_sample}, vel0{vel0} {
			x0 = world.getX0()[0];
			y0 = world.getX0()[1];
			Lx = world.getXd()[0]-x0;
			Ly = world.getXd()[1]-y0;
		}
	void sample(double T) {
		for(int p=0; p<N_sample; p++) {
			double3 pos{x0 + rnd()*Lx, y0 + rnd()*Lx, 0};
			// set up random variables to sample vth in the positive and negative
			double a;
			double b;
			double dir = rnd();
			
			if (dir > 0.5) { a = -1; }
			else { a = 1; } 
			
			dir = rnd();
			if (dir > 0.5) { b = -1; }
			else { b = 1; }
			
			double3 vel =  {a * sp.sampleVth(T), b * sp.sampleVth(T), vel0};
			sp.addParticle(pos, vel);
		}
		
	}
protected:
	Species &sp;
	World &world;
	int N_sample;
	double vel0;		// CC:vel0 should be a double3
	double x0, y0;		// origin
	double Lx, Ly;		// domain length
	double T;			// particle temperature K
	
};
#endif
