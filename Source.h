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
	void sample() {
		for (int p = 0; p < N_sample; p++) {
			double3 pos { x0 + rnd() * Lx,  y0 + rnd() * Ly,  0};  // random x, y, z=0
			double3 vel {0, 0, vel0}; // this is specifically where vel0 can be updated
			sp.addParticle(pos,vel);
		}
	}
protected:
	Species &sp;
	World &world;
	int N_sample;
	double vel0;		// CC:vel0 should be a double3
	double x0, y0;		// origin
	double Lx, Ly;		// domain length
};

class WarmBeamSource{

public:
	WarmBeamSource(Species &sp, World &world, int N_sample, double vel0, int T0) :
		sp{sp}, world{world}, N_sample{N_sample}, vel0{vel0} {
			x0 = world.getX0()[0];
			y0 = world.getX0()[1];
			Lx = world.getXd()[0]-x0;
			Ly = world.getXd()[1]-y0;
		}
	void sample() {
		for(int p=0; p<N_sample; p++) {
			double3 pos{x0 + 0.99*rnd()*Lx, y0 + 0.99*rnd()*Lx, 0};
			double3 vd = {0, 0, vel0}; // drifting velocity 
			//double3 vth = {0,0,0};
			
			double c = rnd();
			double d;
			double f = rnd();
			double g;
			
			if (c > 0.5) { d = -1; }
			else if (c <= 0.5) {d = 1;} 
			
			if (f > 0.5) { g = -1; }
			else if (f <= 0.5) {g = 1;}
			
			double3 vth =  {g * sp.sampleVth(T0), d * sp.sampleVth(T0), 0};
			
			//double3 vth = {d * 1000, g * 1000, 0};
			double3 vel = vd + vth;
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
	int T0;				// particle temperature K
	
};


#endif
