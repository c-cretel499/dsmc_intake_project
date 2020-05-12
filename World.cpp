#include <iostream>
#include "World.h"
#include "Vec.h"

/*** Diffuse reflection of a particle off a wall ***/
double3 World::wallDiffuseVector(const double3& x)
{
	// pick angles, theta=off normal, psi = azimuthal rotation
	// https://www.particleincell.com/2015/cosine-distribution/
	double sin_theta = rnd();		
	double cos_theta = sqrt(1 - sin_theta * sin_theta); // cos*cos + sin*sin = 1
	double psi = 2 * Const::PI * rnd();			// random in plane angle 		
	
	// set n depending on what wall you're on
	double3 n;
	double3 t1;
	double near_x = dh[0] * 0.1;
	double near_y = dh[1] * 0.1;
	if ( x[0] - x0[0] <= near_x){ 
		//std::cout << "Hit x=0 wall" << std::endl;
		n = {1,0,0};
		t1 = {0, rnd(), rnd()};
		t1 = unit(t1);
	}
	else if (x[1] - x0[1] <= near_y) { 
		//std::cout << "Hit y=0 wall" << std::endl;
		n = {0,1,0};
		t1 = {rnd(),0,rnd()};
		t1 = unit(t1);
	}
	else if (xd[0]-x[0] <= near_x ) {
		//std::cout << "Hit x=xd wall" << std::endl;
		n = {-1,0,0};
		t1 = {0,rnd(),rnd()};
		t1 = unit(t1);
	}
	else if (xd[1]-x[1] <= near_y ) { 
		//std::cout << "Hit y=yd wall" << std::endl;
		n = {0,-1,0};
		t1 = {rnd(),0,rnd()};
		t1 = unit(t1);
	}
	
	double3 t2 = cross(n, t1); 			//second tangent
	
	return sin_theta * cos(psi) * t1 + sin_theta * sin(psi) * t2 + cos_theta * n;
}



/* Calculates the time step from the current position to the wall */
double World::lineWallIntersect(const double3 &x1, const double3 &x2){
	// check which plane the particle will intersect, determine the plane's normal vector
	// https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
	double3 n;
	double3 p0;
	if (x2[0] <= 0) {
		n  = {1, 0, 0};
		p0 = x0;
		}
	else if (x2[0] >= xd[0]) {
		n  = {-1, 0, 0};
		p0 = xd;
		}
	else if (x2[1] <= 0) {
		n  = {0, 1, 0};
		p0 = x0;
		}
	else if (x2[1] >= xd[1]) {
		n  = {0, -1, 0};
		p0 = xd;
		}
	
	// calculate the line connecting x1 and x2
	double3 B = x2 - x1;
	double3 pw;
	double d;
	// check that B is not parallel to n
	if (dot(B, n) == 0){
		return 0.5; // return the midpoint if the particle is moving along the wall
	}
	else{
		d = dot(p0-x2, n)/dot(B,n);
		pw = x2 + d*B;
	}
	double b = dot(B,B);
	double a = dot(pw-x1, pw-x1);
	double tp = a / b;
	return tp;
}

/* Count how often a particle hits each wall
	hitWall = [+x, -x, +y, -y, +z, -z]*/
	/*
void World::wallColCount(double3 x){
	if (x[0] <= 0) {
		wallHit[0] ++;
		}
	else if (x[0] >= xd[0]) {
		wallHit[1] ++;
		}
	else if (x[1] <= 0) {
		wallHit[2] ++;
		}
	else if (x[1] >= xd[1]) {
		wallHit[3] ++;
		}
	}
*/

/*returns true if point x is inside or on the sphere*/
bool World::inSphere(const double3 &x){
	double3 r = x-sphere_x0;	//ray to x

    double r_mag2 = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if (r_mag2<=(sphere_r*sphere_r)) return true;
    else return false;
}

// not used in this project, left over from previous
double World::lineSphereIntersect(const double3 &x1, const double3 &x2){

	double3 B = x2-x1; 			// trajectory of the particle
	double3 A = x1-sphere_x0;	// vector from sphere center to particle position
	double a = dot(B,B);		// dot(B,B) becomes the law of cosines
	double b = 2*dot(A,B);		// 2 times the dot product of A to B
	double c = dot(A,A)-(sphere_r*sphere_r);
	double det = b*b-4*a*c;     // this is the discriminant
	if (det<0) return 0.5;

	double tp = (-b + sqrt(det))/(2*a);
	if (tp<0 || tp>1.0)
	{
		tp = (-b - sqrt(det))/(2*a);
		if (tp<0 || tp>1.0)
		{
			std::cerr<<"Failed to find a line-sphere intersection!"<<std::endl;
			tp=0.5;	//set as midpoint
		}
	}
	return tp;
}

/*returns random diffuse direction sampled at a given surface point*/
double3 World::sphereDiffuseVector(const double3 &x)
{
	//pick angles, theta=off normal, psi = azimuthal rotation
	double sin_theta = rnd();
	double cos_theta = sqrt(1-sin_theta*sin_theta);
	double psi = 2*Const::PI*rnd();

	double3 n = unit(x-sphere_x0);	// normal vector
	double3 t1;						// create the first tangent
	if (dot(n,{1,0,0})!=0) t1 = cross(n,{1,0,0});
	else t1 = cross(n,{0,1,0});
	double3 t2 = cross(n,t1); //second tangent

	return sin_theta*cos(psi)*t1+sin_theta*sin(psi)*t2+cos_theta*n; // <- check in on this
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}


