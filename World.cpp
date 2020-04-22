#include <iostream>
#include "World.h"


/*returns true if point x is inside or on the sphere*/
bool World::inSphere(const double3 &x)
{
	double3 r = x-sphere_x0;	//ray to x

    double r_mag2 = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if (r_mag2<=(sphere_r*sphere_r)) return true;
    else return false;
}

/*returns true if point x is at the wall*/
bool World::onEdge(const double3& x) {
	if (x[1] == xd[1] || x[1] == 0) { return true; }
	if (x[2] == xd[2] || x[2] == 0) { return true; }
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

	return sin_theta*cos(psi)*t1+sin_theta*sin(psi)*t2+cos_theta*n;
}

/* 
	wallDiffuseVector is a duplicate of sphereDiffuseVector
	but its position is at the x,y extema
*/
double3 World::wallDiffuseVector(const double3& x)
{
	//pick angles, theta=off normal, psi = azimuthal rotation
	double sin_theta = rnd();
	double cos_theta = sqrt(1 - sin_theta * sin_theta);
	double psi = 2 * Const::PI * rnd();

	double3 n = unit(x - sphere_x0);	// normal vector
	double3 t1;						// create the first tangent
	if (dot(n, { 1,0,0 }) != 0) t1 = cross(n, { 1,0,0 });
	else t1 = cross(n, { 0,1,0 });
	double3 t2 = cross(n, t1); //second tangent

	return sin_theta * cos(psi) * t1 + sin_theta * sin(psi) * t2 + cos_theta * n;
}

/* Can I get an explanation of what this does?*/
double World::lineSphereIntersect(const double3 &x1, const double3 &x2){
	// What are x1 and x2?
	double3 B = x2-x1;
	double3 A = x1-sphere_x0;
	double a = dot(B,B);
	double b = 2*dot(A,B);
	double c = dot(A,A)-(sphere_r*sphere_r);
	double det = b*b-4*a*c;
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

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}


