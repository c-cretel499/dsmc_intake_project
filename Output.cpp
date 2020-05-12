#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "Output.h"
#include "Species.h"

using namespace std;

/*saves fields in VTK image data format*/
void Output::saveVTI(World &world, Species &sp)
{
	/*update gas macroscopic properties*/
	sp.computeGasProperties();

	stringstream name;
	name << "results/fields_" << setfill('0') << setw(5) << world.getTs() << ".vti";

    /*open output file*/
    ofstream out(name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	double3 x0 = world.getX0();
	double3 dh = world.getDh();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<"\" ";
	out<<"Spacing=\""<<dh[0]<<" "<<dh[1]<<" "<<dh[2]<<"\" ";
	out<<"WholeExtent=\"0 "<<world.ni-1<<" 0 "<<world.nj-1<<" 0 "<<world.nk-1<<"\">\n";

	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*species number densities*/
	out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<sp.den;
	out<<"</DataArray>\n";

	/*species stream velocity*/
	out<<"<DataArray Name=\"vel."<<sp.name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<sp.vel;
	out<<"</DataArray>\n";

	/*close out tags*/
	out<<"</PointData>\n";

	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
}
