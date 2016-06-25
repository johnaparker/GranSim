#ifndef GUARD_system_h
#define GUARD_system_h

#include <vector>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <unordered_map>
#include <mutex>
#include "vec.h"
#include "molecule.h"
#include "schedule.h"
#include "parallel.h"

#define CORES     8

typedef std::pair<sphere,sphere> contact;
typedef struct
{
  size_t operator() (const contact &k) const { 
  	return (*k.first.mol).num; //*(*k.second.mol).num + (*k.first.mol).num; 
  }
} pairHash;
 
typedef struct
{
  bool operator() (const contact &x, const contact &y) const { 
  	return ((*x.first.mol).num == (*y.first.mol).num) && ((*x.second.mol).num == (*y.second.mol).num); 
  }
} pairEquals;
typedef std::unordered_map<contact,Vector3d,pairHash,pairEquals> map;

struct dimensions {
	std::vector<double> dim, vdim;
	dimensions(std::vector<double> dim, std::vector<double> vdim = std::vector<double>()):dim(dim),vdim(vdim) {};
};

struct IO {
	const char* dIn; 
	const char* dOut;
	const char* info;
	IO(const char* dIn, const char* dOut = NULL):dIn(dIn),dOut(dOut) {};
};

struct parameters {
	double k_static, k_n, k_t, cutoff;
	int num_molecules, num_spheres;
	double freq;
	std::vector<double> dim;
	double time;

	void write_data(const char* filename) {
		double dx = dim[1] - dim[0], dy = dim[3] -dim[2], dz = dim[5] - dim[4];
		std::ofstream myfile(filename);
		myfile << "k_n:\t\t" << k_n << "\n";
		myfile << "k_t:\t\t" << k_t << "\n"; 
		myfile << "k_static:\t\t" << k_static << "\n"; 
		myfile << "Cutoff:\t\t" << cutoff << "\n\n"; 
		myfile << "num_molecules:\t\t" << num_molecules << "\n";
		myfile << "num_spheres:\t\t" << num_molecules << "\n\n";
		myfile << "Dimensions:\t\t" << dx << "x" << dy << "x" << dz << "\n";
		myfile << "Frequency:\t\t" << freq << "\n";
		myfile << "Angle Range:\t\t" << freq*time << "\n";
		myfile.close();
	}
};

class msystem {
public:
	std::vector<molecule>* molecules;
	map contacts;
	std::mutex contacts_mutex;
	parameters params;
	std::vector<double> dim;
	std::vector<double> vdim;
	IO files;
	const char* dataIn;
	const char* dataOut;
	const char* dataOutDir;
	std::ofstream dataFile;
	std::vector<molecule>::size_type size;
	double time;
	timeline agenda;
	std::vector<std::vector<sphere*> > grid;
	double dl, timeStep;
	int nx,ny,nz,grid_sz;
	Vector3d gravity;

	int sz_th, sz_main;
	pthread_t threads[CORES-1];
	struct parallel_data td[CORES-1];

public:
	msystem(std::vector<molecule> &,parameters,IO, dimensions, timeline);
	molecule operator[] (int i) {molecules[i];}
	double update_dl();
	void createGrid();
	void updateGrid(molecule&);
	void interaction(molecule&);
	void exportData();
	void locateOutput();
	void load();
	void update(double dt);
	void run(double dt);

	int toCell(int i, int j, int k);
	int toCell(sphere sp);


};
#endif