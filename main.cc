#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <numeric>
#include "vec.h"
#include "molecule.h"
#include "system.h"
#include "schedule.h"

using namespace std;

int main() {
	Vector3d x(0,0,0); Vector3d y(0,0,2); Vector3d z(0,2,0);
	sphere sphere1(x,1,1);
	sphere sphere2(y,1,1);
	sphere sphere3(z,1,1);

	sphere arr[] = {sphere1,sphere2,sphere3};
	vector<sphere> spheres(arr,arr+3);
	moleculeType triatomic(spheres);

	molecule mol1(triatomic,Vector3d(5,5,3));
	molecule mol2(triatomic,Vector3d(5,2,3));
	vector<molecule> molecules = {mol1, mol2};
	


	double dim_vals[] = {0,5,0,5,0,5};
	vector<double> dim(dim_vals,dim_vals+6);

	double vdim_vals[] = {-5,10,-5,10,-5,10};
	vector<double> vdim(vdim_vals,vdim_vals+6);

	dimensions sysDim(dim,vdim);

	IO files("smaller.dat");

	parameters param;
	param.k_n = 2e8; param.k_t = 1500; param.k_static = 1800, param.cutoff = 0.1;
	param.freq = 10;

	timeline myTL(.5);
	msystem mysys(molecules, param, files, sysDim, myTL);

	mysys.run(2.5e-6);
}
