#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "vec.h"
#include "molecule.h"
#include "schedule.h"
#include "interactions.h"
#include "system.h"
#include "parallel.h"
#include "functions.h"
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Core>
#include <unistd.h>
#include <cstdlib>
#include <pthread.h>
#include <string>

using namespace std;
using std::vector;
using std::pow;

//TURN ROTATIONS INTO EVENTS

msystem::msystem(vector<molecule>& mols,parameters params, IO files, 
				dimensions sysDim, timeline agenda):agenda(agenda),files(files),params(params){
	
	molecules = &mols;
	dataIn = files.dIn;
	load();
	size = molecules->size();

	params.num_molecules = size;
	params.num_spheres = size;
	params.dim = sysDim.dim;
	params.time = agenda.tf;
	locateOutput();
	params.write_data(this->files.info);

	time = agenda.ti;
	update_dl();
	dim = sysDim.dim;  vdim = sysDim.vdim;
	createGrid();
	gravity.set(0,0,-981);

	//CREATE A VOID POINTER IN AGENDA WHICH WILL BE ASSIGNED TO THE SYSTEM HERE

	sz_th = size/CORES;
	sz_main = size - sz_th*(CORES-1);
	for (int i = 0; i != CORES-1; ++i) {
		td[i].sys = this;
		td[i].num = sz_main + sz_th*i;
		td[i].size = sz_th;
	}
}

int msystem::toCell(int i, int j, int k){
	return k*(nx*ny) + j*(nx) + i;
}
int msystem::toCell(sphere sp) {
	int i = (sp[0]-vdim[0])/dl, j = (sp[1]-vdim[2])/dl, k = (sp[2]-vdim[4])/dl;
	return toCell(i,j,k);
}

double msystem::update_dl() {
	double rmax = 0;
	for (vector<molecule>::const_iterator mol = (*molecules).begin(); mol != (*molecules).end(); ++mol){
		for (vector<sphere>::const_iterator sp = (mol->spheres).begin(); sp != (mol->spheres).end(); ++sp) {
			rmax = max(rmax,sp->r);
		}
	}
	dl = rmax*2;
	return dl;
}

void msystem::createGrid() {
	if (vdim.size() == 0) {
		vdim.push_back(dim[0]-2*dl); vdim.push_back(dim[1]+2*dl); vdim.push_back(dim[2]-2*dl); 
		vdim.push_back(dim[3]+2*dl); vdim.push_back(dim[4]-2*dl); vdim.push_back(dim[5]+2*dl);
	}

	nx = (vdim[1]-vdim[0])/dl + 1;
	ny = (vdim[3]-vdim[2])/dl + 1;
	nz = (vdim[5]-vdim[4])/dl + 1;
	grid_sz = nx*ny*nz;

	for (int i = 0; i != grid_sz; ++i) {
		vector<sphere*> x;
		grid.push_back(x);
	}

	int count = 0;
	for (vector<molecule>::iterator mol = (*molecules).begin(); mol != (*molecules).end(); ++mol){
		for (vector<sphere>::iterator sp = (mol->spheres).begin(); sp != (mol->spheres).end(); ++sp) {
			int cell_num = toCell(*sp);
			sp->cell_num = cell_num;
			sphere* sp_p = &(*sp);
			grid[cell_num].push_back(sp_p);
			sp-> mol = &(*mol);
		}
		mol->num = count;
		count ++;
	}
}
void msystem::updateGrid(molecule& mol) {
	for (vector<sphere>::iterator sp = (mol.spheres).begin(); sp != (mol.spheres).end(); ++sp) {
		int cell_num = toCell(*sp);
		if (cell_num == sp->cell_num) continue;
		else {
			sphere* sp_p = &(*sp);
			grid[cell_num].push_back(sp_p);
			for (vector<sphere*>::iterator sp_p_o = grid[sp->cell_num].begin(); sp_p_o != grid[sp->cell_num].end(); ++sp_p_o) {
				if(*sp_p_o == sp_p) {
					grid[sp->cell_num].erase(sp_p_o);
					sp->cell_num = cell_num;
					break;
				}
			}
		}	
	}
}

void msystem::interaction(molecule& mol) {
	for (vector<sphere>::iterator sp = (mol.spheres).begin(); sp != (mol.spheres).end(); ++sp) {
		int x = (sp->pos[0]-vdim[0])/dl, y = (sp->pos[1]-vdim[2])/dl, z = (sp->pos[2]-vdim[4])/dl;
		for (int i = -1; i != 2; ++i) {
			for (int j = -1; j != 2; ++j){
				int cell_num = toCell(x+i,y+j,z+1);
				for (vector<sphere*>::iterator sp_p = grid[cell_num].begin(); sp_p != grid[cell_num].end(); ++sp_p) {
					interact(*sp, **sp_p,true,true,this);
				}
			}
		}
		for (int i = 0; i != 2; ++i) {
			for (int j = -1; j != 2; ++j){
				int cell_num = toCell(x+i,y+j,z);
				if (i==0 && j ==-1) continue;
				for (vector<sphere*>::iterator sp_p = grid[cell_num].begin(); sp_p != grid[cell_num].end(); ++sp_p) {
					interact(*sp, **sp_p,true,true,this);
				}
			}
		}
		mol.push(gravity*sp->m, sp->pos);
		wallInteract(*sp, dim,this);
	}
}

void msystem::load() {
	ifstream inFile(dataIn, ios::binary);
	float x;
	inFile.read((char*)&x, sizeof(float));
	vector<float> temp;
	while (inFile) {
		inFile.read((char*)&x, sizeof(float));
		temp.push_back(x);
		if (temp.size() == 4) {
			double r = temp[0], x = temp[1], y = temp[2], z = temp[3];
			double m = pow(r/(.321/2.0),3)*.03555;
			Vector3d pos(x,y,z);
				Vector3d zeros(0,0,0);
				sphere sphereX(zeros,r, m);
				vector<sphere> spheresX; spheresX.push_back(sphereX);
				moleculeType* mType = new moleculeType(spheresX);
				molecule mol(*mType, pos);
			// sphere sp(pos,r,m);
			// molecule mol = makeComposite(sp);
			molecules->push_back(mol);
			temp.clear();
		}
	}
	inFile.close();
	cout << molecules->size() << " molecules loaded" << endl;
}

void msystem::locateOutput() {
	if (files.dOut) {
		dataOut = files.dOut;
		dataOutDir = "./";
	}
	else {
		string date = getMonthYear();

		string dataOutDir_str = date;
		char* a1 = new char[dataOutDir_str.length() + 1];
		strcpy(a1,dataOutDir_str.c_str());

		string systemCall = ("mkdir " + dataOutDir_str + " 2> /dev/null");
		char* a2 = new char[systemCall.length() + 1];
		strcpy(a2,systemCall.c_str());
		system(a2);

		string fileName = getFileName(a1);
		string dataOut_str = dataOutDir_str + fileName + fileName + ".dat";
		char* a3 = new char[dataOut_str.length() + 1];
		strcpy(a3,dataOut_str.c_str());
		dataOut = a3;

		string info_str = dataOutDir_str + fileName + "/info.txt";
		char* a4 = new char[info_str.length() + 1];
		strcpy(a4,info_str.c_str());
		files.info = a4;

		dataOutDir_str += fileName;
		char* a5 = new char[dataOutDir_str.length() + 1];
		strcpy(a5,dataOutDir_str.c_str());
		dataOutDir = a5;

		delete [] a1, a2;
	}
	dataFile.open(dataOut, ios::binary);
	
}

void msystem::exportData() {
	for (vector<molecule>::iterator mol = molecules->begin(); mol != molecules->end(); ++mol) {
		for (vector<sphere>::iterator sp = (*mol).spheres.begin(); sp != (*mol).spheres.end(); ++sp) {
			float r = sp->r;
			float x = sp->pos[0], y = sp->pos[1], z = sp->pos[2];
			dataFile.write(reinterpret_cast<const char *>(&r), sizeof(float));
			dataFile.write(reinterpret_cast<const char *>(&x), sizeof(float));
			dataFile.write(reinterpret_cast<const char *>(&y), sizeof(float));
			dataFile.write(reinterpret_cast<const char *>(&z), sizeof(float));
		}
	}
	cout << "Exporting Data...\t" << time << endl; 
}

void msystem::update(double dt) {
	void *status;

	for (int i = 0; i != CORES-1; ++i) {
		pthread_create(&threads[i],NULL,execThread1,(void*)&td[i]);
	}
	for (vector<molecule>::iterator mol = molecules->begin(); mol != molecules->begin() + sz_main; ++mol) {
		mol->predict(dt);
		updateGrid(*mol);
	}
	for (int i = 0; i != CORES-1; ++i) {
		pthread_join(threads[i], &status);
	}

	for (int i = 0; i != CORES-1; ++i) {
		pthread_create(&threads[i],NULL,execThread2,(void*)&td[i]);
	}
	for (vector<molecule>::iterator mol = molecules->begin(); mol != molecules->begin() + sz_main; ++mol) {
		interaction(*mol);
	}
	for (int i = 0; i != CORES-1; ++i) {
		pthread_join(threads[i], &status);
	}

	for (int i = 0; i != CORES-1; ++i) {
		pthread_create(&threads[i],NULL,execThread3,(void*)&td[i]);
	}
	for (vector<molecule>::iterator mol = molecules->begin(); mol != molecules->begin() + sz_main; ++mol) {
		mol->correct(dt);
		mol->clearForces();
	}
	for (int i = 0; i != CORES-1; ++i) {
		pthread_join(threads[i], &status);
	}

	for (map::iterator i = contacts.begin(); i != contacts.end(); ++i) {
		const sphere& sp1 = (i->first).first; const sphere& sp2 = (i->first).second;
		Vector3d dr = sp1.pos - sp2.pos;
		double overlap = pow(dr.dot(dr),.5) - (sp1.r + sp2.r);
		if (overlap > 0) {
			i = contacts.erase(i);
			//fix i pointer --i
		}
	}
}


void msystem::run(double dt) {
	double save_time = .01-dt/2.0; double lap = save_time + 1;
	double omega = params.freq*M_PI/180.0;
	for (int i = 0; i != CORES-1; ++i) {
		td[i].timeStep = dt;
	}

	timeStep = dt;
	std::ofstream myfile(files.info, ios::app);
	myfile << "dt:\t\t" << dt << "\n";
	myfile.close();

	while (time < agenda.tf) {
		//cout << "Time:   " << time << endl;
		if (lap > save_time) {
			exportData();
			lap = 0;
		}
		agenda.step(dt);
		update(dt);
		time += dt;
		//gravity.set(0,-981*sin(omega*time),-981*cos(omega*time));
		lap += dt;
	}
	dataFile.close();
}
