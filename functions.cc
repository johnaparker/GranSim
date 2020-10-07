#include "molecule.h"
#include "functions.h"
#include "vec.h"
#include <iomanip>
#include <ctime>
#include <stdlib.h> 
#include <vector>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Core>
#include <string>
#include <sstream>
#include <cstdlib>
#include <fstream>

using namespace std;

molecule makeComposite(const sphere& sp) {
	double r = sp.r; Vector3d pos = sp.pos; double m = sp.m;
	double rout = r/4.0;
	double rin = 2*rout;
	double V = 8*pow(rout,3) + pow(rin,3);
	vector<sphere> spheres;
	
	Vector3d origin(0,0,0);
	spheres.push_back(sphere(origin,rin, m*pow(rin,3)/V));
	double phis[] = {M_PI/2.0 - M_PI/5.0, M_PI/2.0 + M_PI/5.0};
	double thetas[] = {M_PI/4.0, 3.0*M_PI/4.0, 5.0*M_PI/4.0, 7.0*M_PI/4.0};
	for (int i = 0; i != 2; ++i) {
		for (int j = 0; j != 4; ++j) {
			double phi = phis[i], theta = thetas[j];
			double x = (rin + rout)*sin(phi)*cos(theta);
			double y = (rin + rout)*sin(phi)*sin(theta);
			double z = (rin + rout)*cos(phi);
			origin = Vector3d(x,y,z);
			spheres.push_back(sphere(origin,rout, m*pow(rout,3)/V));
		}
	}
	moleculeType* comp = new moleculeType(spheres);
	Eigen::Vector3d rando = Eigen::Vector3d::Random(); Vector3d rando2 = Vector3d(rando);
	molecule mol(*comp, pos, Vector3d(), rando2);
	return mol;
}


double getMax(vector<molecule>& molecules) {
	double zmax = 0;
	for (vector<molecule>::iterator mol = molecules.begin(); mol != molecules.end(); ++mol) {
		for (vector<sphere>::iterator sp = (*mol).spheres.begin(); sp != (*mol).spheres.end(); ++sp) { 
			zmax = max(zmax, sp->pos[2]);
		}
	}
	return zmax;
}

string getMonthYear() {
	time_t t = time(0);   
    struct tm * now = localtime( & t );

    int year = now->tm_year + 1900;  
    int month = now->tm_mon + 1;
    int day = now->tm_mday;

    static const string months[] = {"jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"};
    const string str_month = months[month-1];
    string str_year = to_string(year);


    string monthyear = str_month + str_year;
    return monthyear;
}

string getFileName(const char* dirName) {
	string dirName_str = dirName;
	string systemCall = "ls -1 " + dirName_str + " | wc -l > ./TEMP.txt";
	char* a1 = new char[systemCall.length() + 1];
	strcpy(a1,systemCall.c_str());
	system(a1);

	ifstream myFile("./TEMP.txt");
	int n;
	myFile >> n;
	string n_str = to_string(n);
	system("rm ./TEMP.txt");

	systemCall = ("mkdir " + dirName_str + "/sim" + n_str + " > /dev/null");
	char* a2 = new char[systemCall.length() + 1];
	strcpy(a2,systemCall.c_str());
	system(a2);

	delete [] a1, a2;
	
	string x = "/sim" + n_str;
	return x;
}
