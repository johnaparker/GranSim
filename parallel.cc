#include "molecule.h"
#include "system.h"
#include "parallel.h"
using namespace std;

void* execThread1(void* n) {
	struct parallel_data *my_data;
	my_data = (struct parallel_data *) n;
	double dt = my_data->timeStep;
	for (int i = 0; i != my_data->size; ++i){
		((*(my_data->sys->molecules))[i+my_data->num]).predict(dt);
		my_data->sys->updateGrid(((*(my_data->sys->molecules))[i+my_data->num]));
	}
}
void* execThread2(void* n) {
	struct parallel_data *my_data;
	my_data = (struct parallel_data *) n;
	double dt = my_data->timeStep;
	for (int i = 0; i != my_data->size; ++i){
		my_data->sys->interaction(((*(my_data->sys->molecules))[i+my_data->num]));
	}
}
void* execThread3(void* n) {
	struct parallel_data *my_data;
	my_data = (struct parallel_data *) n;
	double dt = my_data->timeStep;
	for (int i = 0; i != my_data->size; ++i){
		((*(my_data->sys->molecules))[i+my_data->num]).correct(dt);
		((*(my_data->sys->molecules))[i+my_data->num]).clearForces();
	}
}