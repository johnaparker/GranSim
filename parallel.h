#ifndef GUARD_parallel_h
#define GUARD_parallel_h

class msystem;
struct parallel_data {
	msystem* sys;
	int num;
	int size;
	double timeStep;
};

void* execThread1(void*);
void* execThread2(void*);
void* execThread3(void*);



#endif