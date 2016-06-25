#ifndef GUARD_schedule_h
#define GUARD_schedule_h

#include <vector>

class event {
public:
	double ti,tf,duration;
	void (*fp)(void);
public:
	event(void (*fp)(void), double tf, double ti = 0);
	void execute();
};
bool compEvents(event,event);

class timeline {
public:
	double ti,tf,t,duration;
	std::vector<event> future_events, live_events;
public:
	timeline(double tf, double ti = 0);
	timeline(std::vector<event>, double tf, double ti = 0);
	void append(event);
	void update_events();
	void step(double dt);
	void terminate();
	event operator[](int i) {return live_events[i];}
};


#endif