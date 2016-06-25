#include <iostream>
#include <vector>
#include <algorithm>
#include "schedule.h"

using namespace std;
using std::vector;
using std::sort;

//COPY CONSTRUCTOR

event::event(void (*f)(void), double tf, double ti):tf(tf),ti(ti){
	fp = f;
	duration = tf - ti;
}
void event::execute() {
	fp();
}

bool compEvents(event e1,event e2) {
	return (e2.ti > e1.ti);
}

timeline::timeline(vector<event> events, double tf, double ti):ti(ti),tf(tf) {
	t = ti;
	duration = tf - ti;
	sort(events.begin(),events.end(),compEvents);
	future_events =  events;
	update_events();
}
timeline::timeline(double tf, double ti):ti(ti),tf(tf) {
	t = ti;
	duration = tf - ti;
}
void timeline::append(event e){
	future_events.push_back(e);
}
void timeline::update_events(){
	vector<event> toAdd;
	for (vector<event>::iterator ev = future_events.begin(); ev!= future_events.end(); ++ev) {
		if (t >= ev->ti) {
			live_events.push_back(*ev);
			ev = future_events.erase(ev);
		}
		else if (t >= ev->tf) {
			live_events.erase(ev);
		}
	}
}
void timeline::step(double dt){
	t += dt;
	update_events();
	for (vector<event>::iterator ev = live_events.begin(); ev!= live_events.end(); ++ev) {
		ev->execute();
	}
}
