# all:
# 	g++ -std=c++11 molecule.cc schedule.cc system.cc interactions.cc vec.cc parallel.cc functions.cc -lpthread -o gran

gran: functions.o interactions.o molecule.o parallel.o schedule.o system.o vec.o main.o
	g++ -std=c++11 functions.o interactions.o molecule.o parallel.o schedule.o system.o vec.o main.o -lpthread -o gran

molecule.o: molecule.cc molecule.h vec.o
	g++ -std=c++11 -c molecule.cc

functions.o: functions.cc functions.h molecule.o vec.o
	g++ -std=c++11 -c functions.cc

interactions.o: interactions.cc interactions.h vec.o molecule.o system.o
	g++ -std=c++11 -c interactions.cc

parallel.o: parallel.cc parallel.h molecule.o system.o
	g++ -std=c++11 -c parallel.cc

schedule.o: schedule.cc schedule.h
	g++ -std=c++11 -c schedule.cc

system.o: system.cc system.h vec.o molecule.o schedule.o interactions.o parallel.o functions.o
	g++ -std=c++11 -c system.cc -lpthread

vec.o: vec.cc vec.h
	g++ -std=c++11 -c vec.cc

main.o: main.cc vec.o molecule.o system.o schedule.o
	g++ -std=c++11 -c main.cc