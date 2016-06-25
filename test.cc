#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>
#include <sstream>

using namespace std;

string getMonthYear() {
	time_t t = time(0);   
    struct tm * now = localtime( & t );

    int year = now->tm_year + 1900;  
    int month = now->tm_mon + 1;
    int day = now->tm_mday;

    static const string months[] = {"jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"};
    const string str_month = months[month-1];
    string str_year = static_cast<ostringstream*>( &(ostringstream() << year) )->str();


    string monthyear = str_month + str_year;
    return monthyear;
}

int main() {
	string x = getMonthYear();
	cout << "hi" + x << endl;
}