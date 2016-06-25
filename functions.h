#ifndef GUARD_functions_h
#define GUARD_functions_h

#include "molecule.h"
#include <vector>
#include <string>

molecule makeComposite(const sphere& );
double getMax(std::vector<molecule>& );
std::string getMonthYear();
std::string getFileName(const char*);



#endif