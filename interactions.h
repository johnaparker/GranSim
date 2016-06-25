#ifndef GUARD_interactions_h
#define GUARD_interactions_h

#include <vector>
#include <eigen3/Eigen/Geometry>
#include "vec.h"
#include "molecule.h"
#include "system.h"

void interact(sphere&, sphere&, bool tanForce, bool staticFriction, msystem*);
void virtualInteract(sphere&, Vector3d, msystem*);
void wallInteract(sphere&,const std::vector<double>, msystem*);





#endif