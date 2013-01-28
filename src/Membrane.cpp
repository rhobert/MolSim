
#include "Membrane.h"

#define SQRT2 1.414213562
#define SQRT3 1.732050808

using namespace std;

Membrane::Membrane( double stiffnessConstant, double averageBondLength )
{
	this->stiffnessConstant = stiffnessConstant;
	this->averageBondLength = averageBondLength;
	this->averageBondLengthTyped[0] = averageBondLength;
	this->averageBondLengthTyped[1] = averageBondLength * SQRT2;
	this->averageBondLengthTyped[2] = averageBondLength * SQRT3;
}