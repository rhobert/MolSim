
#ifndef MEMBRANE_H_
#define MEMBRANE_H_

using namespace std;

class Membrane
{
	double stiffnessConstant;
	double averageBondLength;
	double averageBondLengthTyped[3];
	
public:

	Membrane( double stiffnessConstant, double averageBondLength );
	
	friend class MembraneParticle;
};

#endif /* MEMBRANE_H_ */


