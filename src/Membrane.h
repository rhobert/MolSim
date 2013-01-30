
#ifndef MEMBRANE_H_
#define MEMBRANE_H_

using namespace std;

/**
 * @brief This class represents a Membrane
 */
class Membrane
{
	/**
	* @brief The stiffness constant of the membrane
	*/
	double stiffnessConstant;
	
	/**
	* @brief The average bond length of the membrane
	*/
	double averageBondLength;
	
	/**
	* @brief The average bond length of the membrane adapted to the type of neighbourhood (direct: 0, diagonal: 1, diagonal over two dimensions: 2)
	*/
	double averageBondLengthTyped[3];
	
public:

	/**
	* @brief Create an instance of Membrane
	* 
	* @param stiffnessConstant The stiffness constant of the membrane
	* @param averageBondLength The average bond length of the membrane
	*/
	Membrane( double stiffnessConstant, double averageBondLength );
	
	friend class MembraneParticle;
};

#endif /* MEMBRANE_H_ */


