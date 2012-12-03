

#ifndef PARTICLECONTAINER_H_
#define PARTICLECONTAINER_H_

#include "Particle.h"
#include <list>
#include <vector>

using namespace std;


/**
* @brief This class offers the possibility to to iterate over single particles or pairs of particles.
**/
class ParticleContainer {

private:
	
	/**
	 * @brief Logger for ParticleContainer class
	 */
	static log4cxx::LoggerPtr logger;

public:
	
	/**
	* @brief Destructor of the class ParticleContainer
	*/
	virtual ~ParticleContainer();
	
	/**
	* @brief Iterate over all sinlge particles
	* 
	* @param singleFunction function to apply on all single particles
	**/
	virtual void applyToSingleParticles( void (*singleFunction)(Particle&) ) = 0;

		/**
	* @brief Iterate over all sinlge particles
	* 
	* @param singleFunction function to apply on all particle pairs
	**/
	virtual void applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) ) = 0;
	
	/**
	* @brief Count of all single particles
	* 
	* @return Particle count
	**/
	virtual int size() = 0;
};


#endif /* PARTICLECONTAINER_H_ */


