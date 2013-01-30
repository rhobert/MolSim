

#ifndef SIMPLEPARTICLECONTAINER_H_
#define SIMPLEPARTICLECONTAINER_H_

#include "../Particle.h"
#include "ParticleContainer.h"
#include <list>
#include <vector>

using namespace std;


/**
* @brief This class offers the possibility to to iterate over single particles or pairs of particles without boundary conditions.
**/
class SimpleParticleContainer: public ParticleContainer {
	
private:
	
	/**
	 * @brief Logger for ParticleContainer class
	 */
	static log4cxx::LoggerPtr logger;
	
	/**
	* @brief All single Particles
	**/
	SingleList singleList;
	
	/**
	* @brief All particles Pairs
	**/
	PairList pairList;
	
	/** 
	* @brief Creates a list of all distinct particles
	* 
	* @param sList List of single particles
	* 
	* @return List of particle pairs
	**/
	PairList createPairs(SingleList& sList);

public:
	
	/**
	 * @brief Create an instance of LinkedCellParticleContainer which contains no particles yet
	 */
	SimpleParticleContainer();
	
	void addParticles( SingleList pList );
	
	void applyToSingleParticles( void (*singleFunction)(Particle&) );
	
	void applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) );
	
	SingleList& getParticles();
	
	int size();
};

#endif /* SIMPLEPARTICLECONTAINER_H_ */


