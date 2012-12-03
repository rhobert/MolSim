

#ifndef SIMPLEPARTICLECONTAINER_H_
#define SIMPLEPARTICLECONTAINER_H_

#include "Particle.h"
#include "ParticleContainer.h"
#include <list>
#include <vector>

using namespace std;


/**
* @brief This class offers the possibility to to iterate over single particles or pairs of particles without boundary conditions.
**/
class SimpleParticleContainer: public ParticleContainer {
	
public:
	/**
	* @brief Type of a list of single particles
	**/
	typedef vector<Particle> SingleList;
	
	/**
	* @brief Type of a list of particles pairs
	**/
	typedef list<pair<Particle*, Particle*> > PairList;
	
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
	* @param list List of single particles
	* 
	* @return List of particle pairs
	**/
	PairList createPairs(SingleList& list);	

public:
	/**
	* @brief Create an instance of SimpleParticleContainer constructed of a list of particles
	* 
	* @param list Particle list
	**/
	SimpleParticleContainer(list<Particle> list);
	
	/**
	* @brief Iterate over all sinlge particles
	* 
	* @param singleFunction function to apply on all single particles
	**/
	void applyToSingleParticles( void (*singleFunction)(Particle&) );

	/**
	* @brief  Iterate over all particle pairs
	* 
	* @param pairFunction function to apply on all particle pairs
	**/
	void applyToParticlePairs( void (*pairFunction)(Particle&, Particle&) );
	
	/**
	* @brief Count of all single particles
	* 
	* @return Particle count
	**/
	int size();
};

#endif /* SIMPLEPARTICLECONTAINER_H_ */


