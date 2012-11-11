

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
	* @brief Create an instance of ParticleContainer constructed of a list of particles
	* 
	* @param list Particle list
	**/
	ParticleContainer(list<Particle> list);
	
	/**
	* @brief Destructor of the class ParticleContainer
	*/
	virtual ~ParticleContainer();
	
	/**
	* @brief Begin iteration over all sinlge particles
	* 
	* @return Iterator which points on first single particle
	**/
	SingleList::iterator beginSingle();
	
	/**
	* @brief End iteration over all sinlge particles
	* 
	* @return Iterator which points after last single particle
	**/
	SingleList::iterator endSingle();

	/**
	* @brief Begin iteration over all particle pairs
	* 
	* @return Iterator which points on first particle pair
	**/
	PairList::iterator beginPair();

	/**
	* @brief End iteration over all particle pairs
	* 
	* @return Iterator which points after last particle pair
	**/
	PairList::iterator endPair();
	
	/**
	* @brief Count of all single particles
	* 
	* @return Particle count
	**/
	int size();
};

#endif /* PARTICLECONTAINER_H_ */


