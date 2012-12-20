/*
 * ParticleContainerTest.h
 *
 *  Created on: 09.11.2012
 *      Author: chris
 */


#include "../particleContainer/ParticleContainer.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestFixture.h>

#ifndef PARTICLECONTAINERTEST_H_
#define PARTICLECONTAINERTEST_H_

#define RANDOM_SIZE_MAX 256
#define RANDOM_SIZE_MIN 16

#define PARTICLE_POSITION_STD 1.
#define PARTICLE_VELOCITY_STD 2.

#define PARTICLE_POSITION_ERR 0.

/**
* @class ParticleContainerTest
*
* @brief Tests for class ParticleContainer
*/
class ParticleContainerTest : public CppUnit::TestFixture 
{
protected:
	
	/**
	* @brief Standard Particle used by the tests
	*/
	Particle particle;
	
	/**
	* @brief ParticleContainer with no particles
	*/
	ParticleContainer * emptyContainer;
	
	/**
	* @brief ParticleContainer with one particle
	*/
	ParticleContainer * oneContainer;
	
	/**
	* @brief ParticleContainer with two particles
	*/
	ParticleContainer * twoContainer;
	
	/**
	* @brief ParticleContainer with a random number of particles
	*/
	ParticleContainer * randomContainer;
		
	/**
	* @brief Count of particles in the randomContainer (between RANDSIZEMIN and RANDSIZEMAX)
	*/
	int randomSize;
	
	/**
	* @brief Set up an empty ParticleContainer
	* 
	* @return Pointer to an empty ParticleContainer
	*/
	virtual ParticleContainer * setUpParticleContainer() = 0;
	
	/**
	* @brief Count particles returned  by ParticleContainer::getParticles()
	* 
	* @param container ParticleContainer which particles are counted
	* 
	* @return Count of particles
	*/
	int countParticles(ParticleContainer * container);
	
	/**
	* @brief Count particle pairs
	* 
	* @param container ParticleContainer which particle pairs are counted
	* 
	* @return Count of particle pairs
	*/
	int countParticlePairs(ParticleContainer * container);
	
	/**
	* @brief Checks if all particles of a ParticleContainer have the right values
	* 
	* @param container ParticleContainer which particles should be checked
	*/
	bool checkParticles(ParticleContainer * container);
	
	/**
	* @brief Checks if all particles of a ParticleContainer were modified
	* 
	* @param container ParticleContainer which particles should been modified
	*/
	bool checkModfication(ParticleContainer * container);
	
	/**
	* @brief Modifies the values of a particle
	* 
	* @param p Particle to modify
	*/
	static void modifyParticle(Particle & p);
	
	/**
	* @brief Modifies the values of not distinct particles
	* 
	* @param p1 First Particle of pair which to modify if not distinct to second Particle
	* @param p2 Second Particle of pair which to modify if not distinct to first Particle
	*/
	static void modifyNotDistinctParticlePair(Particle & p1, Particle & p2);
	
	/**
	* @brief Increments the velocity for particles in a particle pair 
	* 
	* @param p1 First Particle to modify
	* @param p2 Second Particle  to modify
	*/
	static void incrementParticlePair(Particle & p1, Particle & p2);
	
public:
	
	/**
	* @brief Set up a test by creating the ParticleContainers and fill them with particles
	*/	
	void setUp();
	
	/**
	* @brief Delete test-variables
	*/
	void tearDown();
	
	/**
	* @brief Tests if ParticleContainer has the right number of particles
	*/
	void testSize();
	
	/**
	* @brief Tests if ParticleContainer returns the right count of particles
	*/
	void testParticleCount();
	
	/**
	* @brief Tests if particles of ParticleContainer have the right values
	*/
	void testParticleValues();
	
	/**
	 * @brief Tests if functions are applied on all particle
	 */
	void testApplyToSingleParticles();
	
	/**
	 * @brief Tests if all particle pairs are distinct
	 */
	void testDistinctParticlePairs();
	
	/**
	* @brief Tests if ParticleContainer returns the right count of particle pairs
	*/
	void testParticlePairCount();

};


#endif /* PARTICLECONTAINERTEST_H_ */
