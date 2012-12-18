
#include "ParticleContainerTest.h"
#include "particleContainer/LinkedCellParticleContainer.h"

#ifndef LINKEDCELLPARTICLECONTAINERTEST_H_
#define LINKEDCELLPARTICLECONTAINERTEST_H_

#define DOMAIN_SIZE_STD 5.
#define CUTOFF_RADIUS_STD 1.

/**
 * @class LinkedCellParticleContainerTest
 *
 * @brief test class for LinkedCellParticleContainerTest
 */
class LinkedCellParticleContainerTest : public ParticleContainerTest
{
	CPPUNIT_TEST_SUITE(LinkedCellParticleContainerTest);

	CPPUNIT_TEST(testSize);
	CPPUNIT_TEST(testParticleCount);
	CPPUNIT_TEST(testParticleValues);
	CPPUNIT_TEST(testApplyToSingleParticles);
	CPPUNIT_TEST(testDistinctParticlePairs);
	CPPUNIT_TEST(testParticlePairCount);
/*	
	CPPUNIT_TEST(testDeleteHaloParticles);
	CPPUNIT_TEST(testApplyToBoundaryParticles);
	CPPUNIT_TEST(testCutOffRadius);
*/
	CPPUNIT_TEST_SUITE_END();
	
	LinkedCellParticleContainer * container;
	Particle innerParticle, boundaryParticle, haloParticle;
	
	/**
	 * @brief Move the particle to a halo cell
	 * 
	 * @param p Particle to move to halo
	 */
	static void moveToHalo( Particle& p );
	
	/**
	 * @brief  Move the particle pair to a halo cell
	 * 
	 * @param p1 Particle to move to halo
	 * @param p2 Particle to move to halo
	 */
	static void moveToHalo( Particle& p1, Particle& p2 );
	
protected:
	
	LinkedCellParticleContainer * setUpParticleContainer();

public:
	
	void setUp();
	
	/**
	 * @brief Tests if all, and only, halo particle are deleted
	 */
//	void testDeleteHaloParticles();
	
	/**
	 * @brief Tests if all, and only, boundary particles are modified
	 */
//	void testApplyToBoundaryParticles();
	
	/**
	 * @brief Tests if only particles with distance less equal to the the cutoff radius are pairs
	 */
//	void testCutOffRadius();
	
};


#endif /* LINKEDCELLPARTICLECONTAINERTEST_H_ */
