
#include "ParticleContainerTest.h"
#include "particleContainer/LinkedCellParticleContainer.h"

#ifndef LINKEDCELLPARTICLECONTAINERTEST_H_
#define LINKEDCELLPARTICLECONTAINERTEST_H_

#define DOMAIN_SIZE_STD 1.
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
	CPPUNIT_TEST_SUITE_END();

protected:
	
	LinkedCellParticleContainer * setUpParticleContainer();

public:

	
};


#endif /* LINKEDCELLPARTICLECONTAINERTEST_H_ */
