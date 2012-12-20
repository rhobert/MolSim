
#include "ParticleContainerTest.h"
#include "particleContainer/SimpleParticleContainer.h"

#ifndef SIMPLEPARTICLECONTAINERTEST_H_
#define SIMPLEPARTICLECONTAINERTEST_H_

/**
 * @class SimpleParticleContainerTest
 *
 * @brief test class for SimpleParticleContainer
 */
class SimpleParticleContainerTest : public ParticleContainerTest
{
	CPPUNIT_TEST_SUITE(SimpleParticleContainerTest);
	CPPUNIT_TEST(testSize);
	CPPUNIT_TEST(testParticleCount);
	CPPUNIT_TEST(testParticleValues);
	CPPUNIT_TEST(testApplyToSingleParticles);
	CPPUNIT_TEST(testDistinctParticlePairs);
	CPPUNIT_TEST(testParticlePairCount);
	CPPUNIT_TEST_SUITE_END();

protected:
	
	SimpleParticleContainer * setUpParticleContainer();

public:

	
};


#endif /* SIMPLEPARTICLECONTAINERTEST_H_ */
