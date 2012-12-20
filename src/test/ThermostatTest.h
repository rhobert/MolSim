
#include "Thermostat.h"
#include "particleContainer/SimpleParticleContainer.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestFixture.h>

#ifndef THERMOSTATTEST_H_
#define THERMOSTATTEST_H_

#define RANDOM_SIZE_MAX 256
#define RANDOM_SIZE_MIN 16

#define RANDOM_TEMPERATURE_MAX 200
#define RANDOM_TEMPERATURE_MIN 0

#define PARTICLE_VELOCITY_STD 2.
#define MASS 1.0

#define DIMENSION 3

/**
 * @class ThermostatTest
 *
 * @brief test class for Thermostat
 */
class ThermostatTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(ThermostatTest);
	CPPUNIT_TEST(testInitialization);
	CPPUNIT_TEST(testRegulation);
	CPPUNIT_TEST_SUITE_END();

/**
 * @brief Thermostat used by the tests
 */
	Thermostat * thermostat;

/**
 * @brief Simple particle container
 */
	SimpleParticleContainer * container;

/**
 * @brief Standard particle used used by tests
 */
	Particle particle;

/**
 * @brief A random temperature
 */
	double randomInitialT;

/**
 * @brief A random temperature
 */
	double randomTargetT;

/**
* @brief Count of particles in the container (between RANDSIZEMIN and RANDSIZEMAX)
*/
	int randomSize;

protected:

/**
* @brief Set up a ParticleContainer
* 
* @return Pointer to a ParticleContainer
*/
	SimpleParticleContainer * setUpParticleContainer();

public:

/**
* @brief Set up a test by creating a thermostat and particles
*/	
    void setUp();

/**
* @brief Tests the initialization to a given temperature
*/	
    void testInitialization();

/**
* @brief tests the rgulation to a given temperature
*/	
    void testRegulation();

/**
* @brief Delete test-variables
*/
	void tearDown();

};

#endif /* THERMOSTATTEST_H_ */
