
#include <stdlib.h>
#include <time.h>

#include "test/ThermostatTest.h"

SimpleParticleContainer * ThermostatTest::setUpParticleContainer()
{
	return new SimpleParticleContainer();
}

void ThermostatTest::setUp()
{
	srand( time(NULL) );
	randomSize = rand() % (RANDOM_SIZE_MAX - RANDOM_SIZE_MIN + 1) + RANDOM_SIZE_MIN;
	randomInitialT = (double) ( rand() % (RANDOM_TEMPERATURE_MAX - RANDOM_TEMPERATURE_MIN + 1) + RANDOM_TEMPERATURE_MIN );
	randomTargetT = (double) ( rand() % (RANDOM_TEMPERATURE_MAX - RANDOM_TEMPERATURE_MIN + 1) + RANDOM_TEMPERATURE_MIN );

	container = setUpParticleContainer();

	list<Particle> particles;

	utils::Vector<double,3> x (0.0);
	utils::Vector<double,3> v (PARTICLE_VELOCITY_STD);

	particle = Particle( x, v, MASS );
	Particle * testParticles = new Particle[randomSize];

	for ( int i=0; i<randomSize; i++ )
	{
		testParticles[i] = particle;
		particles.push_back(testParticles[i]);
	}

	container->addParticles(particles);

	thermostat = new Thermostat( *container, randomInitialT, DIMENSION );

}

void ThermostatTest::tearDown()
{
	delete container;
	delete thermostat;
}

void ThermostatTest::testInitialization()
{
	ParticleContainer::SingleList pList = container->getParticles();
	double size = (double) pList.size();

	double currentEnergy = 0;

	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = *i;

		currentEnergy += p.getM() / 2.0 * p.getV().innerProduct();
	}

	double currentT = ( currentEnergy * 2.0 ) / ( ((double) DIMENSION) * size * kB );

	CPPUNIT_ASSERT( currentT == randomInitialT );
}

void ThermostatTest::testRegulation()
{
	thermostat->regulateTemperature( randomTargetT );

	ParticleContainer::SingleList pList = container->getParticles();
	double size = (double) pList.size();

	double currentEnergy = 0;

	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = *i;

		currentEnergy += p.getM() / 2.0 * p.getV().innerProduct();
	}

	double currentT = ( currentEnergy * 2.0 ) / ( ((double) DIMENSION) * size * kB );

	CPPUNIT_ASSERT( currentT == randomTargetT );
}
