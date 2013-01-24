
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
	randomSize = rand() % (THERMO_RANDOM_SIZE_MAX - THERMO_RANDOM_SIZE_MIN + 1) + THERMO_RANDOM_SIZE_MIN;
	randomInitialT = (double) ( rand() % (RANDOM_TEMPERATURE_MAX - RANDOM_TEMPERATURE_MIN + 1) + RANDOM_TEMPERATURE_MIN );
	randomTargetT = (double) ( rand() % (RANDOM_TEMPERATURE_MAX - RANDOM_TEMPERATURE_MIN + 1) + RANDOM_TEMPERATURE_MIN );

	container = setUpParticleContainer();

	list<Particle*> particles;

	utils::Vector<double,3> x (0.0);
	utils::Vector<double,3> v (THERMO_PARTICLE_VELOCITY_STD);

	particle = Particle( x, v, THERMO_MASS );

	for ( int i=0; i<randomSize; i++ )
	{
		particles.push_back( new Particle(particle) );
	}

	container->addParticles(particles);

	thermostat = new Thermostat( *container, randomInitialT, THERMO_DIMENSION );

}

void ThermostatTest::tearDown()
{
	delete container;
	delete thermostat;
}

void ThermostatTest::testInitialization()
{	
	CPPUNIT_ASSERT( thermostat->getTemperature() > 0 );
}

void ThermostatTest::testRegulation()
{
	thermostat->regulateTemperature( randomTargetT );
	
	CPPUNIT_ASSERT( abs(thermostat->getTemperature() - randomTargetT) < 1e-3 );
}

void ThermostatTest::testTemperature()
{
	thermostat->regulateTemperature( randomTargetT );

	ParticleContainer::SingleList pList = container->getParticles();
	double size = (double) pList.size();

	double currentEnergy = 0;

	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = **i;

		currentEnergy += p.getM() / 2.0 * p.getV().innerProduct();
	}

	double currentT = ( currentEnergy * 2.0 ) / ( ((double) THERMO_DIMENSION) * size * kB );
	
	CPPUNIT_ASSERT( abs(currentT - randomTargetT) < 1e-3 );
}