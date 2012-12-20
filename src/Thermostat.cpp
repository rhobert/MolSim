
#include "Thermostat.h"
#include "MaxwellBoltzmannDistribution.h"

log4cxx::LoggerPtr Thermostat::logger(log4cxx::Logger::getLogger("Thermostat"));

Thermostat::Thermostat( ParticleContainer& pc, double initialT, int dimensionCount )
{
	this->pc = &pc;
	this->dimensionCount = dimensionCount;
	
	Thermostat::meanVelocity = initializeTemperature( initialT );
	Thermostat::dimensions = dimensionCount;
	
	LOG4CXX_INFO(logger, "Mean velocity for Maxwell Boltzmann Distribution is " << Thermostat::meanVelocity );
	
	this->pc->applyToSingleParticles( Thermostat::applyMaxwellBoltzmannDistribution );
}

double Thermostat::initializeTemperature( double initialT )
{
	return sqrt( kB * initialT / thermostat_mass );
}

void Thermostat::regulateTemperature( double targetT )
{
	ParticleContainer::SingleList pList = pc->getParticles();
	double size = (double) pList.size();
	
	currentEnergy = 0;
	
	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = *i;

		currentEnergy += p.getM() / 2.0 * p.getV().innerProduct();
	}
	
	currentT = ( currentEnergy * 2.0 ) / ( ((double) dimensionCount) * size * kB );
	
	beta = sqrt( ( currentT + diffT ) / currentT );

	for (ParticleContainer::SingleList::iterator i = pList.begin(); i != pList.end(); i++)
	{
		Particle& p = *i;

		p.setV( p.getV() * beta );
	}
}

double Thermostat::getEnergy()
{
	return currentEnergy;
}
	
double Thermostat::getTemperature()
{
	return currentT;
}

double Thermostat::meanVelocity = 0.1;
int Thermostat::dimensions = 3;

void Thermostat::applyMaxwellBoltzmannDistribution( Particle& p )
{
	MaxwellBoltzmannDistribution(p, Thermostat::meanVelocity, Thermostat::dimensions);
}